"""Polygenic expression-specificity tester for the HMBA v0.5 human atlas."""

from __future__ import annotations

from functools import partial
from io import BytesIO
import os
import re
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "hmba_streamlit_matplotlib"),
)

import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
import streamlit as st

from analysis import (
    convert_to_human,
    load_ortholog_map,
    load_rank_matrix,
    parse_gene_list,
    run_auroc_analysis,
)


APP_DIR = Path(__file__).resolve().parent
DATA_DIR = APP_DIR / "data"
DATASETS = {
    "Subclass": ("WB_HMBA_Human_AIT_subclass_rank.csv.gz", "subclass"),
    "Subclass by brain region": (
        "WB_HMBA_Human_AIT_brain_region_subclass_pseudobulk_ranked.csv.gz",
        "brain_region_subclass",
    ),
    "Brain region": ("WB_HMBA_Human_AIT_brain_region_rank.csv.gz", "brain_region"),
    "Class": ("WB_HMBA_Human_AIT_class_rank.csv.gz", "class"),
}
SOURCE_URL = "https://github.com/AllenInstitute/HMBA_WB_Atlas"

DEFAULT_GENES = """APOE
ABCA1
ABCA7
ABI3
ACE
ADAM17
ADAMTS1
ANK3
ANKH
APH1B
APP
BCKDK
BIN1
BLNK
CASS4
CD2AP
CLNK
CLU
COX7C
CR1
CTSB
CTSH
DOC2A
EED
EPDR1
EPHA1
FERMT2
FOXF1
GRN
HLA-DQA1
HS3ST5
ICA1
IDUA
IGHG3
IGHV3-65
IL34
INPP5D
JAZF1
KLF16
LILRB2
MAF
MINDY2
MME
MS4A4A
MYO15A
NCK2
PLCG2
PLEKHA1
PRDM7
PRKD3
PTK2B
RASGEF1C
RBCK1
RHOH
SCIMP
SEC61G
SHARPIN
SIGLEC11
SLC24A4
SLC2A4RG
SNX1
SORL1
SORT1
SPDYE3
SPI1
SPPL2A
TMEM106B
TNIP1
TPCN1
TREM2
TREML2
TSPAN14
TSPOAP1
UMAD1
UNC5CL
USP6NL
WDR12
WDR81
WNT3"""


@st.cache_resource(show_spinner="Loading HMBA human rank matrix...")
def get_rank_matrix(dataset: str, file_version: tuple[int, int]):
    return load_rank_matrix(DATA_DIR / DATASETS[dataset][0])


@st.cache_resource(show_spinner=False)
def get_ortholog_map():
    return load_ortholog_map(DATA_DIR / "homologene_to_human.csv.gz")


def format_subclass(label: str) -> str:
    return re.sub(r"^(\d+)\s+(.+)$", r"\2 (ID:\1)", label)


def format_profile(label: str) -> str:
    region, separator, subclass = label.partition(" | ")
    return f"{region} | {format_subclass(subclass)}" if separator else format_subclass(label)


def format_results_table(table: pd.DataFrame, dataset: str) -> pd.DataFrame:
    formatted = table.copy()
    labels = formatted.pop("Cell type")
    if dataset == "Subclass by brain region":
        parts = labels.str.split(" | ", n=1, expand=True, regex=False)
        formatted.insert(0, "Subclass", parts[1].map(format_subclass))
        formatted.insert(1, "Brain region", parts[0])
    elif dataset == "Brain region":
        formatted.insert(0, "Brain region", labels)
    else:
        formatted.insert(0, dataset, labels.map(format_subclass))
    return formatted


@st.cache_data(show_spinner=False)
def create_heatmap_pdf(heatmap: pd.DataFrame, dataset: str) -> bytes:
    clustered_heatmap = heatmap
    if heatmap.shape[0] > 1:
        row_order = leaves_list(
            linkage(heatmap.to_numpy(), method="complete", metric="euclidean")
        )
        clustered_heatmap = clustered_heatmap.iloc[row_order, :]
    if heatmap.shape[1] > 1:
        column_order = leaves_list(
            linkage(heatmap.to_numpy().T, method="complete", metric="euclidean")
        )
        clustered_heatmap = clustered_heatmap.iloc[:, column_order]

    width = min(40, max(8, 3 + heatmap.shape[1] * 0.15))
    height = max(6, 3 + heatmap.shape[0] * 0.2)
    figure, axis = plt.subplots(figsize=(width, height))
    image = axis.imshow(clustered_heatmap.to_numpy(), aspect="auto", cmap="viridis")
    axis.set_xticks(range(clustered_heatmap.shape[1]))
    axis.set_xticklabels(
        [format_profile(label) for label in clustered_heatmap.columns],
        rotation=90, fontsize=5,
    )
    axis.set_yticks(range(clustered_heatmap.shape[0]))
    axis.set_yticklabels(clustered_heatmap.index, fontsize=6)
    axis.set_xlabel(dataset)
    axis.set_ylabel("Gene")
    axis.set_title(
        f"HMBA v0.5 human — {dataset}: {heatmap.shape[0]} genes\n"
        "hierarchically clustered genes and profiles; "
        "higher ranks are more expression-specific"
    )
    figure.colorbar(image, ax=axis, label="Expression-specificity rank")
    figure.tight_layout()

    output = BytesIO()
    figure.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(figure)
    return output.getvalue()


def format_p_value(value: float) -> str:
    if pd.isna(value):
        return ""
    if 0 < value < 0.001:
        return f"{value:.2e}"
    return f"{value:.3g}"


st.set_page_config(
    page_title="HMBA human polygenic expression specificity tester",
    layout="wide",
)

st.title("HMBA human polygenic expression specificity tester")

with st.sidebar:
    st.markdown(
        "**Test whether a gene set has higher expression-specificity ranks "
        "in human brain classes, subclasses, or regions using the HMBA v0.5 atlas.**"
    )
    dataset = st.selectbox("Profile grouping:", options=list(DATASETS))
    species = st.selectbox(
        "Species of input genes:", options=["Human", "Mouse", "Rhesus Macaque"]
    )
    st.caption(
        "Use the selected species for both target and background lists. "
        "Mouse and rhesus macaque symbols are converted to human orthologs. "
        "The default example uses human symbols."
    )

    gene_text = st.text_area(
        "Input your gene list:",
        value=DEFAULT_GENES,
        height=220,
    )
    background_text = st.text_area(
        "Background gene list (optional):",
        value="",
        height=90,
    )
    st.divider()
    st.markdown("**Source data:**")
    st.markdown(
        f"[Human and Mammalian Brain Atlas (HMBA), v0.5]({SOURCE_URL})  \n"
        "Allen Institute / BRAIN Initiative Cell Atlas Network (BICAN). "
        "Human AIT data."
    )
    st.caption("v0.5 is a working release; atlas annotations may change.")
    st.markdown("**Example gene set:**")
    st.markdown(
        "Default genes are from the [Bellenguez et al. Alzheimer's disease GWAS]"
        "(https://www.nature.com/articles/s41588-022-01024-z) (Table S5)."
    )


saved = None
try:
    parsed_targets = parse_gene_list(gene_text)
    parsed_background = parse_gene_list(background_text)
    if not parsed_targets:
        raise ValueError("Enter at least one target gene.")

    ortholog_map = None if species == "Human" else get_ortholog_map()
    conversion_unmapped = (
        [] if species == "Human" else [
            gene for gene in parsed_targets
            if not ortholog_map.get(species, {}).get(gene)
        ]
    )
    target_genes = convert_to_human(parsed_targets, species, ortholog_map)
    background_genes = (
        convert_to_human(parsed_background, species, ortholog_map)
        if parsed_background else None
    )
    if not target_genes:
        raise ValueError("None of the input genes could be converted to human symbols.")
    matrix_stat = (DATA_DIR / DATASETS[dataset][0]).stat()
    matrix = get_rank_matrix(dataset, (matrix_stat.st_mtime_ns, matrix_stat.st_size))
    result = run_auroc_analysis(
        matrix=matrix,
        target_genes=target_genes,
        background_genes=background_genes,
    )
    saved = {
        "table": format_results_table(result.table, dataset),
        "heatmap": result.heatmap,
        "matched": result.matched_gene_count,
        "input": result.input_gene_count,
        "background": result.background_gene_count,
        "conversion_unmapped": conversion_unmapped,
        "matrix_unmatched": list(result.unmatched_genes),
    }
except (OSError, ValueError, pd.errors.ParserError) as error:
    st.error(str(error))


if saved is not None:
    results_table = saved["table"]

    summary_columns = st.columns([4, 1])
    with summary_columns[0]:
        if species != "Human":
            st.caption("Gene counts below refer to unique converted human symbols.")
        st.write(f"Genes found in data: {saved['matched']} of {saved['input']}")
    unmapped_count = len(saved["conversion_unmapped"]) + len(saved["matrix_unmatched"])
    with summary_columns[1]:
        with st.popover(f"Unmapped genes ({unmapped_count})"):
            if saved["conversion_unmapped"]:
                st.markdown("**Input genes without a human ortholog:**")
                st.code("\n".join(saved["conversion_unmapped"]), language=None)
            if saved["matrix_unmatched"]:
                label = (
                    "Converted human symbols not found in the selected data/background:"
                    if species != "Human"
                    else "Input genes not found in the selected data/background:"
                )
                st.markdown(f"**{label}**")
                st.code("\n".join(saved["matrix_unmatched"]), language=None)
            if unmapped_count == 0:
                st.success("All input genes were mapped and found in the selected data.")
    st.caption(f"HMBA v0.5 human · {dataset} · {len(matrix.profiles):,} profiles")
    st.write(f"Background genes: {saved['background']}")

    table_for_display = results_table.copy()
    for label_column in ("Subclass", "Class"):
        if label_column in table_for_display.columns:
            table_for_display[label_column] = table_for_display[label_column].str.replace(
                r" \(ID:\d+\)$", "", regex=True
            )
    download_prefix = f"hmba_v0.5_human_{DATASETS[dataset][1]}"
    st.html("""
        <style>
        .st-key-result_downloads button {
            color: #0068c9;
            text-decoration: underline;
            padding: 0;
            min-height: 0;
        }
        </style>
    """)
    with st.container(key="result_downloads", horizontal=True):
        st.download_button(
            "Download as CSV",
            data=results_table.to_csv(index=False).encode("utf-8"),
            file_name=f"{download_prefix}_results.csv",
            mime="text/csv",
            type="tertiary",
            on_click="ignore",
        )
        st.download_button(
            "Download heatmap (PDF)",
            data=partial(create_heatmap_pdf, saved["heatmap"], dataset),
            file_name=f"{download_prefix}_heatmap.pdf",
            mime="application/pdf",
            type="tertiary",
            on_click="ignore",
        )

    display_table = table_for_display.style.format(
        {"AUROC": lambda value: f"{value:.3g}",
         "pValue": format_p_value, "FDR": format_p_value}
    )
    st.dataframe(
        display_table,
        hide_index=True,
        width="stretch",
        height=600,
    )
    st.caption(
        "AUROC compares target genes against non-target background genes using "
        "within-profile rank sums. Two-sided p-values use the AUROC-based "
        "normal approximation with continuity correction and no tie correction. "
        "FDR values use the Benjamini-Hochberg correction across all selected profiles applied to "
        "full-precision p-values."
    )
