# HMBA human polygenic expression specificity app

This Streamlit app compares a human gene set's expression-specificity scores
across profiles from the **Human and Mammalian Brain Atlas (HMBA), v0.5**.
Choose Subclass (the default), Subclass by brain region, Brain region, or Class in the sidebar.

## Data and attribution

Source: [Allen Institute HMBA Whole-Brain Atlas](https://github.com/AllenInstitute/HMBA_WB_Atlas),
part of the BRAIN Initiative Cell Atlas Network (BICAN). This app uses the **human AIT data**.
HMBA v0.5 is a working release, so annotations may change. The source repository
documents the atlas, metadata, and original data downloads.

Keep these files in `data/` (no decompression is needed):

| Grouping | File | Genes | Profiles |
| --- | --- | ---: | ---: |
| Class | `WB_HMBA_Human_AIT_class_rank.csv.gz` | 22,738 | 42 |
| Subclass | `WB_HMBA_Human_AIT_subclass_rank.csv.gz` | 22,738 | 387 |
| Brain region | `WB_HMBA_Human_AIT_brain_region_rank.csv.gz` | 22,738 | 16 |
| Subclass by brain region | `WB_HMBA_Human_AIT_brain_region_subclass_pseudobulk_ranked.csv.gz` | 22,737 | 1,357 |

Each file starts with `gene_symbol`, followed by numeric profile columns.
Regional matrix columns use `brain region | subclass` labels. Results and CSV
downloads split these into Brain region and Subclass columns. The displayed table omits subclass and class IDs. CSV downloads and heatmaps retain
them after the name as `Microglia_NN (ID:408)`. Brain-region-only results display a Brain region column.
The app formats labels in memory without modifying the source CSVs. The updated
matrices contain rank scores ranging from 1 to the number of genes in each file.
Replacing a rank file invalidates its cached matrix using its modification time
and size.

## Run locally

```bash
cd /Users/leon/projects/HMBA_polygenic_app
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

Enter gene symbols separated by whitespace, commas, or semicolons, and select
Human, Mouse, or Rhesus Macaque. Both target and background lists use the selected
species. Mouse and rhesus macaque inputs are converted to human symbols using
the bundled `data/homologene_to_human.csv.gz` file (columns: `species`,
`source_symbol`, `human_symbol`). All mapped human orthologs are retained and
deduplicated. The atlas itself remains human, and the default example uses human
symbols. The unmapped-gene popover distinguishes failed target conversions from
human symbols absent from the selected data/background. Reported analysis gene
counts refer to unique human symbols after conversion. An optional
background list restricts the comparison universe; it should include the target
genes. Without it, the selected matrix's full gene universe is used. Target
genes absent from that universe are reported in the unmapped-gene popover.
The default example is from the [Bellenguez et al. Alzheimer's disease GWAS,
Table S5](https://www.nature.com/articles/s41588-022-01024-z).

## Analysis and downloads

For each profile, AUROC compares the supplied specificity scores of matched
target genes against the non-target genes in the background. The two group sizes are matched
target genes and non-target background genes. Benjamini–Hochberg FDR
is calculated from full-precision p-values across all profiles in the selected
matrix. AUROC above 0.5 indicates higher specificity scores in the target set.

The analysis runs immediately and updates when inputs change. Download the
results as CSV or a PDF heatmap with complete-linkage Euclidean clustering of
genes and profiles. Clustering and PDF generation run only when the heatmap
download is clicked; generated PDFs are cached for identical inputs. The heatmap retains the supplied specificity scores even
with a custom background. Download filenames identify the atlas and grouping.
The regional heatmap includes all 1,357 profiles and is best inspected by zooming.
