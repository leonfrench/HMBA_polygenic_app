"""Core functions for HMBA human expression-specificity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


@dataclass(frozen=True)
class RankMatrix:
    genes: np.ndarray
    profiles: tuple[str, ...]
    ranks: np.ndarray


@dataclass(frozen=True)
class AnalysisResult:
    table: pd.DataFrame
    heatmap: pd.DataFrame
    unmatched_genes: tuple[str, ...]
    input_gene_count: int
    matched_gene_count: int
    background_gene_count: int


def unique_in_order(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def parse_gene_list(text: str) -> list[str]:
    """Split gene symbols on whitespace, commas, or semicolons."""
    return unique_in_order(
        [value for value in re.split(r"[\s,;]+", text.strip()) if value]
    )


def load_rank_matrix(path: Path) -> RankMatrix:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    if not header or header[0] != "gene_symbol":
        raise ValueError(f"{path.name} does not start with a gene_symbol column.")

    # The rank file serializes values with a decimal suffix (for example,
    # 9089.0), so load them as floating-point values.
    dtypes = {column: np.float64 for column in header[1:]}
    frame = pd.read_csv(path, dtype=dtypes)
    genes = frame["gene_symbol"].astype(str).to_numpy()
    if pd.Index(genes).has_duplicates:
        raise ValueError(f"{path.name} contains duplicate gene symbols.")

    ranks = frame.iloc[:, 1:].to_numpy(dtype=np.float64, copy=False)
    if not np.isfinite(ranks).all():
        raise ValueError(f"{path.name} contains missing or non-finite rank values.")
    genes.setflags(write=False)
    ranks.setflags(write=False)
    return RankMatrix(genes=genes, profiles=tuple(header[1:]), ranks=ranks)


def load_ortholog_map(path: Path) -> dict[str, dict[str, tuple[str, ...]]]:
    frame = pd.read_csv(path, dtype=str)
    required = {"species", "source_symbol", "human_symbol"}
    if not required.issubset(frame.columns):
        raise ValueError(f"{path.name} is missing required ortholog columns.")

    output: dict[str, dict[str, tuple[str, ...]]] = {}
    for (species, source_symbol), group in frame.groupby(
        ["species", "source_symbol"], sort=False
    ):
        human_symbols = tuple(unique_in_order(group["human_symbol"].tolist()))
        output.setdefault(species, {})[source_symbol] = human_symbols
    return output


def convert_to_human(
    genes: list[str],
    species: str,
    ortholog_map: dict[str, dict[str, tuple[str, ...]]] | None = None,
) -> list[str]:
    if species == "Human":
        return unique_in_order(genes)
    if ortholog_map is None:
        raise ValueError("The ortholog map is required for non-human input genes.")

    species_map = ortholog_map.get(species, {})
    converted: list[str] = []
    for gene in genes:
        converted.extend(species_map.get(gene, ()))
    return unique_in_order(converted)


def _significant(value: float, digits: int = 3) -> float:
    if value == 0 or not np.isfinite(value):
        return float(value)
    return float(format(value, f".{digits}g"))


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment using full-precision p-values."""
    number_of_tests = len(p_values)
    order = np.argsort(p_values)
    ordered_p_values = p_values[order]
    adjusted_ordered = ordered_p_values * number_of_tests / np.arange(
        1, number_of_tests + 1
    )
    adjusted_ordered = np.minimum.accumulate(adjusted_ordered[::-1])[::-1]
    adjusted = np.empty(number_of_tests, dtype=np.float64)
    adjusted[order] = np.clip(adjusted_ordered, 0.0, 1.0)
    return adjusted


def run_auroc_analysis(
    matrix: RankMatrix,
    target_genes: list[str],
    background_genes: list[str] | None,
) -> AnalysisResult:
    input_gene_count = len(target_genes)
    matrix_gene_set = set(matrix.genes)

    if background_genes is None:
        selected_genes = matrix.genes
        ranks = matrix.ranks
    else:
        background_set = set(background_genes).intersection(matrix_gene_set)
        selected_mask = np.fromiter(
            (gene in background_set for gene in matrix.genes),
            dtype=bool,
            count=len(matrix.genes),
        )
        selected_genes = matrix.genes[selected_mask]
        if len(selected_genes) == 0:
            raise ValueError("None of the background genes were found in the matrix.")
        ranks = matrix.ranks[selected_mask, :]

    selected_gene_set = set(selected_genes)
    unmatched_targets = tuple(
        gene for gene in target_genes if gene not in selected_gene_set
    )
    matched_targets = selected_gene_set.intersection(target_genes)
    target_mask = np.fromiter(
        (gene in matched_targets for gene in selected_genes),
        dtype=bool,
        count=len(selected_genes),
    )
    n_target = int(target_mask.sum())
    n_background = int(len(selected_genes))
    n_negative = n_background - n_target

    if n_target == 0:
        raise ValueError(
            "None of the submitted genes were found in the selected rank-matrix background."
        )
    if n_negative == 0:
        raise ValueError("The background must contain at least one non-target gene.")

    # Mann–Whitney ranks the supplied scores within each profile and handles ties.
    # Work one profile at a time to bound memory for the regional matrix.
    aucs = np.empty(len(matrix.profiles), dtype=np.float64)
    p_values = np.empty_like(aucs)
    for column in range(len(matrix.profiles)):
        if np.ptp(ranks[:, column]) == 0:
            aucs[column] = 0.5
            p_values[column] = 1.0
            continue
        test = mannwhitneyu(
            ranks[target_mask, column],
            ranks[~target_mask, column],
            alternative="two-sided",
            method="asymptotic",
            use_continuity=True,
        )
        aucs[column] = test.statistic / (n_target * n_negative)
        p_values[column] = test.pvalue

    rounded_aucs = np.asarray([_significant(value) for value in aucs])
    rounded_p_values = np.asarray([_significant(value) for value in p_values])
    fdr_values = _benjamini_hochberg(p_values)
    rounded_fdr = np.asarray([_significant(value) for value in fdr_values])

    table = pd.DataFrame({"Cell type": matrix.profiles})
    table["AUROC"] = rounded_aucs
    table["pValue"] = rounded_p_values
    table["FDR"] = rounded_fdr
    table = table.sort_values("AUROC", ascending=False, kind="mergesort").reset_index(
        drop=True
    )
    heatmap = pd.DataFrame(
        ranks[target_mask, :],
        index=selected_genes[target_mask],
        columns=matrix.profiles,
    )
    heatmap.index.name = "gene_symbol"

    return AnalysisResult(
        table=table,
        heatmap=heatmap,
        unmatched_genes=unmatched_targets,
        input_gene_count=input_gene_count,
        matched_gene_count=n_target,
        background_gene_count=n_background,
    )
