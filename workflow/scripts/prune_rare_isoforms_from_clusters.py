#!/usr/bin/env python3
"""
prune_rare_isoforms_from_clusters.py

Clusters transcripts by splice junction overlap and filters rare isoforms using
per-sample support logic (relative expression, fraction of locus expression, or
cumulative expression retention).
"""

import argparse
import gzip
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional, Any, TextIO

import numpy as np
import pandas as pd


def ts() -> str:
    """
    Return current wall-clock timestamp formatted as HH:MM:SS.

    Returns
    -------
    str
        Formatted timestamp string.
    """
    return datetime.now().strftime("%H:%M:%S")


def smart_open(path: Path) -> TextIO:
    """
    Open plain-text or gzip-compressed files transparently.

    Parameters
    ----------
    path : Path
        Path to the file to open.

    Returns
    -------
    TextIO
        An open file stream in read-text mode.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_gtf(
    gtf_path: Path,
    transcript_attr: str = "transcript_id",
    gene_attr: str = "gene_id",
) -> Tuple[
    Dict[str, List[Tuple[int, int]]],   # tx2chain (junction list)
    Dict[str, str],                     # tx2gene
    Dict[str, str],                     # tx2locus "chrom|strand"
]:
    """
    Parse exon features from a GTF file to construct transcript junction chains.

    Parameters
    ----------
    gtf_path : Path
        Path to input GTF file (can be gzip-compressed).
    transcript_attr : str, default="transcript_id"
        GTF attribute key for transcript identifier.
    gene_attr : str, default="gene_id"
        GTF attribute key for gene identifier.

    Returns
    -------
    tx2chain : Dict[str, List[Tuple[int, int]]]
        Map of transcript ID to list of splice junctions `(donor, acceptor)`.
    tx2gene : Dict[str, str]
        Map of transcript ID to gene ID.
    tx2locus : Dict[str, str]
        Map of transcript ID to locus string formatted as `"chrom|strand"`.
    """
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx2meta: Dict[str, Tuple[str, str, str]] = {}

    with smart_open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue

            chrom, _src, _feat, start, end, _score, strand, _frame, attrs = parts

            attr_dict: Dict[str, str] = {}
            for field in attrs.split(";"):
                field = field.strip()
                if not field:
                    continue
                if " " in field:
                    k, v = field.split(" ", 1)
                elif "=" in field:
                    k, v = field.split("=", 1)
                else:
                    continue
                attr_dict[k] = v.strip().strip('"').strip("'")

            tx = attr_dict.get(transcript_attr)
            gene = attr_dict.get(gene_attr)
            if not tx or not gene:
                continue

            exons_by_tx[tx].append((int(start), int(end)))
            if tx not in tx2meta:
                tx2meta[tx] = (gene, chrom, strand)

    tx2chain: Dict[str, List[Tuple[int, int]]] = {}
    tx2gene: Dict[str, str] = {}
    tx2locus: Dict[str, str] = {}

    for tx, exons in exons_by_tx.items():
        gene, chrom, strand = tx2meta[tx]
        exons_sorted = sorted(exons, key=lambda x: x[0])

        juncs: List[Tuple[int, int]] = []
        for i in range(len(exons_sorted) - 1):
            donor = exons_sorted[i][1]        # exon end
            acceptor = exons_sorted[i + 1][0] # next exon start
            if donor < acceptor:
                juncs.append((donor, acceptor))

        tx2chain[tx] = juncs
        tx2gene[tx] = gene
        tx2locus[tx] = f"{chrom}|{strand}"

    return tx2chain, tx2gene, tx2locus


def load_expression_table(path: Path, tx_col: str, expr_col: str) -> pd.DataFrame:
    """
    Load expression matrix from CSV/TSV table into DataFrame.

    Parameters
    ----------
    path : Path
        Path to expression table file.
    tx_col : str
        Name of column containing transcript IDs.
    expr_col : str
        Expression column specification. Accepts single column name, comma-separated
        list of names, or `'*'` / `'ALL'` to use all non-transcript ID columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized `'transcript_id'` column and numeric expression columns.

    Raises
    ------
    ValueError
        If the transcript or expression columns are missing, or if transcript IDs are
        not unique. Duplicate IDs are rejected here rather than tolerated downstream:
        they would be double-counted in cluster totals, silently shifting every
        cumulative expression cutoff in the affected cluster.
    """
    df = pd.read_csv(path, sep=None, engine="python")

    cols_lower = {c.lower(): c for c in df.columns}
    tx_c = cols_lower.get(tx_col.lower(), tx_col)
    if tx_c not in df.columns:
        raise ValueError(f"Missing transcript column '{tx_col}'. Available: {list(df.columns)}")

    expr_col_norm = expr_col.strip().upper()
    if expr_col_norm in ["*", "ALL"]:
        expr_cols = [c for c in df.columns if c != tx_c]
    else:
        expr_cols = [c.strip() for c in expr_col.split(",") if c.strip()]
        if not expr_cols:
            raise ValueError("--expr-col resolved to empty list.")
        missing = [c for c in expr_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Missing expression column(s): {missing}. Available: {list(df.columns)}")

    out = df[[tx_c] + expr_cols].rename(columns={tx_c: "transcript_id"}).copy()
    out["transcript_id"] = out["transcript_id"].astype(str)

    dups = out["transcript_id"][out["transcript_id"].duplicated()].unique()
    if len(dups):
        preview = ", ".join(dups[:5])
        suffix = ", ..." if len(dups) > 5 else ""
        raise ValueError(
            f"Duplicate transcript IDs in expression table ({len(dups)} distinct): {preview}{suffix}. "
            "Deduplicate the input before running; duplicates would be double-counted in cluster totals."
        )

    for c in expr_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def jaccard(a: Set[Tuple[int, int]], b: Set[Tuple[int, int]]) -> float:
    """
    Calculate Jaccard similarity index between two junction sets.

    Parameters
    ----------
    a : Set[Tuple[int, int]]
        First set of splice junctions.
    b : Set[Tuple[int, int]]
        Second set of splice junctions.

    Returns
    -------
    float
        Jaccard similarity coefficient in range [0.0, 1.0].
    """
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    uni = len(a | b)
    return inter / uni if uni else 0.0


def match_ok(
    A: Set[Tuple[int, int]],
    B: Set[Tuple[int, int]],
    mode: str,
    min_shared: int = 1,
    min_jaccard: float = 0.5,
    bridge_strong_min_shared: int = 2,
    bridge_weak_min_jaccard: float = 0.15,
) -> bool:
    """
    Determine if two junction sets satisfy specified structural matching criteria.

    Parameters
    ----------
    A : Set[Tuple[int, int]]
        Junction set for transcript A.
    B : Set[Tuple[int, int]]
        Junction set for transcript B.
    mode : str
        Matching criterion (`'bridge_safe'`, `'any_shared'`, `'subset'`, or `'jaccard'`).
    min_shared : int, default=1
        Minimum number of shared junctions required.
    min_jaccard : float, default=0.5
        Minimum Jaccard index for `'jaccard'` mode.
    bridge_strong_min_shared : int, default=2
        Strong overlap threshold for `'bridge_safe'` mode.
    bridge_weak_min_jaccard : float, default=0.15
        Weak Jaccard threshold required when strong junction overlap is not met in `'bridge_safe'`.

    Returns
    -------
    bool
        True if connection criteria are met, False otherwise.
    """
    inter = len(A & B)

    if mode == "any_shared":
        return inter >= min_shared
    if mode == "subset":
        return inter >= min_shared and (A.issubset(B) or B.issubset(A))
    if mode == "jaccard":
        return inter >= min_shared and jaccard(A, B) >= min_jaccard
    if mode == "bridge_safe":
        if inter >= bridge_strong_min_shared:
            return True
        return inter >= 1 and jaccard(A, B) >= bridge_weak_min_jaccard

    raise ValueError(f"Unknown match mode: {mode}")


def connected_components(nodes: List[str], edges: Dict[str, List[str]]) -> List[List[str]]:
    """
    Compute connected components for an undirected graph via iterative depth-first traversal.

    Parameters
    ----------
    nodes : List[str]
        List of graph node IDs.
    edges : Dict[str, List[str]]
        Adjacency list mapping node ID to connected neighbor IDs.

    Returns
    -------
    List[List[str]]
        List of connected components, where each component is a list of node IDs.
    """
    seen: Set[str] = set()
    comps: List[List[str]] = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        seen.add(n)
        comp: List[str] = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in edges.get(cur, []):
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        comps.append(comp)
    return comps


def build_clusters_for_group(
    tx_list: List[str],
    tx2chain: Dict[str, List[Tuple[int, int]]],
    match_mode: str,
    min_shared: int,
    min_jaccard: float,
    bridge_strong_min_shared: int,
    bridge_weak_min_jaccard: float,
) -> List[List[str]]:
    """
    Cluster transcripts within a gene/locus group into connected structural components.

    Parameters
    ----------
    tx_list : List[str]
        List of transcript IDs belonging to the locus.
    tx2chain : Dict[str, List[Tuple[int, int]]]
        Global junction chain lookup.
    match_mode : str
        Clustering matching strategy (`'bridge_safe'`, `'any_shared'`, `'subset'`, or `'jaccard'`).
    min_shared : int
        Minimum shared junction requirement (validated in `main`).
    min_jaccard : float
        Jaccard threshold for overlap modes.
    bridge_strong_min_shared : int
        Strong overlap threshold for bridge-safe clustering.
    bridge_weak_min_jaccard : float
        Weak Jaccard threshold for bridge-safe clustering.

    Returns
    -------
    List[List[str]]
        List of isoform clusters (each cluster is a list of transcript IDs).
    """
    if len(tx_list) <= 1:
        return [tx_list]

    junc_to_tx: Dict[Tuple[int, int], List[str]] = defaultdict(list)
    tx_sets: Dict[str, Set[Tuple[int, int]]] = {}

    for t in tx_list:
        s = set(tx2chain.get(t, []))
        tx_sets[t] = s
        for j in s:
            junc_to_tx[j].append(t)

    possible_pairs: Set[Tuple[str, str]] = set()
    for shared_txs in junc_to_tx.values():
        n = len(shared_txs)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = shared_txs[i], shared_txs[j]
                possible_pairs.add((a, b) if a <= b else (b, a))

    edges: Dict[str, List[str]] = defaultdict(list)
    for a, b in possible_pairs:
        if match_ok(
            tx_sets[a], tx_sets[b],
            mode=match_mode,
            min_shared=min_shared,
            min_jaccard=min_jaccard,
            bridge_strong_min_shared=bridge_strong_min_shared,
            bridge_weak_min_jaccard=bridge_weak_min_jaccard,
        ):
            edges[a].append(b)
            edges[b].append(a)

    return connected_components(tx_list, edges)


def _keepers_by_cumexpr_fraction(
    values: pd.Series,
    retain_frac: float,
    tie_policy: str = "include_all_ties",
    tie_abs_tol: float = 1e-9,
) -> Set[str]:
    """
    Select transcripts in a sample achieving cumulative expression retention target.

    Parameters
    ----------
    values : pd.Series
        Non-negative expression values for one sample within a cluster, indexed by transcript ID.
    retain_frac : float
        Target fraction of total cluster expression mass to retain (e.g. 0.98 for 98%).
    tie_policy : str, default="include_all_ties"
        Policy for handling boundary ties:
        - `'include_all_ties'`: Retain all non-zero isoforms tied at cutoff value.
        - `'exclude_all_ties'`: Exclude boundary ties entirely (yielding no support if all tied).
        - `'truncate_ties'`: Keep minimal prefix hitting cumulative target, breaking ties by transcript ID.
    tie_abs_tol : float, default=1e-9
        Absolute tolerance for float comparison around the boundary cutoff value.

    Returns
    -------
    Set[str]
        Set of supported transcript IDs for the given sample.
    """
    v = values.clip(lower=0.0)
    total = float(v.sum())
    if total <= 0.0:
        return set()

    # Sort deterministically: value descending, then transcript_id ascending
    sorted_df = (
        pd.DataFrame({"v": v, "tx": v.index.astype(str)})
        .sort_values(by=["v", "tx"], ascending=[False, True], kind="mergesort")
    )

    v_sorted = sorted_df["v"]
    csum = v_sorted.cumsum()
    cutoff = retain_frac * total

    pass_mask = csum >= cutoff
    if not pass_mask.any():
        # Cutoff unreachable through float error (e.g. retain_frac == 1.0): keep all
        # non-zero isoforms rather than collapsing to the single top isoform.
        return set(sorted_df[sorted_df["v"] > 0.0]["tx"].tolist())

    cutoff_idx = int(pass_mask.to_numpy().argmax())
    cutoff_val = float(v_sorted.iloc[cutoff_idx])

    if tie_policy == "truncate_ties":
        keepers = sorted_df.iloc[: cutoff_idx + 1]
        keepers = keepers[keepers["v"] > 0.0]["tx"].tolist()

    elif tie_policy == "exclude_all_ties":
        # Strictly greater than cutoff value using absolute tolerance (rtol=0.0)
        above_mask = (v_sorted > cutoff_val) & ~np.isclose(v_sorted, cutoff_val, rtol=0.0, atol=tie_abs_tol)
        keepers = sorted_df[above_mask & (v_sorted > 0.0)]["tx"].tolist()

    else:  # "include_all_ties"
        # Greater than or close to cutoff value using absolute tolerance (rtol=0.0)
        ge_mask = (v_sorted >= cutoff_val) | np.isclose(v_sorted, cutoff_val, rtol=0.0, atol=tie_abs_tol)
        keepers = sorted_df[ge_mask & (v_sorted > 0.0)]["tx"].tolist()

    return set(keepers)


def keep_by_sample_support(
    comp: List[str],
    sub: pd.DataFrame,
    expr_cols: List[str],
    sample_min_rel_expr: float,
    sample_min_frac_expr: float,
    retain_locus_expr_pct_per_sample: Optional[float],
    min_abs_expr: float,
    min_cluster_expr: float,
    min_support_samples: int,
    min_keep: int,
    tie_policy: str = "include_all_ties",
    tie_abs_tol: float = 1e-9,
) -> Tuple[Set[str], List[str], Dict[str, Any]]:
    """
    Filter transcripts per cluster based on multi-sample support criteria.

    Parameters
    ----------
    comp : List[str]
        Cluster transcript IDs.
    sub : pd.DataFrame
        Locus expression table.
    expr_cols : List[str]
        Sample expression column names.
    sample_min_rel_expr : float
        Fraction of within-cluster max expression required per sample.
    sample_min_frac_expr : float
        Fraction of within-cluster total expression required per sample.
    retain_locus_expr_pct_per_sample : float, optional
        Percentage of cumulative locus expression mass to retain per sample.
    min_abs_expr : float
        Absolute minimum expression floor required per sample.
    min_cluster_expr : float
        Absolute minimum cluster total expression floor required per sample.
    min_support_samples : int
        Minimum number of samples in which an isoform must pass support filter.
    min_keep : int
        Minimum transcripts to retain per cluster (0 disables fallback retention).
    tie_policy : str, default="include_all_ties"
        Tie-handling policy at cumulative expression boundary.
    tie_abs_tol : float, default=1e-9
        Absolute floating point tolerance for boundary ties.

    Returns
    -------
    keep_tx : Set[str]
        Set of kept transcript IDs.
    drop_tx : List[str]
        List of dropped transcript IDs.
    stats : Dict[str, Any]
        Dictionary summarizing filtering threshold metrics, per-transcript drop reasons,
        and transcripts retained only by the `min_keep` fallback.
    """
    comp_df = (
        sub[sub["transcript_id"].isin(comp)]
        .set_index("transcript_id")[expr_cols]
        .astype(float)
    )

    passes = pd.DataFrame(False, index=comp_df.index, columns=comp_df.columns)
    max_per_sample = comp_df.max(axis=0)
    total_per_sample = comp_df.sum(axis=0)

    # 1. Relative-to-max support
    if sample_min_rel_expr > 0:
        thresh = max_per_sample * float(sample_min_rel_expr)
        passes |= comp_df.ge(thresh, axis=1) & (max_per_sample > 0.0)

    # 2. Relative-to-total support
    if sample_min_frac_expr > 0:
        thresh = total_per_sample * float(sample_min_frac_expr)
        passes |= comp_df.ge(thresh, axis=1) & (total_per_sample > 0.0)

    # 3. Cumulative expression retention per sample
    if retain_locus_expr_pct_per_sample is not None:
        retain_frac = float(retain_locus_expr_pct_per_sample) / 100.0
        for sample in expr_cols:
            keepers = _keepers_by_cumexpr_fraction(
                comp_df[sample],
                retain_frac=retain_frac,
                tie_policy=tie_policy,
                tie_abs_tol=tie_abs_tol,
            )
            if keepers:
                passes.loc[list(keepers), sample] = True

    # Track pre-noise-floor passes to distinguish filtering causes
    pass_before_floors = passes.copy()

    # Enforce noise floors (minimum absolute expression & cluster total)
    if min_abs_expr > 0.0:
        passes &= comp_df.ge(min_abs_expr, axis=1)
    if min_cluster_expr > 0.0:
        passes &= (total_per_sample >= min_cluster_expr)

    pass_counts = passes.sum(axis=1)
    supported_tx = set(pass_counts[pass_counts >= int(min_support_samples)].index.astype(str))

    # Keep all supported isoforms
    keep_tx = set(supported_tx)
    rescued_tx: Set[str] = set()

    # Top up shortfall using fallback ranking if keep_tx count is below min_keep.
    # This only ADDS transcripts: a supported isoform is never displaced by the fallback.
    if len(keep_tx) < int(min_keep):
        n_needed = int(min_keep) - len(keep_tx)
        sum_scores = comp_df.sum(axis=1)
        fallback_df = pd.DataFrame({"score": sum_scores, "tx": comp_df.index.astype(str)})
        fallback_df.sort_values(by=["score", "tx"], ascending=[False, True], inplace=True)

        rescued = fallback_df[~fallback_df["tx"].isin(keep_tx)].head(n_needed)
        rescued_tx = set(rescued["tx"].tolist())
        keep_tx |= rescued_tx

    drop_tx = [t for t in comp_df.index.astype(str).tolist() if t not in keep_tx]

    # Annotate drop reasons: separate "the filter ranked this out" from
    # "this never cleared the noise floor".
    pass_counts_before_floors = pass_before_floors.sum(axis=1)
    drop_reasons: Dict[str, str] = {}
    for t in drop_tx:
        if int(pass_counts_before_floors.loc[t]) >= int(min_support_samples):
            drop_reasons[t] = "below_noise_floor"
        else:
            drop_reasons[t] = "failed_sample_support"

    stats = {
        "sample_min_rel_expr": float(sample_min_rel_expr),
        "sample_min_frac_expr": float(sample_min_frac_expr),
        "retain_locus_expr_pct_per_sample": float(retain_locus_expr_pct_per_sample) if retain_locus_expr_pct_per_sample is not None else None,
        "min_abs_expr": float(min_abs_expr),
        "min_cluster_expr": float(min_cluster_expr),
        "min_support_samples": int(min_support_samples),
        "pass_count_min": int(pass_counts.min()) if len(pass_counts) else 0,
        "pass_count_max": int(pass_counts.max()) if len(pass_counts) else 0,
        "rescued_by_min_keep": rescued_tx,
        "drop_reasons": drop_reasons,
    }
    return keep_tx, drop_tx, stats


def main():
    """Main CLI execution entry point."""
    ap = argparse.ArgumentParser(description="Cluster isoforms by splice junction overlap and filter rare isoforms per-sample.")

    # Input Options
    io_group = ap.add_argument_group("Input/Output Arguments")
    io_group.add_argument("--gtf", required=True, type=Path, help="Path to input GTF file.")
    io_group.add_argument("--expr", required=True, type=Path, help="Path to input expression matrix.")
    io_group.add_argument("--tx-col", default="transcript_id", help="Transcript ID column name.")
    io_group.add_argument("--expr-col", default="*", help="Comma-separated expression columns or '*' for all numeric columns.")
    io_group.add_argument("--out", required=True, type=Path, help="Output path for filtered expression matrix.")
    io_group.add_argument("--clusters-out", type=Path, default=None, help="Optional output path for transcript cluster assignments.")
    io_group.add_argument("--dropped-out", type=Path, default=None, help="Optional output path for dropped isoform metrics.")

    # Clustering Options
    cluster_group = ap.add_argument_group("Clustering Options")
    cluster_group.add_argument("--match-mode", choices=["bridge_safe", "any_shared", "subset", "jaccard"], default="bridge_safe",
                               help="Clustering mode. Defaults to 'bridge_safe' to prevent cluster collapsing.")
    cluster_group.add_argument("--min-shared", type=int, default=1)
    cluster_group.add_argument("--min-jaccard", type=float, default=0.5)
    cluster_group.add_argument("--bridge-strong-min-shared", type=int, default=2)
    cluster_group.add_argument("--bridge-weak-min-jaccard", type=float, default=0.15)

    # Per-Sample Filtering Options
    filter_group = ap.add_argument_group("Per-Sample Filtering Options")
    filter_group.add_argument("--retain-locus-expr-pct-per-sample", type=float, default=None,
                               help="Keep isoforms that cumulatively account for top X%% of cluster expression in a sample.")
    filter_group.add_argument("--sample-min-rel-expr", type=float, default=0.0,
                               help="Keep isoform if expr >= this * (max cluster expr in sample). 0 disables.")
    filter_group.add_argument("--sample-min-frac-expr", type=float, default=0.0,
                               help="Keep isoform if expr >= this * (total cluster expr in sample). 0 disables.")
    filter_group.add_argument("--min-abs-expr", type=float, default=0.0,
                               help="Noise floor: minimum absolute expression required in a sample to pass support.")
    filter_group.add_argument("--min-cluster-expr", type=float, default=0.0,
                               help="Noise floor: minimum total cluster expression required in a sample to pass support.")
    filter_group.add_argument("--tie-policy", choices=["include_all_ties", "exclude_all_ties", "truncate_ties"], default="include_all_ties",
                               help="Policy for handling expression ties at cumulative boundary.")
    filter_group.add_argument("--tie-abs-tol", type=float, default=1e-9,
                               help="Absolute floating-point tolerance for considering expression values equal during tie checks.")
    filter_group.add_argument("--min-support-samples", type=int, default=1,
                               help="Keep isoform if supported in at least this many samples.")
    filter_group.add_argument("--min-keep", type=int, default=1,
                               help="Always keep at least this many isoforms per cluster (0 disables fallback retention).")

    args = ap.parse_args()

    # Validation
    # Note: min_keep >= 0 allows setting min_keep=0 to completely disable fallback retention
    if args.min_keep < 0:
        raise ValueError("--min-keep must be >= 0.")
    if args.sample_min_rel_expr < 0:
        raise ValueError("--sample-min-rel-expr must be >= 0.")
    if not (0.0 <= args.sample_min_frac_expr <= 1.0):
        raise ValueError("--sample-min-frac-expr must be in [0.0, 1.0].")
    if args.min_abs_expr < 0:
        raise ValueError("--min-abs-expr must be >= 0.")
    if args.min_cluster_expr < 0:
        raise ValueError("--min-cluster-expr must be >= 0.")
    if args.tie_abs_tol < 0:
        raise ValueError("--tie-abs-tol must be >= 0.")
    if args.min_shared < 1:
        raise ValueError("--min-shared must be >= 1 for overlap-based clustering.")
    if not (0.0 <= args.min_jaccard <= 1.0):
        raise ValueError("--min-jaccard must be in [0.0, 1.0].")
    if not (0.0 <= args.bridge_weak_min_jaccard <= 1.0):
        raise ValueError("--bridge-weak-min-jaccard must be in [0.0, 1.0].")

    # Note: min_support_samples >= 1 is required as 0 sample support is semantically meaningless
    if args.min_support_samples < 1:
        raise ValueError("--min-support-samples must be >= 1.")

    if (args.sample_min_rel_expr == 0.0) and (args.sample_min_frac_expr == 0.0) and (args.retain_locus_expr_pct_per_sample is None):
        raise ValueError("At least one sample support threshold (--sample-min-rel-expr, --sample-min-frac-expr, or --retain-locus-expr-pct-per-sample) must be specified.")

    if args.retain_locus_expr_pct_per_sample is not None:
        if not (0.0 < args.retain_locus_expr_pct_per_sample <= 100.0):
            raise ValueError("--retain-locus-expr-pct-per-sample must be in (0, 100].")

    print(f"[{ts()}] Loading GTF...")
    tx2chain, tx2gene, tx2locus = parse_gtf(args.gtf)

    print(f"[{ts()}] Loading expression table...")
    expr_df = load_expression_table(args.expr, args.tx_col, args.expr_col)

    valid_txs = set(tx2chain.keys())
    expr_df = expr_df[expr_df["transcript_id"].isin(valid_txs)].copy()

    expr_df["__gene__"] = expr_df["transcript_id"].map(tx2gene)
    expr_df["__locus__"] = expr_df["transcript_id"].map(tx2locus)
    expr_df = expr_df[~expr_df["__gene__"].isna() & ~expr_df["__locus__"].isna()].copy()

    expr_cols = [c for c in expr_df.columns if c not in ["transcript_id", "__gene__", "__locus__"]]

    keep_set: Set[str] = set()
    rescued_set: Set[str] = set()
    cluster_rows: List[Dict[str, Any]] = []
    dropped_rows: List[Dict[str, Any]] = []
    cluster_id_counter = 0

    print(f"[{ts()}] Clustering + filtering (match-mode={args.match_mode}, tie-policy={args.tie_policy})...")
    for (gene, locus), sub in expr_df.groupby(["__gene__", "__locus__"], sort=False):
        tx_list = sub["transcript_id"].astype(str).tolist()
        clusters = build_clusters_for_group(
            tx_list=tx_list,
            tx2chain=tx2chain,
            match_mode=args.match_mode,
            min_shared=args.min_shared,
            min_jaccard=args.min_jaccard,
            bridge_strong_min_shared=args.bridge_strong_min_shared,
            bridge_weak_min_jaccard=args.bridge_weak_min_jaccard,
        )

        for comp in clusters:
            cid = f"{gene}|{locus}|{cluster_id_counter}"
            cluster_id_counter += 1

            kt, dt, st = keep_by_sample_support(
                comp=comp,
                sub=sub,
                expr_cols=expr_cols,
                sample_min_rel_expr=args.sample_min_rel_expr,
                sample_min_frac_expr=args.sample_min_frac_expr,
                retain_locus_expr_pct_per_sample=args.retain_locus_expr_pct_per_sample,
                min_abs_expr=args.min_abs_expr,
                min_cluster_expr=args.min_cluster_expr,
                min_support_samples=args.min_support_samples,
                min_keep=args.min_keep,
                tie_policy=args.tie_policy,
                tie_abs_tol=args.tie_abs_tol,
            )
            keep_set |= kt
            rescued_set |= st["rescued_by_min_keep"]

            for t in comp:
                cluster_rows.append({
                    "transcript_id": t,
                    "gene": gene,
                    "locus": locus,
                    "cluster_id": cid,
                    "cluster_size": len(comp),
                    "kept": int(t in kt),
                    "rescued_by_min_keep": int(t in st["rescued_by_min_keep"]),
                })

            if args.dropped_out and dt:
                comp_df = sub[sub["transcript_id"].isin(comp)].set_index("transcript_id")[expr_cols].astype(float)
                for t in dt:
                    vals = comp_df.loc[t].to_numpy(dtype=float)
                    dropped_rows.append({
                        "transcript_id": t,
                        "gene": gene,
                        "locus": locus,
                        "cluster_id": cid,
                        "cluster_size": len(comp),
                        "reason": st["drop_reasons"].get(t, "failed_sample_support"),
                        "sample_min_rel_expr": st["sample_min_rel_expr"],
                        "sample_min_frac_expr": st["sample_min_frac_expr"],
                        "retain_locus_expr_pct_per_sample": st["retain_locus_expr_pct_per_sample"],
                        "min_abs_expr": st["min_abs_expr"],
                        "min_cluster_expr": st["min_cluster_expr"],
                        "min_support_samples": st["min_support_samples"],
                        "expr_sum": float(np.sum(vals)),
                        "expr_max": float(np.max(vals)),
                        "expr_mean": float(np.mean(vals)),
                    })

    # Output generation
    args.out.parent.mkdir(parents=True, exist_ok=True)
    filtered = expr_df[expr_df["transcript_id"].isin(keep_set)][["transcript_id"] + expr_cols].copy()
    filtered.sort_values("transcript_id", inplace=True)
    filtered.rename(columns={"transcript_id": "#TranscriptID"}, inplace=True)
    filtered.to_csv(args.out, sep="\t", index=False)

    if args.clusters_out:
        args.clusters_out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            cluster_rows,
            columns=["transcript_id", "gene", "locus", "cluster_id", "cluster_size", "kept", "rescued_by_min_keep"]
        ).to_csv(args.clusters_out, sep="\t", index=False)

    if args.dropped_out:
        args.dropped_out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(dropped_rows).to_csv(args.dropped_out, sep="\t", index=False)

    print(f"[{ts()}] Done. Kept {len(keep_set)} transcripts across {cluster_id_counter} clusters.")
    print(f"[{ts()}]   {len(rescued_set)} of these were retained only by the --min-keep fallback "
          f"(no support criterion passed).")
    print(f"[{ts()}] Wrote: {args.out}")
    if args.clusters_out:
        print(f"[{ts()}] Wrote: {args.clusters_out}")
    if args.dropped_out:
        print(f"[{ts()}] Wrote: {args.dropped_out}")


if __name__ == "__main__":
    main()
