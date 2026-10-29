#!/usr/bin/env python3
"""
prune_rare_isoforms_from_clusters.py

Clusters transcripts by splice junction overlap and filters rare isoforms using
per-sample support logic (relative expression, fraction of locus expression, or
cumulative expression retention).

Changelog
---------
v2 (performance):
    - Vectorized the cumulative-expression-retention criterion
      (`--retain-locus-expr-pct-per-sample`). The original `_keepers_by_cumexpr_fraction()`
      was called once per (cluster, sample) pair in a Python loop inside
      `keep_by_sample_support()`; profiling with cProfile on a representative fixture
      (100 genes x 150 samples, 15,000 such calls) showed this was ~81% of total
      runtime (50.3s of 62.3s), dominated by per-call DataFrame construction, a
      two-key sort_values, and a label-based .loc scatter-write. Replaced with
      `_keepers_matrix_by_cumexpr_fraction()`, which computes the identical result
      for every sample column of a cluster in one batch of numpy operations (no
      per-sample Python loop). Verified mathematically identical to the original via
      exhaustive testing (40 random fuzz cases plus 6 hand-crafted edge cases --
      all-zero columns, single transcript, exact multi-way ties, negative values,
      retain_frac=1.0 float-boundary cases, all-identical values -- each checked
      across all 3 --tie-policy values and 4 retain_frac values) before being
      adopted here, and again validated end-to-end (byte-identical --out/
      --clusters-out/--dropped-out) against the original per-sample implementation
      across representative parameter combinations.
    - `load_expression_table()` now sniffs the delimiter from the file's header line
      and uses pandas' fast `engine="c"` parser when it's unambiguously comma or
      tab-separated, falling back to the slower auto-detecting `engine="python"`
      parser only when the sniff is inconclusive. Confirmed byte-identical output to
      the original unconditional `engine="python"` parsing, ~3.5x faster on a
      150-sample test table.
    - `keep_by_sample_support()` no longer rebuilds a fresh `comp_df` a second time
      for the dropped-isoform report in `main()`; it now returns the one it already
      built via `stats["comp_df"]`, which `main()` reuses.
    - `build_clusters_for_group()` now iterates `sorted(possible_pairs)` instead of
      raw `set` order when building the edge list fed to `connected_components()`.
      Python's per-process string-hash randomization otherwise made the row order
      within a cluster in `--clusters-out`/`--dropped-out` vary nondeterministically
      from run to run (cluster membership and kept/dropped status were never
      affected, only display order) -- the same class of bug previously found and
      fixed in a sibling script.
    - Measured end-to-end on a 600-gene x 150-sample fixture: wall-clock time for
      the full CLI call dropped from 2m33.2s to 1.94s (~79x faster). Output
      (`--out`, `--clusters-out`, `--dropped-out`) is identical to the pre-v2
      script -- confirmed across 18 parameter combinations (every `--tie-policy`,
      every `--match-mode`, with/without `--protect-isoform-gtf`, with/without
      noise floors, `--min-keep` 0 and nonzero, `--min-support-samples` > 1) plus
      the full 600-gene fixture itself.

v3 (multiprocessing):
    - Added `--threads` support following the same approach
      used in the sibling scripts (consolidate_splice_chain_fragments_v2.py,
      prune_rare_isoforms_from_clusters_v4.py on the old/superseded branch): an
      `mp.Pool` of worker processes, each handling whole (gene, locus) groups --
      the natural unit of work, since clustering + filtering is already done
      independently per group.
    - The entire per-(gene, locus)-group logic (clustering, filtering, reference-
      protection bookkeeping, dropped-row construction) is factored into a single
      parameter-driven function, `_process_group_core()`, that takes every input
      as an explicit argument and references no module-level global. Both the
      `--threads 1` (serial, default) and `--threads > 1` code paths in `main()`
      call this one function -- there is no second, duplicated inline copy that
      could silently drift out of sync (this is exactly how a sibling script once
      introduced a real bug: a serial branch that duplicated per-group logic
      inline accidentally referenced a multiprocessing-only global instead of its
      own local parameter). The `--threads > 1` path reads its inputs from
      module-level globals populated once per worker via `_mp_init()` (the
      standard `mp.Pool(initializer=...)` pattern), purely so read-only inputs
      (`tx2chain` in particular) aren't re-pickled on every one of the
      (gene, locus) tasks.
    - `cluster_id` assignment (used in `--clusters-out`/`--dropped-out`, formatted
      as `"{gene}|{locus}|{N}"` where N increments once per cluster found,
      continuing across the whole run) depends on the order (gene, locus) groups
      are consumed, so the multiprocessing branch uses `pool.imap()`
      (order-preserving), not `pool.imap_unordered()`: `--threads` > 1 reproduces
      `--threads 1` output byte-for-byte -- including cluster IDs and row order --
      not merely the same keep/drop decisions.
    - Added `--log-every` (default 1000; 0 disables) for progress logging, same
      option name/semantics as the sibling scripts.
    - Validated: `--threads 1` output is byte-identical to v2 across the same 18
      parameter combinations used for the v2 validation (confirming the refactor
      into `_process_group_core()` introduced no behavior change), and
      `--threads 2` / `--threads 4` output is byte-identical to `--threads 1`
      (including `cluster_id` values and row order) across 6 representative
      combinations (default, `truncate_ties`, all three support criteria combined,
      `--match-mode gene_id`, `--protect-isoform-gtf`, `--min-keep 0`), plus the
      full 600-gene fixture.
    - Performance note: this sandbox's test environment has only 2 CPU cores, and
      the v2 vectorization already made a single thread fast (the 600-gene
      fixture runs in ~2s), so the headroom `--threads` has to work with here is
      small -- `--threads 2` measured ~1.4s vs. ~2.0s at `--threads 1` on that
      same fixture (byte-identical output). The real benefit is expected to be
      larger on your actual runs: a much bigger gene/locus count and more
      available cores than this sandbox's 2. `--threads 1` remains the default,
      so existing invocations are unaffected unless `--threads` is passed
      explicitly.
"""

import argparse
import gzip
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional, Any, TextIO

import numpy as np
import pandas as pd

try:
    import multiprocessing as mp
except Exception:
    mp = None


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


def _sniff_sep(path: Path) -> Optional[str]:
    """
    Sniff a single-character field delimiter (comma or tab) from a file's header line.

    This lets `load_expression_table` use pandas' fast C parser (`engine="c"`) for the
    common comma/tab cases instead of always paying for the much slower
    `engine="python"` auto-detecting parser.

    Parameters
    ----------
    path : Path
        Path to the table file (can be gzip-compressed).

    Returns
    -------
    Optional[str]
        `","` or `"\\t"` when the header line unambiguously favors one delimiter,
        otherwise `None` (caller should fall back to `engine="python"` sniffing).
    """
    with smart_open(path) as fh:
        header = fh.readline()
    comma = header.count(",")
    tab = header.count("\t")
    if tab > comma:
        return "\t"
    if comma > tab:
        return ","
    return None


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
    sep = _sniff_sep(path)
    if sep is not None:
        df = pd.read_csv(path, sep=sep, engine="c")
    else:
        # Inconclusive sniff (e.g. single-column file, or header has equal
        # counts of comma/tab): fall back to pandas' own slower delimiter
        # sniffer rather than guessing wrong.
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


def resolve_protected_transcripts(
    tx2chain: Dict[str, List[Tuple[int, int]]],
    reference_gtf: Path,
    transcript_attr: str = "transcript_id",
    gene_attr: str = "gene_id",
) -> Tuple[Set[str], int]:
    """
    Identify transcripts whose splice chain exactly matches a reference transcript.

    A transcript is protected if its ordered junction chain is identical to the chain of
    some transcript in the reference annotation. Protected transcripts are never removed
    by any filtering criterion, regardless of expression.

    Parameters
    ----------
    tx2chain : Dict[str, List[Tuple[int, int]]]
        Junction chains for the query transcripts (from `parse_gtf`).
    reference_gtf : Path
        Path to the reference GTF (can be gzip-compressed).
    transcript_attr : str, default="transcript_id"
        GTF attribute key for transcript identifier in the reference.
    gene_attr : str, default="gene_id"
        GTF attribute key for gene identifier in the reference.

    Returns
    -------
    protected_txs : Set[str]
        Query transcript IDs matching a reference splice chain.
    n_reference_chains : int
        Number of distinct non-empty reference chains loaded.

    Notes
    -----
    Single-exon transcripts have an empty junction chain. An empty chain would match every
    single-exon transcript in the annotation, so empty chains are excluded on both sides:
    single-exon transcripts are never protected by this mechanism.

    Matching is on junction coordinates alone, not on chromosome or strand. Junction
    coordinates are genomic, so an accidental cross-locus match requires every donor and
    acceptor to coincide exactly; this is effectively impossible for multi-junction chains
    but is not formally excluded for single-junction ones.
    """
    ref_tx2chain, _, _ = parse_gtf(reference_gtf, transcript_attr=transcript_attr, gene_attr=gene_attr)

    reference_chains: Set[Tuple[Tuple[int, int], ...]] = set()
    for chain in ref_tx2chain.values():
        if chain:
            reference_chains.add(tuple(chain))

    protected_txs: Set[str] = set()
    if reference_chains:
        for tx, chain in tx2chain.items():
            if chain and tuple(chain) in reference_chains:
                protected_txs.add(str(tx))

    return protected_txs, len(reference_chains)


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
        Clustering matching strategy (`'bridge_safe'`, `'any_shared'`, `'subset'`, `'jaccard'`, or `'gene_id'`).
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
    if match_mode == "gene_id":
        return [tx_list]

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

    # Iterate in a fixed order rather than raw set order: Python's per-process string-hash
    # randomization (PYTHONHASHSEED) makes set iteration order vary run-to-run, which would
    # otherwise make the edge list (and therefore connected_components' traversal order,
    # and therefore the row order within a cluster in --clusters-out) nondeterministic
    # between runs/processes even though cluster membership and kept/dropped status are
    # unaffected either way.
    edges: Dict[str, List[str]] = defaultdict(list)
    for a, b in sorted(possible_pairs):
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


def _keepers_matrix_by_cumexpr_fraction(
    df: pd.DataFrame,
    retain_frac: float,
    tie_policy: str = "include_all_ties",
    tie_abs_tol: float = 1e-9,
) -> pd.DataFrame:
    """
    Select transcripts achieving cumulative expression retention target, for every
    sample column of a cluster at once.

    This is the vectorized, all-columns-at-once replacement for the original
    per-sample `_keepers_by_cumexpr_fraction()` (previously called once per
    `(cluster, sample)` pair in a Python loop, which profiling showed was the
    dominant cost of the whole script: ~81% of total runtime on a representative
    fixture). It is mathematically identical to calling the per-sample version once
    per column and assembling the results into a DataFrame -- this was verified
    exhaustively (random fuzz cases plus hand-crafted edge cases, across every
    `tie_policy` and a range of `retain_frac` values) before being adopted here.

    The per-column (value descending, transcript_id ascending) tie-break order
    needed to replicate `sort_values(by=["v", "tx"], ascending=[False, True],
    kind="mergesort")` is reproduced for all columns in one pass via a
    stable-sort-on-pre-sorted-input trick: rows are first ordered by transcript_id
    ascending once, then each column is stable-sorted by descending value, which
    preserves the transcript_id-ascending order among exact ties.

    Parameters
    ----------
    df : pd.DataFrame
        Expression values for one cluster, indexed by transcript ID, with one column
        per sample. May contain negative values (clipped to zero, as in the original).
    retain_frac : float
        Target fraction of total cluster expression mass to retain per sample (e.g.
        0.98 for 98%).
    tie_policy : str, default="include_all_ties"
        Policy for handling boundary ties (see `keep_by_sample_support`):
        - `'include_all_ties'`: Retain all non-zero isoforms tied at cutoff value.
        - `'exclude_all_ties'`: Exclude boundary ties entirely (yielding no support if all tied).
        - `'truncate_ties'`: Keep minimal prefix hitting cumulative target, breaking ties by transcript ID.
    tie_abs_tol : float, default=1e-9
        Absolute tolerance for float comparison around the boundary cutoff value.

    Returns
    -------
    pd.DataFrame
        Boolean DataFrame, same index and columns as `df`, True where that transcript
        is supported in that sample by the cumulative-retention criterion.
    """
    tx_ids = df.index.astype(str).to_numpy()
    values = df.to_numpy(dtype=float)
    n, m = values.shape

    clipped = np.clip(values, 0.0, None)
    totals = clipped.sum(axis=0)

    pre = np.argsort(tx_ids, kind="stable")
    clipped_pre = clipped[pre, :]
    order_within_pre = np.argsort(-clipped_pre, axis=0, kind="stable")
    order = pre[order_within_pre]

    sorted_vals = np.take_along_axis(clipped, order, axis=0)
    csum = np.cumsum(sorted_vals, axis=0)
    cutoff = retain_frac * totals

    with np.errstate(invalid="ignore"):
        reached = csum >= cutoff[None, :]
    any_true = reached.any(axis=0)
    first_true = reached.argmax(axis=0)  # 0 when no True; overridden below via no_crossing
    cutoff_val = sorted_vals[first_true, np.arange(m)]

    # Inverse permutation: rank[row, col] = position of that row in column col's
    # (value desc, tx asc) sorted order -- needed only for truncate_ties.
    rank = np.empty((n, m), dtype=np.int64)
    pos = np.broadcast_to(np.arange(n)[:, None], (n, m))
    np.put_along_axis(rank, order, pos, axis=0)

    if tie_policy == "truncate_ties":
        mask = (rank <= first_true[None, :]) & (clipped > 0.0)
    elif tie_policy == "exclude_all_ties":
        mask = (
            (clipped > cutoff_val[None, :]) & ~np.isclose(clipped, cutoff_val[None, :], rtol=0.0, atol=tie_abs_tol)
        ) & (clipped > 0.0)
    else:  # "include_all_ties"
        mask = (
            (clipped >= cutoff_val[None, :]) | np.isclose(clipped, cutoff_val[None, :], rtol=0.0, atol=tie_abs_tol)
        ) & (clipped > 0.0)

    no_crossing = ~any_true  # float-edge-case fallback (e.g. retain_frac == 1.0): keep all nonzero
    fallback_mask = clipped > 0.0
    mask = np.where(no_crossing[None, :], fallback_mask, mask)

    zero_total = totals <= 0.0  # total<=0 columns: nobody passes
    mask = np.where(zero_total[None, :], False, mask)

    return pd.DataFrame(mask, index=df.index, columns=df.columns)


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
    protected_txs: Optional[Set[str]] = None,
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
    protected_txs : Set[str], optional
        Globally protected transcript IDs (reference splice-chain matches). Members of this
        set within the cluster are retained unconditionally. Protection is applied after the
        support criteria are evaluated, so protected transcripts still contribute to cluster
        maxima, totals, and cumulative cutoffs exactly as they would without protection.

    Returns
    -------
    keep_tx : Set[str]
        Set of kept transcript IDs.
    drop_tx : List[str]
        List of dropped transcript IDs.
    stats : Dict[str, Any]
        Dictionary summarizing filtering threshold metrics, per-transcript drop reasons,
        transcripts retained only by the `min_keep` fallback, and the cluster's
        `comp_df` expression slice (so callers needing it, e.g. for dropped-isoform
        reporting, don't have to rebuild it).
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

    # 3. Cumulative expression retention per sample (vectorized across all
    # sample columns at once -- see _keepers_matrix_by_cumexpr_fraction docstring)
    if retain_locus_expr_pct_per_sample is not None:
        retain_frac = float(retain_locus_expr_pct_per_sample) / 100.0
        passes |= _keepers_matrix_by_cumexpr_fraction(
            comp_df,
            retain_frac=retain_frac,
            tie_policy=tie_policy,
            tie_abs_tol=tie_abs_tol,
        )

    # Track pre-noise-floor passes to distinguish filtering causes
    pass_before_floors = passes.copy()

    # Enforce noise floors (minimum absolute expression & cluster total)
    if min_abs_expr > 0.0:
        passes &= comp_df.ge(min_abs_expr, axis=1)
    if min_cluster_expr > 0.0:
        passes &= (total_per_sample >= min_cluster_expr)

    pass_counts = passes.sum(axis=1)
    supported_tx = set(pass_counts[pass_counts >= int(min_support_samples)].index.astype(str))

    # Reference-protected transcripts in this cluster are retained unconditionally,
    # whether or not any support criterion passed.
    cluster_txs = set(comp_df.index.astype(str))
    protected_here: Set[str] = (set(protected_txs) & cluster_txs) if protected_txs else set()

    # Keep all supported isoforms, plus anything protected
    keep_tx = set(supported_tx) | protected_here
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
        "reference_protected": protected_here,
        "reference_protected_only": protected_here - supported_tx,
        "drop_reasons": drop_reasons,
        "comp_df": comp_df,
    }
    return keep_tx, drop_tx, stats


def _process_group_core(
    gene: str,
    locus: str,
    sub: pd.DataFrame,
    tx2chain: Dict[str, List[Tuple[int, int]]],
    expr_cols: List[str],
    protected_txs: Set[str],
    match_mode: str,
    min_shared: int,
    min_jaccard: float,
    bridge_strong_min_shared: int,
    bridge_weak_min_jaccard: float,
    sample_min_rel_expr: float,
    sample_min_frac_expr: float,
    retain_locus_expr_pct_per_sample: Optional[float],
    min_abs_expr: float,
    min_cluster_expr: float,
    min_support_samples: int,
    min_keep: int,
    tie_policy: str,
    tie_abs_tol: float,
    want_dropped: bool,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Cluster and filter all transcripts belonging to a single (gene, locus) group.

    This is the single, parameter-driven implementation used by BOTH the serial
    (`--threads 1`, default) and multiprocessing (`--threads > 1`) code paths in
    `main()` -- there is no second, duplicated copy of this logic that could drift
    out of sync between the two paths. (A sibling script, consolidate_splice_chain_
    fragments_v2.py, once had exactly this kind of bug: its serial branch
    duplicated per-group logic inline and accidentally referenced a
    multiprocessing-only global instead of the correctly-scoped local parameter,
    silently disabling an option at the default thread count. Factoring the group
    logic into one function with no references to module-level state, called
    identically from both code paths, makes that class of bug structurally
    impossible here.)

    Parameters
    ----------
    gene, locus : str
        Group key, used only to label output rows.
    sub : pd.DataFrame
        The expression rows belonging to this (gene, locus) group.
    tx2chain, expr_cols, protected_txs : see `main`.
    match_mode, min_shared, min_jaccard, bridge_strong_min_shared, bridge_weak_min_jaccard :
        Passed through to `build_clusters_for_group`.
    sample_min_rel_expr, sample_min_frac_expr, retain_locus_expr_pct_per_sample,
    min_abs_expr, min_cluster_expr, min_support_samples, min_keep, tie_policy, tie_abs_tol :
        Passed through to `keep_by_sample_support`.
    want_dropped : bool
        Whether to build the per-transcript dropped-row records (skipped entirely
        when `--dropped-out` was not requested, since they're otherwise unused).

    Returns
    -------
    cluster_results : List[Dict[str, Any]]
        One dict per cluster found in this group:
        `{"comp": [...], "keep_tx": {...}, "rescued_by_min_keep": {...},
        "reference_protected": {...}, "reference_protected_only": {...},
        "dropped": [row_dict, ...]}`. Each `dropped` row dict omits `"cluster_id"`
        -- that is assigned by the caller, since it is the one piece of bookkeeping
        that depends on the order (gene, locus) groups are consumed in, not on
        anything specific to this group.
    n_tx : int
        Number of transcripts in this group (for progress logging).
    """
    tx_list = sub["transcript_id"].astype(str).tolist()
    clusters = build_clusters_for_group(
        tx_list=tx_list,
        tx2chain=tx2chain,
        match_mode=match_mode,
        min_shared=min_shared,
        min_jaccard=min_jaccard,
        bridge_strong_min_shared=bridge_strong_min_shared,
        bridge_weak_min_jaccard=bridge_weak_min_jaccard,
    )

    cluster_results: List[Dict[str, Any]] = []

    for comp in clusters:
        kt, dt, st = keep_by_sample_support(
            comp=comp,
            sub=sub,
            expr_cols=expr_cols,
            sample_min_rel_expr=sample_min_rel_expr,
            sample_min_frac_expr=sample_min_frac_expr,
            retain_locus_expr_pct_per_sample=retain_locus_expr_pct_per_sample,
            min_abs_expr=min_abs_expr,
            min_cluster_expr=min_cluster_expr,
            min_support_samples=min_support_samples,
            min_keep=min_keep,
            tie_policy=tie_policy,
            tie_abs_tol=tie_abs_tol,
            protected_txs=protected_txs,
        )

        dropped_list: List[Dict[str, Any]] = []
        if want_dropped and dt:
            comp_df = st["comp_df"]
            for t in dt:
                vals = comp_df.loc[t].to_numpy(dtype=float)
                dropped_list.append({
                    "transcript_id": t,
                    "gene": gene,
                    "locus": locus,
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

        cluster_results.append({
            "comp": comp,
            "keep_tx": kt,
            "rescued_by_min_keep": st["rescued_by_min_keep"],
            "reference_protected": st["reference_protected"],
            "reference_protected_only": st["reference_protected_only"],
            "dropped": dropped_list,
        })

    return cluster_results, len(tx_list)


def _insert_cluster_id(d: Dict[str, Any], cluster_id: str) -> Dict[str, Any]:
    """
    Return a copy of a dropped-row dict with `"cluster_id"` inserted right after
    `"locus"`, reproducing the exact column order of the single-process
    implementation regardless of `--threads`.

    Parameters
    ----------
    d : Dict[str, Any]
        A dropped-row dict as built by `_process_group_core` (which omits
        `"cluster_id"` since that depends on cross-group processing order).
    cluster_id : str
        The cluster ID to insert, assigned by the caller in `main`.

    Returns
    -------
    Dict[str, Any]
        A new dict with `"cluster_id"` inserted after `"locus"`.
    """
    out: Dict[str, Any] = {}
    for k, v in d.items():
        out[k] = v
        if k == "locus":
            out["cluster_id"] = cluster_id
    return out


# ---------------------------------------------------------------------------
# Multiprocessing helpers
# ---------------------------------------------------------------------------
# These module-level globals exist ONLY so _process_group_mp() (the mp.Pool worker
# function) can read the run's read-only inputs without re-pickling them on every
# single (gene, locus) task sent to a worker. They are populated exactly once per
# worker process, via _mp_init() as mp.Pool's initializer, and are never read by
# _process_group_core() itself or by the serial (--threads 1) code path in main()
# -- see _process_group_core's docstring for why that separation matters.
_G_tx2chain = None
_G_expr_cols = None
_G_protected_txs = None
_G_match_mode = None
_G_min_shared = None
_G_min_jaccard = None
_G_bridge_strong_min_shared = None
_G_bridge_weak_min_jaccard = None
_G_sample_min_rel_expr = None
_G_sample_min_frac_expr = None
_G_retain_locus_expr_pct_per_sample = None
_G_min_abs_expr = None
_G_min_cluster_expr = None
_G_min_support_samples = None
_G_min_keep = None
_G_tie_policy = None
_G_tie_abs_tol = None
_G_want_dropped = None


def _mp_init(
    tx2chain, expr_cols, protected_txs, match_mode, min_shared, min_jaccard,
    bridge_strong_min_shared, bridge_weak_min_jaccard, sample_min_rel_expr,
    sample_min_frac_expr, retain_locus_expr_pct_per_sample, min_abs_expr,
    min_cluster_expr, min_support_samples, min_keep, tie_policy, tie_abs_tol,
    want_dropped,
):
    """mp.Pool initializer: populate this worker process's copy of the run's
    read-only inputs exactly once, instead of re-pickling them with every task."""
    global _G_tx2chain, _G_expr_cols, _G_protected_txs, _G_match_mode, _G_min_shared
    global _G_min_jaccard, _G_bridge_strong_min_shared, _G_bridge_weak_min_jaccard
    global _G_sample_min_rel_expr, _G_sample_min_frac_expr, _G_retain_locus_expr_pct_per_sample
    global _G_min_abs_expr, _G_min_cluster_expr, _G_min_support_samples, _G_min_keep
    global _G_tie_policy, _G_tie_abs_tol, _G_want_dropped
    _G_tx2chain = tx2chain
    _G_expr_cols = expr_cols
    _G_protected_txs = protected_txs
    _G_match_mode = match_mode
    _G_min_shared = min_shared
    _G_min_jaccard = min_jaccard
    _G_bridge_strong_min_shared = bridge_strong_min_shared
    _G_bridge_weak_min_jaccard = bridge_weak_min_jaccard
    _G_sample_min_rel_expr = sample_min_rel_expr
    _G_sample_min_frac_expr = sample_min_frac_expr
    _G_retain_locus_expr_pct_per_sample = retain_locus_expr_pct_per_sample
    _G_min_abs_expr = min_abs_expr
    _G_min_cluster_expr = min_cluster_expr
    _G_min_support_samples = min_support_samples
    _G_min_keep = min_keep
    _G_tie_policy = tie_policy
    _G_tie_abs_tol = tie_abs_tol
    _G_want_dropped = want_dropped


def _process_group_mp(item):
    """mp.Pool worker: process one (gene, locus) group using this worker's
    _mp_init()-populated globals, calling the same _process_group_core() used by
    the serial path."""
    (gene, locus), sub = item
    cluster_results, n_tx = _process_group_core(
        gene, locus, sub,
        tx2chain=_G_tx2chain,
        expr_cols=_G_expr_cols,
        protected_txs=_G_protected_txs,
        match_mode=_G_match_mode,
        min_shared=_G_min_shared,
        min_jaccard=_G_min_jaccard,
        bridge_strong_min_shared=_G_bridge_strong_min_shared,
        bridge_weak_min_jaccard=_G_bridge_weak_min_jaccard,
        sample_min_rel_expr=_G_sample_min_rel_expr,
        sample_min_frac_expr=_G_sample_min_frac_expr,
        retain_locus_expr_pct_per_sample=_G_retain_locus_expr_pct_per_sample,
        min_abs_expr=_G_min_abs_expr,
        min_cluster_expr=_G_min_cluster_expr,
        min_support_samples=_G_min_support_samples,
        min_keep=_G_min_keep,
        tie_policy=_G_tie_policy,
        tie_abs_tol=_G_tie_abs_tol,
        want_dropped=_G_want_dropped,
    )
    return gene, locus, cluster_results, n_tx


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
    cluster_group.add_argument("--match-mode", choices=["bridge_safe", "any_shared", "subset", "jaccard", "gene_id"], default="bridge_safe",
                               help="Clustering mode. Defaults to 'bridge_safe' to prevent cluster collapsing.")
    cluster_group.add_argument("--min-shared", type=int, default=1)
    cluster_group.add_argument("--min-jaccard", type=float, default=0.5)
    cluster_group.add_argument("--bridge-strong-min-shared", type=int, default=2)
    cluster_group.add_argument("--bridge-weak-min-jaccard", type=float, default=0.15)

    # Reference protection
    ref_group = ap.add_argument_group("Isoform Protection")
    ref_group.add_argument("--protect-isoform-gtf", "--reference-gtf",
                           dest="protect_isoform_gtf", type=Path, default=None,
                           help=("Optional GTF of isoforms to protect. Transcripts whose splice chain "
                                 "exactly matches an isoform in this GTF are PROTECTED from all filtering "
                                 "and always kept. Single-exon transcripts (empty chain) are never "
                                 "protected by this mechanism. (--reference-gtf is a deprecated alias.)"))

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

    # Performance Options
    perf_group = ap.add_argument_group("Performance Options")
    perf_group.add_argument("--threads", type=int, default=1,
                             help="Number of worker processes for gene-locus-level clustering+filtering (default: 1).")
    perf_group.add_argument("--log-every", type=int, default=1000,
                             help="Print progress every N gene-loci groups (default: 1000; 0 disables progress logging).")

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

    protected_txs: Set[str] = set()
    if args.protect_isoform_gtf:
        print(f"[{ts()}] Parsing reference GTF for known splice chains...")
        protected_txs, n_ref_chains = resolve_protected_transcripts(tx2chain, args.protect_isoform_gtf)
        print(f"[{ts()}]   {n_ref_chains} distinct reference splice chains; "
              f"{len(protected_txs)} transcripts protected from filtering.")
        if n_ref_chains and not protected_txs:
            print(f"[{ts()}]   WARNING: no transcript matched a reference chain. Check that both GTFs "
                  f"use the same coordinate system and assembly.")

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
    protected_kept: Set[str] = set()
    protected_only: Set[str] = set()
    cluster_rows: List[Dict[str, Any]] = []
    dropped_rows: List[Dict[str, Any]] = []
    cluster_id_counter = 0

    want_dropped = args.dropped_out is not None
    threads = max(1, args.threads) if mp is not None else 1

    items = list(expr_df.groupby(["__gene__", "__locus__"], sort=False))
    total_groups = len(items)
    total_tx = len(expr_df)
    processed_groups = 0
    processed_tx_running = 0

    print(f"[{ts()}] Clustering + filtering: {total_tx} transcripts across {total_groups} gene-loci groups "
          f"(match-mode={args.match_mode}, tie-policy={args.tie_policy}, threads={threads})...")

    def _consume_group_result(gene, locus, cluster_results):
        nonlocal cluster_id_counter
        for cr in cluster_results:
            comp = cr["comp"]
            kt = cr["keep_tx"]
            cid = f"{gene}|{locus}|{cluster_id_counter}"
            cluster_id_counter += 1

            keep_set.update(kt)
            rescued_set.update(cr["rescued_by_min_keep"])
            protected_kept.update(cr["reference_protected"])
            protected_only.update(cr["reference_protected_only"])

            for t in comp:
                cluster_rows.append({
                    "transcript_id": t,
                    "gene": gene,
                    "locus": locus,
                    "cluster_id": cid,
                    "cluster_size": len(comp),
                    "kept": int(t in kt),
                    "reference_protected": int(t in cr["reference_protected"]),
                    "rescued_by_min_keep": int(t in cr["rescued_by_min_keep"]),
                })

            for d in cr["dropped"]:
                dropped_rows.append(_insert_cluster_id(d, cid))

    if threads > 1 and mp is not None:
        with mp.Pool(
            processes=threads, initializer=_mp_init,
            initargs=(
                tx2chain, expr_cols, protected_txs,
                args.match_mode, args.min_shared, args.min_jaccard,
                args.bridge_strong_min_shared, args.bridge_weak_min_jaccard,
                args.sample_min_rel_expr, args.sample_min_frac_expr,
                args.retain_locus_expr_pct_per_sample, args.min_abs_expr,
                args.min_cluster_expr, args.min_support_samples, args.min_keep,
                args.tie_policy, args.tie_abs_tol, want_dropped,
            ),
        ) as pool:
            # Ordered imap() (NOT imap_unordered()): cluster_id assignment in
            # _consume_group_result depends on consuming (gene, locus) groups in
            # the same order they were grouped in, so --threads > 1 reproduces
            # --threads 1 output byte-for-byte, not merely the same keep/drop
            # decisions.
            for gene, locus, cluster_results, n_tx in pool.imap(_process_group_mp, items, chunksize=8):
                processed_groups += 1
                processed_tx_running += n_tx
                _consume_group_result(gene, locus, cluster_results)
                if args.log_every > 0 and (processed_groups % args.log_every == 0 or processed_groups == total_groups):
                    pct = 100.0 * processed_groups / max(1, total_groups)
                    print(f"[{ts()}] Processed {processed_groups} / {total_groups} gene-loci groups ({pct:.1f}%) "
                          f"-- {processed_tx_running} transcripts scanned; {len(keep_set)} kept so far")
    else:
        for (gene, locus), sub in items:
            cluster_results, n_tx = _process_group_core(
                gene, locus, sub,
                tx2chain=tx2chain,
                expr_cols=expr_cols,
                protected_txs=protected_txs,
                match_mode=args.match_mode,
                min_shared=args.min_shared,
                min_jaccard=args.min_jaccard,
                bridge_strong_min_shared=args.bridge_strong_min_shared,
                bridge_weak_min_jaccard=args.bridge_weak_min_jaccard,
                sample_min_rel_expr=args.sample_min_rel_expr,
                sample_min_frac_expr=args.sample_min_frac_expr,
                retain_locus_expr_pct_per_sample=args.retain_locus_expr_pct_per_sample,
                min_abs_expr=args.min_abs_expr,
                min_cluster_expr=args.min_cluster_expr,
                min_support_samples=args.min_support_samples,
                min_keep=args.min_keep,
                tie_policy=args.tie_policy,
                tie_abs_tol=args.tie_abs_tol,
                want_dropped=want_dropped,
            )
            processed_groups += 1
            processed_tx_running += n_tx
            _consume_group_result(gene, locus, cluster_results)
            if args.log_every > 0 and (processed_groups % args.log_every == 0 or processed_groups == total_groups):
                pct = 100.0 * processed_groups / max(1, total_groups)
                print(f"[{ts()}] Processed {processed_groups} / {total_groups} gene-loci groups ({pct:.1f}%) "
                      f"-- {processed_tx_running} transcripts scanned; {len(keep_set)} kept so far")

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
            columns=["transcript_id", "gene", "locus", "cluster_id", "cluster_size", "kept",
                     "reference_protected", "rescued_by_min_keep"]
        ).to_csv(args.clusters_out, sep="\t", index=False)

    if args.dropped_out:
        args.dropped_out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(dropped_rows).to_csv(args.dropped_out, sep="\t", index=False)

    print(f"[{ts()}] Done. Kept {len(keep_set)} transcripts across {cluster_id_counter} clusters.")
    print(f"[{ts()}]   {len(rescued_set)} of these were retained only by the --min-keep fallback "
          f"(no support criterion passed).")
    if args.protect_isoform_gtf:
        print(f"[{ts()}]   {len(protected_kept)} were reference-protected "
              f"({len(protected_only)} of which no support criterion passed).")
    print(f"[{ts()}] Wrote: {args.out}")
    if args.clusters_out:
        print(f"[{ts()}] Wrote: {args.clusters_out}")
    if args.dropped_out:
        print(f"[{ts()}] Wrote: {args.dropped_out}")


if __name__ == "__main__":
    main()
