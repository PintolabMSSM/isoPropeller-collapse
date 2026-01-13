#!/usr/bin/env python3
"""
filter_isoforms_by_junction_cluster.py
Optimized + fixed + per-sample cumulative-expression retention

Clustering:
- Group transcripts within (gene_id, chrom|strand)
- Build splice-junction sets per transcript: junction = (exon_end, next_exon_start)
- Create clusters as connected components based on match-mode.

Filtering modes:

A) Global-per-cluster isoform-count percentile:
   - Score per transcript across samples (sum/max/mean)
   - Keep top X% transcripts per cluster by that score

B) Per-sample support filter:
   - For each sample, within each cluster:
       (1) Relative-to-max support (optional):
           expr(t,s) >= sample_min_rel_expr * max_expr_in_cluster(s)
       (2) Cumulative-expression retention:
           Sort isoforms by expr desc; keep smallest prefix whose cumulative sum
           reaches retain_locus_expr_pct_per_sample% of total cluster expr in that sample.
           (This keeps isoforms that together make up the top fraction of expression.)
   - An isoform is "supported" in a sample if it passes ANY enabled criterion.
   - Keep isoforms supported in >= min_support_samples samples.
   - Always keep at least min_keep per cluster (fallback by global score).

"""

import argparse
import gzip
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional, Any

import numpy as np
import pandas as pd


def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def smart_open(path: Path):
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
    exons_by_tx: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    tx2meta: Dict[str, Tuple[str, str, str]] = {}  # tx -> (gene, chrom, strand)

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
    for c in expr_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def jaccard(a: Set[Tuple[int, int]], b: Set[Tuple[int, int]]) -> float:
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
    inter = len(A & B)

    if mode == "any_shared":
        return inter >= min_shared
    if mode == "exact":
        return A == B
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
    if len(tx_list) <= 1:
        return [tx_list]

    if match_mode == "exact":
        groups: Dict[frozenset, List[str]] = defaultdict(list)
        for t in tx_list:
            groups[frozenset(tx2chain.get(t, []))].append(t)
        return list(groups.values())

    if min_shared < 1:
        raise ValueError("For optimized overlap-based clustering, --min-shared must be >= 1.")

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


def compute_score(expr_sub: pd.DataFrame, expr_cols: List[str], score_mode: str) -> pd.Series:
    m = expr_sub[expr_cols].to_numpy(dtype=float)
    idx = expr_sub["transcript_id"].astype(str).to_numpy()
    if score_mode == "sum":
        return pd.Series(np.sum(m, axis=1), index=idx)
    if score_mode == "max":
        return pd.Series(np.max(m, axis=1), index=idx)
    if score_mode == "mean":
        return pd.Series(np.mean(m, axis=1), index=idx)
    raise ValueError(f"Unknown score mode: {score_mode}")


def keep_by_global_percentile(
    comp: List[str],
    sub: pd.DataFrame,
    expr_cols: List[str],
    retain_top_pct: float,
    score_mode: str,
    min_keep: int,
) -> Tuple[Set[str], List[str], float]:
    comp_df = sub[sub["transcript_id"].isin(comp)].copy()
    comp_scores = compute_score(comp_df, expr_cols, score_mode).astype(float)
    comp_sorted = comp_scores.sort_values(ascending=False)

    n = len(comp_sorted)
    n_keep_target = int(np.ceil((retain_top_pct / 100.0) * n))
    n_keep = max(int(min_keep), min(n, n_keep_target))

    keep_tx = set(comp_sorted.head(n_keep).index.astype(str).tolist())
    drop_tx = [t for t in comp_sorted.index.astype(str).tolist() if t not in keep_tx]
    cutoff = float(comp_sorted.iloc[n_keep - 1]) if n_keep > 0 else float("nan")
    return keep_tx, drop_tx, cutoff


def _keepers_by_cumexpr_fraction(values: pd.Series, retain_frac: float) -> Set[str]:
    """
    values: index=transcript_id, non-negative expression values for ONE sample within ONE cluster.
    retain_frac: e.g. 0.98 to retain isoforms that cumulatively account for 98% of expression.

    Returns transcript_ids to keep for that sample.

    Behavior:
    - If total expression is 0, keep empty set (caller will enforce min_keep / fallback).
    - Otherwise, sort desc; keep minimal prefix achieving cum >= retain_frac*total.
      (Always keeps at least 1 isoform when total>0.)
    """
    v = values.clip(lower=0.0)
    total = float(v.sum())
    if total <= 0.0:
        return set()

    v_sorted = v.sort_values(ascending=False)
    csum = v_sorted.cumsum()
    cutoff = retain_frac * total

    # minimal prefix with cum >= cutoff
    k = int((csum >= cutoff).to_numpy().argmax()) + 1  # argmax of first True
    k = max(1, min(k, len(v_sorted)))

    return set(v_sorted.iloc[:k].index.astype(str).tolist())


def keep_by_sample_support(
    comp: List[str],
    sub: pd.DataFrame,
    expr_cols: List[str],
    sample_min_rel_expr: float,
    retain_locus_expr_pct_per_sample: Optional[float],
    min_support_samples: int,
    fallback_score_mode: str,
    min_keep: int,
) -> Tuple[Set[str], List[str], Dict[str, Any]]:
    """
    Per-sample support:
    - Criterion A (optional): relative-to-max within cluster per sample
    - Criterion B (optional): cumulative-expression retention within cluster per sample
                              keep isoforms that together explain top X% of expression mass
    Supported in sample if A OR B passes for that sample.
    Keep if supported in >= min_support_samples samples.
    """
    comp_df = (
        sub[sub["transcript_id"].isin(comp)]
        .set_index("transcript_id")[expr_cols]
        .astype(float)
    )

    if (sample_min_rel_expr <= 0.0) and (retain_locus_expr_pct_per_sample is None):
        raise ValueError("Sample-support filtering requires sample_min_rel_expr>0 and/or retain_locus_expr_pct_per_sample set.")

    passes = pd.DataFrame(False, index=comp_df.index, columns=comp_df.columns)

    # A) Relative-to-max support
    if sample_min_rel_expr and sample_min_rel_expr > 0:
        max_per_sample = comp_df.max(axis=0)
        thresh = max_per_sample * float(sample_min_rel_expr)
        passes |= comp_df.ge(thresh, axis=1)

    # B) Cumulative-expression retention per sample (your requested "top percentile of locus expression")
    if retain_locus_expr_pct_per_sample is not None:
        pct = float(retain_locus_expr_pct_per_sample)
        if not (0.0 < pct <= 100.0):
            raise ValueError("--retain-locus-expr-pct-per-sample must be in (0,100].")
        retain_frac = pct / 100.0

        for sample in expr_cols:
            keepers = _keepers_by_cumexpr_fraction(comp_df[sample], retain_frac=retain_frac)
            if keepers:
                passes.loc[list(keepers), sample] = True

    pass_counts = passes.sum(axis=1)
    keep_tx = set(pass_counts[pass_counts >= int(min_support_samples)].index.astype(str))

    # Ensure at least min_keep per cluster
    if len(keep_tx) < int(min_keep):
        fallback_sub = sub[sub["transcript_id"].isin(comp)].copy()
        scores = compute_score(fallback_sub, expr_cols, fallback_score_mode).astype(float)
        keep_tx = set(scores.sort_values(ascending=False).head(int(min_keep)).index.astype(str).tolist())

    drop_tx = [t for t in comp_df.index.astype(str).tolist() if t not in keep_tx]

    stats = {
        "sample_min_rel_expr": float(sample_min_rel_expr),
        "retain_locus_expr_pct_per_sample": float(retain_locus_expr_pct_per_sample) if retain_locus_expr_pct_per_sample is not None else None,
        "min_support_samples": int(min_support_samples),
        "pass_count_min": int(pass_counts.min()) if len(pass_counts) else 0,
        "pass_count_max": int(pass_counts.max()) if len(pass_counts) else 0,
    }
    return keep_tx, drop_tx, stats


def main():
    ap = argparse.ArgumentParser(description="Cluster isoforms by splice junction overlap and filter within clusters.")

    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--expr", required=True, type=Path)
    ap.add_argument("--tx-col", default="transcript_id")
    ap.add_argument("--expr-col", default="expression")

    ap.add_argument("--match-mode", choices=["any_shared", "exact", "subset", "jaccard", "bridge_safe"],
                    default="any_shared")
    ap.add_argument("--min-shared", type=int, default=1)
    ap.add_argument("--min-jaccard", type=float, default=0.5)
    ap.add_argument("--bridge-strong-min-shared", type=int, default=2)
    ap.add_argument("--bridge-weak-min-jaccard", type=float, default=0.15)

    # Global filtering (isoform COUNT percentile)
    ap.add_argument("--retain-top-pct", type=float, default=98.0,
                    help="Global mode: retain top X percent of isoforms (by score) per cluster.")
    ap.add_argument("--score", choices=["sum", "max", "mean"], default="sum")

    # Sample-support filtering (recommended)
    ap.add_argument("--sample-support-filter", action="store_true",
                    help="Enable per-sample support logic and require support in >= min-support-samples samples.")

    ap.add_argument("--sample-min-rel-expr", type=float, default=0.0,
                    help="Per sample, isoform is supported if expr >= this * (cluster max in sample). 0 disables.")

    ap.add_argument("--retain-locus-expr-pct-per-sample", type=float, default=None,
                    help=("Per sample, keep isoforms that cumulatively account for the top X%% of cluster expression "
                          "(e.g. 98 keeps isoforms whose cumulative sum reaches 98%% of total)."))

    ap.add_argument("--min-support-samples", type=int, default=1,
                    help="Keep isoforms supported in at least this many samples (sample-support mode).")

    ap.add_argument("--min-keep", type=int, default=1,
                    help="Always keep at least this many isoforms per cluster (fallback).")

    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--clusters-out", type=Path, default=None)
    ap.add_argument("--dropped-out", type=Path, default=None)

    args = ap.parse_args()

    if not (0.0 < args.retain_top_pct <= 100.0):
        raise ValueError("--retain-top-pct must be in (0,100].")
    if args.min_keep < 1:
        raise ValueError("--min-keep must be >= 1.")
    if args.sample_min_rel_expr < 0:
        raise ValueError("--sample-min-rel-expr must be >= 0.")
    if args.min_support_samples < 1:
        raise ValueError("--min-support-samples must be >= 1.")

    if args.sample_support_filter:
        if (args.sample_min_rel_expr == 0.0) and (args.retain_locus_expr_pct_per_sample is None):
            raise ValueError(
                "--sample-support-filter requires --sample-min-rel-expr > 0 and/or "
                "--retain-locus-expr-pct-per-sample set."
            )

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
    cluster_rows: List[Tuple[str, str, str, str, int]] = []
    dropped_rows: List[Dict[str, Any]] = []
    cluster_id_counter = 0

    print(f"[{ts()}] Clustering + filtering (match-mode={args.match_mode})...")
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
            for t in comp:
                cluster_rows.append((t, gene, locus, cid, len(comp)))

            if args.sample_support_filter:
                kt, dt, st = keep_by_sample_support(
                    comp=comp,
                    sub=sub,
                    expr_cols=expr_cols,
                    sample_min_rel_expr=args.sample_min_rel_expr,
                    retain_locus_expr_pct_per_sample=args.retain_locus_expr_pct_per_sample,
                    min_support_samples=args.min_support_samples,
                    fallback_score_mode=args.score,
                    min_keep=args.min_keep,
                )
                keep_set |= kt

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
                            "reason": "failed_sample_support",
                            "sample_min_rel_expr": st["sample_min_rel_expr"],
                            "retain_locus_expr_pct_per_sample": st["retain_locus_expr_pct_per_sample"],
                            "min_support_samples": st["min_support_samples"],
                            "expr_sum": float(np.sum(vals)),
                            "expr_max": float(np.max(vals)),
                            "expr_mean": float(np.mean(vals)),
                        })
            else:
                kt, dt, cutoff = keep_by_global_percentile(
                    comp=comp,
                    sub=sub,
                    expr_cols=expr_cols,
                    retain_top_pct=args.retain_top_pct,
                    score_mode=args.score,
                    min_keep=args.min_keep,
                )
                keep_set |= kt

                if args.dropped_out and dt:
                    comp_df = sub[sub["transcript_id"].isin(comp)].set_index("transcript_id")[expr_cols].astype(float)
                    scores = compute_score(sub[sub["transcript_id"].isin(comp)], expr_cols, args.score).astype(float)
                    for t in dt:
                        vals = comp_df.loc[t].to_numpy(dtype=float)
                        dropped_rows.append({
                            "transcript_id": t,
                            "gene": gene,
                            "locus": locus,
                            "cluster_id": cid,
                            "cluster_size": len(comp),
                            "reason": "below_cluster_cutoff",
                            "retain_top_pct": float(args.retain_top_pct),
                            "score_mode": args.score,
                            "score": float(scores.loc[t]),
                            "cutoff_score": float(cutoff),
                            "expr_sum": float(np.sum(vals)),
                            "expr_max": float(np.max(vals)),
                            "expr_mean": float(np.mean(vals)),
                        })

    # Write outputs
    args.out.parent.mkdir(parents=True, exist_ok=True)
    filtered = expr_df[expr_df["transcript_id"].isin(keep_set)][["transcript_id"] + expr_cols].copy()
    filtered.sort_values("transcript_id", inplace=True)
    filtered.rename(columns={"transcript_id": "#TranscriptID"}, inplace=True)  # Rename the first columns to "#TranscriptID"
    filtered.to_csv(args.out, sep="\t", index=False)

    if args.clusters_out:
        args.clusters_out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            cluster_rows,
            columns=["transcript_id", "gene", "locus", "cluster_id", "cluster_size"]
        ).to_csv(args.clusters_out, sep="\t", index=False)

    if args.dropped_out:
        args.dropped_out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(dropped_rows).to_csv(args.dropped_out, sep="\t", index=False)

    print(f"[{ts()}] Done. Kept {len(keep_set)} transcripts.")
    print(f"[{ts()}] Wrote: {args.out}")
    if args.clusters_out:
        print(f"[{ts()}] Wrote: {args.clusters_out}")
    if args.dropped_out:
        print(f"[{ts()}] Wrote: {args.dropped_out}")


if __name__ == "__main__":
    main()
