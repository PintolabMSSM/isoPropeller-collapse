#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional

import numpy as np
import pandas as pd


def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


def smart_open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gtf(
    gtf_path: Path,
    transcript_attr: str = "transcript_id",
    gene_attr: str = "gene_id",
) -> Tuple[
    Dict[str, Tuple[Tuple[int, int], ...]],  # tx2chain (tuple of junctions for hashing)
    Dict[str, Set[Tuple[int, int]]],         # tx2junction_set (pre-computed sets)
    Dict[str, str],                          # tx2gene
    Dict[str, str],                          # tx2locus "chrom|strand"
]:
    exons_by_tx: Dict[str, List[Tuple[str, str, int, int]]] = defaultdict(list)
    tx2gene: Dict[str, str] = {}

    with smart_open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != "exon":
                continue

            # Robust attribute parsing
            attr_dict = {}
            for field in attrs.split(";"):
                field = field.strip()
                if not field:
                    continue
                if " " in field:
                    key, val = field.split(" ", 1)
                elif "=" in field:
                    key, val = field.split("=", 1)
                else:
                    continue
                val = val.strip().strip('"').strip("'")
                attr_dict[key] = val

            tx = attr_dict.get(transcript_attr)
            gene = attr_dict.get(gene_attr)
            if tx is None or gene is None:
                continue

            s, e = int(start), int(end)
            exons_by_tx[tx].append((chrom, strand, s, e))
            tx2gene[tx] = gene

    tx2chain: Dict[str, Tuple[Tuple[int, int], ...]] = {}
    tx2set: Dict[str, Set[Tuple[int, int]]] = {}
    tx2locus: Dict[str, str] = {}

    for tx, exons in exons_by_tx.items():
        # Ensure all exons are on same chrom/strand
        chroms = {c for c, _, _, _ in exons}
        strands = {st for _, st, _, _ in exons}
        if len(chroms) != 1 or len(strands) != 1:
            continue
        chrom = next(iter(chroms))
        strand = next(iter(strands))

        # Sort by genomic coordinate to identify introns
        exons_sorted = sorted(exons, key=lambda x: x[2])
        juncs: List[Tuple[int, int]] = []
        for i in range(len(exons_sorted) - 1):
            _, _, s_i, e_i = exons_sorted[i]
            _, _, s_j, e_j = exons_sorted[i + 1]
            # Intron is defined as (End_Exon_1, Start_Exon_2)
            donor = e_i
            acceptor = s_j
            if donor < acceptor:
                juncs.append((donor, acceptor))

        # Convert to tuple for storage/hashing and set for fast comparison
        tx2chain[tx] = tuple(juncs)
        tx2set[tx] = set(juncs)
        tx2locus[tx] = f"{chrom}|{strand}"

    return tx2chain, tx2set, tx2gene, tx2locus


def load_expression_table(path: Path, tx_col: str, expr_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    cols_lower = {c.lower(): c for c in df.columns}

    tx_c = cols_lower.get(tx_col.lower(), tx_col)
    if tx_c not in df.columns:
        raise ValueError(f"Could not find transcript id column '{tx_col}' in {list(df.columns)}")

    if expr_col == "*" or expr_col.upper() == "ALL":
        # Select all columns except ID, but filter for numeric types only
        potential_cols = [c for c in df.columns if c != tx_c]
        numeric_df = df[potential_cols].select_dtypes(include=[np.number])
        expr_cols = numeric_df.columns.tolist()
        if len(expr_cols) < len(potential_cols):
            print(f"[WARN] Dropped {len(potential_cols) - len(expr_cols)} non-numeric columns from expression consideration.")
    elif "," in expr_col:
        expr_cols = [c.strip() for c in expr_col.split(",") if c.strip()]
        for c in expr_cols:
            if c not in df.columns:
                raise ValueError(f"Requested expr column '{c}' not in columns {list(df.columns)}")
    else:
        expr_c = cols_lower.get(expr_col.lower(), expr_col)
        if expr_c not in df.columns:
            raise ValueError(f"Could not find expression column '{expr_col}' in {list(df.columns)}")
        expr_cols = [expr_c]

    out = df[[tx_c] + expr_cols].rename(columns={tx_c: "transcript_id"}).copy()
    
    # Force numeric and fill NA
    for c in expr_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
        
    return out


# ---------------------------
# Junction matching strategies
# ---------------------------

def check_overlap(
    set_a: Set[Tuple[int, int]], 
    set_b: Set[Tuple[int, int]], 
    mode: str, 
    min_shared: int = 1, 
    min_jaccard: float = 0.5
) -> bool:
    
    # Quick short circuit for empty sets (single exon transcripts)
    if not set_a or not set_b:
        # If both are empty (single exon), do they overlap? 
        # Without coords of exons, we assume NO match unless exact strategy on single exons is handled elsewhere.
        # Here we only match based on junctions. Single exon transcripts have 0 junctions.
        return False

    intersection_size = len(set_a.intersection(set_b))

    if mode == "any_shared":
        return intersection_size >= min_shared

    if mode == "exact":
        return set_a == set_b

    if mode == "subset":
        return intersection_size >= min_shared and (set_a.issubset(set_b) or set_b.issubset(set_a))

    if mode == "jaccard":
        if intersection_size < min_shared:
            return False
        union_size = len(set_a.union(set_b))
        j_score = intersection_size / union_size if union_size > 0 else 0.0
        return j_score >= min_jaccard

    return False


# ---------------------------
# Clustering (connected comps)
# ---------------------------

def connected_components(nodes: List[str], edges: Dict[str, List[str]]) -> List[List[str]]:
    seen: Set[str] = set()
    comps: List[List[str]] = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        seen.add(n)
        comp = []
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
    tx2set: Dict[str, Set[Tuple[int, int]]],
    match_mode: str,
    min_shared: int,
    min_jaccard: float,
) -> List[List[str]]:
    
    # Pre-fetch sets to avoid dict lookups in the loop
    # Filter out transcripts that have NO junctions if min_shared > 0
    # (Single-exon transcripts cannot cluster by splice junctions)
    valid_tx = [t for t in tx_list if t in tx2set]
    
    # Optimization: Map index to transcript to use integer array indexing if needed, 
    # but strictly Python lists are fine for <10k items per group.
    
    edges: Dict[str, List[str]] = defaultdict(list)
    n = len(valid_tx)
    
    # O(N^2) comparison
    for i in range(n):
        u = valid_tx[i]
        set_u = tx2set[u]
        if not set_u: continue 

        for j in range(i + 1, n):
            v = valid_tx[j]
            set_v = tx2set[v]
            
            if check_overlap(set_u, set_v, match_mode, min_shared, min_jaccard):
                edges[u].append(v)
                edges[v].append(u)

    # Add back single-exon transcripts as isolated components?
    # Current logic: If they have no junctions, they won't match anything, 
    # so they will be returned as singleton clusters by connected_components.
    return connected_components(tx_list, edges)


def compute_score(expr_sub: pd.DataFrame, expr_cols: List[str], score_mode: str) -> pd.Series:
    m = expr_sub[expr_cols].to_numpy(dtype=float)
    if score_mode == "sum":
        vals = np.sum(m, axis=1)
    elif score_mode == "max":
        vals = np.max(m, axis=1)
    elif score_mode == "mean":
        vals = np.mean(m, axis=1)
    else:
        raise ValueError(f"Unknown score mode: {score_mode}")
    
    return pd.Series(vals, index=expr_sub["transcript_id"].values)


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Cluster isoforms by splice-junction matching within (gene,locus), "
            "then remove low-expressed isoforms per cluster so the top X%% are retained."
        )
    )
    ap.add_argument("--gtf", required=True, type=Path, help="Input GTF (optionally .gz)")
    ap.add_argument("--expr", required=True, type=Path, help="Expression table (CSV/TSV)")
    ap.add_argument("--tx-col", default="transcript_id", help="Transcript ID column in expression table.")
    ap.add_argument("--expr-col", default="expression",
                    help="Expression columns: single, comma-separated, or '*' for all non-ID numeric columns.")

    ap.add_argument("--match-mode", choices=["any_shared", "exact", "subset", "jaccard"],
                    default="any_shared",
                    help="How to connect isoforms into clusters (default: any_shared).")
    ap.add_argument("--min-shared", type=int, default=1,
                    help="Minimum number of shared junctions to consider a match (default: 1).")
    ap.add_argument("--min-jaccard", type=float, default=0.5,
                    help="Minimum Jaccard similarity for match-mode=jaccard (default: 0.5).")

    ap.add_argument("--retain-top-pct", type=float, default=98.0,
                    help="Retain the top X percent of isoforms by score within each cluster (default: 98).")
    ap.add_argument("--score", choices=["sum", "max", "mean"], default="sum",
                    help="Score metric across expression columns (default: sum).")
    ap.add_argument("--min-keep", type=int, default=1,
                    help="Always keep at least this many isoforms per cluster (default: 1).")

    ap.add_argument("--out", required=True, type=Path, help="Output filtered expression table (TSV).")
    ap.add_argument("--clusters-out", type=Path, default=None,
                    help="Optional TSV: transcript -> cluster_id stats.")
    ap.add_argument("--dropped-out", type=Path, default=None,
                    help="Optional TSV: dropped transcripts stats.")

    args = ap.parse_args()

    if not (0.0 < args.retain_top_pct <= 100.0):
        raise ValueError("--retain-top-pct must be in (0, 100].")

    print(f"[INFO {ts()}] Loading GTF + expression...")
    tx2chain, tx2set, tx2gene, tx2locus = parse_gtf(args.gtf)
    expr_df = load_expression_table(args.expr, args.tx_col, args.expr_col)

    expr_cols = [c for c in expr_df.columns if c != "transcript_id"]
    print(f"[INFO {ts()}] Parsed {len(tx2chain)} transcripts from GTF; {len(expr_cols)} numeric expression column(s)")

    # restrict to transcripts present in expr and in GTF parse
    # Use set intersection for speed
    common_ids = set(tx2chain.keys()).intersection(expr_df["transcript_id"])
    expr_df = expr_df[expr_df["transcript_id"].isin(common_ids)].copy()
    print(f"[INFO {ts()}] Using {len(expr_df)} transcripts present in both expression and GTF")

    # Map gene/locus efficiently
    # Note: map(dict) is faster than apply(lambda)
    expr_df["__gene__"] = expr_df["transcript_id"].map(tx2gene)
    expr_df["__locus__"] = expr_df["transcript_id"].map(tx2locus)
    
    # Drop rows where gene or locus mapping failed (though intersection above should handle this)
    expr_df = expr_df.dropna(subset=["__gene__", "__locus__"]).copy()

    groups = expr_df.groupby(["__gene__", "__locus__"], sort=False)
    print(f"[INFO {ts()}] Clustering within {len(groups)} (gene,locus) groups (mode={args.match_mode})...")

    keep_set: Set[str] = set()
    cluster_rows = []
    dropped_rows = []

    cluster_id_counter = 0

    for (gene, locus), sub in groups:
        tx_list = sub["transcript_id"].astype(str).tolist()
        
        # Trivial case: 1 transcript
        if len(tx_list) == 1:
            t = tx_list[0]
            keep_set.add(t)
            if args.clusters_out:
                cluster_rows.append((t, gene, locus, f"{gene}|{locus}|{cluster_id_counter}", 1))
            cluster_id_counter += 1
            continue

        clusters = build_clusters_for_group(
            tx_list=tx_list,
            tx2set=tx2set,
            match_mode=args.match_mode,
            min_shared=args.min_shared,
            min_jaccard=args.min_jaccard,
        )

        score_s = compute_score(sub, expr_cols, args.score)

        for comp in clusters:
            cluster_id = f"{gene}|{locus}|{cluster_id_counter}"
            cluster_id_counter += 1

            n = len(comp)
            
            # Extract scores for this cluster
            # Create a DataFrame to allow sorting by Score THEN ID (deterministic tie-breaking)
            comp_df = pd.DataFrame({"score": score_s.loc[comp].astype(float), "id": comp})
            comp_df.set_index("id", inplace=True)
            
            # Sort descending by score, then ascending by ID (for ties)
            comp_sorted = comp_df.sort_values(by=["score", "id"], ascending=[False, True])

            # Calculate cutoff
            n_keep_target = int(np.ceil((args.retain_top_pct / 100.0) * n))
            n_keep = max(int(args.min_keep), min(n, n_keep_target))

            keep_tx_ids = comp_sorted.head(n_keep).index.tolist()
            drop_tx_ids = comp_sorted.tail(n - n_keep).index.tolist()

            keep_set.update(keep_tx_ids)

            # Logging
            if args.clusters_out:
                for t in comp:
                    cluster_rows.append((t, gene, locus, cluster_id, n))

            if args.dropped_out and drop_tx_ids:
                cutoff_score = comp_sorted.iloc[n_keep - 1]["score"]
                for t in drop_tx_ids:
                    dropped_rows.append({
                        "transcript_id": t,
                        "gene": gene,
                        "locus": locus,
                        "cluster_id": cluster_id,
                        "cluster_size": n,
                        "score": comp_df.loc[t, "score"],
                        "cutoff_score": cutoff_score,
                        "reason": "below_cluster_cutoff"
                    })

    # Final filtering
    filtered = expr_df[expr_df["transcript_id"].isin(keep_set)].copy()
    filtered = filtered[["transcript_id"] + expr_cols].sort_values("transcript_id").reset_index(drop=True)
    filtered.rename(columns={"transcript_id": "#TranscriptID"}, inplace=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.out, sep="\t", index=False)
    print(f"[INFO {ts()}] Kept {len(filtered)} transcripts; wrote: {args.out}")

    if args.clusters_out:
        args.clusters_out.parent.mkdir(parents=True, exist_ok=True)
        cl_df = pd.DataFrame(cluster_rows, columns=["transcript_id", "gene", "locus", "cluster_id", "cluster_size"])
        cl_df.to_csv(args.clusters_out, sep="\t", index=False)
        print(f"[INFO {ts()}] Wrote cluster assignments: {args.clusters_out}")

    if args.dropped_out:
        args.dropped_out.parent.mkdir(parents=True, exist_ok=True)
        dr_df = pd.DataFrame(dropped_rows)
        if not dr_df.empty:
            dr_df.to_csv(args.dropped_out, sep="\t", index=False)
        else:
            # write empty file with header
            with open(args.dropped_out, "w") as fh:
                fh.write("transcript_id\tgene\tlocus\tcluster_id\tcluster_size\tscore\tcutoff_score\treason\n")
        print(f"[INFO {ts()}] Wrote dropped list: {args.dropped_out}")

    # Summary
    total = len(expr_df)
    kept = len(filtered)
    dropped = total - kept
    print(f"[INFO {ts()}] Summary: total={total}, kept={kept}, dropped={dropped}, "
          f"retain-top-pct={args.retain_top_pct}")


if __name__ == "__main__":
    main()
