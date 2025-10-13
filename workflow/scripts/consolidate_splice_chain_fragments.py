#!/usr/bin/env python3
import argparse
import gzip
import time
from datetime import datetime
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

try:
    import multiprocessing as mp
except Exception:
    mp = None


def ts() -> str:
    return datetime.now().strftime("%H:%M:%S")


@dataclass(frozen=True)
class Junction:
    donor: int
    acceptor: int

    def as_tuple(self) -> Tuple[int, int]:
        return (self.donor, self.acceptor)


def smart_open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_bed_regions(bed_path: Optional[Path]) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    regions: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    if bed_path is None:
        return regions
    with smart_open(bed_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            strand = parts[5] if len(parts) > 5 else "."
            if strand not in ["+", "-", "."]:
                strand = "."
            regions[(chrom, strand)].append((start, end))
    return regions


def point_within_window_of_any_interval(chrom: str, strand: str, pos: int,
                                        regions: Dict[Tuple[str, str], List[Tuple[int, int]]],
                                        window: int) -> bool:
    for key in [(chrom, strand), (chrom, ".")]:
        for (s, e) in regions.get(key, []):
            if s - window <= pos <= e + window:
                return True
    return False


def parse_gtf(gtf_path: Path,
              transcript_attr: str = "transcript_id",
              gene_attr: str = "gene_id") -> Tuple[
                  Dict[str, List[Tuple[int, int]]],
                  Dict[str, str],
                  Dict[str, str],
                  Dict[str, Tuple[str, str, int, int]],
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

    tx2chain: Dict[str, List[Tuple[int, int]]] = {}
    tx2locus: Dict[str, str] = {}
    tx2_tss_tes: Dict[str, Tuple[str, str, int, int]] = {}

    for tx, exons in exons_by_tx.items():
        chroms = {c for c, _, _, _ in exons}
        strands = {st for _, st, _, _ in exons}
        if len(chroms) != 1 or len(strands) != 1:
            continue
        chrom = next(iter(chroms))
        strand = next(iter(strands))

        exons_sorted = sorted(exons, key=lambda x: x[2])
        juncs: List[Tuple[int, int]] = []
        starts, ends = [], []
        for i in range(len(exons_sorted) - 1):
            _, _, s_i, e_i = exons_sorted[i]
            _, _, s_j, e_j = exons_sorted[i + 1]
            donor = e_i; acceptor = s_j
            if donor < acceptor:
                juncs.append((donor, acceptor))
        for _, _, s, e in exons_sorted:
            starts.append(s); ends.append(e)

        if strand == "+":
            tss = min(starts); tes = max(ends)
        else:
            tss = max(ends); tes = min(starts)

        tx2chain[tx] = juncs
        tx2locus[tx] = f"{chrom}|{strand}"
        tx2_tss_tes[tx] = (chrom, strand, tss, tes)

    return tx2chain, tx2gene, tx2locus, tx2_tss_tes


def is_contiguous_subsequence(subseq: List[Tuple[int, int]], seq: List[Tuple[int, int]]) -> bool:
    if not subseq:
        return True
    n, m = len(subseq), len(seq)
    for i in range(m - n + 1):
        if seq[i:i+n] == subseq:
            return True
    return False


# ---- Multiprocessing helpers ----
_G_tx2chain = None
_G_tx2_tss_tes = None
_G_terminal_only = None
_G_minimal_superset = None
_G_tss_regions = None
_G_tes_regions = None
_G_protect_window = None


def _mp_init(tx2chain, tx2_tss_tes, terminal_only, minimal_superset,
             tss_regions, tes_regions, protect_window):
    global _G_tx2chain, _G_tx2_tss_tes, _G_terminal_only, _G_minimal_superset
    global _G_tss_regions, _G_tes_regions, _G_protect_window
    _G_tx2chain = tx2chain
    _G_tx2_tss_tes = tx2_tss_tes
    _G_terminal_only = terminal_only
    _G_minimal_superset = minimal_superset
    _G_tss_regions = tss_regions
    _G_tes_regions = tes_regions
    _G_protect_window = protect_window


def _process_gene(args) -> Tuple[Dict[str, List[str]], int, int]:
    (gene, locus), tlist = args
    tx2chain = _G_tx2chain
    tx2_tss_tes = _G_tx2_tss_tes
    terminal_only = _G_terminal_only
    minimal_superset = _G_minimal_superset
    tss_regions = _G_tss_regions
    tes_regions = _G_tes_regions
    protect_window = _G_protect_window

    chain_sets = {tx: frozenset(tx2chain.get(tx, [])) for tx in tlist}
    chain_lists = {tx: tx2chain.get(tx, []) for tx in tlist}

    partial: Dict[str, List[str]] = {}
    contained = 0

    for a in tlist:
        Aset = chain_sets[a]
        Alist = chain_lists[a]
        chrom, strand, a_tss, a_tes = tx2_tss_tes.get(a, ("?", "?", 0, 0))

        cand: List[Tuple[str, int]] = []
        for b in tlist:
            if a == b:
                continue
            Bset = chain_sets[b]
            Blist = chain_lists[b]

            if not Aset:
                if not Bset:
                    continue
            else:
                if not (Aset.issubset(Bset) and len(Bset) > len(Aset)):
                    continue

            if terminal_only and not is_contiguous_subsequence(Alist, Blist):
                continue

            b_chrom, b_strand, b_tss, b_tes = tx2_tss_tes.get(b, (chrom, strand, 0, 0))

            a_near_tss = bool(tss_regions) and point_within_window_of_any_interval(chrom, strand, a_tss, tss_regions, protect_window)
            b_near_tss = bool(tss_regions) and point_within_window_of_any_interval(b_chrom, b_strand, b_tss, tss_regions, protect_window)
            a_near_tes = bool(tes_regions) and point_within_window_of_any_interval(chrom, strand, a_tes, tes_regions, protect_window)
            b_near_tes = bool(tes_regions) and point_within_window_of_any_interval(b_chrom, b_strand, b_tes, tes_regions, protect_window)

            if tss_regions and (a_near_tss or b_near_tss):
                if abs(a_tss - b_tss) > protect_window:
                    continue
            if tes_regions and (a_near_tes or b_near_tes):
                if abs(a_tes - b_tes) > protect_window:
                    continue

            extra = max(0, len(Blist) - len(Alist))
            cand.append((b, extra))

        if minimal_superset and cand:
            min_extra = min(extra for _, extra in cand)
            cand = [(tx, extra) for tx, extra in cand if extra == min_extra]

        supers = sorted(set([tx for tx, _ in cand]))
        partial[a] = supers
        if supers:
            contained += 1

    return partial, contained, len(tlist)


def find_supersets_with_filters_and_logging(
    tx2chain: Dict[str, List[Tuple[int, int]]],
    tx2gene: Dict[str, str],
    tx2locus: Dict[str, str],
    tx2_tss_tes: Dict[str, Tuple[str, str, int, int]],
    terminal_only: bool,
    minimal_superset: bool,
    tss_regions: Dict[Tuple[str, str], List[Tuple[int, int]]],
    tes_regions: Dict[Tuple[str, str], List[Tuple[int, int]]],
    protect_window: int,
    log_every: int = 1000,
    threads: int = 1,
) -> Tuple[Dict[str, List[str]], int, int]:
    gene_groups: Dict[Tuple[str, str], List[str]] = defaultdict(list)
    for tx, gene in tx2gene.items():
        locus = tx2locus.get(tx, "?|?")
        gene_groups[(gene, locus)].append(tx)

    total_genes = len(gene_groups)
    total_tx = len(tx2chain)
    print(f"[INFO {ts()}] Starting containment scan: {total_tx} transcripts across {total_genes} gene-loci (threads={threads})")

    tx_supersets: Dict[str, List[str]] = {tx: [] for tx in tx2chain.keys()}

    processed_genes = 0
    processed_tx_running = 0
    contained_tx_count = 0

    items = list(gene_groups.items())

    if threads > 1 and mp is not None:
        with mp.Pool(processes=threads, initializer=_mp_init,
                     initargs=(tx2chain, tx2_tss_tes, terminal_only, minimal_superset,
                               tss_regions, tes_regions, protect_window)) as pool:
            for partial, contained, n_tx in pool.imap_unordered(_process_gene, items, chunksize=32):
                processed_genes += 1
                processed_tx_running += n_tx
                contained_tx_count += contained
                tx_supersets.update(partial)

                if log_every > 0 and (processed_genes % log_every == 0 or processed_genes == total_genes):
                    pct = 100.0 * processed_genes / max(1, total_genes)
                    print(f"[INFO {ts()}] Processed {processed_genes} / {total_genes} genes ({pct:.1f}%) — {processed_tx_running} transcripts scanned; {contained_tx_count} contained so far")
    else:
        for (gene, locus), tlist in items:
            processed_genes += 1
            processed_tx_running += len(tlist)

            chain_sets = {tx: frozenset(tx2chain.get(tx, [])) for tx in tlist}
            chain_lists = {tx: tx2chain.get(tx, []) for tx in tlist}

            for a in tlist:
                Aset = chain_sets[a]
                Alist = chain_lists[a]
                chrom, strand, a_tss, a_tes = tx2_tss_tes.get(a, ("?", "?", 0, 0))

                cand: List[Tuple[str, int]] = []
                for b in tlist:
                    if a == b:
                        continue
                    Bset = chain_sets[b]
                    Blist = chain_lists[b]

                    if not Aset:
                        if not Bset:
                            continue
                    else:
                        if not (Aset.issubset(Bset) and len(Bset) > len(Aset)):
                            continue

                    if terminal_only and not is_contiguous_subsequence(Alist, Blist):
                        continue

                    b_chrom, b_strand, b_tss, b_tes = tx2_tss_tes.get(b, (chrom, strand, 0, 0))

                    a_near_tss = bool(tss_regions) and point_within_window_of_any_interval(chrom, strand, a_tss, tss_regions, protect_window)
                    b_near_tss = bool(tss_regions) and point_within_window_of_any_interval(b_chrom, b_strand, b_tss, tss_regions, protect_window)
                    a_near_tes = bool(tes_regions) and point_within_window_of_any_interval(chrom, strand, a_tes, tes_regions, protect_window)
                    b_near_tes = bool(tes_regions) and point_within_window_of_any_interval(b_chrom, b_strand, b_tes, tes_regions, protect_window)

                    if tss_regions and (a_near_tss or b_near_tss):
                        if abs(a_tss - b_tss) > protect_window:
                            continue
                    if tes_regions and (a_near_tes or b_near_tes):
                        if abs(a_tes - b_tes) > protect_window:
                            continue

                    extra = max(0, len(Blist) - len(Alist))
                    cand.append((b, extra))

                if minimal_superset and cand:
                    min_extra = min(extra for _, extra in cand)
                    cand = [(tx, extra) for tx, extra in cand if extra == min_extra]

                supers = sorted(set([tx for tx, _ in cand]))
                tx_supersets[a] = supers
                if supers:
                    contained_tx_count += 1

            if log_every > 0 and (processed_genes % log_every == 0 or processed_genes == total_genes):
                pct = 100.0 * processed_genes / max(1, total_genes)
                print(f"[INFO {ts()}] Processed {processed_genes} / {total_genes} genes ({pct:.1f}%) — {processed_tx_running} transcripts scanned; {contained_tx_count} contained so far")

    return tx_supersets, total_genes, contained_tx_count


def load_expression_table(path: Path,
                          tx_col: str,
                          expr_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    cols_lower = {c.lower(): c for c in df.columns}
    tx_c = cols_lower.get(tx_col.lower(), tx_col)
    if tx_c not in df.columns:
        raise ValueError(f"Could not find transcript id column '{tx_col}' in {list(df.columns)}")

    if expr_col == "*" or expr_col.upper() == "ALL":
        expr_cols = [c for c in df.columns if c != tx_c]
    elif "," in expr_col:
        expr_cols = [c.strip() for c in expr_col.split(",") if c.strip()]
        for c in expr_cols:
            if c not in df.columns:
                raise ValueError(f"Requested expr column '{c}' not in columns {list(df.columns)}")
    else:
        if expr_col not in df.columns:
            if expr_col.lower() in cols_lower:
                expr_col = cols_lower[expr_col.lower()]
            else:
                raise ValueError(f"Could not find expression column '{expr_col}' in {list(df.columns)}")
        expr_cols = [expr_col]

    out = df[[tx_c] + expr_cols].rename(columns={tx_c: "transcript_id"}).copy()
    for c in expr_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0.0)
    return out


def _cascade_per_column(col_values: Dict[str, float],
                        tx2supersets: Dict[str, List[str]],
                        tx2chain: Dict[str, List[Tuple[int, int]]],
                        proportional: bool) -> Dict[str, float]:
    """
    Propagate mass upward from contained transcripts to their supersets in a DAG,
    processing from smallest junction-chain to largest.
    """
    mass = dict(col_values)
    nodes = [t for t, supers in tx2supersets.items() if supers]

    def size(t): return len(tx2chain.get(t, []))
    nodes_sorted = sorted(nodes, key=size)

    for t in nodes_sorted:
        amt = mass.get(t, 0.0)
        supers = tx2supersets.get(t, [])
        if not supers or amt <= 0:
            continue

        if proportional:
            weights = [mass.get(s, 0.0) for s in supers]
            total_w = sum(w for w in weights if w > 0)
            if total_w > 0:
                shares = [w / total_w for w in weights]
            else:
                shares = [1.0 / len(supers) for _ in supers]
        else:
            shares = [1.0 / len(supers) for _ in supers]

        for s, sh in zip(supers, shares):
            mass[s] = mass.get(s, 0.0) + amt * sh
        mass[t] = 0.0

    return mass


def redistribute_expression(
    expr_df: pd.DataFrame,
    tx2supersets: Dict[str, List[str]],
    tx2chain: Dict[str, List[Tuple[int, int]]],
    drop_contained: bool = True,
    proportional: bool = True,
) -> pd.DataFrame:
    """
    Cascade redistribution is the default (and only) behavior:
    for each expression column, propagate mass upward from contained transcripts
    to their minimal supersets in ascending junction-chain order.
    """
    all_tx = set(expr_df["transcript_id"])
    expr_cols = [c for c in expr_df.columns if c != "transcript_id"]
    base = expr_df.set_index("transcript_id")

    adjusted_frames = []
    for col in expr_cols:
        col_base = base[col].to_dict()
        mass = _cascade_per_column(col_base, tx2supersets, tx2chain, proportional)
        series = pd.Series(mass, name=col)
        adjusted_frames.append(series)

    out = pd.concat(adjusted_frames, axis=1).reset_index().rename(columns={"index": "transcript_id"})

    if drop_contained:
        contained_ids = {t for t, supers in tx2supersets.items() if supers and t in all_tx}
        out = out[~out["transcript_id"].isin(contained_ids)].reset_index(drop=True)

    out = out.sort_values("transcript_id").reset_index(drop=True)
    return out


def round_preserve_sum(values: np.ndarray, target_sum: Optional[float] = None, method: str = "nearest") -> np.ndarray:
    if target_sum is None:
        target_sum = float(np.nansum(values))
    x = np.asarray(values, dtype=float)
    x = np.where(np.isnan(x), 0.0, x)
    x = np.where(x < 0, 0.0, x)

    if method == "floor":
        base_int = np.floor(x)
    elif method == "ceil":
        base_int = np.ceil(x)
    else:
        base_int = np.floor(x)

    diff = int(round(target_sum - base_int.sum()))
    frac = x - base_int
    order = np.argsort(-frac)
    out = base_int.astype(int)

    if diff >= 0:
        out[order[:diff]] += 1
    else:
        order_rev = np.argsort(frac)
        out[order_rev[:abs(diff)]] = np.maximum(0, out[order_rev[:abs(diff)]] - 1)

    return out


def apply_count_rounding(adjusted_df: pd.DataFrame, method: str) -> pd.DataFrame:
    out = adjusted_df.copy()
    expr_cols = [c for c in out.columns if c != "transcript_id"]
    for c in expr_cols:
        original_total = float(out[c].sum())
        out[c] = round_preserve_sum(out[c].to_numpy(), target_sum=original_total, method=method)
    return out


def write_containment_map(path: Path,
                          tx2supersets: Dict[str, List[str]]):
    with open(path, "w") as fh:
        fh.write("contained_transcript\tsuperset_transcripts\n")
        for t, supers in sorted(tx2supersets.items()):
            if supers:
                fh.write(f"{t}\t{','.join(supers)}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Merge contained transcripts by redistributing expression to isoforms that contain their splice chain (cascade redistribution), with long-read safeguards, logging, multiprocessing, multi-sample support, and counts rounding."
    )
    ap.add_argument("--gtf", required=True, type=Path, help="Input GTF (optionally .gz) with exon features.")
    ap.add_argument("--expr", required=True, type=Path, help="Expression table (CSV/TSV). Must contain transcript IDs.")
    ap.add_argument("--tx-col", default="transcript_id", help="Column name in expression table for transcript IDs (default: transcript_id).")
    ap.add_argument("--expr-col", default="expression", help="Expression columns: single name, comma-separated list, or '*' for all non-ID columns (default: 'expression').")
    ap.add_argument("--drop-contained", action="store_true", help="Drop contained transcripts from output (default: keep zeroed unless flag set).")
    ap.add_argument("--proportional", action="store_true", help="Redistribute proportional to supersets' expression (default: equal split unless flag set).")

    # Long-read specific controls
    ap.add_argument("--terminal-only", action="store_true", help="Require that contained vs. parent differ only at transcript ends (no internal novel junctions).")
    ap.add_argument("--minimal-superset", action="store_true", help="Route only to the closest (fewest extra junctions) supersets.")
    ap.add_argument("--protect-tss-bed", type=Path, default=None, help="BED file with known CAGE/TSS regions; enforces TSS proximity within --protect-window.")
    ap.add_argument("--protect-tes-bed", type=Path, default=None, help="BED file with known PAS/TES regions; enforces TES proximity within --protect-window.")
    ap.add_argument("--protect-window", type=int, default=50, help="Window (bp) to consider TSS/TES proximity (default: 50bp).")

    # Logging & parallelism
    ap.add_argument("--log-every", type=int, default=1000, help="Print progress every N genes (default: 1000).")
    ap.add_argument("--threads", type=int, default=1, help="Number of worker processes to use for gene-level scanning (default: 1).")

    # Raw counts
    ap.add_argument("--counts", action="store_true", help="Input values are raw counts; enable integer rounding in output.")
    ap.add_argument("--round-counts", choices=["none", "nearest", "floor", "ceil"], default=None,
                    help="Rounding mode for counts (default: 'nearest' if --counts, else 'none'). Column sums preserved.")

    ap.add_argument("--out", required=True, type=Path, help="Output CSV with adjusted expression.")
    ap.add_argument("--map-out", type=Path, default=None, help="Optional TSV mapping each contained transcript to its supersets.")

    args = ap.parse_args()

    start_time = time.time()

    print(f"[INFO {ts()}] Loading GTF and expression...")
    tx2chain, tx2gene, tx2locus, tx2_tss_tes = parse_gtf(args.gtf)
    expr_df = load_expression_table(args.expr, args.tx_col, args.expr_col)
    n_cols = len([c for c in expr_df.columns if c != "transcript_id"])
    print(f"[INFO {ts()}] Parsed {len(tx2chain)} transcripts from GTF; expression table has {n_cols} expression column(s)")

    # Load protected regions
    tss_regions = parse_bed_regions(args.protect_tss_bed) if args.protect_tss_bed else {}
    tes_regions = parse_bed_regions(args.protect_tes_bed) if args.protect_tes_bed else {}

    if args.protect_tss_bed:
        print(f"[INFO {ts()}] Loaded TSS protection regions from {args.protect_tss_bed}")
    if args.protect_tes_bed:
        print(f"[INFO {ts()}] Loaded TES protection regions from {args.protect_tes_bed}")

    tx2supersets, total_genes, contained_tx_count = find_supersets_with_filters_and_logging(
        tx2chain=tx2chain,
        tx2gene=tx2gene,
        tx2locus=tx2locus,
        tx2_tss_tes=tx2_tss_tes,
        terminal_only=args.terminal_only,
        minimal_superset=args.minimal_superset,
        tss_regions=tss_regions,
        tes_regions=tes_regions,
        protect_window=args.protect_window,
        log_every=args.log_every,
        threads=max(1, args.threads) if mp is not None else 1,
    )

    print(f"[INFO {ts()}] Redistributing expression across {n_cols} column(s) with cascade...")
    adjusted_df = redistribute_expression(
        expr_df,
        tx2supersets,
        tx2chain,
        drop_contained=args.drop_contained,
        proportional=args.proportional,
    )

    round_mode = args.round_counts
    if round_mode is None:
        round_mode = "nearest" if args.counts else "none"
    if round_mode != "none":
        print(f"[INFO {ts()}] Rounding counts with mode='{round_mode}' while preserving per-column totals")
        adjusted_df = apply_count_rounding(adjusted_df, method=round_mode)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adjusted_df.to_csv(args.out, index=False, sep="\t")

    if args.map_out:
        args.map_out.parent.mkdir(parents=True, exist_ok=True)
        write_containment_map(args.map_out, tx2supersets)

    elapsed = time.time() - start_time
    mins, secs = divmod(int(elapsed), 60)
    print(f"[INFO {ts()}] Completed in {mins}m {secs}s — {len(tx2chain)} transcripts, {total_genes} genes; {contained_tx_count} contained transcripts detected")
    print(f"[INFO {ts()}] Wrote adjusted expression to: {args.out}")
    if args.map_out:
        print(f"[INFO {ts()}] Wrote containment map to: {args.map_out}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        print("""\
Usage examples:

# Multi-sample counts (cascade redistribution is always on):
  python redistribute_splice_chain_expression.py \
    --gtf transcripts.gtf \
    --expr counts.tsv \
    --tx-col '#TranscriptID' \
    --expr-col '*' \
    --terminal-only --minimal-superset \
    --proportional --drop-contained \
    --counts --round-counts nearest \
    --threads 8 --log-every 1000 \
    --out adjusted_counts.csv \
    --map-out containment_map.tsv

# Abundances (e.g., TPM) across selected columns:
  python redistribute_splice_chain_expression.py \
    --gtf transcripts.gtf \
    --expr abundances.tsv \
    --tx-col transcript_id \
    --expr-col TPM,FPKM \
    --protect-tss-bed known_TSS_clusters.bed --protect-window 75 \
    --proportional \
    --out adjusted_expression.csv
""")
    else:
        main()
