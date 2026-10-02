#!/usr/bin/env python3
"""
scripts/liftover_hm3.py
-----------------------
Give every HM3 SNP in the PRS-CSx reference (snpinfo_mult_1kg_hm3, GRCh37)
its GRCh38 position, for two uses:

  * matching GRCh38 GWAS files by position (scripts/format_gwas_sumstats.py)
  * finding the HM3 SNPs in the GRCh38 All of Us genotype data
    (aou/00_extract_hm3.sh, aou/01_match_variants.py)

SNPs that fail to lift, land on another chromosome, or map to the minus strand
are dropped (a minus-strand hit would need complemented alleles; there are very
few among HM3 SNPs). Two SNPs landing on one GRCh38 position are both dropped.

Output (tab-separated, header):  SNP CHR BP37 BP38 A1 A2

Usage:
    python scripts/liftover_hm3.py --snpinfo data/ld_ref/snpinfo_mult_1kg_hm3 \\
        --chain data/ld_ref/hg19ToHg38.over.chain.gz --out data/ld_ref/hm3_hg38.tsv
"""

import argparse
from collections import Counter

from pyliftover import LiftOver


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--snpinfo", required=True)
    ap.add_argument("--chain", required=True, help="hg19ToHg38.over.chain(.gz)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lo = LiftOver(args.chain)
    rows, n_in = [], 0
    fail = Counter()
    with open(args.snpinfo) as fh:
        next(fh)
        for line in fh:
            chrom, snp, bp, a1, a2 = line.split()[:5]
            n_in += 1
            hits = lo.convert_coordinate("chr" + chrom, int(bp) - 1)   # 0-based
            if not hits:
                fail["unmapped"] += 1
                continue
            new_chrom, new_pos0, strand = hits[0][0], hits[0][1], hits[0][2]
            if new_chrom != "chr" + chrom:
                fail["other_chromosome"] += 1
                continue
            if strand != "+":
                fail["minus_strand"] += 1
                continue
            rows.append((snp, chrom, int(bp), new_pos0 + 1, a1, a2))

    pos_count = Counter((r[1], r[3]) for r in rows)
    kept = [r for r in rows if pos_count[(r[1], r[3])] == 1]
    fail["shared_grch38_position"] = len(rows) - len(kept)

    with open(args.out, "w") as fh:
        fh.write("SNP\tCHR\tBP37\tBP38\tA1\tA2\n")
        for r in kept:
            fh.write("%s\t%s\t%d\t%d\t%s\t%s\n" % r)

    print("HM3 SNPs read: %d; lifted to GRCh38: %d" % (n_in, len(kept)))
    for reason, n in sorted(fail.items()):
        print("  dropped (%s): %d" % (reason, n))


if __name__ == "__main__":
    main()
