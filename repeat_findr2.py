#!/usr/bin/env python3
import sys
import subprocess
import os
from Bio import SeqIO
import argparse

"""
repeat_findr.py: Identify simple sequence repeats in a single-gene FASTA using Tandem Repeats Finder.
Usage: repeat_findr.py [-h] [--min-period MIN_PERIOD] [--max-period MAX_PERIOD]
                       [--min-copies MIN_COPIES] gene_fasta

By default, searches for repeats of period 1–7 bp with at least 3 copies.
"""

def run_trf(fasta, trf_params):
    cmd = ["trf", fasta] + trf_params + ["-d", "-h"]
    result = subprocess.run(cmd)
    if result.returncode not in (0, 1):
        result.check_returncode()
    base = os.path.basename(fasta)
    dat = f"{base}.{'-'.join(trf_params)}.dat"
    # fallback: remove d and h flags
    if not os.path.exists(dat):
        # original naming convention for compatibility
        dat = f"{base}.2.7.7.80.10.50.500.dat"
    if os.path.exists(dat):
        return dat
    raise FileNotFoundError(f"TRF .dat not found: {dat}")


def parse_trf(dat, min_period, max_period, min_copies):
    with open(dat) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            cols = line.split()
            if len(cols) < 14:
                continue
            seq_id = cols[0]
            start  = int(cols[1])
            end    = int(cols[2])
            period = int(cols[3])
            copies = float(cols[4])
            motif  = cols[13]
            # filter by user thresholds
            if period < min_period or period > max_period:
                continue
            if copies < min_copies:
                continue
            yield seq_id, start, end, period, copies, motif


def main():
    parser = argparse.ArgumentParser(description="Find simple sequence repeats via TRF.")
    parser.add_argument('fasta', help='Input single-sequence FASTA file')
    parser.add_argument('--min-period', type=int, default=1,
                        help='Minimum repeat unit length (bp)')
    parser.add_argument('--max-period', type=int, default=7,
                        help='Maximum repeat unit length (bp)')
    parser.add_argument('--min-copies', type=float, default=3.0,
                        help='Minimum number of repeat copies')
    parser.add_argument('--trf-params', nargs=7, metavar=('MATCH','MISMATCH','DELTA','PM','PI','MIN_SCORE','MAX_PERIOD'),
                        default=["2","7","7","80","10","50","500"],
                        help='Seven TRF parameters: match, mismatch, delta, PM, PI, min_score, max_period')
    args = parser.parse_args()

    # read and validate FASTA
    recs = list(SeqIO.parse(args.fasta, 'fasta'))
    if len(recs) != 1:
        sys.exit('Error: supply a single-sequence FASTA.')

    # run TRF
    dat = run_trf(args.fasta, args.trf_params)

    # parse and filter
    repeats = list(parse_trf(dat, args.min_period, args.max_period, args.min_copies))

    # output TSV
    print('seq_id\tstart\tend\tperiod\tcopies\tmotif')
    for rpt in repeats:
        print('\t'.join(map(str, rpt)))

if __name__ == '__main__':
    main()
# Usage examples:
# 1. Basic (default thresholds 1–7 bp, ≥3 copies):
#    python3 repeat_findr.py geneseqs/PA1003_pqsR_CDS.fasta \
#      > pqsR_repeats.tsv
#
# 2. Relaxed thresholds (2–10 bp, ≥2 copies):
#    python3 repeat_findr.py geneseqs/PA4547_pilR_CDS.fasta \
#      --min-period 2 --max-period 10 --min-copies 2 \
#      > pilR_relaxed_repeats.tsv
#
# 3. Custom TRF parameters:
#    python3 repeat_findr.py geneseqs/PA1003_pqsR_CDS.fasta \
#      --trf-params 2 5 5 80 10 20 1000 \
#      > pqsR_custom_trf.tsv
#
# For help: python3 repeat_findr.py -h
