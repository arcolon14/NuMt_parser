#!/usr/bin/env python3

# Working script for Numt Finder (pending name)
# Identify NUMTs in a set of raw reads
# Identify NUMTs in the genome
# (c) 2024 Angel G. Rivera-Colon & Alida de Flamingh

import gzip, argparse, sys
from os import path
from datetime import datetime

DATE = datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = '''Use a set of nuclear and mitochondrial reference sequences
to identify NUMTs in a set of sequencing reads.'''

def parse_args() -> argparse.ArgumentParser:
    # 
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-m', '--mt-fasta', required=True,
                   help='(str) Path to mitochondrial reference FASTA')
    p.add_argument('-n', '--nuc-fasta', required=True,
                   help='(str) Path to the nuclear genome reference FASTA')
    p.add_argument('-r', '--reads-fastq', required=True,
                   help='(str) Path to sequencing reads FASTQ')
    # TODO: what should default K be?
    k_len = 17
    p.add_argument('-k', '--kmer-len', required=False, type=int, default=k_len,
                   help=f'(int) Length of k for generating query and reference k-mers  [default={k_len}]')
    # TODO: what should min kmers be?
    min_k = 5
    p.add_argument('-i', '--min-kmers', required=False, type=str, default=min_k,
                   help=f'(int) Minimum of k-mers required to match two sequences  [default={min_k}]')
    p.add_argument('-o', '--outdir', required=False, type=str, default='.',
                   help='(str) Path to output directory')
    args = p.parse_args()
    # Check inputs
    if not path.exists(args.mt_fasta):
        sys.exit(f'Error: Mt FASTA ({args.mt_fasta}) does not exist.')
    if not path.exists(args.nuc_fasta):
        sys.exit(f'Error: Nuclear FASTA ({args.nuc_fasta}) does not exist.')
    if not path.exists(args.outdir):
        sys.exit(f'Error: Output directory ({args.outdir}) does not exist.')
    if not args.kmer_len > 0:
        sys.exit(f'Error: length of k ({args.kmer_len}) must be > 0.')
    if not args.min_kmers >= 1:
        sys.exit(f'Error: length of k ({args.min_kmers}) must be >= 1.')
    return args


def now() -> str:
    '''Print date and time'''
    return f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def main():
    print(f'{PROG} started on {now()}\n')
    # Parse Arguments
    args = parse_args()
    parse_args()
    # Finish
    print(f'\n{PROG} finished on {now()}')


if __name__ == '__main__':
    main()
