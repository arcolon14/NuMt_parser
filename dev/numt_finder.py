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
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-m', '--mt-fasta', required=True,
                   help='(str) Path to mitochondrial reference FASTA')
    p.add_argument('-n', '--nuc-fasta', required=True,
                   help='(str) Path to the nuclear genome reference FASTA')
    p.add_argument('-r', '--reads-fastq', required=True,
                   help='(str) Path to sequencing reads FASTQ')
    # TODO: what should default K be?
    p.add_argument('-k', '--kmer-len', required=False, type=int, default=17,
                   help='(int) Length of k for generating query and reference k-mers  [default=17]')
    # TODO: what should min kmers be?
    p.add_argument('-i', '--min-kmers', required=False, type=str, default=5,
                   help='(int) Minimum of k-mers required to match two sequences  [default=5]')
    p.add_argument('-o', '--outdir', required=False, type=str, default='.',
                   help='(str) Path to output directory')
    args = p.parse_args()

    if not path.exists(args.mt_fasta):
        sys.exit(f'Error: `{args.mt_fasta}`: Mt FASTA does not exist.')
    if not path.exists(args.nuc_fasta):
        sys.exit(f'Error: `{args.nuc_fasta}`: Nuclear FASTA does not exist.')
    if not path.exists(args.outdir):
        sys.exit(f'Error: `{args.outdir}`: Output directory does not exist.')
    # TODO: Checks for kmers

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
