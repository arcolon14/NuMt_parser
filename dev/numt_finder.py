#!/usr/bin/env python3

# Working script for Numt Finder (pending name)
# Identify NUMTs in a set of raw reads
# Identify NUMTs in the genome
# (c) 2024 Angel G. Rivera-Colon & Alida de Flamingh

import gzip, argparse, sys
from os import path
from datetime import datetime

# Globals
DATE = datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = '''Use a set of nuclear and mitochondrial reference sequences
to identify NUMTs in a set of sequencing reads.'''
K_LEN = 17 # Default length of kmers.  TODO: what should this be?
MIN_K = 5  # Default min number of kmers to have a match.  TODO: what should this be?
MAX_N = 3  # Default max number of uncalled bases (Ns) seen before k-mers is discarded. TODO: what should this be?

def parse_args() -> argparse.ArgumentParser:
    # 
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-m', '--mt-fasta', required=True,
                   help='(str) Path to mitochondrial reference FASTA')
    p.add_argument('-n', '--nuc-fasta', required=True,
                   help='(str) Path to the nuclear genome reference FASTA')
    p.add_argument('-r', '--reads-fastq', required=True,
                   help='(str) Path to sequencing reads FASTQ')
    p.add_argument('-o', '--outdir', required=False, type=str, default='.',
                   help='(str) Path to output directory')
    p.add_argument('-k', '--kmer-len', required=False, type=int, default=K_LEN,
                   help=f'(int) Length of k for generating query and reference k-mers  [default={K_LEN}]')
    p.add_argument('-i', '--min-kmers', required=False, type=str, default=MIN_K,
                   help=f'(int) Minimum of k-mers required to match two sequences  [default={MIN_K}]')
    p.add_argument('-a', '--max-uncalled', required=False, type=int, default=MAX_N,
                   help=f'(int) Maximum number of uncalled bases (Ns) seen before a k-mer is discarded.  [default={MAX_N}]')

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
    if not args.max_uncalled > 0:
        sys.exit(f'Error: max number of uncalled bases ({args.max_uncalled}) in a k-mer must be > 0.')
    if args.max_uncalled >= args.kmer_len:
        sys.exit(f'Error: max uncalled bases ({args.max_uncalled}) must be smaller than length of k ({args.kmer_len})')
    return args

#
# Classes
#

class KmerCoordinate:
    '''Stores the coordinate of a k-mer in a sequence'''
    def __init__(self, seq_id: str, seq_bp: int):
        assert seq_bp >= 0
        self.id = seq_id
        self.bp = seq_bp
    def __str__(self):
        return f'{self.id} {self.bp}'
    def extract_kmer_seq(self, sequence: str, k_len: int=K_LEN) -> str:
        '''Extract this k-mer from the source sequence'''
        assert len(sequence) > 0
        kmer = sequence[self.bp:(self.bp+k_len)]
        return kmer

#
# Functions
#

def now() -> str:
    '''Print date and time'''
    return f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def kmerize_sequence(sequence: str, seq_id: str, kmers_dict: dict, k_len: int=K_LEN, 
                     max_n: int=MAX_N) -> list:
    '''
    Generate kmers from a given input sequence

    Args:
        sequence: (str) nucleotide sequence
        seq_id: (str) ID of the target sequence
        kmers_dict : (dict) Dictionary of kmers to be populated. Key is k-mer 
            sequence and value are a list of k-mer coordinates. 
        k_len: (int) K, length of k-mers
        max_n: (int) Max number of uncalled bases before a k-mer id discarded

    Returns:
        kmers_dict : (dict) Populated dictionary of kmers. Key is k-mer sequence 
            and value are a list of k-mer coordinates.
    '''
    assert len(sequence) > 0
    assert k_len > 0
    assert len(sequence) >= k_len
    # For log
    unique_kmers = set()
    kept_kmers = 0
    # Determine the total number of k-mers in the sequence
    n_kmers = len(sequence) - k_len + 1
    # Extract the kmers
    for pos in range(n_kmers):
        kmer = sequence[pos:(pos+k_len)]
        # Filter k-mers with too many Ns
        n = kmer.count('N')
        if n >=max_n:
            continue
        # Process the kept kmers
        kept_kmers += 1
        kmer_coord = KmerCoordinate(seq_id, pos)
        # Inititalize output dictionary
        kmers_dict.setdefault(kmer, [])
        kmers_dict[kmer].append(kmer_coord)
        # For record keeping
        unique_kmers.add(kmer)
    # Report to log
    n_unique = len(unique_kmers)
    print(f'''    Processed sequence {seq_id}, of total length {len(sequence):,} bp
        Expected {n_kmers:,} {k_len}-mers
        Retained {kept_kmers:,} {k_len}-mers, {n_unique:,} ({(n_unique/kept_kmers):0.2%}) are unique.''', flush=True)
    return kmers_dict

def kmerize_input_fasta(fasta_f: str, k_len: int=K_LEN, max_n: int=MAX_N) -> dict:
    '''
    Process an input FASTA and kmerize all available sequences

    Args:
        fasta_f: (str) Path to a sequence in FASTA format
        k_len: (int) K, length of k-mers
        max_n: (int) Max number of uncalled bases before a k-mer id discarded

    Returns:
        kmers_dict : (dict) dictionary of kmers. Key is k-mer sequence and value 
            are a list of k-mer coordinates.
    '''
    print(f'Processing sequence {fasta_f}...', flush=True)
    # Check inputs
    assert path.exists(fasta_f), f'Error: {fasta_f} not found'
    # Main output
    kmers_dict = dict()
    # Checks sequences and length
    total_seqs = 0
    total_lens = 0
    # Parse FASTA
    with gzip.open(fasta_f, 'rt') if fasta_f.endswith('.gz') else open(fasta_f) as f:
        name = None
        seq = ''
        for line in f:
            line = line.strip('\n')
            # Ignore empty lines and comments
            if len(line) == 0 or line[0] in ['#', '.']:
                continue
            # Check if its ID or sequence
            if line.startswith('>'):
                if name is not None:
                    # Kmerize sequence
                    kmers_dict = kmerize_sequence(seq, name, kmers_dict, k_len)
                    total_len += len(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                # If it is sequence
                seq += line.upper()
        # Process the last sequence
        kmers_dict = kmerize_sequence(seq, name, kmers_dict, k_len)
        total_len += len(seq)
    # Count the total number of kmers seen
    n_kmers = sum([ len(kmers_dict[k]) for k in kmers_dict ])
    n_unique = len(kmers_dict)
    # Report to log
    print(f'''\n    Total report for {fasta_f}:
        Sequence length:     {total_len:,} bp
        Total k-mers kept:   {n_kmers:,}
        Total unique k-mers: {n_unique:,} ({(n_unique/n_kmers):0.2%})''')
    return kmers_dict


def main():
    print(f'{PROG} started on {now()}\n')
    # Parse Arguments
    args = parse_args()

    # Finish
    print(f'\n{PROG} finished on {now()}')


if __name__ == '__main__':
    main()
