#!/usr/bin/env python3
import gzip, argparse, sys, datetime
from os import path

#
# Script to compare read similarity to Mitochondrial and numt sequences
# (c) 2020 Angel G. Rivera-Colon & Alida de Flamingh
#

DATE = datetime.datetime.now().strftime("%Y%m%d")
PROG = sys.argv[0].split('/')[-1]
DESC = 'Parse a set of mitochondrial reads and compare similarity between mitochondrial and NUMT references.'
ALL_CIGAR_OPS = list('MIDNSHP=X')
VALID_CIGAR_OPS = list('MIDSH')

# -----------
# Input files
# -----------

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('--mt-fasta',   required=True, help='Mitochondrial Sequence FASTA')
    p.add_argument('--numt-fasta', required=True, help='numt Sequence FASTA')
    p.add_argument('--mt-sam',     required=True, help='Read alignments to the mt reference in SAM format')
    p.add_argument('--numt-sam',   required=True, help='Read alignments to the numt reference in SAM format')
    p.add_argument('--outfile',    required=True, help='Path and name to the output TSV file. Example: ./<sample_id>.tsv')
    args = p.parse_args()

    if not path.exists(args.numt_fasta):
        sys.exit(f'Error: `{args.numt_fasta}`: Numt FASTA does not exist.')
    if not path.exists(args.mt_fasta):
        sys.exit(f'Error: `{args.mt_fasta}`: Mt FASTA does not exist.')
    if not path.exists(args.numt_sam):
        sys.exit(f'Error: `{args.numt_sam}`: Numt SAM does not exist.')
    if not path.exists(args.mt_sam):
        sys.exit(f'Error: `{args.mt_sam}`: Mt SAM does not exist.')

    return args

# TODO: Make the script work for several samples at a time
# TODO: Make compatible with BAM files. Check library `pysam`

# -------
# Classes
# -------

#
# Create a RefSequence class
class RefSequence:
    def __init__(self, sequence_id, sequence):
        self.id  = sequence_id
        self.seq = sequence
        self.len = len(self.seq)
    def __str__(self):
        return f'{self.id} : {self.len:,} bp'

#
# Create a ReadAlignment class
class ReadAlignment:
    def __init__(self, read_id, sam_flag, ctg_name, pos_bp, cigar, sequence):
        assert type(sam_flag) is int
        assert type(pos_bp) is int
        self.rid = read_id
        self.ctg = ctg_name
        self.pos = pos_bp
        self.seq = sequence
        # Process SAM flags
        self.sam = sam_flag
        # Determine if reverse complimented
        self.dir = 'forward'
        if (self.sam & 16) == 16:
            self.dir = 'reverse'
        # Determine if read is mapped
        self.map = 'mapped'
        if (self.sam & 4) == 4:
            self.map = 'unmapped'
        # Determine the ALignment TYpe (aka, if primary, secondary, supplementary)
        self.alty = 'primary'
        if (self.sam & 4) == 4:
            self.alty = 'unmapped'
        elif (self.sam & 256) == 256:
            self.alty = 'secondary'
        elif (self.sam & 2048) == 2048:
            self.alty = 'supplementary'
        # Process Paired information
        self.paired = False
        self.read   = 0
        if (self.sam & 1) == 1:
            self.paired = True
            # If the read IS paired, check of read 1 or read 2
            if (self.sam & 64) == 64:
                self.read = 1
            elif (self.sam & 128) == 128:
                self.read = 2
        # Get a read-pair specific read ID
        self.rpid = f'{self.rid}/{self.read}'
        # Process CIGARs
        self.cig = CIGAR(cigar)
    def __str__(self):
        return f'''Read ID:   {self.rid}
SAM flag:  {self.sam}
Ctg Name:  {self.ctg}
Align BP:  {self.pos:,}
Paired:    {self.paired} (read {self.read})
CIGAR:     {self.cig}
Seq:       {self.seq}
Direction: {self.dir}
Status:    {self.map}
Aln Type:  {self.alty}
'''

#
# CIGAR string class
class CIGAR:
    # All "legal" cigar operations: MIDNSHP=X
    def __init__(self, cigar_str):
        self.str = cigar_str
        self.rev_str = ''.join(split_cigar_str(self.str)[::-1])
    def __str__(self):
        return f'{self.str} {self.rev_str}'
    # CIGAR class function to generate a list of per-nucleotide CIGAR elements
    # i.e.: is nucleotide 'n' a 'M'?
    def per_nt_cigar(self, reverse=False):
        cigar = self.str
        if reverse is True:
            cigar = self.rev_str
        nt_cigar_list = list()
        # Loop over the split CIGAR elements
        for cig in split_cigar_str(cigar):
            # Break down CIGAR elements into type and count
            count = int(cig[:-1])
            ctype = cig[-1]
            assert ctype in VALID_CIGAR_OPS
            nt_cigar_list += [ctype] * count
        return nt_cigar_list
    # Function to split cigar into list; same as `split_cigar_str`, but works automatically with the class.
    def split_cigar(self, reverse=False):
        cigar = self.str
        if reverse is True:
            cigar = self.rev_str
        splitcig = []
        prev = 0
        for i, c in enumerate(cigar):
            if c in ALL_CIGAR_OPS:
                splitcig.append(cigar[prev:i+1])
                prev=i+1
                if c not in VALID_CIGAR_OPS:
                    sys.exit(f"Error: \'{c}\' in CIGAR {cigar} is an incompatible opertation. Compatible CIGAR operations are: {','.join(VALID_CIGAR_OPS)}. See README for more info.")
        return splitcig

#
# Create a class to store Alignment Comparisons
class AlignmentComparison:
    def __init__(self, read_id, aln_status, alignment_size, num_mismatches):
        assert aln_status in [True, False]
        self.rid = read_id
        self.status = aln_status
        self.aln_len = None
        self.n_mismatch = None
        self.per_identity = None
        if self.status is True:
            self.aln_len = alignment_size
            self.n_mismatch = num_mismatches
            self.per_identity = 1 - (self.n_mismatch / self.aln_len)
    def __str__(self):
        return f'Read ID: {self.rid}  Alignment Size: {self.aln_len}bp  Number Mismatches: {self.n_mismatch}  % Identity {self.per_identity}'

# ---------
# Functions
# ---------

#
# Print the current date and time.
def now():
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

#
# function to reverse complete sequence. It finds complement and inverts the order
def rev_comp(sequence):
    rev = []
    for nt in sequence.upper():
        if nt == 'A':
            rev.append('T')
        elif nt == 'C':
            rev.append('G')
        elif nt == 'G':
            rev.append('C')
        elif nt == 'T':
            rev.append('A')
        elif nt not in ['A', 'C', 'G', 'T']:
            rev.append('N')
    return ''.join(rev[::-1])

#
# Function to read fasta file and return a dictionary with all sequence objects available
# { RefSequence_id : RefSequence_object }
def read_fasta(fasta_f):
    # Check input file
    assert path.exists(fasta_f) is True
    fasta_dict = dict()
    total_len = 0
    with gzip.open(fasta_f, 'rt') if fasta_f.endswith('.gz') else open(fasta_f) as f:
        name = None
        seq = []
        for line in f:
            line = line.strip('\n')
            # Ignore empty lines
            if len(line) == 0:
                continue
            # Check if its ID or sequence
            if line[0] == '>':
                if name is not None:
                    fasta_dict[name] = RefSequence(name, ''.join(seq))
                    total_len += fasta_dict[name].len
                name = line[1:].split()[0]
                seq = []
            elif line[0] in ['#', '.']:
                # Ignore comment lines
                continue
            else:
                # if it is sequence
                seq.append(line.upper())
    fasta_dict[name] = RefSequence(name, ''.join(seq))
    total_len += fasta_dict[name].len
    # Print to log:
    print(f'    Read file: {fasta_f}\n    Composed of {len(fasta_dict):,} sequence(s) with total length of {total_len:,} bp.')
    # Return genome dictionary
    return fasta_dict

#
# Function to read both mt and numt fasta files and merge into a single final dictionary
# { sequence_1 : RefSequence, sequence_2 : RefSequence, ... }
def extract_ref_sequence_dictionary(mt_fasta, nuMt_fasta):
    # Read mt fasta
    print('Reading mitochondrial reference...')
    mt_seq_dict = read_fasta(mt_fasta)
    # Read numt fasta
    print('\nReading NUMT reference...')
    numt_seq_dict = read_fasta(nuMt_fasta)
    # Merge both into a single reference dictionary
    ref_sequence_dictionary = {**mt_seq_dict, **numt_seq_dict}
    # Return the final dictonary
    return ref_sequence_dictionary

#
# Function to parse SAM file
# Output is a list containing ReadAlignment objects
# [ ReadAlignment, ReadAlignment, ... ]
def read_sam_file(sam_f):
    # Check input file
    assert path.exists(sam_f) is True
    align_list = list()
    total_aln = 0
    kept_aln = 0
    for line in open(sam_f):
        # Skip headers if they haven't been removed
        if line[0] in ['@', '#']:
            continue
        fields = line.strip('\n').split('\t')
        # Sam line structure
        # 0   read name
        # 1   sam flag
        # 2   contig name
        # 3   map pos start
        # 4   MAPQ
        # 5   CIGAR
        # 6   name of mate
        # 7   position of mate
        # 8   template length
        # 9   sequence
        # 10  PHRED
        # 11+ Additional tag info

        # Alignment elements to keep
        read_id  = fields[0]
        sam_flag = int(fields[1])
        ctg_name = fields[2]
        pos_bp   = int(fields[3])
        cigar    = fields[5]
        sequence = fields[9]
        # Generate ReadAlignment object
        alignment = ReadAlignment(read_id, sam_flag, ctg_name, pos_bp, cigar, sequence)
        total_aln += 1
        # Populate Alignment List if read is mapped
        if alignment.map == 'mapped':
            align_list.append(alignment)
            kept_aln += 1
    # Print to log
    print(f'    Read file: {sam_f}\n    Parsed {total_aln:,} total alignments.\n    Retained {kept_aln:,} ({(kept_aln/total_aln):.3%}) alignments.')
    # Return the alignment list
    return align_list

# Split CIGAR string into CIGAR list
def split_cigar_str(cigar):
    splitcig = []
    prev = 0
    for i, c in enumerate(cigar):
        if c in ALL_CIGAR_OPS:
            splitcig.append(cigar[prev:i+1])
            prev=i+1
            if c not in VALID_CIGAR_OPS:
                sys.exit(f"Error: \'{c}\' in CIGAR {cigar} is an incompatible opertation. Compatible CIGAR operations are: {','.join(VALID_CIGAR_OPS)}. See README for more info.")
    return splitcig

#
# Function to generate an Alignment Pair dictionary from alignment list
# Dictionary will hold the mt and numt alignment for each individual read
# { read_1_id : [ mt_ReadAlignment, nuMt_ReadAlignment ], ... }
def generate_alignment_pair(mt_alignments, numt_alignments):
    # Check inputs
    assert type(mt_alignments) is list
    assert type(numt_alignments) is list
    align_pair_dict = dict()
    # Process mt alignments
    for alignment in mt_alignments:
        # Check inputs
        assert isinstance(alignment, ReadAlignment)
        # Dictionary default is [ mt_None, numt_None ]
        align_pair_dict.setdefault(alignment.rpid, [ None, None ] )
        align_pair_dict[alignment.rpid][0] = alignment
    # Process numt alignments
    for alignment in numt_alignments:
        # Check inputs
        assert isinstance(alignment, ReadAlignment)
        # Dictionary default is [ mt_None, numt_None ]
        align_pair_dict.setdefault(alignment.rpid, [ None, None ] )
        align_pair_dict[alignment.rpid][1] = alignment
    # Return Alignment Pair Dictionary
    return align_pair_dict

#
# Function to process SAMs and generate alignment read pair dictionary
def sam_to_alignment_pair(mt_sam_f, numt_sam_f):
    # Process mt SAM file
    print('\nReading mitochondrial alignment SAM...')
    mt_align = read_sam_file(mt_sam_f)
    # Process numt SAM file
    print('\nReading NUMT alignment SAM...')
    numt_align = read_sam_file(numt_sam_f)
    # Generate alignment read pair
    align_pair_dict = generate_alignment_pair(mt_align, numt_align)
    # Return Alignment Pair
    return align_pair_dict

#
# Compare an alignment against the referece sequences
def compare_alignment(alignment, ref_sequence_dictionary):
    assert type(ref_sequence_dictionary) is dict
    # Finish if alignment is empty
    if alignment is None:
        return AlignmentComparison(None, False, None, None)
    # Continue otherwise
    assert isinstance(alignment, ReadAlignment)
    # Output elements
    aln_size   = 0
    mismatches = 0
    # Determine reference sequence to use
    ref_sequence = ref_sequence_dictionary.get(alignment.ctg, None)
    assert ref_sequence is not None, 'Alignment not in Ref Sequences'
    # Generate elements for comparison
    cigar_list = alignment.cig.per_nt_cigar()
    align_seq = alignment.seq
    aln_size = len(align_seq)
    algn_start = alignment.pos-1
    aln_i = 0 # Index for alignment sequence
    ref_i = 0 # Index for reference sequence
    r_dir = 1

    # Loop over the cigar list index
    for idx in range(len(cigar_list)):
        cig = cigar_list[idx]
        aln_nt = None
        ref_nt = None
        # If match
        if cig == 'M':
            aln_nt = align_seq[aln_i]
            ref_nt = ref_sequence.seq[algn_start + ref_i]
            # Compare the sequences
            if aln_nt != ref_nt:
                mismatches += 1
            aln_i += 1
            ref_i += r_dir
        # If deletion
        elif cig == 'D':
            aln_nt = '-'
            ref_nt = ref_sequence.seq[algn_start + ref_i]
            # Count indel as mismatch
            mismatches += 1
            # Move along reference
            ref_i += r_dir
        # If insertion
        elif cig == 'I':
            aln_nt = align_seq[aln_i]
            ref_nt = '-'
            # Count indel as mismatch
            mismatches += 1
            # Move along alignment
            aln_i += 1
        # If soft clipped
        elif cig == 'S':
            aln_nt = '*'
            ref_nt = '*'
            # If sequence is soft clipped, increase both counters and reduce length of compared sequence
            aln_size -= 1
            aln_i += 1
        # If hard clipped
        elif cig == 'H':
            aln_nt = '*'
            ref_nt = '*'
            # When sequences are hard clipped, only the aligned portion of the read is written to the BAM
            # Thus, DO NOT modify the size of the alignment nor the index

    # Return AlignmentComparison object
    return AlignmentComparison(alignment.rpid, True, aln_size, mismatches)

#
# Compare all alignment pair
# Compare each alignment pair against the corresponding reference and obtain percentage identity
def compare_all_alignments(alignment_pair_dictionary, ref_sequence_dictionary):
    assert type(alignment_pair_dictionary) is dict
    assert type(ref_sequence_dictionary) is dict
    print('\nComparing alignments...')
    n_pairs = 0
    # Pre-generate output
    per_identity_dictionary = dict()
    # Loop over read ids in the pair dictionary
    for read_id in sorted(alignment_pair_dictionary.keys()):
        per_identity_dictionary.setdefault(read_id, [ None, None ] )
        # Extract specific alignment pair
        algn_pair = alignment_pair_dictionary[read_id]
        n_pairs += 1
        # Test mt alignment
        mt_identity = compare_alignment(algn_pair[0], ref_sequence_dictionary)
        per_identity_dictionary[read_id][0] = mt_identity
        # Test numt alignment
        numt_identity = compare_alignment(algn_pair[1], ref_sequence_dictionary)
        per_identity_dictionary[read_id][1] = numt_identity
    # Print to log
    print(f'    Compared alignments pairs for {n_pairs:,} reads.')
    # Return the populated dictonary
    return per_identity_dictionary

#
# Generate output TSV
# Following format
# #read_ID<tab>mt_aln_bp<tab>mt_mismatch<tab>mt_identity<tab>numt_aln_bp<tab>numt_mismatch<tab>numt_identity<tab>Candidate
def generate_output_tsv(per_identity_dictionary, output_f):
    out = open(output_f, 'w')
    out.write('#read_ID\tread_in_pair\tmt_aln_bp\tmt_mismatch\tmt_identity\tnumt_aln_bp\tnumt_mismatch\tnumt_identity\tcandidate\n')
    # Get current identity pair
    for read in sorted(per_identity_dictionary):
        mt_identity = per_identity_dictionary[read][0]
        numt_identity = per_identity_dictionary[read][1]
        candidate = None
        mt_perc = None
        numt_perc = None
        # Check types
        assert isinstance(mt_identity, AlignmentComparison)
        assert isinstance(numt_identity, AlignmentComparison)
        # Check for highest identity
        if numt_identity.status is False:
            candidate = 'cymt'
            mt_perc = f'{mt_identity.per_identity:.6f}'
        elif mt_identity.status is False:
            candidate = 'numt'
            numt_perc = f'{numt_identity.per_identity:.6f}'
        elif mt_identity.per_identity > numt_identity.per_identity:
            candidate = 'cymt'
            mt_perc = f'{mt_identity.per_identity:.6f}'
            numt_perc = f'{numt_identity.per_identity:.6f}'
        elif numt_identity.per_identity > mt_identity.per_identity:
            candidate = 'numt'
            mt_perc = f'{mt_identity.per_identity:.6f}'
            numt_perc = f'{numt_identity.per_identity:.6f}'
        elif mt_identity.per_identity == numt_identity.per_identity:
            candidate = 'undetermined'
            mt_perc = f'{mt_identity.per_identity:.6f}'
            numt_perc = f'{numt_identity.per_identity:.6f}'
        
        # Process some read information
        read_id = read[:-2]
        rp = read[-1:]
        assert rp.isnumeric(), 'Error: incorrect read paired int'
        rp = int(rp)

        # Print into file
        out.write(f'{read_id}\t{rp}\t{mt_identity.aln_len}\t{mt_identity.n_mismatch}\t{mt_perc}\t{numt_identity.aln_len}\t{numt_identity.n_mismatch}\t{numt_perc}\t{candidate}\n')

#
# Main function
def main():
    print(f'{PROG} started on {now()}\n')
    # Parse Arguments
    args = parse_args()

    # Reference fasta files
    mt_fa = args.mt_fasta
    numt_fa = args.numt_fasta
    # Alignment files
    mt_sam = args.mt_sam
    numt_sam = args.numt_sam
    # Output file
    output_f = args.outfile

    # Run numt parser

    # 1. Read mt and numt fasta and generate a ref sequence dictionary
    ref_seq_dict = extract_ref_sequence_dictionary(mt_fa, numt_fa)
    # 2. Load SAM/BAM into read dictionary
    alignment_pair_dict = sam_to_alignment_pair(mt_sam, numt_sam)
    # 3. Compare all alignments against reference sequences
    per_identity_dictionary = compare_all_alignments(alignment_pair_dict, ref_seq_dict)
    # 4. Save output
    generate_output_tsv(per_identity_dictionary, output_f)

    # Finish
    print(f'\n{PROG} finished on {now()}')


# --------
# Run Code
# --------

if __name__ == '__main__':
    main()
