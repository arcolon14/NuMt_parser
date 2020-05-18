#!/usr/bin/env python3
import gzip, os.path, argparse
from os import path
#
# Script to compare read similarity to Mitochondrial and Nu-MT sequences
# (c) 2020 Angel G. Rivera-Colon & Alida de Flamingh
#

# -----------
# Input files
# -----------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--mt-fasta', required=True, help='Mitochondrial Sequence FASTA')
    p.add_argument('--numt-fasta', required=True, help='nuMT Sequence FASTA')
    p.add_argument('--mt-sam', required=True, help='Read alignments to the Mt reference in SAM format')
    p.add_argument('--numt-sam', required=True, help='Read alignments to the nuMt reference in SAM format')
    p.add_argument('--outfile', required=True, help='Path and name to the output TSV file. Example: ./<sample_id>.tsv')
    args = p.parse_args()
    return args

args = parse_args()

# Reference fasta files
mt_fa = args.mt_fasta
numt_fa = args.numt_fasta
# Alignment files
mt_sam = args.mt_sam
numt_sam = args.numt_sam
# Output file
output_f = args.outfile

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
        return '{} : {:,} bp'.format(self.id, self.len)

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
        if len(f'{self.sam:b}') >= 5:
            if int(list(f'{self.sam:b}')[-5]) is 1:
                self.dir = 'reverse'
        # Determine if read is mapped
        self.map = 'mapped'
        if len(f'{self.sam:b}') >= 3:
            if int(list(f'{self.sam:b}')[-3]) is 1:
                self.map = 'unmapped'
        # Process CIGARs
        self.cig = CIGAR(cigar)
    def __str__(self):
        return '''Read ID: {}
SAM flag: {}
Ctg Name: {}
Align BP: {:,}
CIGAR: {}
Seq: {}
Direction: {}
Status: {}
'''.format(self.rid, self.sam, self.ctg, self.pos, self.cig, self.seq, self.dir, self.map)

#
# CIGAR string class
class CIGAR:
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
            assert ctype in ('M', 'I', 'D', 'S', 'H')
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
            if c in ('M', 'I', 'D', 'S', 'H'):
                splitcig.append(cigar[prev:i+1])
                prev=i+1
        return splitcig

# ---------
# Functions
# ---------

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
                name = line[1:].split()[0]
                seq = []
            elif line[0] in ['#', '.']:
                # Ignore comment lines
                continue
            else:
                # if it is sequence
                seq.append(line.upper())
    fasta_dict[name] = RefSequence(name, ''.join(seq))
    # Return genome dictionary
    return fasta_dict

#
# Function to read both MT and nuMT fasta files and merge into a single final dictionary
# { sequence_1 : RefSequence, sequence_2 : RefSequence, ... }
def extract_ref_sequence_dictionary(mt_fasta, nuMt_fasta):
    # Read MT fasta
    mt_seq_dict = read_fasta(mt_fasta)
    # Read nuMT fasta
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
        # Populate Alignment List if read is mapped
        if alignment.map is 'mapped':
            align_list.append(alignment)
    # Return the alignment list
    return align_list

# Split CIGAR string into CIGAR list
def split_cigar_str(cigar):
    splitcig = []
    prev = 0
    for i, c in enumerate(cigar):
        if c in ('M', 'I', 'D', 'S', 'H'):
            splitcig.append(cigar[prev:i+1])
            prev=i+1
    return splitcig

#
# Function to generate an Alignment Pair dictionary from alignment list
# Dictionary will hold the Mt and nuMt alignment for each individual read
# { read_1_id : [ mt_ReadAlignment, nuMt_ReadAlignment ], ... }
def generate_alignment_pair(mt_alignments, numt_alignments):
    # Check inputs
    assert type(mt_alignments) is list
    assert type(numt_alignments) is list
    align_pair_dict = dict()
    # Process Mt alignments
    for alignment in mt_alignments:
        # Check inputs
        assert isinstance(alignment, ReadAlignment)
        # Dictionary default is [ mt_None, numt_None ]
        align_pair_dict.setdefault(alignment.rid, [ None, None ] )
        align_pair_dict[alignment.rid][0] = alignment
    # Process nuMt alignments
    for alignment in numt_alignments:
        # Check inputs
        assert isinstance(alignment, ReadAlignment)
        # Dictionary default is [ mt_None, numt_None ]
        align_pair_dict.setdefault(alignment.rid, [ None, None ] )
        align_pair_dict[alignment.rid][1] = alignment
    # Return Alignment Pair Dictionary
    return align_pair_dict

#
# Function to process SAMs and generate alignment read pair dictionary
def sam_to_alignment_pair(mt_sam_f, numt_sam_f):
    # Process Mt SAM file
    mt_align = read_sam_file(mt_sam_f)
    # Process nuMt SAM file
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
        return None
    # Continue otherwise
    assert isinstance(alignment, ReadAlignment)
    # Output elements
    read_size  = 0
    mismatches = 0
    # Determine reference sequence to use
    ref_sequence = ref_sequence_dictionary.get(alignment.ctg, None)
    assert ref_sequence is not None, 'Alignment not in Ref Sequences'
    # Generate elements for comparison
    cigar_list = alignment.cig.per_nt_cigar()
    align_seq = alignment.seq
    read_size = len(align_seq)
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
        if cig is 'M':
            aln_nt = align_seq[aln_i]
            ref_nt = ref_sequence.seq[algn_start + ref_i]
            # Compare the sequences
            if aln_nt != ref_nt:
                mismatches += 1
            aln_i += 1
            ref_i += r_dir
        # If deletion
        if cig is 'D':
            aln_nt = '-'
            ref_nt = ref_sequence.seq[algn_start + ref_i]
            # Count indel as mismatch
            mismatches += 1
            # Move along reference
            ref_i += r_dir
        # If insertion
        if cig is 'I':
            aln_nt = align_seq[aln_i]
            ref_nt = '-'
            # Count indel as mismatch
            mismatches += 1
            # Move along alignment
            aln_i += 1
        # If clipped
        if cig in ['S', 'H']:
            aln_nt = '*'
            ref_nt = '*'
            # If sequence is clipped, increase both counters and reduce length of compared sequence
            read_size -= 1
            aln_i += 1
        # print(aln_i, ref_i, aln_nt, ref_nt, mismatches)

    # Calculate percent identity
    per_identify = 1 - (mismatches/read_size)
    return per_identify

#
# Compare all alignment pair
# Compare each alignment pair against the corresponding reference and obtain percentage identity
def compare_all_alignments(alignment_pair_dictionary, ref_sequence_dictionary):
    assert type(alignment_pair_dictionary) is dict
    assert type(ref_sequence_dictionary) is dict
    # Pre-generate output
    per_identity_dictionary = dict()
    # Loop over read ids in the pair dictionary
    for read_id in sorted(alignment_pair_dict.keys()):
        per_identity_dictionary.setdefault(read_id, [ None, None ] )
        # Extract specific alignment pair
        algn_pair = alignment_pair_dictionary[read_id]
        # Test Mt alignment
        mt_identity = compare_alignment(algn_pair[0], ref_sequence_dictionary)
        per_identity_dictionary[read_id][0] = mt_identity
        # Test nuMt alignment
        numt_identity = compare_alignment(algn_pair[1], ref_sequence_dictionary)
        per_identity_dictionary[read_id][1] = numt_identity
    # Return the populated dictonary
    return per_identity_dictionary

#
# Generate output TSV
# Following format
# #Read_id<tab>Mt_identity<tab>nuMT_identity<tab>Candidate
def generate_output_tsv(per_identity_dictionary, output_f):
    out = open(output_f, 'w')
    out.write('#Read_ID\tMt_identity\tnuMt_identity\tCandidate\n')
    # Get current identity pair
    for read in sorted(per_identity_dictionary):
        mt_identity = per_identity_dictionary[read][0]
        numt_identity = per_identity_dictionary[read][1]
        candidate = None
        # Check for highest identity
        if numt_identity is None:
            candidate = 'Mt'
            mt_identity = f'{mt_identity:.6f}'
        elif mt_identity is None:
            candidate = 'nuMt'
            numt_identity = f'{numt_identity:.6f}'
        elif mt_identity > numt_identity:
            candidate = 'Mt'
            mt_identity = f'{mt_identity:.6f}'
            numt_identity = f'{numt_identity:.6f}'
        elif numt_identity > mt_identity:
            candidate = 'nuMt'
            mt_identity = f'{mt_identity:.6f}'
            numt_identity = f'{numt_identity:.6f}'
        elif mt_identity == numt_identity:
            candidate = 'Unknown'
            mt_identity = f'{mt_identity:.6f}'
            numt_identity = f'{numt_identity:.6f}'
        # Print into file
        out.write(f'{read}\t{mt_identity}\t{numt_identity}\t{candidate}\n')


# --------
# Run Code
# --------

# 1. Read MT and nuMT fasta and generate a ref sequence dictionary
ref_seq_dict = extract_ref_sequence_dictionary(mt_fa, numt_fa)
# 2. Load SAM/BAM into read dictionary
alignment_pair_dict = sam_to_alignment_pair(mt_sam, numt_sam)
# 3. Compare all alignments against reference sequences
per_identity_dictionary = compare_all_alignments(alignment_pair_dict, ref_seq_dict)
# 4. Save output
generate_output_tsv(per_identity_dictionary, output_f)
