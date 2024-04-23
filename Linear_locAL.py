import argparse
from typing import Tuple
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Perform local alignment on two DNA sequences")
    parser.add_argument('seq_files', type=str, help='FASTA file containing two sequences')
    parser.add_argument('-m', '--match', type=float, default=2, help='Match reward')
    parser.add_argument('-s', '--mismatch', type=float, default=-1, help='Mismatch penalty')
    parser.add_argument('-d', '--indel', type=float, default=-1, help='Indel penalty')
    parser.add_argument('-a', '--alignment', action='store_true', help='Display the alignment')
    # From Python Docs:
    # action='store_true'
    # if the specific argument is provided in the command line, 
    # then the corresponding attribute is set to True. If the 
    # argument is not provided, then the attribute is set to 
    # False by default.
    return parser.parse_args()

def read_sequences(filename):
    sequences = list(SeqIO.parse(filename, 'fasta'))
    if len(sequences) != 2:
        raise ValueError("Input file must contain exactly two sequences.")
    return sequences[0].seq, sequences[1].seq