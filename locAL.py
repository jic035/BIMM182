import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Perform local alignment on two DNA sequences. Enter signed scores.")
    parser.add_argument('seq_files', type=str, help='FASTA file containing two sequences')
    parser.add_argument('-m', '--match', type=int, default=2, help='Match reward')
    parser.add_argument('-s', '--mismatch', type=int, default=-1, help='Mismatch penalty')
    parser.add_argument('-d', '--indel', type=int, default=-1, help='Indel penalty')
    parser.add_argument('-a', '--alignment', action='store_true', help='Display the alignment')
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

def initialize_scoring_matrix(rows, cols):
    return [[0] * cols for _ in range(rows)] # Local alignment: all 0's, including first column/row

def calculate_local_alignment(seq1, seq2, match, mismatch, indel): # All scores are with signs
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = initialize_scoring_matrix(rows, cols)
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            match = score_matrix[i-1][j-1] + (match if seq1[j-1] == seq2[i-1] else mismatch)
            delete = score_matrix[i-1][j] + indel
            insert = score_matrix[i][j-1] + indel
            score = max(0, match, delete, insert)  # Adding free-ride for local alignment
            score_matrix[i][j] = score

            if score > max_score:
                max_score = score
                max_pos = (i, j)

    return max_score, max_pos, score_matrix

def traceback(score_matrix, seq1, seq2, start_pos, match, mismatch, indel):
    alignment = []
    i, j = start_pos
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if seq1[i-1] == seq2[j-1]:
            alignment.append((seq1[i-1], seq2[j-1]))
            i -= 1
            j -= 1
        elif score_matrix[i-1][j] + indel == score_matrix[i][j]:
            alignment.append((seq1[i-1], '-'))
            i -= 1
        else:
            alignment.append(('-', seq2[j-1]))
            j -= 1
    alignment.reverse()
    return alignment

def main():
    args = parse_args()
    seq1, seq2 = read_sequences(args.seq_files)
    max_score, max_pos, score_matrix = calculate_local_alignment(seq1, seq2, args.match, args.mismatch, args.indel)

    if args.alignment:
        alignment = traceback(score_matrix, seq1, seq2, max_pos, args.match, args.mismatch, args.indel)
        aligned_seq1, aligned_seq2 = zip(*alignment)
        print("Alignment:")
        print(''.join(aligned_seq1))
        print(''.join(aligned_seq2))

    print(f"Score of the best local alignment: {max_score}")
    print(f"Length of the best local alignment: {len(alignment)}")

if __name__ == "__main__":
    main()
