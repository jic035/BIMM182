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

def LocalAlignment(s: str, t: str, match_reward: float, mismatch_penalty: float, indel_penalty: float) -> Tuple[float, str, str]:
    s_len = len(s)
    t_len = len(t)
    score_matrix = [[0] * (s_len + 1) for _ in range(t_len + 1)]
    backtrack = [[None] * (s_len + 1) for _ in range(t_len + 1)]
    max_score = 0
    max_position = (0, 0)
    # Do not initialize first column/row to indel penalty 
    # Fill the scoring matrix
    for i in range(1, t_len + 1):
        for j in range(1, s_len + 1):
            match = score_matrix[i-1][j-1] + (match_reward if s[j-1] == t[i-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + indel_penalty
            insert = score_matrix[i][j-1] + indel_penalty
            score_matrix[i][j] = max(0, match, delete, insert)  # Adding free-ride
            if score_matrix[i][j] == match:
                backtrack[i][j] = "diag"
            elif score_matrix[i][j] == delete:
                backtrack[i][j] = "up"
            elif score_matrix[i][j] == insert:
                backtrack[i][j] = "left"
            # Update max score
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)

    # Traceback from the best local alignment
    i, j = max_position
    aligned_s, aligned_t = '', ''
    while score_matrix[i][j] != 0:
        if backtrack[i][j] == "diag":
            aligned_s = s[j-1] + aligned_s
            aligned_t = t[i-1] + aligned_t
            i -= 1
            j -= 1
        elif backtrack[i][j] == "up":
            aligned_s = '-' + aligned_s
            aligned_t = t[i-1] + aligned_t
            i -= 1
        elif backtrack[i][j] == "left":
            aligned_s = s[j-1] + aligned_s
            aligned_t = '-' + aligned_t
            j -= 1

    return max_score, aligned_s, aligned_t

def main():
    args = parse_args()
    seq1, seq2 = read_sequences(args.seq_files)
    max_score, aligned_s, aligned_t = LocalAlignment(seq1, seq2, args.match, args.mismatch, args.indel)
    if args.alignment:
        print("Alignment:")
        print(''.join(aligned_s))
        print(''.join(aligned_t))

    print(f"Score of the best local alignment: {max_score}")
    print(f"Length of the best local alignment: {len(aligned_s)}")

if __name__ == "__main__":
    main()