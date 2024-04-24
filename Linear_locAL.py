import argparse
from Bio import SeqIO
import numpy as np
from locAL import LocalAlignment

def parse_args():
    parser = argparse.ArgumentParser(description="Perform local alignment on two DNA sequences")
    parser.add_argument('seq_files', type=str, help='FASTA file containing two sequences')
    parser.add_argument('-m', '--match', type=float, default=2, help='Match reward')
    parser.add_argument('-s', '--mismatch', type=float, default=-1, help='Mismatch penalty')
    parser.add_argument('-d', '--indel', type=float, default=-1, help='Indel penalty')
    parser.add_argument('-a', '--alignment', action='store_true', help='Display the alignment')
    return parser.parse_args()

def read_sequences(filename):
    sequences = list(SeqIO.parse(filename, 'fasta'))
    if len(sequences) != 2:
        raise ValueError("Input file must contain exactly two sequences.")
    return sequences[0].seq, sequences[1].seq

def linear_space_local_alignment(seq1, seq2, match_score, mismatch_penalty, indel_penalty):
    n = len(seq1)
    m = len(seq2)

    # Initialize rows for scoring and length
    prev_score_row = np.zeros(m + 1, dtype=int)
    current_score_row = np.zeros(m + 1, dtype=int)
    prev_length_row = np.zeros(m + 1, dtype=int)
    current_length_row = np.zeros(m + 1, dtype=int)

    max_score = 0
    max_length = 0
    max_pos = (0, 0)

    # Fill scoring matrix and track length
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = prev_score_row[j - 1] + match_score
            else:
                score = prev_score_row[j - 1] + mismatch_penalty

            # Calculating max score and choosing path
            if score > max(prev_score_row[j] + indel_penalty, current_score_row[j - 1] + indel_penalty, 0):
                current_score_row[j] = score
                if score == prev_score_row[j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
                    current_length_row[j] = prev_length_row[j - 1] + 1
                else:
                    current_length_row[j] = 0  # Reset if no continuation of alignment
            else:
                current_score_row[j] = max(prev_score_row[j] + indel_penalty, current_score_row[j - 1] + indel_penalty, 0)
                # Track length for indels
                if current_score_row[j] == prev_score_row[j] + indel_penalty:
                    current_length_row[j] = prev_length_row[j] + 1
                elif current_score_row[j] == current_score_row[j - 1] + indel_penalty:
                    current_length_row[j] = current_length_row[j - 1] + 1
                else:
                    current_length_row[j] = 0  # Reset if alignment does not continue

            # Check if we have a new maximum score
            if current_score_row[j] > max_score:
                max_score = current_score_row[j]
                max_length = current_length_row[j]
                max_pos = (i, j)

        # Swap current and previous rows for next iteration
        prev_score_row, current_score_row = current_score_row, np.zeros(m + 1, dtype=int)
        prev_length_row, current_length_row = current_length_row, np.zeros(m + 1, dtype=int)

    return max_score, max_length, max_pos
        
# seq1 = "AGCCCACAATATAGCTAAACATCGCCACTAGTCCTAAAACCACTCCGCAGAAAAAGAATAAGGCCAAAACACGACTAAAATCGAAAGACATGACAAGTAAACGAGAAAAGAAAAATAAACGACATACACACTTGTAGGAAAAATAAGAAAAAGGGGGAGACGAAGCAAAGAAAGGGCAGCTAACCCTCAAGGAAGAACCAGACAGAATAAGAAAAACCCGAAAGCCACCAAATGAAAGGACAATAACACCTAAGAGCAAAATCAATAAAACACCGATCCTCCGAGGATAACCAAGAGAGACCTAAGAACGACAAGAAACCAATGAAAGAAAAAGAAAATGGACATCAGAACGACTTAGAATGCTGGGAAAAAGAAAAATTATAAACGAAGGATGGGCATAAATTGGACGAAGCCAAGAGATAGGCCGAGATAAAACGGAGAACAATAAGGGAGACCATGGAGAGCAAACCAACCGCAACAAATAAAGGGGGGGACAAAAACAAGACCAACCCAAACTGTCAGACAGGAAGAGCAATAACCAAGACAGAAGAAGAAACAGGAGACAAACAACATAATATAAGAGCACCTAGCTAACAAAAAAGACCAGCAAACGGATTAAGAAGATAAAGAAAACGTAAAGAACAGTCAAGGAACAAGCGATAATAAATGCAGGGAAAAAATGGGGACAGACGAAGGAAACAACCAGAAATAATCTAACGCATCGCAGAAGATGACACTGCGAGAAAATACGAGCCGTATACGACACAAAACCGGGAATAAAGAAAAAAACCATACCCAAAAAGAACAACGCGAAAGATGAAACGCTCCCAACTCGGATGAGCAAAGCCGCCAGGCCAAAAAAGAGAACCAGAGCAGAGCGAAGCTATGGGTAGAAAACACCCTAAGCGCGGGTAGTAGAGACGAAAAATAAAAACAGGCTGACCCGAACATAAGAGCCCACACAAGTAGAAGAACGGAAAGAAAACGAAAAGACCGGAGC"
# seq2 = "AGTAAGGCTGCAGCAAACCAGAAACCATCCTAAAAGAAGGAAAGCAAATAGGAGAACAAAACAATAATCAAAGCGATATAAAGGAGAAACATAACAACCGCTACACCAATCACACCACAAAGGAAAGATAAGCGCAAAGAGAAGTACCCTGCATACCTACAACCAAATAAAAGAGGAGAACTGAGAAACGCCACCAAAACAAACGTGACGATAACTAATGAAACGAATGAAAAAGGAGGATAGACCTCAAATTCAAAGAGATGAACATGACAGCTAAAAGACAACGAGCAAAAATGCTAGGAGACATAACCAAGCTAAAGACCAGGACCCAACGACCGAACGCAAGAATGAAAATTAGCCCCAAAAAACGCGCACGAGAAAGAACTAAAGAGCCACAAAGACACAATAGAAAGTGCTCGACGGACAAAAAAAAAAAAGAAGAGACACAAATAGAAACAAAAACAAAAGAAGCAAACGAATATGACAAAGAAAACCAACTACCAAAAGCAGTACATGACACATCATGCACACAGCGAAAACAACAAAATGAACAACGAAAACAACCAAAAAGACGAATCGAACAGGAGAGAGGATCCCCAAAAAAGAGGGCCCAACTAAGACAATGCAAAGAAACGCGACAAAGCCTCGCCAACAGAATCAACCAAAGCATGAACAGCACTTTTAAAACATGTGGCGCGGCGTCGAGCAGTACGGTTCAAATGCAAAAATTACAAAAAAGACATGCACTGAACCCCGTAAAGAACGAGAAACTTCAAGAAGGAAGAAGACATAAGCCAAAAAACCAAATAATAGACACAGCTCGAGAAAAAGCCCAACAAACGACAGAAAAGAAGGGGAGTAAGGAGAAGCAAAACAGAACAGGGGCGAGAACGACTGCGGAGTAATCGAAAGACATGACAAGTAAACGAGAAAAGATTAATAAACGACATACACAAACAAAAACAATAGAAACGAAAAAAAAAAAAAATTAAGCAAGCAC"
# max_score, max_pos, aligned_1, aligned_2 = linear_space_local_alignment(seq1, seq2, match_score=1, mismatch_penalty=-10, indel_penalty=-10)
# print(f"Maximum alignment score: {max_score}")
# print(aligned_1)
# print(aligned_2)
# print(len(aligned_1))


def main():
    args = parse_args()
    seq1, seq2 = read_sequences(args.seq_files)
    max_score, max_length, max_pos = linear_space_local_alignment(seq1, seq2, args.match, args.mismatch, args.indel)
    if args.alignment:
        i, j= max_pos
        p_seq1 = seq1[i-max_length:i+1]
        p_seq2 = seq2[j-max_length:j+1]
        print('sequences:', p_seq1, p_seq2)
        _, aligned_s, aligned_t = LocalAlignment(p_seq1, p_seq2, args.match, args.mismatch, args.indel)
        print("Alignment:")
        print(''.join(aligned_s))
        print(''.join(aligned_t))

    print(f"Score of the best local alignment: {max_score}")
    print(f"Length of the best local alignment: {max_length}")

if __name__ == "__main__":
    main()