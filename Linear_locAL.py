import argparse
from Bio import SeqIO
import numpy as np
from locAL import LocalAlignment
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Perform local alignment on two DNA sequences")
    parser.add_argument('seq_files', type=str, help='FASTA file containing two sequences')
    parser.add_argument('-m', '--match', type=float, default=2, help='Match reward')
    parser.add_argument('-s', '--mismatch', type=float, default=-1, help='Mismatch penalty')
    parser.add_argument('-d', '--indel', type=float, default=-1, help='Indel penalty')
    parser.add_argument('-a', '--alignment', action='store_true', help='Display the alignment')
    parser.add_argument('-t', '--time', action='store_true', help='Display the running time')
    return parser.parse_args()

def read_sequences(filename):
    sequences = list(SeqIO.parse(filename, 'fasta'))
    if len(sequences) != 2:
        raise ValueError("Input file must contain exactly two sequences.")
    return sequences[0].seq, sequences[1].seq

def linear_space_local_alignment(seq1, seq2, match_reward, mismatch_penalty, indel_penalty):
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

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = prev_score_row[j - 1] + (match_reward if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = prev_score_row[j] + indel_penalty
            insert = current_score_row[j - 1] + indel_penalty
            current_score_row[j] = max(0, match, delete, insert)

            if current_score_row[j] == match:
                current_length_row[j] = prev_length_row[j - 1] + 1
            elif current_score_row[j] == delete:
                current_length_row[j] = prev_length_row[j] + 1
            elif current_score_row[j] == insert:
                current_length_row[j] = current_length_row[j - 1] + 1
            else:
                current_length_row[j] = 0

            if current_score_row[j] > max_score:
                max_score = current_score_row[j]
                max_length = current_length_row[j]
                max_pos = (i, j)

        prev_score_row, current_score_row = current_score_row, np.zeros(m + 1, dtype=int)
        prev_length_row, current_length_row = current_length_row, np.zeros(m + 1, dtype=int)

    return max_score, max_length, max_pos
        
# seq1 = "AGCCCACAATATAGCTAAACATCGCCACTAGTCCTAAAACCACTCCGCAGAAAAAGAATAAGGCCAAAACACGACTAAAATCGAAAGACATGACAAGTAAACGAGAAAAGAAAAATAAACGACATACACACTTGTAGGAAAAATAAGAAAAAGGGGGAGACGAAGCAAAGAAAGGGCAGCTAACCCTCAAGGAAGAACCAGACAGAATAAGAAAAACCCGAAAGCCACCAAATGAAAGGACAATAACACCTAAGAGCAAAATCAATAAAACACCGATCCTCCGAGGATAACCAAGAGAGACCTAAGAACGACAAGAAACCAATGAAAGAAAAAGAAAATGGACATCAGAACGACTTAGAATGCTGGGAAAAAGAAAAATTATAAACGAAGGATGGGCATAAATTGGACGAAGCCAAGAGATAGGCCGAGATAAAACGGAGAACAATAAGGGAGACCATGGAGAGCAAACCAACCGCAACAAATAAAGGGGGGGACAAAAACAAGACCAACCCAAACTGTCAGACAGGAAGAGCAATAACCAAGACAGAAGAAGAAACAGGAGACAAACAACATAATATAAGAGCACCTAGCTAACAAAAAAGACCAGCAAACGGATTAAGAAGATAAAGAAAACGTAAAGAACAGTCAAGGAACAAGCGATAATAAATGCAGGGAAAAAATGGGGACAGACGAAGGAAACAACCAGAAATAATCTAACGCATCGCAGAAGATGACACTGCGAGAAAATACGAGCCGTATACGACACAAAACCGGGAATAAAGAAAAAAACCATACCCAAAAAGAACAACGCGAAAGATGAAACGCTCCCAACTCGGATGAGCAAAGCCGCCAGGCCAAAAAAGAGAACCAGAGCAGAGCGAAGCTATGGGTAGAAAACACCCTAAGCGCGGGTAGTAGAGACGAAAAATAAAAACAGGCTGACCCGAACATAAGAGCCCACACAAGTAGAAGAACGGAAAGAAAACGAAAAGACCGGAGC"
# seq2 = "AGTAAGGCTGCAGCAAACCAGAAACCATCCTAAAAGAAGGAAAGCAAATAGGAGAACAAAACAATAATCAAAGCGATATAAAGGAGAAACATAACAACCGCTACACCAATCACACCACAAAGGAAAGATAAGCGCAAAGAGAAGTACCCTGCATACCTACAACCAAATAAAAGAGGAGAACTGAGAAACGCCACCAAAACAAACGTGACGATAACTAATGAAACGAATGAAAAAGGAGGATAGACCTCAAATTCAAAGAGATGAACATGACAGCTAAAAGACAACGAGCAAAAATGCTAGGAGACATAACCAAGCTAAAGACCAGGACCCAACGACCGAACGCAAGAATGAAAATTAGCCCCAAAAAACGCGCACGAGAAAGAACTAAAGAGCCACAAAGACACAATAGAAAGTGCTCGACGGACAAAAAAAAAAAAGAAGAGACACAAATAGAAACAAAAACAAAAGAAGCAAACGAATATGACAAAGAAAACCAACTACCAAAAGCAGTACATGACACATCATGCACACAGCGAAAACAACAAAATGAACAACGAAAACAACCAAAAAGACGAATCGAACAGGAGAGAGGATCCCCAAAAAAGAGGGCCCAACTAAGACAATGCAAAGAAACGCGACAAAGCCTCGCCAACAGAATCAACCAAAGCATGAACAGCACTTTTAAAACATGTGGCGCGGCGTCGAGCAGTACGGTTCAAATGCAAAAATTACAAAAAAGACATGCACTGAACCCCGTAAAGAACGAGAAACTTCAAGAAGGAAGAAGACATAAGCCAAAAAACCAAATAATAGACACAGCTCGAGAAAAAGCCCAACAAACGACAGAAAAGAAGGGGAGTAAGGAGAAGCAAAACAGAACAGGGGCGAGAACGACTGCGGAGTAATCGAAAGACATGACAAGTAAACGAGAAAAGATTAATAAACGACATACACAAACAAAAACAATAGAAACGAAAAAAAAAAAAAATTAAGCAAGCAC"
# max_score, max_length, max_pos = linear_space_local_alignment(seq1, seq2, match_reward=1, mismatch_penalty=-10, indel_penalty=-1)

# print(f"Maximum alignment score: {max_score}")
# i, j= max_pos
# print('length',max_length)
# if i < max_length:
#     p_seq1 = seq1[:i+1]
#     p_seq2 = seq2[:j+1]
# else:
#     p_seq1 = seq1[i-max_length:i+1]
#     p_seq2 = seq2[j-max_length:j+1]
# print('pos',max_pos)
# print(p_seq1, p_seq2)
# score, aligned_s, aligned_t = LocalAlignment(p_seq1, p_seq2, match_reward=1, mismatch_penalty=-10, indel_penalty=-1)
# print("Alignment:")
# print(''.join(aligned_s))
# print(''.join(aligned_t))
# print("Alignment length:", max_length)
# print("locAL score:", score)

def main():
    start_time = time.time()
    args = parse_args()
    seq1, seq2 = read_sequences(args.seq_files)
    max_score, max_length, max_pos = linear_space_local_alignment(seq1, seq2, args.match, args.mismatch, args.indel)
    if args.alignment:
        i, j= max_pos
        if i < max_length:
            p_seq1 = seq1[:i+1]
            p_seq2 = seq2[:j+1]
        else:
            p_seq1 = seq1[i-max_length:i+1]
            p_seq2 = seq2[j-max_length:j+1]
        print('sequences:', p_seq1, p_seq2)
        _, aligned_s, aligned_t = LocalAlignment(p_seq1, p_seq2, args.match, args.mismatch, args.indel)
        print("Alignment:")
        print(''.join(aligned_s))
        print(''.join(aligned_t))

    print(f"Score of the best local alignment: {max_score}")
    print(f"Length of the best local alignment: {max_length}")
    end_time = time.time()
    if args.time:
        print(f"Execution time: {end_time - start_time} seconds")

if __name__ == "__main__":
    main()