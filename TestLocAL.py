import numpy as np

def local_alignment_with_length(seq1, seq2, match_score, mismatch_penalty, indel_penalty):
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


seq1 = "AGCCCACAATATAGCTAAACATCGCCACTAGTCCTAAAACCACTCCGCAGAAAAAGAATAAGGCCAAAACACGACTAAAATCGAAAGACATGACAAGTAAACGAGAAAAGAAAAATAAACGACATACACACTTGTAGGAAAAATAAGAAAAAGGGGGAGACGAAGCAAAGAAAGGGCAGCTAACCCTCAAGGAAGAACCAGACAGAATAAGAAAAACCCGAAAGCCACCAAATGAAAGGACAATAACACCTAAGAGCAAAATCAATAAAACACCGATCCTCCGAGGATAACCAAGAGAGACCTAAGAACGACAAGAAACCAATGAAAGAAAAAGAAAATGGACATCAGAACGACTTAGAATGCTGGGAAAAAGAAAAATTATAAACGAAGGATGGGCATAAATTGGACGAAGCCAAGAGATAGGCCGAGATAAAACGGAGAACAATAAGGGAGACCATGGAGAGCAAACCAACCGCAACAAATAAAGGGGGGGACAAAAACAAGACCAACCCAAACTGTCAGACAGGAAGAGCAATAACCAAGACAGAAGAAGAAACAGGAGACAAACAACATAATATAAGAGCACCTAGCTAACAAAAAAGACCAGCAAACGGATTAAGAAGATAAAGAAAACGTAAAGAACAGTCAAGGAACAAGCGATAATAAATGCAGGGAAAAAATGGGGACAGACGAAGGAAACAACCAGAAATAATCTAACGCATCGCAGAAGATGACACTGCGAGAAAATACGAGCCGTATACGACACAAAACCGGGAATAAAGAAAAAAACCATACCCAAAAAGAACAACGCGAAAGATGAAACGCTCCCAACTCGGATGAGCAAAGCCGCCAGGCCAAAAAAGAGAACCAGAGCAGAGCGAAGCTATGGGTAGAAAACACCCTAAGCGCGGGTAGTAGAGACGAAAAATAAAAACAGGCTGACCCGAACATAAGAGCCCACACAAGTAGAAGAACGGAAAGAAAACGAAAAGACCGGAGC" 
seq2 = "AGTAAGGCTGCAGCAAACCAGAAACCATCCTAAAAGAAGGAAAGCAAATAGGAGAACAAAACAATAATCAAAGCGATATAAAGGAGAAACATAACAACCGCTACACCAATCACACCACAAAGGAAAGATAAGCGCAAAGAGAAGTACCCTGCATACCTACAACCAAATAAAAGAGGAGAACTGAGAAACGCCACCAAAACAAACGTGACGATAACTAATGAAACGAATGAAAAAGGAGGATAGACCTCAAATTCAAAGAGATGAACATGACAGCTAAAAGACAACGAGCAAAAATGCTAGGAGACATAACCAAGCTAAAGACCAGGACCCAACGACCGAACGCAAGAATGAAAATTAGCCCCAAAAAACGCGCACGAGAAAGAACTAAAGAGCCACAAAGACACAATAGAAAGTGCTCGACGGACAAAAAAAAAAAAGAAGAGACACAAATAGAAACAAAAACAAAAGAAGCAAACGAATATGACAAAGAAAACCAACTACCAAAAGCAGTACATGACACATCATGCACACAGCGAAAACAACAAAATGAACAACGAAAACAACCAAAAAGACGAATCGAACAGGAGAGAGGATCCCCAAAAAAGAGGGCCCAACTAAGACAATGCAAAGAAACGCGACAAAGCCTCGCCAACAGAATCAACCAAAGCATGAACAGCACTTTTAAAACATGTGGCGCGGCGTCGAGCAGTACGGTTCAAATGCAAAAATTACAAAAAAGACATGCACTGAACCCCGTAAAGAACGAGAAACTTCAAGAAGGAAGAAGACATAAGCCAAAAAACCAAATAATAGACACAGCTCGAGAAAAAGCCCAACAAACGACAGAAAAGAAGGGGAGTAAGGAGAAGCAAAACAGAACAGGGGCGAGAACGACTGCGGAGTAATCGAAAGACATGACAAGTAAACGAGAAAAGATTAATAAACGACATACACAAACAAAAACAATAGAAACGAAAAAAAAAAAAAATTAAGCAAGCAC"  
max_score, max_length, max_pos = local_alignment_with_length(seq1, seq2, 1, -10, -2)
print(f"Maximum score: {max_score}, Alignment length: {max_length}, Position: {max_pos}")