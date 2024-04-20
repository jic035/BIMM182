import sys
import random

def random_dna(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))

def generate_pairs(num_pairs, length):
    return [(random_dna(length), random_dna(length)) for _ in range(num_pairs)]

def write_fasta(sequences, filename="randomDNA.fasta"):
    with open(filename, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">Sequence_{i+1}\n{seq}\n")

def calculate_frequencies(sequences):
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    total_length = 0
    for seq in sequences:
        for nucleotide in seq:
            if nucleotide in counts:
                counts[nucleotide] += 1
        total_length += len(seq)
    
    frequencies = {nuc: count / total_length for nuc, count in counts.items()}
    return frequencies

def main():
    if len(sys.argv) != 3:
        print("Incorrect number of arguments.")
        sys.exit(1)

    num_sequences = int(sys.argv[1])
    length_of_sequences = int(sys.argv[2])
    sequences = []
    for _ in range(num_sequences):
        sequences.append(random_dna(length_of_sequences))    
    write_fasta(sequences)

    frequencies = calculate_frequencies(sequences)
    for nucleotide, freq in frequencies.items():
        print(f"{nucleotide}: {freq:.4f}")

if __name__ == "__main__":
    main()