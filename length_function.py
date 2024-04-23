import argparse
import matplotlib.pyplot as plt
from randomDNA import generate_pairs
from locAL import LocalAlignment

def run_alignment_tests(sequence_pairs, penalties):
    results = []
    for penalty in penalties:
        print(f"Running alignments with mismatch/indel penalty: {penalty}")
        total_length = 0
        for seq1, seq2 in sequence_pairs:
            _, aligned_seq1, _ = LocalAlignment(seq1, seq2, 1, penalty, penalty)
            total_length += len(aligned_seq1)
        average_length = total_length / len(sequence_pairs)
        results.append((penalty, average_length))
    return results

def plot_results(results):
    penalties, lengths = zip(*results)  # Unpack results into separate lists
    plt.figure(figsize=(10, 6))
    plt.bar(penalties, lengths, color='blue')
    plt.xlabel('Penalty (Mismatch/Indel)')
    plt.ylabel('Average Alignment Length')
    plt.title('Effect of Mismatch/Indel Penalties on Alignment Length')
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Run alignment tests with varying penalties and plot results.")
    parser.add_argument('penalties', nargs='+', type=float, help="List of mismatch/indel penalties to test")
    args = parser.parse_args()

    # Generate 500 pairs of 1000 length random sequences
    sequence_pairs = generate_pairs(500, 1000)
    
    # Run alignment tests
    results = run_alignment_tests(sequence_pairs, args.penalties)
    
    # Plot results
    plot_results(results)

if __name__ == "__main__":
    main()
    
