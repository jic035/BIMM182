import matplotlib.pyplot as plt
import argparse

def read_lengths_from_file(filename):
    lengths = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()  # Split the line into words
            length = int(parts[-1])  # The length value is the last element
            lengths.append(length)
    return lengths

def plot_histogram(lengths, bin):
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=bin)
    plt.title('Histogram of Local Alignment Lengths')
    plt.xlabel('Length of Alignment')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot a histogram of alignment lengths.")
    parser.add_argument('filename', type=str, help="The filename containing the alignment lengths.")
    parser.add_argument('-b', '--bin', type=int, default=10, help='Histogram bin size.')
    args = parser.parse_args()
    filename = "./Alignment_out/" + args.filename
    lengths = read_lengths_from_file(filename)
    plot_histogram(lengths, args.bin)

if __name__ == "__main__":
    main()
