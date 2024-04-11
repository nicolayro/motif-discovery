import math
import numpy as np
import sys

DEBUG = False
CUTOFF = 0.00001
MAX_ITERATIONS = 100
INIT_WEIGHT = 0.7
PSEUDO_COUNTS = 0.1

alphabet = "ACGT"
letter_to_index = {"A": 0, "C": 1, "G": 2, "T": 3}
letter_counts = [0, 0, 0, 0]

# Example data
sequences = ["GTCAGG", "GAGAGT", "ACGGAG", "CCAGTC"]
L = 6
W = 3  # Length of motif that is being searched for
M = L - W + 1  # Number of possible starting positions


def count_letters():
    counts = []
    for letter in alphabet:
        count = 0
        for sequence in sequences:
            count += sum([c == letter for c in sequence])
        counts.append(count)
    return counts


def print_matrix(m):
    for row in m:
        for element in row:
            print(f"{element:.4f}", end='  ')
        print()


def find_starting_motifs():
    motifs = {}
    for sequence in sequences:
        for i in range(M):
            motif = sequence[i:i+W]
            motifs[motif] = True
    return motifs


def init_pwm(motif):
    """
    Initializes a Probability Weight Matrix (PWM) based on the given motif
    The format of the PWM is as follows:

    Motif           GAG
    Sequence        CCGAGCC
    Total seqs.     4
    Starting pos.   3

        background  |
            0       |   1   2   3
    --------------------------------------
    A     0.25      |   0.1 0.7 0.1
    C     0.25      |   0.1 0.1 0.1
    G     0.25      |   0.7 0.1 0.7
    T     0.25      |   0.1 0.1 0.1
    """
    PWM = np.full((len(alphabet), len(motif) + 1), 1/len(sequences))
    for i, letter in enumerate(alphabet):
        for j, s in enumerate(motif):
            # +1 as we don't set the background values
            PWM[i, j+1] = INIT_WEIGHT if s == letter else (1-INIT_WEIGHT)/W
    return PWM


def P(i, j, PWM):
    sequence = sequences[i]
    before, motif, after = sequence[:j], sequence[j:j+W], sequence[j+W:M+1]
    p = 1
    for letter in before:
        p *= PWM[letter_to_index[letter], 0]
    for i, letter in enumerate(motif):
        p *= PWM[letter_to_index[letter], i+1]
    for letter in after:
        p *= PWM[letter_to_index[letter], 0]
    return p


def nck(c, k, Z):
    """
    expected number appearances of character c in position k (based on Z)
    """
    nck = 0  # Expected number of character c in position k
    for i in range(len(sequences)):
        for j in range(M):
            s = sequences[i][j:j+W]
            if s[k] == c:
                nck += Z[i, j]

    return nck


def pck(nck):
    """
    probability of character c in position k (based on Z)
    """
    return (nck + PSEUDO_COUNTS) / (len(alphabet) * PSEUDO_COUNTS + len(sequences))


def expectation(PWM):
    Z = np.empty((len(sequences), M))
    for i in range(len(sequences)):
        for j in range(M):
            Z[i, j] = P(i, j, PWM)

    # Normalize Z
    for i in range(len(Z)):
        Z[i] /= np.sum(Z[i])

    if DEBUG:
        print("=== Z ===")
        print_matrix(Z)

    return Z


def maximization(Z):
    PWM = np.empty((len(alphabet), W + 1))
    for i, c in enumerate(alphabet):
        probabilities = [nck(c, k, Z) for k in range(W)]

        background = letter_counts[letter_to_index[c]] - sum(probabilities)
        probabilities.insert(0, background)

        probabilities = [pck(n) for n in probabilities]

        PWM[i] = probabilities

    # Normalize background
    PWM[:, 0] /= np.sum(PWM, axis=0)[0]

    if DEBUG:
        print("=== PWM ===")
        print_matrix(PWM)

    return PWM


def log_likelihood(pwm):
    total_likelihood = 0
    for i in range(len(sequences)):
        likelihood = 0
        for j in range(M):
            probability = P(i, j, pwm)
            likelihood += probability
        likelihood /= M
        total_likelihood += math.log(likelihood)
    return total_likelihood


def find_candidate_motif(motifs):
    # Choose which motif to search with
    best = -np.Infinity, "XXX", None, None
    for i, motif in enumerate(motifs):
        # Initializing
        PWM = init_pwm(motif)

        # Run one iteration
        Z = expectation(PWM)
        PWM = maximization(Z)
        likelihood = log_likelihood(PWM)

        if likelihood > best[0]:
            best = likelihood, Z, PWM

        print(f"{i}/{len(motifs)} analysed")

    _, Z, PWM = best
    return Z, PWM


def EM(Z, PWM):
    likelihood = 0
    for i in range(MAX_ITERATIONS):
        Z = expectation(PWM)
        PWM = maximization(Z)

        new_likelihood = log_likelihood(PWM)
        if abs(new_likelihood - likelihood) < CUTOFF:
            print("Converged at iteration", i + 1, new_likelihood)
            break
        likelihood = new_likelihood

    return Z, PWM


def parse_args(argv):
    if len(argv) == 1:
        return None

    if len(argv) != 3:
        print(f"Usage: python3 {argv[1]} <filename> <motif_width>")
        exit(1)

    _, filename, motif_width = argv

    global W
    W = int(motif_width)

    return filename


def parse_fasta(filename):
    with open(filename, 'r') as f:
        entries = f.read().split('>')[1:]

    if entries is None:
        print("Error parsing fasta file", filename)
        exit(1)

    sequences = []
    for entry in entries:
        meta, *sequence = entry.splitlines()
        sequence = "".join(sequence).upper()

        # Check for letters outside of alphabet
        if any(letter not in alphabet for letter in sequence):
            print(f"Unable to parse sequence {sequence}")
            continue

        sequences.append(sequence)

    return sequences


def to_jaspar(filename, name, PWM):
    with open(filename, "w") as f:
        f.write(f">{name}\n")
        for letter, row in zip(alphabet, PWM):
            f.write(letter)
            f.write(" [")
            for x in row[1:]:
                count = f"{x * len(sequences):.2f}"
                f.write(f"{count:>8}")
            f.write(" ]\n")


def visualize_starting_positions():
    print("=== Alignment ===")
    for i, row in enumerate(Z):
        start = np.argmax(row)
        for j in range(M + W - 1):
            if j == start:
                print("[", end='')
            elif j == start+W:
                print("]", end='')
            elif j > start and j < start+W:
                print("=", end='')
            else:
                print("-", end='')
        print()


if __name__ == "__main__":
    filename = parse_args(sys.argv)

    if filename is not None:
        sequences = parse_fasta(filename)
        M = min([len(sequence) for sequence in sequences]) - W + 1

    letter_counts = count_letters()

    if DEBUG:
        for i, sequence in enumerate(sequences):
            if i > 100:
                break
            print(f"{i:2} {sequence}")
            print(f"{L=} {W=} {M=}")
        for letter in alphabet:
            print(letter, letter_counts[letter_to_index[letter]])

    motifs = find_starting_motifs()
    Z, PWM = find_candidate_motif(motifs)
    Z, PWM = EM(Z, PWM)
    visualize_starting_positions()

    print("\n=== PWM ===")
    print_matrix(PWM)

    print("\n=== MOTIF ===")
    motif = [alphabet[i] for i in np.argmax(PWM, axis=0)]
    print("".join(motif[1:]))

    print("\nWriting to jaspar...")
    to_jaspar("output.jaspar", "MOTIF", PWM)
