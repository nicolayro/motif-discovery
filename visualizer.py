from Bio import motifs
import sys


if __name__ == "__main__":
    argv = sys.argv
    if len(argv) != 2:
        print(f"Usage: python {argv[0]} <filename>")
        exit(1)

    filename = argv[1]

    with open(filename, "r") as handle:
        m = motifs.read(handle, "jaspar")
        m.weblogo("weblogo.png")
