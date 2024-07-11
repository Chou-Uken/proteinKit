import argparse
from biofkit.seqKit.convKit import readFasta


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Output the sequences in a fasta file.')
    parser.add_argument("fastaFile", help='Fasta file with at least one sequence.')
    args = parser.parse_args()
    contents: OrderedDict = readFasta(args.fastaFile)
    for num, (key, value) in enumerate(contents.items()):
        print("{0}, {1}:{2}".format(num, key, value))

if __name__ == "__main__":
    main()
