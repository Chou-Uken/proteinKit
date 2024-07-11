import argparse
from biofkit.seqKit.convKit import pairwiseDnaAlign


def main():
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Output pairwise alignment.')
    groupA = parser.add_argument_group()
    groupB = parser.add_argument_group()
    groupA.add_argument("-1", "--seq1", help="the first sequence", required=False)
    groupA.add_argument("-2", "--seq2", help="the second sequence", required=False)
    groupB.add_argument("-f", "--file", help="fasta file with at least 2 sequences", required=False)
    parser.add_argument()
    args = parser.parse_args()
    if ((args.seq1 and args.seq2) and (not args.file)):
        pairwiseDnaAlign()
    elif (args.file and not (args.seq1 or args.seq2)):
        pairwiseDnaAlign(fasta=args.file, )
    else:
        print("sad")


if __name__ == "__main__":
    main()
