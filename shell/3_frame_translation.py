#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def translate_3_frames(input_fasta, output_fasta):
    with open(output_fasta, "w") as out_f:

        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = record.seq.upper()

            for frame in range(3):
                shifted = seq[frame:]

                # trim to full codons
                trimmed = shifted[:len(shifted) - (len(shifted) % 3)]

                protein = str(trimmed.translate(to_stop=False))

                out_f.write(f">{record.id}_frame{frame+1}\n")
                out_f.write(protein + "\n")


def main():
    parser = argparse.ArgumentParser(description="3-frame translation of FASTA sequences")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output amino acid FASTA file")

    args = parser.parse_args()

    translate_3_frames(args.input, args.output)


if __name__ == "__main__":
    main()