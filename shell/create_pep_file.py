import argparse
import pandas as pd



def main():
    parser = argparse.ArgumentParser(description='Create a .pep file for binding prediction tools')
    parser.add_argument('--input', required=True, help='Path to input psm.tsv file (output from the MS peptidomis search engine)')
    parser.add_argument('--output', required=True, help='Path to output .pep file for binding prediction')
    args = parser.parse_args()

    psm_file = pd.read_csv(args.input, sep="\t")
    peptides = psm_file["Peptide"].unique().tolist()
    with open(args.output, 'w') as file:
        for pep in peptides:
            if 7<len(pep)<15:
                file.write(pep+'\n')




if __name__ == "__main__":
    main()