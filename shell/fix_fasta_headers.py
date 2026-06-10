import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Fix fasta headers for ORFfinder')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    args = parser.parse_args()

    input_fasta = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    df_input_fasta = pd.DataFrame.from_dict({'decription': input_fasta.keys(), 'sequence': input_fasta.values()})
    with open(args.output, 'w') as new_db_fa:
        for index, row in df_input_fasta.iterrows():
            new_descr = f'{row['decription'].replace('.', '_').replace('|', '_').replace(' ', '_')}'
            new_record = SeqRecord(row['sequence'].seq, id='', description=new_descr)
            SeqIO.write(new_record, new_db_fa, 'fasta-2line')





if __name__ == "__main__":
    main()