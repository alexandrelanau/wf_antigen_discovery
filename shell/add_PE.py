#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import time

def main():
    parser = argparse.ArgumentParser(description='Add PE=2 tag to the FASTA database for group variable FDR calculation in MS Fragger')
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output FASTA file')
    parser.add_argument('--verbose', action='store_true', default=False, help='Verbose output')
    args = parser.parse_args()

    start_time = time.time()

    #Add PE=2 to the neotranscript database and add the _NEOTRANSCRIPT tag to facilitate peptide class assignment in the Frapipe results downstream analysis
    #PE stands for protein evidence. In this pipeline we assign PE=1 to the uniprot proteins in _data/uniprot and PE=2 to the neotranscript proteins.
    #This is done to then use group FDR calculation in Fragpipe with group tag being protein evidence from fasta file
    #This way FDR is estimated separately for canonical and non-canonical peptides. This is more stringent, but ensures there is no FDR underestimation 
    with open(args.output, 'w') as new_db_fa:
        for record in SeqIO.parse(args.input, 'fasta'):
            new_descr = 'PE=2'
            new_record = SeqRecord(record.seq, id=f'{record.id}_NEOTRANSCRIPT', description=new_descr)
            SeqIO.write(new_record, new_db_fa, 'fasta')
    
    # If verbose, print information
    if args.verbose:
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()