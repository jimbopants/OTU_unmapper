#!bin/bash
#Written by JG 3/16/17

""" ######## Documentation: ##########
Background: After OTU clustering, vsearch creates an OTU map that assigns every sequence in the initial 
dataset to the created OTUs. This file type is documented here: 
    https://drive5.com/usearch/manual/opt_uc.html

Objective:
This script iterates through a uc file to identify sequences associated with specific OTUs and writes the headers to files.
Afterwards, the sequences can be extracted from a fasta file based on the header, using the fasta_

Arguments (passed as an ordered list at run time):
1. OTU list - Text file with OTUID as header and 1 OTU per row
2. OTU_map - the OTU_map.uc created by vsearch used for a particular project
3. Output directory - directory to store OTU sequence lists as text files. Will be created if it does not exist.

"""

# Imports:
import sys
import os
import pandas as pd
from multiprocessing import Pool


# Function definitions
def otu_mapper(otu_map, otu_of_interest, seq_dir):
    with open(otu_map, 'r') as f:
        output_name = seq_dir + otu_of_interest.split(';')[0]+'.txt'
                
        with open( output_name, 'w') as g:
            for line in f:
                otu =  line.split()[-1]
                seq = line.split()[-2]
                if otu == otu_of_interest:
                    g.write(seq+'\n')
    return


if __name__ == '__main__':

# Command line arguments
    OTU_list = sys.argv[1] # - 
    otu_map = sys.argv[2]
    seq_dir = sys.argv[3] # -
    
# Make directory if it doesn't exist:
    if not os.path.exists(seq_dir):
        os.makedirs(seq_dir)

    OTU_df = pd.read_csv(OTU_list, sep='\t', header=0)
    OTU_names = OTU_df['OTUID'].values.tolist()
    pool_args = [(otu_map, i, seq_dir) for i in OTU_names]
    with Pool(processes=4) as p:
        p.starmap(otu_mapper, pool_args)
