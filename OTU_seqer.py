"""
OTU_seqer (canonically pronounced seeker)
Written by JG
v1 3/20/17 for oligotyping
Updated and uploaded to github: 3/12/18 

Outline: During vsearch OTU clustering, sequences are clustered into OTUs and assigned taxonomy. 
Information about microdiversity within OTUs or phylogentically named clusters is lost in favor of parsimony. 
This script retrieves the full sequence data associated with sequence names output from otu_demapper and outputs them to text.

Inputs:
    1. Directory - Directory containing text files with rows as sequence names
    2. combined_seqs - A fasta file to query for sequences whose headers match names listed in files in the cirectory
    3. out_dir - An output directory to store fasta files for individual OTUs.

"""

# Imports
import sys
from multiprocessing import Pool
import os
import itertools
import io
import argparse

# Functions:
def get_line_no(file_path): 
    fp = io.open(file_path,'r', errors='ignore')
    lines = 0
    while fp.readline():
        lines+=1
    return lines

def make_otu_list(otu_files): #
    otu_list = []
    with open(otu_files, 'r') as f:
        for line in f:
            otu_list.append(line.strip())
    return otu_list


def detect_formatting(fasta):
    """ This script checks a fasta file to determine how sequences are formatted.
    Typically, lines are either 80 characters long or the entire sequence is on one line.
    If lines/seq is an integer, this value is used to parse the fasta.
    If lines/seq is variable, a slower method of writing lines using an if comparison is used instead."""
    with open(fasta, 'r') as f:
        headers = 0
        lines = 0
        for line in f:
            
            if line.startswith('>'):
                header_count +=1

    with open(fasta, 'r') as f:
        headers = 0
        lines = 0
        for line in f:
            if line.startswith('>'):
                headers +=1
            lines +=1
        lines_per_seq = (lines / (headers -1))
        if lines_per_seq.is_integer():
            return lines_per_seq
        else:
            print ('Sequence line length does not appear to be consistent. lines/seq = {0}, consider passing var_length=True to search_db.'.format(lines_per_seq))
            return None
            
# any other functions I can split this into?
def search_db(fasta_db, otu_file, out_dir):
    """
    Input:
        fasta_db: A fasta file of 16S sequences typically labeled "combined_seqs.fa" or similar
        key: An OTU targeted for further classification.
        val: A list of sequence identifiers that are classified as OTU labeled "key"
    Overview:
        Pass the results of the ID_collector function to this along with the path to a combinedseqs fasta file.
        This will iterate through fastadb looking for hits. 
        
    This function holds its place in the fasta file after each hit, expecting sequences to be in order. 
    Occasionally order may not agree between the OTU map and the combined seqs file, so it will go through the entire fasta 2X before going to the next OTU.
    
    Output: A results file named after the OTU with all of the sequences collected.

    """
    print('Starting OTU file {0}'.format(otu_file))
    with open(otu_file, 'r', errors='ignore') as f:
        data = f.readlines()
        name = otu_file.split('/')[-1]    # Split filename to create output filename.
        file_out = out_dir+name
        line_max = get_line_no(otu_file)

        lines_ind = 0
        loops = 0
        
        with open(fasta_db, 'r', errors='ignore') as g:
        # Open the new output file for writing.
            with open(file_out, 'w') as h:
                
                while lines_ind < line_max:
                    line = data[lines_ind]
                    line = '>'+line.strip()
                    lines_ind+=1

                    for seq in g:
                        if seq.startswith(">"):                    
                            seq_header = seq.split()[0]     # Just take first part so it is comparable to above.
                            if seq_header == line:               # If match, write header and lines until next header.
                                h.write(line+'\n')
                                while(True): # It took way too long to realize I can't iterate & compare againt next(g) because it advances twice.
                                    b = next(g)
                                    if b.startswith('>'):
                                        break
                                    else:
                                        h.write(b)                             
                                break
                    # Only enter this else loop if sequene name isn't found.
                    # Basically on first full loop, restart sequence and on 2nd
                    # loop, just skip completely.
                    else:
                        g.seek(0)
                        if loops == 0:
                            lines_ind = lines_ind - 1
                            loops =1
                        else:
                            loops = 0
    print('File {0} Done!'.format(otu_file))
    return

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--directory", help="Directory containing text files with rows as sequence names")
    parser.add_argument("--seqs", help="A fasta file to query for sequences whose headers match names listed in files in the directory")
    parser.add_argument("--outdir", help="An output directory to store fasta files for individual OTUs.")

    try:
        result = parser.parse_args()
        return result
    except Exception as e: 
        parser.print_help()
        print(e)
        sys.exit(0)

    
if __name__ == '__main__':
    # passed arguments
    args = parse_arguments()
    directory = args.directory
    combined_seqs = args.seqs
    out_dir = args.outdir
    

    print('\nOTU Seeker Starting!')
    print('Version 1.1, updated 3/12/18\n')

    # Set up argument lists to pass to subfunction.
    paths = [os.path.join(directory,fn) for fn in next(os.walk(directory))[2]]
    paths = [x for x in paths if '.txt' in x]
    mapped_inputs = list(itertools.product([combined_seqs], paths, [out_dir]))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

# Parallel process OTUs    
    with Pool() as pool:
        pool.starmap(search_db, mapped_inputs)
