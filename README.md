# OTU_unmapper
Two short scripts for extracting raw amplicon sequences associated with specific OTUs.

OTU_mapper:
After OTU clustering, vsearch creates an OTU map that assigns every sequence in the initial 
dataset to the created OTUs. This file type is documented here: 
  https://drive5.com/usearch/manual/opt_uc.html

This script iterates through a .uc file to identify sequences associated with specific OTUs and writes the headers to files.
Afterwards, the sequences can be extracted from a fasta file based on the header, using the OTU_seqer.py script


OTU_seqer:
This script retrieves the full sequence data associated with sequence names output from otu_mapper and outputs them to text files.

See documentation and arguments at the top of scripts for more details.



