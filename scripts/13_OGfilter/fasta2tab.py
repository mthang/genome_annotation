#!/usr/bin/env python3
import os
import re
import sys



#directory="shell_genes"
#output_filepath = 'shell_genes.tsv'
directory=sys.argv[1]
output_filepath=sys.argv[1] + '.tsv'
#print(directory)
#print(output_filepath)
#quit()
output_file = open(output_filepath,'w+')
header_list = ['Ortholog_Filename','Taxonomy','Gene','Sequence_Name','Description','Preferred_Name','Protein_Sequence']
header = "\t".join(header_list)+"\n"
output_file.write(header)


def parse_fasta_file(input_file):
    """Return a dict of {id:gene_seq} pairs based on the sequences in the input FASTA file
    input_file -- a file handle for an input fasta file
    """
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []

    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if curr_seq_id is not None:
                parsed_seqs[curr_seq_id] = ''.join(curr_seq)

            curr_seq_id = line[1:]
            curr_seq = []
            continue

        curr_seq.append(line)

    #Add the final sequence to the dict
    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
    return parsed_seqs

#code reference : https://colab.research.google.com/github/zaneveld/full_spectrum_bioinformatics/blob/master/content/06_biological_sequences/reading_and_writing_fasta_files.ipynb#scrollTo=cnjzDxUQ-Qy1

for filename in os.listdir(directory):
    ogfile = os.path.join(directory,filename)
    if os.path.isfile(ogfile):
        print(ogfile)
#Normally this would be determined
#by user input via argparse. Hard-coded for now
#input_file = 'OG0050310.fa'
#input_file = 'OG0019577.fa'

        f = open(ogfile)
        parsed_seqs = parse_fasta_file(ogfile)
        print(parsed_seqs)


#        output_filepath = 'cloud_genes.txt'
#        output_file = open(output_filepath,'w+')
        for seq_id,seq in parsed_seqs.items():
            first_split=re.split('\\|',seq_id)
            second_split=re.split('^(\d+)\.',first_split[0])
            join_list=second_split + first_split[1:]
            output_line = "\t".join([filename,"\t".join(join_list),seq])+"\n"
            print(output_line)
            output_file.write(output_line)
