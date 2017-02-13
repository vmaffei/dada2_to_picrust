#!/usr/bin/env python
import argparse
import pandas
from Bio import SeqIO

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

__authors__ = "Gene Blanchard, Vince Maffei"

'''
dada2_to_picrust conversion using @vmaffei's workflow
'''

def subset_img(seqtab, ko, ref, img):
    img_ko = pd.DataFrame.from_csv(ko, sep='\t')
    id_img = pd.DataFrame.from_csv(img, sep='\t')\
    id_img_sub = id_img.loc[id_img["img_genome_id"].isin(list(set(id_img["img_genome_id"]).intersection(img_ko.index)))]
    ref_dict = SeqIO.to_dict(SeqIO.parse(ref, "fasta"))
    # Subset the gg refrence  fastq
    ref_sub = {}
    with open("gg_13_5_study_db.fasta", "w") as handle:
        for key in map(str, list(id_img_sub.index)):
            try:
                SeqIO.write(ref_dict[key], handle, "fasta")
            except KeyError:
                pass
    
    

def main():
    # Argument Parser
    parser = argparse.ArgumentParser(description='Get your dada2 data ready for picrust')

    # Seqtab
    parser.add_argument('-s', '--seqtab', dest='seqtab', required=True, help='dada2 seqtab.nochim rds file')
    # IMG_ko_counts.tab
    parser.add_argument('-k', '--ko', dest='ko', required=True, help='File: IMG_ko_counts.tab, see:https://github.com/picrust/picrust/blob/master/tutorials/picrust_starting_files.zip')
    # GG ref
    parser.add_argument('-r', '--ref', dest='ref', required=True, help='File: gg_13_5.fasta, see http://greengenes.secondgenome.com/downloads/database/13_5')
    # GG img 
    parser.add_argument('-i', '--img', dest='img', required=True, help='File: gg_13_5_img.txt, see http://greengenes.secondgenome.com/downloads/database/13_5')
    
    # Output file
    parser.add_argument('-o', '--output', dest='output', help='The output file')

    # Parse arguments
    args = parser.parse_args()
    seqtab = args.seqtab
    # ko = args.ko
    # ref = args.ref
    # img = args.img
    ko = "data/IMG_ko_counts.tab"
    ref = "data/gg_13_5.fasta"
    img = "data/gg_13_5_img.txt"

    # Do stuff here

if __name__ == '__main__':
    main()
