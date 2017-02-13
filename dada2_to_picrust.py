#!/usr/bin/env python
import argparse
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

__authors__ = "Gene Blanchard, Vince Maffei"

'''
dada2_to_picrust conversion using @vmaffei's workflow
'''
def R_code(seqtab, ko, ref, img):
    ShortRead = importr("ShortRead")
    Biom = importr("biom")
    DADA = importr("dada2")
    Base = importr("base")
    ref = ShortRead.readFasta
    

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
    ko = args.ko
    ref = args.ref
    img = args.img

    # Do stuff here

if __name__ == '__main__':
    main()
