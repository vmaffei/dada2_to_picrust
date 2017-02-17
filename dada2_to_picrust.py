#!/usr/bin/env python
import argparse
import qiime_default_reference as qdr
import numpy as np
import pandas as pd
import itertools
from skbio.parse.sequences import parse_fasta
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import DNASequence
from pynast.util import pynast_seqs, pairwise_alignment_methods
from pynast.logger import NastLogger
import subprocess

__authors__ = "Gene Blanchard, Vince Maffei"

'''
dada2_to_picrust conversion using @vmaffei's workflow
'''

def remove_quotes(line):
    line = line.replace('"', '')
    return line


def dict_to_repset(repsetfile, otudict):
    with open(repsetfile, 'w') as outhandle:
        for key in otudict:
            outhandle.write(">{}\n{}\n".format(key, otudict[key]))


def format_seqtab_to_repset_table(seqtab):
    otutable = seqtab + ".otutable"
    rep_set = seqtab + ".repset"
    parsed_lines = []
    with open(seqtab) as seqtab_handle, open(otutable, 'w') as otuhandle:
        for index, line in enumerate(seqtab_handle):
            if index == 0:
                otus = map(remove_quotes, line.rstrip('\n').split(','))[1::]
                otu_dict = {"node{}".format(n): seq for n, seq in enumerate(otus)}
                dict_to_repset(rep_set, otu_dict)
                parsed_lines.append(["#OTU ID"]+["node"+n for n in (map(str, range(0, len(otus))))])
            else:
                parsed_lines.append(remove_quotes(line).rstrip().split(','))
        otuhandle.write('\n'.join(['\t'.join(line) for line in list(np.array(parsed_lines).T)]))
    return otu_dict
    


def subset_img(ko, ref, img):
    # Subset IMG & KO
    img_ko = pd.DataFrame.from_csv(ko, sep='\t')
    id_img = pd.DataFrame.from_csv(img, sep='\t')
    id_img_sub = id_img.loc[id_img["img_genome_id"].isin(list(set(id_img["img_genome_id"]).intersection(img_ko.index)))]

    # Create refrence dict
    ref_dict = {}
    key_list = map(str, list(id_img_sub.index))
    with open(ref, 'r') as fasta_h:
        for key, seq in itertools.izip_longest(*[fasta_h] * 2):
            key = key.lstrip('>').rstrip('\n')
            seq = seq.rstrip('\n')
            if key in  key_list:
                ref_dict[key] = seq
    return ref_dict


def pynasty(fasta, min_pct=0.75, min_len=90):
    # load template sequences
    template_alignment = []
    template_alignment_fp = qdr.get_template_alignment()
    for seq_id, seq in parse_fasta(open(template_alignment_fp)):
        # replace '.' characters with '-' characters
        template_alignment.append((seq_id, seq.replace('.', '-').upper()))
    template_alignment = Alignment.from_fasta_records(template_alignment, DNASequence, validate=True)
    # load candidate sequences
    with open(fasta, 'U') as seq_file:
        candidate_sequences = parse_fasta(seq_file)
        pynast_aligned, pynast_failed = pynast_seqs(candidate_sequences, template_alignment, min_pct=min_pct, min_len=min_len)

    failed = []
    for i, seq in enumerate(pynast_failed):
        skb_seq = DNASequence(str(seq), id=seq.Name)
        failed.append(skb_seq)
    pynast_failed = SequenceCollection(failed)
    aligned = []
    for i, seq in enumerate(pynast_aligned):
        skb_seq = DNASequence(str(seq), id=seq.Name)
        aligned.append(skb_seq)
    pynast_aligned = Alignment(aligned)

    with open("pynast_failed.fa", 'w') as failed_h:
        failed_h.write(pynast_failed.to_fasta())

    with open("pynast_aligned.fa", 'w') as pynast_h:
        pynast_h.write(pynast_aligned.to_fasta())


def fasttree():
    args = "./bin/FastTree -nt -gamma -fastest -no2nd -spr 4 pynast_aligned.fa".split(' ')
    with open("gg_13_5_study_db.tree", 'w') as tree:
        popen = subprocess.Popen(args, stdout=tree)
        popen.wait()


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
    # Template
    parser.add_argument('-t', '--template', dest='template', default=qdr.get_template_alignmet(), help='Pynast Alignmet Template - Default:85_otus.pynast.fasta')
    
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

    # Qiime files from seqtab
    seqtab_dict = format_seqtab_to_repset_table(seqtab)
    ref_dict = subset_img(ko, ref, img)
    
    # Combine dictionaries
    # Check duplicate keys
    if set(seqtab_dict.keys()).isdisjoint(ref_dict.keys()):
        with open("gg_13_5_dada_db.fasta", 'w') as fasta_h:
            for key in ref_dict:
                fasta_h.write(">{}\n{}\n".format(key, ref_dict[key]))
            for key in seqtab_dict:
                fasta_h.write(">{}\n{}\n".format(key, seqtab_dict[key]))
    else:
        print "Somehow your samples have the same names as GG ids"
    
    
    
    
    
    

if __name__ == '__main__':
    main()
