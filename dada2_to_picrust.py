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
from cogent.parse.tree import DndParser
from picrust.format_tree_and_trait_table import *
from picrust.util import make_output_dir, PicrustNode
from picrust.parse import parse_trait_table
from os.path import join,splitext, basename

__author__ = "Gene Blanchard"
__credits__ = ["Vince Maffei", "Gene Blanchard"]


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

def format_tree_and_traits(tree, traits, mapping, output):
    with open(tree, 'r') as tree_h:
        tree = DndParser(tree_h)
    with open(traits, 'U') as traits_h:
        traits_lines = traits_h.readlines()
    with open(mapping, 'U') as mapping_h:
        trait_to_tree_mapping = make_id_mapping_dict(parse_id_mapping_file(mapping_h))
    new_reference_tree, not_useful_trait_table_lines =\
      reformat_tree_and_trait_table(\
      tree=tree,\
      trait_table_lines = [],\
      trait_to_tree_mapping = None,\
      input_trait_table_delimiter= None,\
      output_trait_table_delimiter= None,\
      filter_table_by_tree_tips=False,\
      convert_trait_floats_to_ints=False,\
      filter_tree_by_table_entries=False,\
      convert_to_bifurcating=True,\
      add_branch_length_to_root=False,\
      name_unnamed_nodes=True,\
      min_branch_length=0.0001,\
      verbose=True)
    new_reference_tree_copy=new_reference_tree.deepcopy()
    new_tree, new_trait_table_lines = \
      reformat_tree_and_trait_table(tree=new_reference_tree_copy,\
      trait_table_lines = traits_lines,\
      trait_to_tree_mapping = trait_to_tree_mapping,\
      input_trait_table_delimiter= 'tab',\
      output_trait_table_delimiter='tab',\
      filter_table_by_tree_tips=True,\
      convert_trait_floats_to_ints=False,\
      filter_tree_by_table_entries=True,\
      convert_to_bifurcating=False,\
      add_branch_length_to_root=False,\
      name_unnamed_nodes=False,\
      min_branch_length=0.0001,\
      verbose=True)
    #Set output base file names
    trait_table_base = 'trait_table.tab'
    pruned_tree_base = 'pruned_tree.newick'
    reference_tree_base = 'reference_tree.newick'
    output_dir = make_output_dir(output,strict=False)
    basefile = splitext(basename(tree))[0]
    output_table_fp = join(output,"{}_{}".format(basefile, trait_table_base))
    output_tree_fp = join(output,"{}_{}".format(basefile, pruned_tree_base))
    output_reference_tree_fp = join(output,"{}_{}".format(basefile, reference_tree_base))
    output_trait_table_file = open(output_table_fp,"w+")
    output_tree_file  = open(output_tree_fp,"w+")
    output_reference_tree_file  = open(output_reference_tree_fp,"w+")
    
    # format 16S copy number data
    # format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i ./genome_prediction/picrust_starting_files/IMG_16S_counts.tab -m ./genome_prediction/gg_13_5_img.txt -o ./genome_prediction/format/16S/
    # format kegg IMG data
    # format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i ./genome_prediction/picrust_starting_files/IMG_ko_counts.tab -m ./genome_prediction/gg_13_5_img.txt -o ./genome_prediction/format/KEGG/
    # perform ancestral state reconstruction


def main():
    # Argument Parser
    parser = argparse.ArgumentParser(description='Get your dada2 data ready for picrust')
    # Default refrence
    def_ref = qdr.get_template_alignment()

    # Seqtab
    parser.add_argument('-s', '--seqtab', dest='seqtab', help='dada2 seqtab.nochim rds file')
    # IMG_ko_counts.tab
    parser.add_argument('-k', '--ko', dest='ko', default="data/IMG_ko_counts.tab", help='File: IMG_ko_counts.tab, see:https://github.com/picrust/picrust/blob/master/tutorials/picrust_starting_files.zip')
    # GG ref
    parser.add_argument('-r', '--ref', dest='ref', default="data/gg_13_5.fasta", help='File: gg_13_5.fasta, see http://greengenes.secondgenome.com/downloads/database/13_5')
    # GG img 
    parser.add_argument('-i', '--img', dest='img', default="data/gg_13_5_img.txt", help='File: gg_13_5_img.txt, see http://greengenes.secondgenome.com/downloads/database/13_5')
    # Template
    parser.add_argument('-t', '--template', dest='template', default=def_ref, help='Pynast Alignmet Template - Default:85_otus.pynast.fasta')
    
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
    pynasty("gg_13_5_dada_db.fasta")
    fasttree()
    format_tree_and_traits("gg_13_5_study_db.tree", "data/IMG_16S_counts.tab", "data/gg_13_5_img.txt", "genome_prediction")
    
    
    
    
    
    

if __name__ == '__main__':
    main()
