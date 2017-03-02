#!/usr/bin/env python
import argparse
import qiime_default_reference as qdr
import re
from collections import Counter
import numpy as np
import pandas as pd
import itertools

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

def main():
    # Argument Parser
    parser = argparse.ArgumentParser(description='Get your dada2 data ready for picrust')
    # Default refrence
    def_ref = qdr.get_template_alignment()

    # Seqtab
    parser.add_argument('-s', '--seqtab', dest='seqtab', help='DADA2 seqtab.nochim.csv file')
    # IMG_ko_counts.tab
    parser.add_argument('-k', '--ko', dest='ko', default="data/IMG_ko_counts.tab", help='File: IMG_ko_counts.tab, see:https://github.com/picrust/picrust/blob/master/tutorials/picrust_starting_files.zip')
    # GG ref
    parser.add_argument('-r', '--ref', dest='ref', default="data/gg_13_5.fasta", help='File: gg_13_5.fasta, see http://greengenes.secondgenome.com/downloads/database/13_5')
    # GG img
    parser.add_argument('-i', '--img', dest='img', default="data/gg_13_5_img.txt", help='File: gg_13_5_img.txt, see http://greengenes.secondgenome.com/downloads/database/13_5')
    # Template
    parser.add_argument('-t', '--template', dest='template', default=def_ref, help='Pynast Alignmet Template - Default:85_otus.pynast.fasta')
    # Output file
    parser.add_argument('-o', '--output', dest='output', default="dada2_to_picrust", help='The output basename')

    # Parse arguments
    args = parser.parse_args()
    seqtab = args.seqtab
    # ko = args.ko
    ko = "data/IMG_ko_counts.tab"
    # ref = args.ref
    ref = "data/gg_13_5.fasta"
    # img = args.img
    img = "data/gg_13_5_img.txt"
    template = args.template
    output = args.output
    workflow = "{}_workflow.sh".format(output)
    

    # Qiime files from seqtab csv
    seqtab_dict = format_seqtab_to_repset_table(seqtab)

    # Create a refrence dictionary
    ref_dict = subset_img(ko, ref, img)

    # Combine dictionaries
    # Check duplicate keys
    dada_fasta = "{}.fasta".format(output)
    if set(seqtab_dict.keys()).isdisjoint(ref_dict.keys()):
        with open(dada_fasta, 'w') as fasta_h:
            for key in ref_dict:
                fasta_h.write(">{}\n{}\n".format(key, ref_dict[key]))
            for key in seqtab_dict:
                fasta_h.write(">{}\n{}\n".format(key, seqtab_dict[key]))
    else:
        print "Somehow your samples have the same names as GG ids"
    with open(workflow, 'w') as work_h:
        # VSearch
        work_h.write("vsearch --usearch_global {} -db data/gg_13_5.fasta --matched {}.filtered.fasta --id 0.5\n".format(dada_fasta, output))
        # Filter Vsearch failures and make biom
        work_h.write("head -n1 {0}.otutable > {0}.filtered.otutable\n".format(seqtab))
        work_h.write("cat {}.filtered.fasta | grep '^>node' | cut -c 2- | xargs -I % sh -c 'cat {}.otutable | grep -P \"^%\\t\"' >> {}.filtered.otutable\n".format(output, seqtab, seqtab))
        work_h.write("biom convert -i {0}.filtered.otutable --to-json -o {0}.filtered.biom --table-type 'OTU table'\n".format(seqtab))
        # Pynast
        work_h.write("pynast -i {}.filtered.fasta -l 90 -p 0.1 -t {}\n".format(output, template))
        # Gap Trim
        work_h.write("gap_trimmer.py -i {0}.filtered_pynast_aligned.fasta -o {0}.filtered_pynast_aligned.trimmed.fasta -v\n".format(output))
        # FastTree
        work_h.write("FastTree -nt -gamma -fastest -no2nd -spr 4 {0}.filtered_pynast_aligned.trimmed.fasta > {0}_pynast_aligned.trimmed.tree\n".format(output))
        # Format Tree and Trait Table
        work_h.write("format_tree_and_trait_table.py -t {0}_pynast_aligned.trimmed.tree -i data/IMG_16S_counts.tab -m data/gg_13_5_img.txt -o {0}_genome_prediction_16S\n".format(output))
        work_h.write("format_tree_and_trait_table.py -t {0}_pynast_aligned.trimmed.tree -i data/IMG_ko_counts.tab -m data/gg_13_5_img.txt -o {0}_genome_prediction_kegg\n".format(output))
        # ASR
        work_h.write("ancestral_state_reconstruction.py -i {0}_genome_prediction_16S/trait_table.tab -t {0}_genome_prediction_16S/pruned_tree.newick -o {0}_16S_asr_counts.tab\n".format(output))
        work_h.write("ancestral_state_reconstruction.py -i {0}_genome_prediction_kegg/trait_table.tab -t {0}_genome_prediction_kegg/pruned_tree.newick -o {0}_kegg_asr_counts.tab\n".format(output))
        # Predict Traits
        work_h.write("predict_traits.py -i {0}_genome_prediction_16S/trait_table.tab -t {0}_genome_prediction_16S/reference_tree.newick -r {0}_16S_asr_counts.tab -o {0}_genome_prediction_16S/16S_precalculated.tab\n".format(output))
        work_h.write("predict_traits.py -i {0}_genome_prediction_kegg/trait_table.tab -t {0}_genome_prediction_kegg/reference_tree.newick -r {0}_kegg_asr_counts.tab -o {0}_genome_prediction_kegg/kegg_precalculated.tab\n".format(output))
        # Normalize by copy number
        work_h.write("normalize_by_copy_number.py -i {0}.filtered.biom -c {1}_genome_prediction_16S/16S_precalculated.tab -o {1}_normalized_counts.biom\n".format(seqtab, output))
        # Predict Metagenomes
        work_h.write("predict_metagenomes.py -i {0}_normalized_counts.biom -o {0}_meta_counts.biom -c {0}_genome_prediction_kegg/kegg_precalculated.tab\n".format(output))

if __name__ == '__main__':
    main()
