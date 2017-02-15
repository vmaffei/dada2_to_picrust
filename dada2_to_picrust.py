#!/usr/bin/env python
import argparse
import qiime_default_reference as qdr
import pandas
from Bio import SeqIO
from skbio.parse.sequences import parse_fasta
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import DNASequence
from pynast.util import pynast_seqs, pairwise_alignment_methods
from pynast.logger import NastLogger

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

def pynasty(fasta, template_file, min_pct=0.75, min_len=90):
    # Prepare template
    with open(template_file, 'r') as template_h:
        template = [(seqid, seq.replace('.', '-').upper()) for seqid, seq in parse_fasta(template_h)]
    template = Alignment.from_fasta_records(template, DNASequence, validate=True)
    # Prepare logger
    logger = NastLogger("pynast.log")
    # Align
    # Parse fasta
    with open(fasta, 'r') as fasta_h:
        inseqs = parse_fasta(fasta_h)
        pynast_aligned, pynast_failed = pynast_seqs(inseqs, template, min_pct=min_pct, min_len=min_len, align_unaligned_seqs_f='blast', logger=logger)
    for i, seq in enumerate(pynast_failed):
        seq_record = DNASequence(str(seq), id=seq.Name)
        pynast_failed[i] = seq_record
    pynast_failed = SequenceCollection(pynast_failed)

    for i, seq in enumerate(pynast_aligned):
        seq_record = DNASequence(str(seq), id=seq.Name)
        pynast_aligned[i] = seq_record
    pynast_aligned = Alignment(pynast_aligned)

    with open("pynast_failed.fa", 'w') as failed_h:
        failed_h.write(pynast_failed.to_fasta())
    with open("pynast_aligned.fa", 'w') as pynast_h:
        pynast_h.write(pynast_aligned.to_fasta())

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
    format_seqtab_to_repset_table(seqtab)
    

if __name__ == '__main__':
    main()
