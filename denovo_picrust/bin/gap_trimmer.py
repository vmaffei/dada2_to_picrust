#!/usr/bin/env python
import argparse
import re
from collections import Counter
import operator

__author__ = "Gene Blanchard"
__credits__ = ["Vince Maffei", "Gene Blanchard"]

'''
Trim the gaps semi-intelligently from greengenes alignments
'''

def gap_trimmer(infa, outfa, prefix, verbose):
    fwd_counts = Counter()
    rev_counts = Counter()
    with open(infa, 'r') as fa_h:
        for line in fa_h:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line
                if line.startswith(">{}".format(prefix)):
                    NODE = True
                else:
                    NODE = False
            else:
                if NODE:
                    fwd_count = len(re.split('[A,C,G,T]', line, maxsplit=1)[0])
                    rev_count = len(re.split('[A,C,G,T]', line[::-1], maxsplit=1)[0])
                    fwd_counts[len(re.split('[A,C,G,T]', line, maxsplit=1)[0])] += 1
                    rev_counts[len(re.split('[A,C,G,T]', line[::-1], maxsplit=1)[0])] += 1
    fwd_trim = max(fwd_counts.iteritems(), key=operator.itemgetter(1))[0]
    rev_trim = max(rev_counts.iteritems(), key=operator.itemgetter(1))[0]
    if verbose:
        print "Forward Counts:\n\tLength\tCount\n"
        for key in sorted(fwd_counts.keys()):
            print "\t{}\t{}\n".format(key, fwd_counts[key])
        print "Forward Trim:\t{}\n".format(fwd_trim) 
        print "Reverse Counts:\n\tLength\tCount\n"
        for key in sorted(rev_counts.keys()):
            print "\t{}\t{}\n".format(key, rev_counts[key])
        print "Reverse Trim:\t{}\n".format(rev_trim) 
    with open(infa, 'r') as fa_h, open(outfa, 'w') as out_h:
        for line in fa_h:
            line = line.rstrip('\n')
            if line.startswith('>'):
                out_h.write(line+"\n")
            else:
                line = line[fwd_trim::][::-1][rev_trim::][::-1]
                out_h.write(line+"\n")

def main():
    # Argument Parser
    parser = argparse.ArgumentParser(description='<This is what the script does>')

    # Input file
    parser.add_argument('-i', '--input', dest='input', required=True, help='The input file')
    # Output file
    parser.add_argument('-o', '--output', dest='output', required=True, help='The output file')
    # Prefix
    parser.add_argument('-p', '--prefix', dest='prefix', default="node", help='The prefix on all of your sequences, default is node')
    # Verbose
    parser.add_argument('-v', '--verbose', dest='verbose', action="store_true", help='Print helpful diagnostic information')
    
    # Parse arguments
    args = parser.parse_args()
    infile = args.input
    outfile = args.output
    prefix = args.prefix
    verbose = args.verbose

    gap_trimmer(infile, outfile, prefix, verbose)

if __name__ == '__main__':
    main()
