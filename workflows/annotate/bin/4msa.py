#!/usr/bin/env python


'''
This script takes the alignment result from "blasting" query genomes against
the reference using mmseqs (aln.m8) and turns this result into a fasta file for
MSA using eg mafft. Sequences are translated as well.

Notes:

TTG can be a start codon

> TTG serves as an initiation codon for the ribosomal protein MvaS7 from the archaeon Methanococcus vannielii. -- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC177430/
'''


import argparse
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
import pandas as pd


def translate(x):
    '''
    http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec25
    '''
    return Seq(x).translate().__str__().replace('*', '')


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--aln', required=True, help='Alignment result')
parser.add_argument(
    '--outdir', default='cluster', help='Output directory')
args = parser.parse_args()


# wc -l work/5a/4140d6ff544b24ae7d7523a3352013/aln.m8
# 2103
# fp = 'work/07/8c4ab2899f8dc93928e0929cc3e4c0/aln.m8'
# outdir = 'foo'
# !mkdir foo
outdir = Path(args.outdir)
outdir.mkdir()

'''
> By default (--format-mode 0), alnRes.tab will contain alignment result in a BLAST tabular result (comparable to -m 8 -outfmt 6) with 12 columns: (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score. -- https://github.com/soedinglab/mmseqs2/wiki#alignment-format
'''
# header = 'query target seq_id aln_len n_mismatches n_gaps query_start query_end target_start target_end e_value bit_score'.split(' ')
header = 'query,target,qseq,tseq'.split(',')
df = pd.read_csv(args.aln, sep='\t', names=header)


matches = defaultdict(list)
seqs = {}
for i in df.itertuples():
    matches[i.target].append(i.query)
    seqs[i.query] = i.qseq
    seqs[i.target] = i.tseq


# For each cluster, export fna and faa sequences and remove stop symbol *
for r, qs in matches.items():
    with open(f'{outdir}/{r}.faa', 'w+') as out_faa, open(f'{outdir}/{r}.fna', 'w+') as out_fna:
        r_na = seqs[r]
        r_aa = translate(r_na)
    
        # The reference sequence always comes at the top of the file.
        out_fna.write(f'>{r}\n{r_na}\n')
        out_faa.write(f'>{r}\n{r_aa}\n')
        
        # Following are the corresponding homologue proteins of query genomes.
        for q in qs:
            q_na = seqs[q]
            q_aa = translate(q_na)

            out_fna.write(f'>{q}\n{q_na}\n')
            out_faa.write(f'>{q}\n{q_aa}\n')

