#!/usr/bin/env python


'''
Below the hyphy output is parsed and mapped to the reference sequence. This
means that all gaps in the reference sequence in the MSA will be ignored, as
will be selection values corresponding to these gaps (and which are present in
other proteins in the MSA).
'''


import argparse
import json
import pdb

import screed

from foldvis.io import load_conserved


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--msa', required=True, help='Multiple (amino acid) sequence alignment')
parser.add_argument(
    '--selection', required=True, nargs='+', help='Selection analysis from hyphy')
parser.add_argument(
    '--out', default='result.json', help='Results')
args = parser.parse_args()


with screed.open(args.msa) as file:
    # The reference sequence is always at the top of the files.
    nx = next(file)
    seq = nx.sequence
    seq_ = seq.replace('-', '')  # original state before alignment


gaps = []
for ix, i in enumerate(seq):
    if i == '-':
        gaps.append(1)
    else:
        gaps.append(0)


results = {}
results['sequence'] = seq_


for fp in args.selection:
    # For example <name>__busted
    method = fp.split('__')[-1]
    
    with open(fp, 'r') as file:
        d = json.load(file)
        
    if method == 'fubar':
        # posterior probability
        scores = [i[3:5] for i in d['MLE']['content']['0']]
        neg, pos = zip(*[p for p, gap in zip(scores, gaps) if not gap])
        assert len(seq_) == len(neg)
        e = {'positive': {}, 'negative': {}}
        e['negative']['scores'] = neg
        e['negative']['sites'] = [ix for ix, i in enumerate(neg) if i >= 0.95]
        e['positive']['scores'] = pos
        e['positive']['sites'] = [ix for ix, i in enumerate(pos) if i >= 0.95]

    elif method == 'slac':
        # p-value
        scores = [i[8:10] for i in d['MLE']['content']['0']['by-site']['AVERAGED']]
        pos, neg = zip(*[p for p, gap in zip(scores, gaps) if not gap])
        assert len(seq_) == len(neg)
        e = {'positive': {}, 'negative': {}}
        e['negative']['scores'] = neg
        e['negative']['sites'] = [ix for ix, i in enumerate(neg) if i <= 0.05]
        e['positive']['scores'] = pos
        e['positive']['sites'] = [ix for ix, i in enumerate(pos) if i <= 0.05]

    elif method == 'meme':
        # p-value
        scores = [i[6] for i in d['MLE']['content']['0']]
        pos = [p for p, gap in zip(scores, gaps) if not gap]
        assert len(seq_) == len(neg)
        e = {'positive':{}}
        e['positive']['scores'] = pos
        e['positive']['sites'] = [ix for ix, i in enumerate(pos) if i <= 0.05]

    elif method == 'busted':
        e = {'positive':{}}
        e['positive']['p'] = round(d['test results']['p-value'], 4)
        
    else:
        raise ValueError('Method not implemented, exit.')

    results[method] = e

'''
BUSTED

 "test results":{
   "LRT":1.693664598489704,
   "p-value":0.2143855005993182
  },
'''



# p = 1 - p  # "reverse p-value"

# FUBAR .. 3, 4
# MEME .. 6


# First sequence is reference
results['conserved'] = load_conserved(args.msa)


with open(args.out, 'w+') as out:
    json.dump(results, out, indent=4)


'''
import json


# https://stevenweaver.github.io/hyphy-site/getting-started/#characterizing-selective-pressures
# fubar
fp = 'ref__118122e78336cb41a7ca4c6851be70f8_27.hyphy.json'
# slac
fp = '/Users/phi/tmp/evobiochem/primates/hyphy.json'
# meme

# cd /Users/phi/tmp/evobiochem/workflow/work/d2/591f0cb577ea5d6c0c96e9eef3090e/
fp = 'ref__118122e78336cb41a7ca4c6851be70f8_27.hyphy.json'
# fel
# ...
# busted

with open(fp, 'r') as file:
    d = json.load(file)

d['MLE']['headers'] 

FUBAR

[['alpha', 'Mean posterior synonymous substitution rate at a site'],
 ['beta', 'Mean posterior non-synonymous substitution rate at a site'],
 ['beta-alpha', 'Mean posterior beta-alpha'],
 ['Prob[alpha>beta]', 'Posterior probability of negative selection at a site'],
 ['Prob[alpha<beta]', 'Posterior probability of positive selection at a site'],
 ['BayesFactor[alpha<beta]',
  'Empiricial Bayes Factor for positive selection at a site']]

SLAC

[['ES', 'Expected synonymous sites'],
 ['EN', 'Expected non-synonymous sites'],
 ['S', 'Inferred synonymous substitutions'],
 ['N', 'Inferred non-synonymous substitutions'],
 ['P[S]', 'Expected proportion of synonymous sites'],
 ['dS', 'Inferred synonymous susbsitution rate'],
 ['dN', 'Inferred non-synonymous susbsitution rate'],
 ['dN-dS', 'Scaled by the length of the tested branches'],
 ['P [dN/dS > 1]',
  'Binomial probability that S is no greater than the observed value, with P<sub>s</sub> probability of success'],
 ['P [dN/dS < 1]',
  'Binomial probability that S is no less than the observed value, with P<sub>s</sub> probability of success'],
 ['Total branch length',
  'The total length of branches contributing to inference at this site, and used to scale dN-dS']]

MEME

[['alpha;', 'Synonymous substitution rate at a site'],
 ['&beta;<sup>-</sup>',
  'Non-synonymous substitution rate at a site for the negative/neutral evolution component'],
 ['p<sup>-</sup>',
  'Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection'],
 ['&beta;<sup>+</sup>',
  'Non-synonymous substitution rate at a site for the positive/neutral evolution component'],
 ['p<sup>+</sup>',
  'Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection'],
 ['LRT',
  'Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;'],
 ['p-value',
  'Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;'],
 ['# branches under selection',
  'The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate'],
 ['Total branch length',
  'The total length of branches contributing to inference at this site, and used to scale dN-dS'],
 ['MEME LogL', 'Site Log-likelihood under the MEME model'],
 ['FEL LogL', 'Site Log-likelihood under the FEL model']]

'''


