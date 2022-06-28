#!/usr/bin/env python


'''
Color the protein structure and run some exploratory stats
'''


import argparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from foldvis.io import load_conserved
from foldvis.models import Fold, Complex, AlphaFold, Binding
from foldvis.vis import *  # yes I know, whatever


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--model', required=True, help='Path to AlphaFold v2 results')
parser.add_argument(
    '--pfam', required=True, help='Path to folder with Pfam database')
parser.add_argument(
    '--out', default='domains.csv', help='Results')
args = parser.parse_args()


model = Fold(args.model)

b = Binding(model, 'representable')
'''
- 'confident': 'InteracDome_v0.3-confident.tsv'
- 'representable': 'InteracDome_v0.3-representable.tsv'
- 'non-redundant': 'InteracDome_v0.3-representableNR.tsv'
'''
pfam = Path(args.pfam) / 'Pfam-A.hmm'

b.predict_binding_(pfam)
print(b.domains)

b.domains.to_csv(args.out, index=False)

# with open('domains.csv', 'w+') as out:
#     out.write(model.sequence)


for dom in b.domains['acc'].unique(): pass
    # might be empty
