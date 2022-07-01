#!/usr/bin/env python


'''
Color the protein structure and run some exploratory stats
'''


import argparse
import hashlib
import warnings
warnings.filterwarnings('ignore')
import sys

from foldvis.models import Fold


parser = argparse.ArgumentParser(description='Remove short contigs')
parser.add_argument(
    '--model', required=True, help='Path to AlphaFold v2 results')
args = parser.parse_args()


model = Fold(args.model)
seq = model.sequence

hsh = hashlib.md5(seq.encode('utf-8')).hexdigest()
sys.stdout.write(hsh)