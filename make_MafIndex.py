"""
  Make MafIndex file from maf file 
  Pre-requisite: 
  Usage: python make_MafIndex.py chr22
"""

import sys
sys.path.insert(1, "/home/rcf-40/jqu/panfs/tools/python_packages/alignio-maf")

import Bio
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
import glob

chrom=sys.argv[1]
maf_file=chrom+".7sp.maf" #change the maf file suffix if needed
idx = MafIndex(chrom+".mafindex", maf_file, "hg19."+chrom)  #change the reference name if necessary
