#!/usr/bin/env python
import sys
import os
import re
import argparse

""" net.axt file format
https://genome.ucsc.edu/goldenPath/help/axt.html

1. Summary line
 0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500

The summary line contains chromosomal position and size information about the alignment. It consists of 9 required fields:

    Alignment number -- The alignment numbering starts with 0 and increments by 1, i.e. the first alignment in a file is numbered 0, the next 1, etc.
    Chromosome (primary organism)
    Alignment start (primary organism) -- The first base is numbered 1.
    Alignment end (primary organism) -- The end base is included.
    Chromosome (aligning organism)
    Alignment start (aligning organism)
    Alignment end (aligning organism)
    Strand (aligning organism) -- If the strand value is "-", the values of the aligning organism's start and end fields are relative to the reverse-complemented coordinates of its chromosome.
    Blastz score -- Different blastz scoring matrices are used for different organisms. See the README.txt file in the alignments directory for scoring information specific to a pair of alignments. 
"""

def is_header(line) :
  return(line[0:1] == "#")

def is_start_block(line):
  fields=line.split()
  return( len(fields) == 9)


def read_block(NetAxt):
  seq1 = NetAxt.readline().replace('\n', '')
  seq2 = NetAxt.readline().replace('\n', '')
  na = NetAxt.readline()
  assert  na == "\n"
  return[seq1, seq2]


def process_block(NetAxt, line, achromsizes, out):
  [seq1, seq2] =  read_block(NetAxt)
  [AlignNum, PChr, PStart, PEnd, AChr, AStart, AEnd, AStrand, Bscore] = line.split()

  achrsize = achromsizes[AChr]
  pstart = int(PStart)
  pend = int(PEnd)
  if AStrand == '+':
    astart = int(AStart)
    aend = int(AEnd)
  else : 
    astart = achrsize - int(AEnd) +1   # BED file tradition coordinates for di-nucleotide
    aend = achrsize - int(AStart)  +1   


  L = len(seq1)
  nmatch = 0
  nts = 0
  ntv = 0 
  gap1 = 0
  gap2 = 0
  nother = 0

  match = set(['AA', 'TT', 'GG', 'CC'])
  transition = set([ 'AG', 'GA', 'CT', 'TC'])
  transversion = set(['AC', 'CA', 'AT', 'TA', 'CG', 'GC', 'TG', 'GT'])

  for i in range(L):
    if seq1[i] =='-' :
      gap1 = gap1 +1
    elif seq2[i] == '-':
      gap2 = gap2 +1
    else :
      S = (seq1[i] + seq2[i]).upper()
      if S in match :
        nmatch = nmatch +1
      elif S in transition :
        nts = nts +1
      elif S in transversion :
        ntv = ntv +1
      else :
        nother = nother +1
 
  out.write(PChr + '\t' + str(pstart) + '\t' + str(pend+1) +'\t' + \
            AChr + ':' + str(astart) + ':'+ str(aend+1) + ':' + AStrand + '\t'+ \
            Bscore + '\t'+ '+' + '\t' + \
            str(nmatch) + '\t' + str(nts) + '\t' + str(ntv) + '\t' + \
            str(gap1) + '\t' + str(gap2) + '\t' + str(nother) + '\n')

  return




def main():
    
  parser =  argparse.ArgumentParser(description="Process net.axt file to report #of matches and mismatches."
                                    "Additional fields: #MATCH #TRANSITION #TRNSVERSION #GAP_TARGET #GAP_QUERY #UNKNOWN")
  parser.add_argument("-o", "--outfile", help="name of output bed format file")
  parser.add_argument("-s", "--chromsize", help="chrom size file for the aligning species")
  parser.add_argument('AxtNet', help='<PSpecies>.<ASpecies>.net.axt file')
  args = parser.parse_args()

  achromsizes = {}
  with open(args.chromsize) as sizefile:
    for sizeline in sizefile:
      (key, val) = sizeline.split()
      achromsizes[key] = int(val)

  out = open(args.outfile, 'w')
  axtnet = open(args.AxtNet, 'r')
    
  line = axtnet.readline()
  while is_header(line):
    line = axtnet.readline()

  while is_start_block(line) :
    process_block(axtnet, line, achromsizes, out)
    line = axtnet.readline()

  out.close()
  
if __name__ == "__main__":
  main() 
