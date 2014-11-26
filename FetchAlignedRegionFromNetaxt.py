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

output is in Aspecies coordinates
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
  assert na == "\n"
  return[seq1, seq2]


def process_block(NetAxt, line, out, achromsizes):
  [seq1, seq2] = read_block(NetAxt)
  [AlignNum, PChr, PStart, PEnd, AChr, AStart, AEnd, AStrand, Bscore] = line.split()
  pstart = int(PStart)
  pend = int(PEnd)
  astart = int(AStart)
  aend = int(AEnd)
  L = len(seq1)
  achrsize = achromsizes[AChr]
  i = 0

  while i < L :
    match = re.search('-+', seq1[i:])
    
    if match :
      s = match.start() + i
      e = match.end() + i
      seg_astart = astart + i - seq2[0:i].count('-')
      seg_aend = astart + s - seq2[0:s].count('-')
      if AStrand == '-' :
        seg_astart = achrsize - (astart + s - seq2[0:s].count('-')) + 1  
        seg_aend = achrsize - (astart + i - seq2[0:i].count('-')) + 1
      i = e
    else: 
      seg_astart = astart + i - seq2[0:i].count('-')
      seg_aend = aend
      if AStrand == '-' :
        seg_astart = achrsize - aend + 1  
        seg_aend = achrsize - (astart + i - seq2[0:i].count('-')) + 1
      i = L
      
    out.write(AChr + '\t' + str(seg_astart) + '\t'+ str(seg_aend)+ '\tX\t0\t+\n')
  
  
return


def main():
  parser = argparse.ArgumentParser(description='Get regions in Primary species that is mappable to Aligning species')
  parser.add_argument("-s", "--chromsize", help="chrom size file for the aligning species")
  parser.add_argument("-o", "--outfile", help="name of output file (BED format in A coordinates) ")
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
    process_block(axtnet, line, out, achromsizes)
    line = axtnet.readline()

  out.close()

if __name__ == "__main__":
  main()
