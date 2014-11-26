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


def process_block(NetAxt, line, out):
  [seq1, seq2] =  read_block(NetAxt)
  [AlignNum, PChr, PStart, PEnd, AChr, AStart, AEnd, AStrand, Bscore] = line.split()
  pstart = int(PStart)
  pend = int(PEnd)
  astart = int(AStart)
  aend = int(AEnd)
  L = len(seq1)
  i = 0
  while i < L :
    match = re.search('-+', seq2[i:])
    if match : 
      s = match.start() + i
      e = match.end() + i
      gap = seq1[0:i].count('-')
      seg_pstart = pstart + i - 1 -gap
      seg_pend = pstart + s - 1 - gap
      i = e
    else: 
      seg_pstart = pstart + i - gap -1
      seg_pend = pend  
      i = L
    out.write(PChr + '\t' + str(seg_pstart) + '\t'+ str(seg_pend)+ '\tX\t0\t+\n')
  return 




def main():
    
  parser =  argparse.ArgumentParser(description='Get regions in Primary species that is mappable to Aligning species')
  parser.add_argument("-o", "--outfile", help="name of output file (BED format)")
  parser.add_argument('AxtNet', help='<PSpecies>.<ASpecies>.net.axt file')
  args = parser.parse_args()

  out = open(args.outfile, 'w')
  axtnet = open(args.AxtNet, 'r')
    
  line = axtnet.readline()
  while is_header(line):
    line = axtnet.readline()

  while is_start_block(line) :
    process_block(axtnet, line, out)
    line = axtnet.readline()

  out.close()
  
if __name__ == "__main__":
  main()     
