#!/usr/bin/env python
import sys
import os
import re
import argparse
"""
Create bin-size lift-over map
"""

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


def process_block(NetAxt, line, achromsizes, binsize, out):
  [seq1, seq2] =  read_block(NetAxt)
  [AlignNum, PChr, PStart, PEnd, AChr, AStart, AEnd, AStrand, Bscore] = line.split()
  achrsize = achromsizes[AChr]
  pstart = int(PStart)
  pend = int(PEnd)
  astart = int(AStart)
  aend = int(AEnd)
  L = len(seq1)

  shift = (binsize - (pstart % binsize))%binsize
  if (pend-pstart  < shift + binsize) :
    return
  i = 0
  pcount = 0
  acount = 0
  while (pcount < shift) :
    if (seq1[i] != '-') : 
      pcount +=1
    if (seq2[i] != '-') :
      acount +=1  
    i += 1
  prev_pcount = pcount
  prev_acount = acount
  while i < L and pcount + binsize < pend-pstart: 
    while (prev_pcount + binsize < pend-pstart and pcount < prev_pcount + binsize ) :
      if (seq1[i] != '-') : 
        pcount +=1
      if (seq2[i] != '-') :
        acount +=1  
      i += 1
    if (pcount - prev_pcount == binsize) : 
      if AStrand == '+' : 
        out.write(PChr +"\t" + str(pstart+prev_pcount) + "\t" + \
                  str(pstart + pcount) +"\t" + AChr + "\t"+ str(astart + prev_acount) + \
                  "\t" + str( astart + acount) + "\t+\n" )
      else :
        out.write(PChr +"\t" + str(pstart+prev_pcount) + "\t" + \
                  str(pstart + pcount) +"\t" + AChr + "\t"+ str(achrsize - (astart + acount)) + \
                  "\t" + str(achrsize-(astart + prev_acount)) + "\t+\n" )
    prev_pcount = pcount
    prev_acount = acount
    if not (prev_pcount + binsize <= pend-pstart and pcount < prev_pcount + binsize) :
      break
  return 




def main():
    
  parser =  argparse.ArgumentParser(description='Process net.axt file to map CpGs.')
  parser.add_argument("-o", "--outfile", help="name of output file")
  parser.add_argument("-s", "--chromsize", help="chrom size file for the aligning species")
  parser.add_argument("-b", "--binsize", help="binsize")
  parser.add_argument('--AxtNet', help='<PSpecies>.<ASpecies>.net.axt file')
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

  binsize = int(args.binsize)
  while is_start_block(line) :  
    process_block(axtnet, line, achromsizes, binsize, out)
    line = axtnet.readline()

  out.close()
  
if __name__ == "__main__":
  main()     

