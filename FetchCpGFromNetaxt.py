#!/usr/bin/env python
import gzip
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


def process_block(NetAxt, line, achromsizes, out, dinuc):
  [seq1, seq2] =  read_block(NetAxt)
  [AlignNum, PChr, PStart, PEnd, AChr, AStart, AEnd, AStrand, Bscore] = line.split()
  achrsize = achromsizes[AChr]
  pstart = int(PStart)
  pend = int(PEnd)
  astart = int(AStart)
  aend = int(AEnd)
  L = len(seq1)
  gap1s = gap1e = 0
  gap2s = gap2e = 0
  i = 0
  while i < L :
    match = re.search('c-*g', seq1[i:], re.IGNORECASE)
    if match : 
      s = match.start() + i
      e = match.end() + i
      gap1s = seq1[0:s].count('-')
      gap1e = seq1[0:e].count('-')
      gap2s = seq2[0:s].count('-')
      gap2e = seq2[0:e].count('-')
      cg_pstart = pstart + s - 1 -gap1s
      alignmatch =  re.search('cg', seq2[s]+ seq2[e-1], re.IGNORECASE)
      if AStrand == '-' : 
        cg_astart = achrsize - (astart + e - gap2e) +1   # BED file tradition coordinates for di-nucleotide
        cg_aend = achrsize - (astart + s - gap2s) +1   
      else :
        cg_astart = astart + s - 1 -gap2s
        cg_aend = astart + e - 1 -gap2e

      assert(cg_astart >= 0)
      if( cg_aend > cg_astart):
        pstrandout = '+'
        astrandout = AStrand
        if(dinuc):
          pstartout = str(cg_pstart)
          pendout = str(cg_pstart+2)
          pseqout = seq1[s] + seq1[e-1]
          astartout = str(cg_astart)
          aendout = str(cg_aend)
          aseqout =  seq2[s] + seq2[e-1] 
        else: 
          pstartout = str(cg_pstart)
          pendout = str(cg_pstart+1)
          astrandout = '+'  # astart are coordinates for the C nucleotide in + strand when possible
          if(AStrand == '+' and seq2[s]!='-'): 
            pseqout = seq1[s] 
            astartout = str(cg_astart)
            aendout = str(cg_astart + 1)
            aseqout =  seq2[s] 
          elif(AStrand == '+' and seq2[s]=='-'): 
            pseqout = seq1[e-1] 
            astartout = str(cg_aend-1)
            aendout = str(cg_aend)
            aseqout = seq2[e-1] 
          elif(AStrand == '-' and seq2[e-1]!='-'): 
            pseqout = seq1[e-1] 
            astartout = str(cg_astart)
            aendout = str(cg_astart+1)
            aseqout = seq2[e-1] 
          elif(AStrand == '-' and seq2[e-1]=='-'): 
            pseqout = seq1[s] 
            astartout = str(cg_aend-1)
            aendout = str(cg_aend)
            aseqout = seq2[s]   

        out.write(Bscore + '\t' + PChr + '\t' + pstartout + '\t'+ pendout +'\t' + \
                  pseqout + '\t+\t' + \
                  AChr + '\t' + astartout + '\t'+ aendout + '\t' + \
                  aseqout + '\t'+ astrandout + '\n')
      i = e
    else: return
  return 




def main():
    
  parser =  argparse.ArgumentParser(description='Process net.axt file to map CpGs.')
  parser.add_argument("-o", "--outfile", help="name of output file")
  parser.add_argument("-s", "--chromsize", help="chrom size file for the aligning species")
  parser.add_argument("--dinuc", help="output dinucleotide (CG) mapping location instead of single nucleotide (C)",
                      action="store_true")
  parser.add_argument('AxtNet', help='<PSpecies>.<ASpecies>.net.axt file')
  args = parser.parse_args()

  achromsizes = {}
  with open(args.chromsize) as sizefile:
    for sizeline in sizefile:
      (key, val) = sizeline.split()
      achromsizes[key] = int(val)

  axtnet = gzip.open(args.AxtNet, 'r') if args.AxtNet[-2:]=="gz" else open(args.AxtNet, 'r')
  out = open(args.outfile, 'w')
    
  line = axtnet.readline()
  while is_header(line):
    line = axtnet.readline()

  while is_start_block(line) :
    process_block(axtnet, line, achromsizes, out, args.dinuc)
    line = axtnet.readline()

  out.close()
  
if __name__ == "__main__":
  main()     
