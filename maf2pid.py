#!/home/rcf-40/jqu/panfs/tools/Python-2.7.8/python
""" 
  maf2pid
  Compute PID (percentage identity )for give genomic intervals (in hg19)
    @param species (species list)
    @param interval (bed file name)
    @param index (MafIndex built by Bio.Align.MafIO.MafIndex)
    @param maf (corresponding maf file)
    @param outfile (output file name)
    @return boolean result
"""

import sys
sys.path.insert(1, "/home/rcf-40/jqu/panfs/tools/python_packages/alignio-maf")
import argparse
import Bio
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex

###############################################################################
def parse_bed(bedfile):
  """
  Returns a dictionary with chromosome name as keys, 
  and interval starts and ends in two lists as values
  """
  f = open(bedfile, 'r')
  intervals_by_chrom = {}
  l = f.readline()
  while (len(l)) :
    [chrom, start, end, name] = l.split()[0:4]
    start = int(start)
    end = int(end)
    if chrom in intervals_by_chrom :
      starts = intervals_by_chrom[chrom][0] + [start]
      ends = intervals_by_chrom[chrom][1] + [end]
      names = intervals_by_chrom[chrom][2] + [name]
      del intervals_by_chrom[chrom]
      intervals_by_chrom[chrom] = [starts, ends, names]
    else :
      intervals_by_chrom[chrom] = [[start], [end], [name]]
    l = f.readline()
  return intervals_by_chrom

def read_species(splst) :
  species = []
  f = open(splst, 'r')
  for line in f:
    species.append(line.replace("\n", ""))
  f.close()
  return species

def make_pair_names(splst) :
  """
  Read in species list file and make paired names 
  Returns a list: [species, pairs]
  """
  species = read_species(splst)
  pairwise = []   
  for i in range(len(species)-1):
    for j in range(i+1, len(species)) :
      pairwise.append(species[i]+ "_"+ species[j])
  return [species, pairwise]


def get_PID(seq1, seq2):
  """
  PID = #Identical/#Aligned
  """
  assert len(seq1) == len(seq2)
  SEQ1 = seq1.upper().replace("N", "-")
  SEQ2 = seq2.upper().replace("N", "-")
  identical = 0
  aligned = 0
  for i in range(len(seq1)) :
    if SEQ1[i] != "-" and SEQ2[i] != "-" :
      aligned += 1
      if SEQ1[i] == SEQ2[i] :
        identical += 1
  if aligned :
    return float(identical)/aligned
  else :
    return -1

###############################################################################

def main():
  parser = argparse.ArgumentParser(description='report pair-wise PID in intervals', 
                                   prog='maf2pid')
  parser.add_argument('--splst', dest='splst', help='species list')
  parser.add_argument('--bed', dest='interval', help='bed intervals')
  parser.add_argument('--output', dest='output', help='output file')
  if len(sys.argv) < 4:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()

  intervals_by_chrom = parse_bed(args.interval)
  [species, pairwise] = make_pair_names(args.splst)

  fout = open(args.output, 'w')
  header = "#chrom\tstart\tend\tpme"
  for i in pairwise :
    header += "\t"+i
  fout.write( header + "\n")
  
  for chrom in intervals_by_chrom.keys() :
    sqlite_file = chrom + ".mafindex"
    maf_file = chrom + ".7sp.maf"
    target_seqname = "hg19." + chrom
    idx = MafIndex(sqlite_file, maf_file, target_seqname)
    # MafIndex.get_spliced(starts, ends, strand='+1')
    starts = intervals_by_chrom[chrom][0]
    ends = intervals_by_chrom[chrom][1]
    names = intervals_by_chrom[chrom][2]
    for k in range(len(starts)) :
      result = idx.get_spliced([starts[k]], [ends[k]], strand='+1')
      recsp = [rec.id.split(".")[0] for rec in result]
      PIDscores = {}
      for i in range(len(recsp)-1) :
        for j in range(i+1, len(recsp)) :
          pw = recsp[i]+ "_"+ recsp[j]
          if pw not in pairwise :
            pw = recsp[j]+ "_"+ recsp[i]
          if (pw in pairwise) and (not PIDscores.has_key(pw)): # exclude paralogs
            pid = get_PID(str(result[i].seq), str(result[j].seq))
            PIDscores[pw] = pid
      for pw in pairwise :
        if pw not in PIDscores :
          PIDscores[pw] = -1
      fout.write(chrom + "\t" + str(starts[k]) + "\t" + str(ends[k]) + "\tpme_" + names[k])
      for pw in pairwise :
        fout.write("\t"+ "{:.5f}".format(PIDscores[pw]))
      fout.write("\n")
  fout.close()
  return 

if __name__ == "__main__":
  main()
