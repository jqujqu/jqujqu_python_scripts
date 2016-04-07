#!/usr/bin/env python
import sys
import os
import re
import argparse
import gzip 

def peek_line(f):
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line

def is_SEQ_line(line):
  return  len(line) > 3 and line[:3]== "SEQ"

def is_TREE_line(line):
  return line[:4] == "TREE"

def skip_comments(file_obj, tag) :
  line = file_obj.readline()
  while line[:len(tag)]==tag :
    line = file_obj.readline()
  return line

def get_SEQ_chunk(file_obj) :
  header_lines = []
  line = file_obj.readline()
  while len(line) > 0 and not is_SEQ_line(line) : 
    line = file_obj.readline()
  if len(line) > 0 :
    while is_SEQ_line(line) :
      #print line.rstrip()
      header_lines.append(line.rstrip())
      line = file_obj.readline()
    #print line, len(line)
    if len(header_lines) > 0 :
      assert is_TREE_line(line)
      line = file_obj.readline()
      assert line.rstrip() == "DATA"
      Species = []; Chrom = []; Start = []; End = []; Strand = []
      for hl in header_lines :
        (jnk1, sp, ch, st, en, sd, jnk2) = hl.split() 
        Species.append(sp)
        Chrom.append(ch)
        Start.append(st)
        End.append(en)
        Strand.append(sd)
      seq_header = dict([("species", Species), ("chrom", Chrom), 
                       ("start", Start), ("end", End), ("strand",Strand )] )
  else :
    seq_header = dict()
  return (seq_header)


def pick_ancestor(ancestors) :
  ind = []
  for anc in ancestors : 
    ind.append(int(anc.split('_')[2]))
  smst = 0 
  for i in range(1,len(ind)) :
    if ind[i] < ind[smst] :
      smst = i
  return ancestors[smst]

def label_selection(seq_header, target_species) :
  ind = []
  for i in seq_header["species"] :
    if i == "ancestral_sequences" : 
      ind.append(0)
    elif target_species.count(i) == 1 :
      ind.append(1)
    else :
      ind.append(-1)
  last_valid = len(ind) -1
  while ind[last_valid] != 1 :
    last_valid = last_valid -1
  fist_valid = 0
  while ind[fist_valid] != 1 :
    fist_valid = fist_valid +1
  #
  ancestors_chrom_to_keep = []
  pos = fist_valid
  while pos < last_valid :
    assert ind[pos] ==1
    pos = pos +1 
    ancestors = []
    while ind[pos] != 1 :
      if ind[pos] == 0:
        ancestors.append(seq_header["chrom"][pos])
      pos = pos + 1
    ancestors_chrom_to_keep.append(pick_ancestor(ancestors))
  #
  for pos in range(len(ind)) :
    if ind[pos] == 0:
      if ancestors_chrom_to_keep.count(seq_header["chrom"][pos])==0 :
        ind[pos] = -1
      else :
        ind[pos] = 1
  #
  return ind

def select_species(seq_header, target_species) :
  """
    return a list of indices where target_species occur in the seq_header records
  """
  counter = 0
  for i in target_species :
    if seq_header["species"].count(i) == 1 :
      counter = counter +1
  labels = []
  if counter == len(target_species) :
    #make sure right order
    forward = 0
    for i in target_species :
      if i in seq_header["species"][forward:] :
        forward = seq_header["species"].index(i, forward)+1
      else :
        print "wrong order"
        return labels
    ind = label_selection(seq_header, target_species)
    for i in range(len(ind)) :
      if ind[i] == 1:
        labels.append(i)
  else :
    print "Incomplete or duplicated"
  return labels 


def skip_block(file_obj) :
  line = file_obj.readline()
  while line.rstrip()!="//" and len(line): 
    line = file_obj.readline()
  return line


def parent_index(seq_headers, labels) :
  """ get names of the useful subset"""
  parent_labels = [-1]*len(labels)
  ancestors = []
  for i in labels : 
    if seq_headers["species"][i] == "ancestral_sequences" :
      ancestors.append(int(seq_headers["chrom"][i].split('_')[2]))
  ancestor_sorted_idx = list(reversed(sorted(range(len(ancestors)), 
                                             key=lambda k: ancestors[k]))) # descending
  for i in range(len(ancestor_sorted_idx)) :
    anc_id_in_labels = ancestor_sorted_idx[i]*2 + 1
    if (parent_labels[anc_id_in_labels-1] == -1) : 
      parent_labels[anc_id_in_labels-1] = labels[anc_id_in_labels]
    else :
      node = (ancestor_sorted_idx[i]-1)*2+1
      while parent_labels[node] != -1 : 
        node = labels.index(parent_labels[node])
      parent_labels[node] = labels[anc_id_in_labels]
    
    if ( parent_labels[anc_id_in_labels+1]== -1 ) :
      parent_labels[anc_id_in_labels+1] = labels[anc_id_in_labels]
    else : 
      node = (ancestor_sorted_idx[i]+1)*2+1
      while parent_labels[node] != -1 : 
        node = labels.index(parent_labels[node])
      parent_labels[node] = labels[anc_id_in_labels]
  for i in range(len(parent_labels)) :
    if parent_labels[i] == -1 :
      print "root is ", seq_headers["species"][labels[i]]
      parent_labels[i] = labels[i]
  for i in range(len(labels)):
    print str(i)+"-"+str(labels[i]) + "-" + str(parent_labels[i]), \
          seq_headers["species"][labels[i]], seq_headers["chrom"][labels[i]], \
          "has parent",  seq_headers["species"][parent_labels[i]], \
          seq_headers["chrom"][parent_labels[i]]
  return parent_labels





def no_gap(seq, labels) :
  for i in labels :
    if seq[i] == "-" or seq[i] == 'N':
      return False
  return True

def all_identical(seq, labels):
  s = seq[labels[0]]
  for i in labels :
    if s != seq[i] :
      return False
  return True


def print_mutation(mutout, labels, parent_labels,
                    chrom, start, end, strandname,
                    prevline, line, nextline) : 
  if len(nextline) < len(line) :
    nextline = "N"*len(line)
  for i in range(len(labels)) :
    if line[labels[i]] != line[parent_labels[i]] :
      mutout.write("chr" + chrom + "\t" + str(start) + "\t" + str(end) +
                   "\tnode_" + str(i)+ ":" + prevline[parent_labels[i]] +
                   line[parent_labels[i]] + nextline[parent_labels[i]] + ":" +
                   prevline[labels[i]] + line[labels[i]] + nextline[labels[i]] +
                   "\t0\t" +  strandname + "\n")
  return

"""
Assume:
- ref sequence is the first leaf species 
"""
def process_alignment(file_obj, seq_headers, labels, parent_labels, mutout, alnout) :
  ref_bp_end = 0
  ref_bp_start = 0
  prevline = "N"*len(seq_headers["species"])
  ############################
  # read line code block 
  line = file_obj.readline().rstrip().upper()
  if line[labels[0]] in set(["A", "T", "G","C", "N"]) :
    ref_bp_end = ref_bp_end+1
  NOGAP = no_gap(line, labels)
  IDENTICAL = all_identical(line, labels)
  #############################
  # find first nongap line
  while (not NOGAP) and (line != "//"):
    prevline = line
    line = file_obj.readline().rstrip().upper()
    if (line[labels[0]] in set(["A", "T", "G","C", "N"])) :
      ref_bp_end = ref_bp_end+1
    NOGAP = no_gap(line, labels)
    IDENTICAL = all_identical(line, labels)
  # end while
  ref_bp_start = ref_bp_end
  chrom = seq_headers["chrom"][0]
  chrom_start = int(seq_headers["start"][0])
  chrom_end = int(seq_headers["end"][0])
  #
  INGAP = False
  #############################
  # loop alternate between
  # nongap block and gap block
  while line != "//" :
    if not INGAP :
      while NOGAP :
        if not IDENTICAL:
          nextline = peek_line(file_obj).rstrip().upper()
          # print mutations
          if seq_headers["strand"][0] == '1' :
            pos = chrom_start + ref_bp_end -2
            strandname = "+"
          else :
            pos = chrom_end - ref_bp_end
            strandname = "-"
          print_mutation(mutout, labels, parent_labels,
                         chrom, pos, pos+1, strandname,
                         prevline, line, nextline)
        # maintain conserved block
        prevline = line 
        line = file_obj.readline().rstrip().upper()
        if line[labels[0]] in set(["A", "T", "G","C", "N"]):
          ref_bp_end = ref_bp_end+1
        if line != "//" :
          NOGAP = no_gap(line, labels)
          IDENTICAL = all_identical(line, labels)
        else :
          NOGAP = False
      #while broke
      # mark gap and print previous conserved region
      INGAP = True   
      if seq_headers["strand"][0] == '1' :
        alnout.write("chr" + chrom + "\t" + str(chrom_start + ref_bp_start-2) + "\t" + 
                     str(chrom_start + ref_bp_end-2) + "\tALN" +
                     str(ref_bp_end-ref_bp_start) + "\t+\n")
      else :
        alnout.write("chr"+ chrom + "\t" + str(chrom_end - ref_bp_end) + "\t"+ 
                     str(chrom_end - ref_bp_start) + "\tALN" +
                     str(ref_bp_end-ref_bp_start) + "\t-\n")
    else : # INGAP is True
      while not NOGAP and line != "//" :  
        prevline = line
        line = file_obj.readline().rstrip().upper()
        if (line[labels[0]] in set(["A", "T", "G","C", "N"])) :
          ref_bp_end = ref_bp_end+1
        if  line != "//" :
          NOGAP = no_gap(line, labels)
          IDENTICAL = all_identical(line, labels)
      #end while
      INGAP = False
      ref_bp_start = ref_bp_end
    # end else INGAP is TRUE
  assert line == "//"
  return 


def main():
  """
  output node index is ordered by in-order traversal
  """
  parser =  argparse.ArgumentParser(description='Mutations from EPO compara emf file')
  parser.add_argument("-o", "--outfile", help="name of output file (BED format)")
  parser.add_argument('emf', help='emf file (gzipped ok)')
  args = parser.parse_args()
  alnout = open(args.outfile+".aligned.bed", 'w')
  mutout = open(args.outfile+".mutation.bed", 'w')
  if args.emf[-3:] == ".gz" :
    emf = gzip.open(args.emf, 'r')
  else :
    emf = open(args.emf, 'r')
  #emf = gzip.open("Compara.15_eutherian_mammals_EPO.chr10_10.emf.gz", 'r')
  target_species = [ "homo_sapiens", "pan_troglodytes", "gorilla_gorilla",
                     "macaca_mulatta", "mus_musculus",
                     "rattus_norvegicus","canis_familiaris"]
  line = skip_comments(emf, "#")
  while len(line) :
    seq_header = get_SEQ_chunk(emf)
    if len(seq_header) :
      labels = select_species(seq_header, target_species)  
      if len(labels) :
        parent_labels = parent_index(seq_header, labels)
        process_alignment(emf, seq_header, labels, parent_labels, mutout, alnout) 
        print "+\n"
    else :
      line = skip_block(emf)

if __name__ == "__main__":
  main()     

