#!/home/rcf-40/jqu/panfs/tools/Python-2.7.8/python
import gzip 
import sys
import argparse

"""
parse EPO internal node names
"""

def parse_epo_internal(internal):
  id = []
  for seq in internal :
    parts = seq.split(".")
    assert parts[0] == "ancestral_sequences"
    id.append( int(parts[1].split("_")[2]) )
  m = min (id)
  return [ i-m for i in id]

def tree_build(internal, leaves):
  if len(internal) == 0 and len(leaves) == 1 :
    return {leaves[0]:[]}
  ord = parse_epo_internal(internal)
  rootid = ord.index(0)
  internal_l = internal[:rootid]
  internal_r = internal[(rootid+1):]
  leaves_l = leaves[:(rootid+1)]
  leaves_r = leaves[(rootid+1):]
  lsize = len(leaves_l)
  rsize = len(leaves_r)
  lt = tree_build(internal_l, leaves_l)
  rt = tree_build(internal_r, leaves_r)
  t = {internal[rootid]: [lt, rt]}
  #t = tree_trim(t)
  return t

def in_order_traverse(tree) :
  if len(tree.values()[0])==0 :
    return tree.keys()
  root = tree.keys()[0]
  left = tree.values()[0][0]
  right = tree.values()[0][1]
  return in_order_traverse(left) + [root] + in_order_traverse(right)


###############################################################################
###############################################################################
class PhyloTree:
  def is_empty(self):
    return len(self.tree) == 0
  #
  def get_subtrees(self):
    nchild = len(self.tree.values()[0])
    if  nchild > 0 :
      return [ PhyloTree(self.tree.values()[0][i]) for i in range(nchild)]
    else :
      return [PhyloTree({})]
  #
  def trim(self):
    null_leaf = {'':[]}
    root = self.tree.keys()[0]
    if len(self.tree.values()[0]) == 0 :
      return PhyloTree(self.tree)
    subtrees = self.get_subtrees()
    left = subtrees[0]
    right = subtrees[1]
    if left.tree == null_leaf and right.tree == null_leaf :
      return PhyloTree({})
    elif left.tree == null_leaf :
      return right.trim()
    elif right.tree == null_leaf:
      return left.trim()
    else :
      l = left.trim()
      r = right.trim()
      if l.is_empty() :
        return r
      elif r.is_empty() :
        return l
      else :
        return PhyloTree({root:[l.tree, r.tree]})
  #
  def inorder(self):
    if len(self.tree.values()[0])==0 :
      return self.tree.keys()
    root = self.tree.keys()[0]
    subtrees = self.get_subtrees()
    left = subtrees[0]
    right = subtrees[1]
    return left.inorder() + [root] + right.inorder()
  #
  def inorder_order(self):
    spinorder = self.inorder()
    io = {}
    for i in range(len(spinorder)) :
      io[spinorder[i]] = i
    return io
  #
  def set_childidx(self):
    inordermap = self.inorder_order()
    self.internal_set_childidx(self.tree, inordermap)
  #
  def internal_set_childidx(self, tree, inordermap):
    root = tree.keys()[0]
    size = len(tree.values()[0]) 
    if size == 0 :
      self.childidx[inordermap[root]] = [-1] #means it's a leaf
      return
    else :
      children = [ subtree.keys()[0] for subtree in tree.values()[0]]
      self.childidx[inordermap[root]] = [ inordermap[child] for child in children]
      for subtree in tree.values()[0] :
        self.internal_set_childidx(subtree, inordermap)
      return
  #
  def set_parentidx(self):
    if len(self.childidx) == 0 :
      self.set_childidx()
    for i in self.childidx :
      for j in self.childidx[i] :
        if j != -1 :
          self.parentidx[j] = i
  #
  def get_leaves(self):
    if len(self.childidx) == 0 :
      self.set_childidx()
    return [i for i in self.childidx if self.childidx[i] == [-1]]
  #
  def get_root_id(self):
    for i in self.childidx :
      if i not in self.parentidx :
        return i
  #
  def __init__(self, tree):
    self.tree = tree
    self.childidx = {}
    self.parentidx = {}


###############################################################################
"""
maf parsing
"""
def is_BlockStart(s):
  return s.split()[0] == 'a'

def is_Sequence(s) :
  return s.split()[0]=='s'

def read_block_maf(f) :
  line = f.readline()
  while len(line) and not is_BlockStart(line):
    line = f.readline()
  seqlines =[]
  line = f.readline()
  while (line != '\n') :
    if is_Sequence(line):
      seqlines.append(line)
    line = f.readline()
  return seqlines

class Seq_Rec:
  """
  s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
  """
  def __init__(self, seq):
    s = seq.split()
    self.id = s[1]
    self.strand = s[4]
    self.srcSize = int(s[5])
    self.seq = s[6].upper()
    self.size = int(s[3])
    self.start = int(s[2])
    self.stdstart = 0
    if self.strand[0] == '+' :
      self.stdstart = int(s[2])
    else :
      self.stdstart = self.srcSize - self.start -self.size +1 

def filter_alignment(alignment, target_species) :
  nodes = []
  for i in alignment :
    nodes.append(i.id)
  #endfor
  assert (len(nodes)+1)%2 == 0
  n = (len(nodes)+1)/2
  internal = [ nodes[i] for i in range(1,len(nodes), 2) ]
  leaves = [ nodes[i] for i in range(0,len(nodes), 2) ]
  #mask non-target leaf species name
  for  i in range(n):
    if leaves[i].split(".")[0] not in species :
      leaves[i] = ""
  #endfor
  #get all species in the tree
  T = PhyloTree(tree_build(internal, leaves))
  T = T.trim()
  filtered_species = T.inorder()
  filtered_alignment = []
  for i in alignment :
    if i.id in filtered_species :
      filtered_alignment.append(i)
  return filtered_alignment

def get_MA_column(alignment, i) :
  s = ''
  L = len(alignment[0].seq)
  if i >= L :
    return s
  else :
    for seqrec in alignment :
      s += seqrec.seq[i]
    return s

def get_MA_nodes(alignment): 
  nodes = []
  for seqrec in alignment :
    nodes.append(seqrec.id)
  return nodes


def peek_line(f):
  pos = f.tell()
  line = f.readline()
  f.seek(pos)
  return line


def get_mutation(s, entryid, cidx) :
  m = []
  if cidx[entryid] == [-1] :
    return []
  else :
    for id in cidx[entryid]:
      if s[entryid] != s[id]:
        m.append(id)  
      m += get_mutation(s, id, cidx)
    return m



species = ["homo_sapiens", "pan_troglodytes", "gorilla_gorilla", 
           "mus_musculus", "rattus_norvegicus", "canis_familiaris", 
           "macaca_mulatta"]

def main():
  parser = argparse.ArgumentParser(description='report Mutation sites in aligned positions', 
                                   prog='EPOmaf2Mut')
  parser.add_argument('--mafgz', dest='maffilegz', help='gzipped epo maf file')
  parser.add_argument('--output', dest='output', help='output file')
  if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()

  outfile=args.output # "EPO.chr1_1.nongap_mutations.bed"
  maffilegz=args.maffilegz  #"Compara.17_eutherian_mammals_EPO.chr1_1.maf.gz"

  f = gzip.open(maffilegz, 'r')
  line =  peek_line(f)
  alignments = []
  while len(line):
    alignment =[]
    seq_lines = read_block_maf(f) 
    for seq in seq_lines :
      alignment.append(Seq_Rec(seq))
    naln = filter_alignment(alignment, species)
    if len(naln) == 2*len(species)-1 :
      alignments.append( naln )
    line = peek_line(f)

  #build tree and get childidx
  nodes = get_MA_nodes(alignments[0])
  internal = [ nodes[i] for i in range(1,len(nodes), 2) ]
  leaves = [ nodes[i] for i in range(0,len(nodes), 2) ]
  T = PhyloTree(tree_build(internal, leaves))
  T.set_childidx()
  T.set_parentidx()
  inorder = T.inorder()
  rootid = T.get_root_id()
  cidx = T.childidx
  pidx = T.parentidx

  block_mut_rec = {}
  for i in range(len(species)*2-1) :
    block_mut_rec[i] = 0

  block_mut_rec["identical"] = 0
  block_mut_rec["multiple_hits"] = 0
  block_mut_rec["contain_gap"] = 0

  i = 0
  count = 0
  MA = alignments[0]
  ALNsize = len(MA[0].seq)

  out = open(outfile, 'w')
  while i < ALNsize:
    s = get_MA_column(MA, i)
    if s[0] != "-" :
      count += 1
      if s.count('-') == 0 and s.count('N') == 0  :
        mut = get_mutation(s, rootid, cidx)
        if len(mut) == 0:
          block_mut_rec["identical"] = block_mut_rec["identical"]+1
        else :
          spid = mut[0]
          block_mut_rec[spid] = block_mut_rec[spid]+1
          ref = MA[0].id
          start = MA[0].start + count
          end = start +1
          strand = MA[0].strand[:1]
          if strand == '-' :
            [end, start] = [MA[0].srcSize - start+1, MA[0].srcSize - end +1 ]
          else :
            [start, end] = [start -1, end -1]
          if len(mut) ==1 : 
            name = s[pidx[spid]]+"="+s[spid] 
            score = spid       
            out.write( "%s\t%s\t%s\t%s\t%s\t%s"%(ref, start, end, name, score, strand) + '\n')
          else :
            block_mut_rec["multiple_hits"] = block_mut_rec["multiple_hits"]+1
            for spid in mut :
              name = s[pidx[spid]]+"="+s[spid] 
              score = spid
              out.write( "%s\t%s\t%s\t%s\t%s\t%s"%(ref, start, end, name, score, strand) + '\n')
      else :
        block_mut_rec["contain_gap"] = block_mut_rec["contain_gap"]+1
    i += 1

  f.close()
  out.close()


if __name__ == "__main__":
  main()
