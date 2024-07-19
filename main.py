from Bio import SeqIO
from argparse import ArgumentParser
import numpy as np

class Genome:
  def __init__(self, taxon:str, seq:str):
    self.taxon = taxon
    self.seq = seq

class Contig:
  def __init__(self, name:str, seq:str, target:Genome):
    self.name = name
    self.seq = seq
    self.target = target
  
  def __str__(self):
    return "[" + self.name + ", " + self.target.taxon + "]"

"""
Description:
A helper method that returns all the substrings of a string.
Returns a set containing all subsets of the input string. If a substring is
present multiple times, it is only counted once.

Parameters:
  s | str | the input string

Complexity: O(1) Average
Note:
len(s) = 6 is fixed since we only look at substrings for hexamers.
"""
def substrings(s: str):
  res = set()
  for i in range(len(s)):
    for j in range(i+1, len(s)):
      res.add(s[i:j])
  return res

def todict(G: list):
  res = dict()
  for rule in G:
    term = rule[0].split("->")[1].split("S")[0]
    res[term] = rule[1]
  return res

def set_bases(G: dict):
  G['A'] = 1.0
  G['G'] = 1.0
  G['T'] = 1.0
  G['C'] = 1.0
  return G

def gen_gram(S: str):
  d = []
  f = ""
  for i in range(len(S)):
    f += S[i]
    if f not in S[0:i]:
      d.append(str(f))
      f=""
  return d

def gen_sto_gram(prob: dict, S: str, cuts=0):
  d = gen_gram(S)
  stod = dict()
  for rule in d:
    term = rule
    try:
      p = prob[term]
    except:
      p=0
    stod[term] = p
  e = list(stod.keys())
  e = sorted(e, key=len)
  i = 0
  while (i < int(cuts) and i < len(e)):
    stoe = dict(stod)
    if (len(e[i]) < 2):
      i+=1
    else:
      stoe.pop(e[i],None)
      ptree = parse(stoe, S, get_max(e))
      p = pval(stoe, ptree, S)
      if (p < 1):
        stod.pop(e[i], None)
      i+=1
  stod = set_bases(stod)
  return stod

def get_max(input: list):
  maximum = 0
  for item in input:
    if (len(item) > maximum):
      maximum = len(item)
  return maximum

def kmer_counts(kmers: dict, contig: Contig):
  hexamer_counts = dict()
  sequence = str(contig.seq)
  for i in range(len(sequence) - 5):
    hexamer = sequence[i:i+6]
    if (hexamer in hexamer_counts.keys()):
      hexamer_counts[hexamer] += 1
    else:
      hexamer_counts[hexamer] = 1
  return hexamer_counts

"""
Description:
Outputs a dictionary that contains the probability of randomly choosing a kmer
from the set of contigs that contains any given substring

Parameters:
  kmers | dict | The dictionary of kmer frequencies

Complexity: O(len(kmers.keys()))
"""
def get_kmer_prob(kmers: dict, min_prob:float):
  substring_count = dict()
  for kmer in kmers.keys():
    for string in substrings(kmer):
      try:
        substring_count[string] += 1
      except:
        substring_count[string] = 1

  n = len(kmers.keys())
  for substring in substring_count.keys():
    substring_count[substring] = np.max([float(substring_count[substring]/n), float(min_prob)])
  substring_count = set_bases(substring_count)
  return substring_count

"""
Description:
Parses the input string using the provided probabilistic grammar as a dictionary.

Parameters:
  G | dict | The grammar that is used to parse the input
  input | str | input string
  lookback | int | The size of the buffer used to determine available rules.

Complexity: O(len(input)*len(max_rule)), where max_rule is the rule in grammar
G with the longest production. Note that the following code block runs in O(1)

  n = len(buffer)
      for j in range(n):
        if (buffer[n-j-1:n] in G.keys()):
          if (replaceable[n-j-1]==1):
            max_rule = buffer[n-j-1:n]

since the buffer size is restricted to being less than or equal to the
lookback.

Note: Ideally, one should make sure that the lookback is longer than the
largest rule.
"""
def parse(G: dict, input: str, lookback: int):
  parsetree = []
  buffer = ""
  replaceable = []
  for st in range(len(input)):
    c = input[st]
    if (len(buffer) < lookback):
      buffer+=c
      replaceable.append(1)
    else:
      buffer = buffer[1:len(buffer)]
      buffer+=c
      replaceable.pop(0)
      replaceable.append(1)
    max_rule = c
    n = len(buffer)
    for j in range(n):
      if (buffer[n-j-1:n] in G.keys()):
        if (replaceable[n-j-1]==1):
          max_rule = buffer[n-j-1:n]

    replaced_rules = []
    if (len(parsetree) == 0):
      parsetree.append(max_rule)
    else:
      sum=0
      it=0
      while (sum < len(max_rule)-1):
        replaced_rules.append(parsetree[len(parsetree)-it-1])
        sum+=len(parsetree[len(parsetree)-it-1])
        it+=1

      p1 = 1
      for rule in replaced_rules:
        try:
          p1*=G[rule]
        except:
          p1 = float('inf')
      p2 = 1
      try:
        p2 = G[max_rule]
      except:
        p2 = float('inf')
      thres = 1
      try:
        thres = p1/(p1+p2)
      except:
        pass
      if (np.random.random() < thres):
        for _ in replaced_rules:
          parsetree.pop()
        parsetree.append(max_rule)

        for k in range(len(max_rule)-1):
          replaceable[len(replaceable)-k-1]=0
        replaceable[len(replaceable)-len(max_rule)]=1
      else:
        parsetree.append(c)
  return parsetree

def get_sig(G1: dict, s2: str):
  parsetree = parse(G1,s2,get_max(G1.keys()))
  return pval(G1, parsetree, s2)

def pval(G: dict, tree: list, s: str):
  res = ""
  for rule in tree:
    res += rule

  p=1
  for rule in tree:
    if rule in G.keys():
      p*=G[rule]
  return p

def bin(contigs: list[Contig], sig: float, cut_rules=20, min_prob=0.05):
  clusters = []
  contig_list = list(contigs)
  kmers = dict()

  for contig in contigs:
    kmers = kmer_counts(kmers, contig)
  prob = get_kmer_prob(kmers, min_prob)

  i = 0
  while (i < len(contig_list)):
    icontig = contig_list[i]
    curcluster = [icontig] 
    G1 = gen_sto_gram(prob, icontig.seq, cut_rules)
    G1 = set_bases(G1)
    j=i+1
    while (j < len(contig_list)):
      jcontig = contig_list[j]
      p_val = get_sig(G1, jcontig.seq)
      if (p_val <= sig):
        curcluster.append(jcontig)
        contig_list.pop(j)
      else:
        j+=1
    clusters.append(curcluster)
    i+=1
  return clusters

def disp_bins(bins: list):
  res = ""
  for bi in range(len(bins)):
    res += "Bin " + str(bi) + "\n"
    for contig in bins[bi]:
      res += "\t" + str(contig) + "\n"
  return res

def disp_bin_stats(bins: list):
  res = ""
  for bi in range(len(bins)):
    taxa_counts = dict()

    tot = 0
    for contig in bins[bi]:
      try:
        taxa_counts[contig.target.taxon]+=1
      except:
        taxa_counts[contig.target.taxon]=1
      tot+=1
    
    taxons = taxa_counts.keys()
    for taxa in taxons:
      res += str(np.round(100*taxa_counts[taxa]/tot)) + "% of Bin " + str(bi) + " are " + taxa + "\n"
  return res

def main():
  parser = ArgumentParser()
  parser.add_argument("--file-1", dest="file1", help="fasta file 1")
  parser.add_argument("--taxon-1", dest="taxon1", help="taxon 1")
  parser.add_argument("--file-2", dest="file2", help="fasta file 2")
  parser.add_argument("--taxon-2", dest="taxon2", help="taxon 2")
  parser.add_argument("--sig", dest="sig", help="significance level")
  parser.add_argument("--cuts", dest="cuts", help="Number of grammar cuts")
  parser.add_argument("--min-prob", dest="min_prob", help="Minimum probability allowed")
  args = parser.parse_args()
  contigs = []
  
  it=0
  for record in SeqIO.parse(args.file1, "fasta"):
    contigs.append(Contig("contig_"+str(it+1), str(record.seq),args.taxon1))
    it+=1
  
  it=0
  for record in SeqIO.parse(args.file2, "fasta"):
    contigs.append(Contig("contig_"+str(it+1), str(record.seq),args.taxon2))
    it+=1
  
  with open("bin_results", "w") as output_file:
    bins = bin(contigs,args.sig,args.cuts,args.min_prob)
    output_file.write(disp_bin_stats(bins))
  return 0    

if __name__ == "__main__":
  main()
