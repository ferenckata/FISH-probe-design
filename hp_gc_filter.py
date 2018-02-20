
# coding: utf-8

# In[5]:

### this script is for filtering out oligos with homopolymer stretch > 4nt
### and oligos with suboptimal GC content
import argparse
import re


# In[7]:

# Add script description
parser = argparse.ArgumentParser(
	description = 'Characterize k-mers from fasta file without headers.'
)
# Add mandatory arguments
parser.add_argument('fastaInput', type = str, nargs = 1,
	help = 'Path to input fasta file.')
parser.add_argument('output', type = str, nargs = 1,
	help = 'Path to output fasta file.')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
fin = args.fastaInput[0]
fout = args.output[0]


# In[9]:

p = re.compile('>')
m = 0
with open(fin,'r') as inp:
    with open(fout,'a') as oup:
        for line in inp:
            if p.match(line):
                l = line.split("\t")[0]
                if not int(l.split(':')[4]) == 1 and float(l.split(':')[2])<=0.6 and float(l.split(':')[2])>=0.4:
                    oup.write(line)
                    m = 1
            elif m == 1:
                oup.write(line)
                m = 0

