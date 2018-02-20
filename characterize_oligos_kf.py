#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''

170712 - Gabriele Girelli
modified: Kata Ferenc (31/10/2017)
Project: COSMIC oligos

Aim:
	Generate and characterize oligos: GC content, melting temperature/

'''

# DEPENDENCIES =================================================================

import argparse
import math
import re

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(
	description = 'Characterize k-mers from fasta file without headers.'
)

# Add mandatory arguments
parser.add_argument('fastaInput', type = str, nargs = 1,
	help = 'Path to input fasta file.')
parser.add_argument('output', type = str, nargs = 1,
	help = 'Path to output fasta file.')

# Add arguments with default value
parser.add_argument('-k', '--oligolen', type = int, nargs = 1,
    metavar = 'oligoConc', help = """
	Oligo length in nt. Default: 30 nt
	""", default = [30])
parser.add_argument('-o', '--oligoconc', type = int, nargs = 1,
	metavar = 'oligoConc', help = """
	Oligo molar concentration. Default: 0.25e-6
	""", default = [0.25e-6])
parser.add_argument('-l', '--hplen', type = int, nargs = 1,
	metavar = 'hplen', help = """
	Homopolymer stretch length in nt. Default: 4 nt
	""", default = [4])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
k = args.oligolen[0]
fain = args.fastaInput[0]
out = args.output[0]
oligo_conc = args.oligoconc[0]
hp_len = args.hplen[0]

# FUNCTIONS ====================================================================

# create kmers out of input fasta file
def kmermaker(wseq, k):
    inum = 1
    for s in range(len(wseq)):
        yield str(inum)
        for i in range([len(j) for j in wseq[s]][0]-k+1):
            yield wseq[s][0][i:i+k]
        inum += 1
        
def rc(na, t):
	'''
	Args:
		na (string): nucleic acid sequence.
		t (string): nucleic acid type, either 'dna' or 'rna'.

	Return:
		string: reverse complement of na.
	'''

	# Identify type
	t = t.lower()

	# Select alphabet
	if t == 'dna':
		ab = ["ATCG", "TAGC"]
	elif t == 'rna':
		ab = ["AUCG", "UAGC"]
	else:
		print('ERROR: unknown na type.')
		return()

	rab = ab[1].strip().lower()
	ab = ab[0].strip().lower()

	# Check provided string
	na = na.lower()

	for c in na:
		if not c in ab:
			print('ERROR: provided string conflicts with selected alphabet.')
			return()

	# Calculate reverse
	r = na[::-1]

	# Calculate reverse complement
	rc = []
	for c in r:
		rc.append(rab[ab.index(c)])

	rc=''.join([str(c) for c in rc]).upper()

	return(rc)

def has_hp(seq, th):
	''''''
	c = 0
	for i in range(len(seq) - 1):
		if seq[i] == seq[i + 1]:
			c += 1
		else:
			c = 0

		if c >= th:
			return(True)
	return(False)


def characterize(seq, oligo_conc, hp_len):
	# DEFAULT ------------------------------------------------------------------
	
	# Gas constant
	R = 1.987 / 1000	# kcal / (K mol)

	# Table from Allawi&Santalucia, Biochemistry(36), 1997 - in 1 M NaCl [DNA]
	# dH0: kcal / mol
	# dS0: eu = cal / (K mol)
	# dG0: kcal / mol
	tt = {
		#				 dH0	 dS0 	 dG0
		'AA'		:	(-7.9,	-22.2,	-1.0),
		'TT'		:	(-7.9,	-22.2,	-1.0),
		'AT'		:	(-7.2,	-20.4,	-0.88),
		'TA'		:	(-7.2,	-21.3,	-0.58),
		'CA'		:	(-8.5,	-22.7,	-1.45),
		'TG'		:	(-8.5,	-22.7,	-1.45),
		'GT'		:	(-8.4,	-22.4,	-1.44),
		'AC'		:	(-8.4,	-22.4,	-1.44),
		'CT'		:	(-7.8,	-21.0,	-1.28),
		'AG'		:	(-7.8,	-21.0,	-1.28),
		'GA'		:	(-8.2,	-22.2,	-1.3),
		'TC'		:	(-8.2,	-22.2,	-1.3),
		'CG'		:	(-10.6,	-27.2,	-2.17),
		'GC'		:	(-9.8,	-24.4,	-2.24),
		'GG'		:	(-8.0,	-19.9,	-1.84),
		'CC'		:	(-8.0,	-19.9,	-1.84),
		'endG'		:	(.1,	-2.8,	.98),
		'endC'		:	(.1,	-2.8,	.98),
		'endA'		:	(2.3,	4.1,	1.03),
		'endT'		:	(2.3,	4.1,	1.03),
		'has_end'	:	True,
		'sym'		:	(2.3,	4.1,	1.03),
		'has_sym'	:	True
	}

	# GC content ---------------------------------------------------------------

	# Calculate GC content
	fgc = (seq.count('G') + seq.count('C')) / float(len(seq))
	
	# Melting temperature [1 M Na+] --------------------------------------------

	# Make NN couples
	couples = [seq[i:(i+2)] for i in range(len(seq) - 1)]

	# Calculate dH0(N-N)
	h = sum([tt[c][0] for c in couples])

	# Add initiation enthalpy
	if tt['has_end']:
		h += tt['end%s' % (seq[0],)][0]
		h += tt['end%s' % (seq[-1],)][0]

	# Correct enthalpy for symmetry
	if tt['has_sym'] and seq == rc(seq, 'dna'):
		h += tt['sym'][0]

	# Calculate dS0(N-N) in kcal / (K mol)
	s = sum([tt[c][1] for c in couples])

	# Add initiation enthalpy
	if tt['has_end']:
		s += tt['end%s' % (seq[0],)][1]
		s += tt['end%s' % (seq[-1],)][1]

	# Correct enthalpy for symmetry
	if tt['has_sym'] and seq == rc(seq, 'dna'):
		s += tt['sym'][1]

	s /= 1e3

	# Calculate melting temperature in Celsius
	Tm1 = h / (s + R * math.log(oligo_conc))

	# Homopolymer --------------------------------------------------------------
	hp = int(has_hp(seq, hp_len))

	# Output -------------------------------------------------------------------
	return((fgc, Tm1 - 273.15, hp))

# RUN ==========================================================================
# get the sequence from the fasta file
p = re.compile('>')
seq = []
with open(fain) as seq_in:
    for line in seq_in:
        if p.match(line):
            gnm = line.split('_')[0]
        else:
            seq.append([line.strip(),])


fout = open(out, 'w+')
i = 0
j = 0
onum = 0
inum = 0
for s in kmermaker(seq, k):
    if s.isnumeric():
        inum += 1
        onum = 0
    else:
        onum += 1
        i += 1
        if 0 != s.count('N'):
            j += 1
            continue
        (fgc, tm, hp) = characterize(s, oligo_conc, hp_len)
        fout.write("%s_%d:%d:%f:%f:%d\n%s\n" % (gnm,inum,onum,fgc, tm, hp,s))
fout.close()


print("Skipped %d (out of %d) sequences containing Ns." % (j, i+j,))

# END ==========================================================================

################################################################################