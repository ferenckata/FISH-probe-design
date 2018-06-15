#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.2.0
# Date: 20170628
# Project: RNA FISH oligo design
# Description:	filter BLASTN output based on:
# 					- homology (#PM/k),
# 					- number of OT, and
# 					- common (saturated) OTs.
# 
# Notes:
# 		PM: Perfect Match.
# 		OT: Off Target.
# 		The script should be run on single-gene outputs.
# 
# Changelog:
# 		1.0.0: first implementation.
# 		1.1.0: added argparser support.
# 		1.2.0: changed saturation OT filter.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(
	description = 'Filter BLASTN output.'
)

# Add mandatory arguments
parser.add_argument('blastInput', type = str, nargs = 1,
	help = 'Path to BLAST input fasta file.')
parser.add_argument('blastOutput', type = str, nargs = 1,
	help = 'Path to BLAST output with outfmt 6.')
parser.add_argument('geneTranscriptTable', type = str, nargs = 1,
	help = """Path to table with TRANSCRIPT_ID:GENE_SYMBOL
	tabulation-separated columns.""", default = [30])
parser.add_argument('output', type = str, nargs = 1,
	help = 'Path to output fasta file, after filtering.')

# Add arguments with default value
parser.add_argument('-k', type = int, nargs = 1,
	metavar = 'k', help = """Oligonucleotide length in nt.
	Default: 30.""", default = [30])
parser.add_argument('-t', '--homology-thr', type = float, nargs = 1,
	metavar = 'ht', help = """Threshold on maximum homology, as fraction of k.
	Accepts float values from 0 to 1.
	Default: .85""", default = [.85])
parser.add_argument('-e', '--gene-thr', type = int, nargs = 1,
	metavar = 'gt', help = """Threshold on the number of off-targets gene,
	for a single oligo.
	Default: 20""", default = [20])
parser.add_argument('-s', '--saturation-thr', type = int, nargs = 1,
	metavar = 'st', help = """
	Threshold on the number of oligos off-targeting a gene,
	for the selections of 'saturated' off-target genes. Oligos targeting
	a saturated off-target are filtered out.
	Default: 5""", default = [5])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
blast_input = args.blastInput[0]
blast_output = args.blastOutput[0]
gene_transcript_table = args.geneTranscriptTable[0]
output_file = args.output[0]
k = args.k[0]
homThr = args.homology_thr[0]
gene_ot_thr = args.gene_thr[0]
oligo_ot_thr = args.saturation_thr[0]

# Log to screen the settings
print("""
Settings:
              BLASTN input : %s
             BLASTN output : %s
     Gene-Transcript table : %s
               Output file : %s
                         K : %d
        Homology threshold : %f
       Gene Off-Target thr : %d
 Off-Target saturation thr : %d

""" % (blast_input, blast_output, gene_transcript_table, output_file,
	k, homThr, gene_ot_thr, oligo_ot_thr))


# FUNCTIONS ====================================================================

# RUN ==========================================================================
print("Run:")

# Build TRANSCRIPT_ID:GENE_SYMBOL dictionary -----------------------------------
print(" · Building TRANSCRIPT_ID:GENE_SYMBOL dictionary...")

# Initialize empty TRANSCRIPT_ID:GENE_SYMBOL dictionary
trn_gene_dict = {}

# Read table line by line
with open(gene_transcript_table) as gttf:
	for line in gttf:

		# Split every line in two fields
		tmp = line.strip().split('\t')

		# Store key,value couples in the dictionary
		trn_gene_dict[tmp[0]] = tmp[1]

# Work on BLASTN output --------------------------------------------------------

# Filter based on max homology
# ----------------------------
print(" · Filtering based on maximum homology...")

# Initalize empty maximum homology dictionary
max_homology = {}

# Read table line by line
with open(blast_output) as bof:
	for line in bof:

		# Split every line by column
		tmp = line.strip().split('\t')

		# Identify oligomer ID
		OID = tmp[0].split(':')[1]

		# Identify target gene
		target = tmp[0].split('_')[0]

		# Identify transcript ID
		transcript = tmp[1].split('.')[0]

		# Calculate homology
		homology = (int(tmp[3]) - int(tmp[4])) / float(k)

		# If it is an off-target
		if target.upper() != trn_gene_dict[transcript].upper():
			if not OID in max_homology.keys():
				max_homology[OID] = homology
			else:
				# With homology higher than the threshold
				if homology > max_homology[OID]:
					max_homology[OID] = homology
		elif  not OID in max_homology.keys():
			max_homology[OID] = 0

# Identify oligos that pass the threshold
pass_homology = []
for OID in max_homology.keys():
	if max_homology[OID] < homThr:
		pass_homology.append(OID)

# Log
print(" >>> %d oligos do not have any off-targets." % (len(pass_homology),))
print(" >>> %d oligos have off-targets." % (len(max_homology)-len(pass_homology),))
print(" >>> Saving off-target free oligos. Analyzing further the rest.")

# Check off target location
# -------------------------
print(" · Identifying off-target locations...")

# Initialize empty off-targets dictionary
# OligoID:OTgene:OTtranscript
ot_dict = {}

# Read table line by line
with open(blast_output) as bof:
	for line in bof:

		# Split every line by column
		tmp = line.strip().split('\t')

		# Identify oligomer ID
		OID = tmp[0].split(':')[1]

		# Work only on oligos with off-targets
		# that do not pass the homology filter
		if OID in pass_homology:
			continue

		# Identify target gene
		target = tmp[0].split('_')[0]

		# Identify transcript ID
		ot_trans = tmp[1].split('.')[0]

		# Calculate homology
		homology = (int(tmp[3]) - int(tmp[4])) / float(k)

		# Identify the off-target
		ot_gene = trn_gene_dict[ot_trans]

		# If it is an off-target
		if target.upper() != ot_gene.upper():

			# With homology higher than the threshold
			if homology >= homThr:

				if not OID in ot_dict.keys():
					ot_dict[OID] = {}
					ot_dict[OID][ot_gene] = {}
					ot_dict[OID][ot_gene][ot_trans] = 1
				else:
					if not ot_gene in ot_dict[OID].keys():
						ot_dict[OID][ot_gene] = {}
						ot_dict[OID][ot_gene][ot_trans] = 1
					else:
						if not ot_trans in ot_dict[OID][ot_gene].keys():
							ot_dict[OID][ot_gene][ot_trans] = 1
						else:
							ot_dict[OID][ot_gene][ot_trans] += 1

# Calculate number of off-targets
# -------------------------------
print(" · Filtering based on number of off-target genes...")

ot_gene_count = {}
for OID in ot_dict.keys():
	ot_gene_count[OID] = len(ot_dict[OID].keys())
if len(ot_gene_count.values()) != 0:
	print(""" >>> OT counts summary:
        	     min : %f
      		1st Quart. : %f
          		median : %f
            		mean : %f
      		2nd Quart. : %f
             		max : %f
 	>>> Current threshold at the %d-ith percentile.""" % (
		np.percentile(ot_gene_count.values(), 0),
		np.percentile(ot_gene_count.values(), 25),
		np.percentile(ot_gene_count.values(), 50),
		np.mean(ot_gene_count.values()),
		np.percentile(ot_gene_count.values(), 75),
		np.percentile(ot_gene_count.values(), 100),
		int(sum(np.array(ot_gene_count.values()) < gene_ot_thr)
			/ float(len(ot_gene_count.values())) * 100)
	))

	# Filter based on number of OTs per oligo
	pass_oligo_ot_count = []
	for OID in ot_gene_count.keys():
		if ot_gene_count[OID] < gene_ot_thr:
			pass_oligo_ot_count.append(OID)
	print(" >>> %d oligos pass the OT count filter" % (len(pass_oligo_ot_count),))

	# Calculate number of common off-targets
	# --------------------------------------

	pass_gene_ot_count = pass_oligo_ot_count

	if 1 == len(pass_oligo_ot_count):
		print(" · Skipping saturated off-target transcript filter...")
	else:
		print(" · Filtering based on saturated off-target transcripts...")

		gene_ot_count = {}
		for OID in pass_oligo_ot_count:
			for gene in ot_dict[OID].keys():
				if not gene in gene_ot_count.keys():
					gene_ot_count[gene] = ot_dict[OID][gene]
				else:
					for trans in ot_dict[OID][gene].keys():
						if not trans in gene_ot_count[gene].keys():
							gene_ot_count[gene][trans] = ot_dict[OID][gene][trans]
						else:
							gene_ot_count[gene][trans] += ot_dict[OID][gene][trans]

		# Filter based on number of OTs per transcript
		for gene in gene_ot_count.keys():
			for trans in gene_ot_count[gene].keys():
				# Discard oligos that off-target a transcript
				# shared by too many oligos
				if gene_ot_count[gene][trans] >= oligo_ot_thr:
					for OID in pass_gene_ot_count:
						if gene in ot_dict[OID].keys():
							if trans in ot_dict[OID][gene].keys():
								pass_gene_ot_count.remove(OID)

		print(" >>> %d oligos pass the saturation OT filter"
			% (len(pass_gene_ot_count),))
else:
	pass_oligo_ot_count = []
	pass_gene_ot_count = []
# Merge list of filtered oligos
# -----------------------------

output_list = pass_homology
output_list.extend(pass_gene_ot_count)
output_list.sort()

# Log
print(" · %d oligos can be used for further screening." % (len(output_list),))

# Prepare output ---------------------------------------------------------------
print(" · Generating FASTA output...")

# Prepare output string
s = ''

# Variable to keep the non-header lines
keep = False

# Read Fasta line by line
with open(blast_input) as bif:
	for line in bif:

		# If header line
		if '>' == line[0]:
			# Check if the sequence should be kept
			if line[1:].split(':')[1] in output_list:
				keep = True
				s += line
			else:
				keep = False
		# Otherwise keep it if it passed the filters
		elif keep:
			s += line

# Write output
f = open(output_file, 'w')
f.write(s)
f.close()

# END ==========================================================================

print("""
DONE!
""")

################################################################################
