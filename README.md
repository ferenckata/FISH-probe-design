# FISH-probe-design
For RNA or DNA FISH probes of any gene/transcript of interest


### Retrieving the sequences
For cDNA or CDS
- The script is the `cds_cdna.r`
- Retrieving the CDS or CDNA sequences from mm10 or hg38
- The input is a CSV file with the Gene Symbol, Ensembl Gene ID or Transcript ID (header must include these names starting with capital letter)
- Output is a fasta file

For intronic seq
- The script is the `intron_by_coord.r`
- Currently it's getting the coordinates of the exons from either Mus musculus mm10 or Homo sapiens hg38 versions (the user defines the species: `h` for human, `m` for mouse)
- It calculates the shortest distance between consecutive exons irrespectively of which transcript the exons are from
- The output is a fasta file with the sequences. The header contains the Ensembl gene ID and the position of the sequence
- Additionally it produces a bed file with the coordinates that is compatible with UCSC Genome Browser

### Creating all possible k-mers and characterize them
- The code is called `characterize_oligos.py` (in python 3, modified from [Gabriele Girelli](https://github.com/ggirelli)'s code)
- The mandatory inputs are the fasta file output of the first step and the name of a fasta output. The optional inputs are in the help of the script.
- The header of the output contains the gene name, the intron number, the oligo number in the intron, the GC content, the melting temperature and the presence of homopolymers


### Plot and filter based on GC and HP
Plotter
- The plotter is an R script that is runnable from the terminal: `oligo_stat.R`
- The input file or folder containing the input files should be defined. (not implemented yet)
- The output is two graphs: one about GC content, the other is about melting temperature. If there is only one gene of interest, the plots are in oligo resolution, otherwise they show the median for each gene.
- The log shows some statistics about the GC content and the melting temperature.
Filter
- The filter is a python script: `hp_gc_filter.py`
- The input is the same fasta file as for the plotteR and a name for the output.
- The output is also a fasta file containing only those oligos that passed the hp and GC filters.

### BLAST and homology filter
BLAST
- We are using local blast from the terminal, called with `blastn` command.
- The parameters are the following:
  - `db` (Followed by the indexed database made with `makeblastdb` command. In our case it is a joined database of ncRNA and cDNA.)
  - `query $file` (The input fasta file)
  - `outfmt 6` (This is the default, further information available online (eg [here](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6)))
  - `out` (The name of the output tsv file)
  - `word_size 10` (This is the length of the searching window, the smaller the following number is, the more precise and slow the search is, however it shouldn't be lower than 7bp. We set to 10 to be sensitive enough for our 30nt long sequences)
  - `num_threads 4` (The number of cores to use. The more we use, the faster the search is. On p700 we should use no more than 4, on iMac one can use up to 6 if nothing else is running.)
  - `penalty -2 -reward 3` (These are the scores for mismatch and perfect match respectively. The BLAST score is calculated by summing the scores. The results are sorted in descending order by the BLAST score.)
  - `gapopen 100 -gapextend 100` (These are scores for opening a gap and extending it. We set both quite big, so any gap in the alignment will have a big negative effect on the BLAST score of the oligo.)
  - `evalue 1000` (The expected value is set here. It shows how likely is that our sequence is found by chance rather than true matching. Thus it is the threshold on how many entries we want to have in our list. In our case it is set high, because we have short oligos that have low statistical significance.)
  - `strand plus` (Our oligos are created from the genomic sequence. At the least step they will be reverse complemented. We want to keep only those that match to a specific location, thus we only compare with one strand.)
		
Homology filter (modified from the original code of [Gabriele Girelli](https://github.com/ggirelli))
- This step processes the result of BLAST
- The script is `blast_filter.py` (currently `blast_filter_for_introns.py` and `blast_filter_for_cdna_cds.py` )
- It requires the blast input fasta file, the blast output tsv file, a gene-transcript table file that contains all transcripts and the corresponding genes of a certain species
- Optional arguments are the thresholds of each filter and the oligo length
- The output is a filtered fasta file
- The filter is done in three steps:
  1. Homology filtering
     - If an oligo is less homologous to any other genes than its origin, than a set threshold it passes the homology filter
     - Default threshold is 0.85, calculated as (# of perfect matches - # of mismatches) / length of the oligo
  2. Off-taget filtering
     - If an oligo has a higher similarity to another gene than the threshold, but the number of such off-target genes is lower than the threshold, it passes the off-target filter
     - Default threshold is 20
  3. Saturation filter
     - If an oligo fails both the homology filter and the off-target number filter, it has one more chance. If the total number of oligos targeting an off-target transcript (not gene!) is lower than the trheshold, these oligos passes the saturation filter.
     - Default threshold is 5

### Selection
- Chose a set of oligos from the overlapping set left after all filtering steps
- The latest version choses those oligos that have the GC content closest to 50% among the overlapping ones, IF it doesn't decrease the total number of oligos
  - Ie. Two oligos of 60-40% GC is chosen over one with 50% GC
  - It could easily be changed to be the other way around, but not implemented yet
- The script is called selection.py
- The input is the folder with only the blastfiltered fasta files, the output is the a "...probes.tsv"  with the properties of groups of oligos, and a "...chosen.fa" with the chosen oligos
- _after all filtering steps there are regions with and without oligos. "group of oligos" is a name for a region with (often overlapping) oligos, the groups are separated by empty regions_

### Barcoding
- The script is the `barcode_appender.r` , not yet runnable from the terminal
- The input is the barcode list, the fasta file with the oligos to be used
- The output is a list of barcoded oligos in the format seen below and a table about the barcode combination

```
Fwbcsn = forward barcode short name, fwbcln = forward barcode long name
Rbcsn= reverse barcode short name, rbcln = reverse barcode long name
	
probeID_fwbcsn_rbcsn <tab> 5'-3' sequence <tab> oligocoordinate_fwbcln_rbcln
```

