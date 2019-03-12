#!/usr/bin/env python3

'''
Author: Mitchell Rezzonico

Description: This script parses FASTA files and plots protein binding
motifs on a gene. Introns are represented as a thin horizantal line and
exons are represented as a horizantal black box on the line.
'''

################################
####### Import Libraries #######
################################

import cairo
import math
import argparse
import re
import random
import numpy as np

################################
###### Argparse Statements #####
################################

def get_arguments():
	parser = argparse.ArgumentParser(description="Instructions to navigate this program")
	parser.add_argument("-f", "--file", help="File path to special FASTA. Format must be 'intron', 'exon', 'intron'.", required=True, type=str)
	parser.add_argument("-m", "--motif", help="File path to motifs of interest.", required=True, type=str)
	return parser.parse_args()

################################
###### Argparse Variables ######
################################

args = get_arguments()
fasta_file = args.file # fasta file
motif_file = args.motif # motif file

################################
######### Open Files ###########
################################

motif_original = [] # original motifs
gene_dict = {} # key: header, value: gene (unedited)

with open(motif_file, "rt") as mot, open(fasta_file, "rt") as fh:
	while True:
		motif = mot.readline().strip()
		if motif == "":
			break

		motif_original.append(motif)

	while True:
		full_fasta = fh.readline().strip()
		if full_fasta == "":
			break
		if full_fasta.startswith(">"):
			header = full_fasta[1:] # get rid of >
			gene_dict[header] = ""
		else:
			gene_dict[header] += full_fasta

################################
########## Functions ###########
################################

def ColorGenerator() -> (float, float, float):
	'''
	Returns 3 random integers between 0-1 to be used for the
	rgb color pallete. Will be called by DrawFig.
	'''
	# generate random numbers between 0-1
	np.random.random(1)[0]
	#r,g,b = np.random.random((3,1))
	r = np.random.random((1,1))
	g = np.random.random((1,1))
	b = np.random.random((1,1))

	# make sure colors are always distinguishable
	diff_rg = abs(r-g)
	diff_rb = abs(r-b)
	diff_gb = abs(g-b)

	if diff_rg <= 0.3:
		r = np.random.random((1,1))
	if diff_rb <= 0.3:
		b = np.random.random((1,1))
	if diff_gb <= 0.3:
		g = np.random.random((1,1))

	return(float(r),float(g),float(b))


def ProcessMotifs(motifs:list, genes:dict) -> dict:
	'''
	motifs: List of raw motif sequences
	genes: Dictionary with keys as headers and values as gene sequences

	Returns dictionary with key as gene header stripped of ">" and value
	as another dictionary with the key as the motif and the value as a list
	of all start sites for that motif.

	Ambiguity in nucleotides:
	R.................A or G
	Y.................C or T
	S.................G or C
	W.................A or T
	K.................G or T
	M.................A or C
	B.................C or G or T
	D.................A or G or T
	H.................A or C or T
	V.................A or C or G
	N.................any base
	U.................T
	'''
	motif_dict = {} # key: original motif, value: spans of motifs
	for motif in motifs:
		# iterate through raw motif sequences
		motif = motif.upper() # convert motif to uppercase to search
		motif_dict[motif] = ""

		# single search string that, if required, would contain every ambiguous element for that NT
		motif_temp = motif.replace("R","[AG]")
		motif_temp = motif_temp.replace("Y","[CT]")
		motif_temp = motif_temp.replace("S","[GC]")
		motif_temp = motif_temp.replace("W","[AT]")
		motif_temp = motif_temp.replace("K","[GT]")
		motif_temp = motif_temp.replace("M","[AC]")
		motif_temp = motif_temp.replace("B","[CGT]")
		motif_temp = motif_temp.replace("D","[AGT]")
		motif_temp = motif_temp.replace("H","[ACT]")
		motif_temp = motif_temp.replace("V","[ACG]")
		motif_temp = motif_temp.replace("N","[ATCG]")
		motif_temp = motif_temp.replace("U","[T]")

		motif_dict[motif] = motif_temp

	pos_dict_gene = {} # key: gene header, value: dictionary with key as motif and value as a list of all start sites for that motif
	for gene in genes:
		# iterate through dictionary to grab gene sequences
		pos_dict_gene[gene] = {} # key: header, value: dict with keys=motifs, values=span of motifs
		for m in motif_dict:
			# iterate through uppercase motifs
			pos_dict_gene[gene][m] = [] # if the motif isn't in the gene the list will be empty
			search_genes = genes[gene].upper() # uppercase genes for regex to search

			for match in re.finditer(motif_dict[m], search_genes): # "thing being searched", "thing you're searching"
				# search for motifs in the genes
				start = match.start() # capture start position of motif
				pos_dict_gene[gene][m].append(start)

	return(pos_dict_gene) # feed into DrawFig


def FindExons(fasta_gene:dict) -> dict:
	'''
	fasta_gene: Dictionary with gene header as key and gene sequence as value.

	Exons are distinguished by capitalization. Identifies the start position
	and length of the exon in the FASTA file. Returns a dictionary with the
	gene header as the key the position of the exon as the value.
	'''
	exon_format = "[A-Z]+"

	exons = {} # key: gene header, value: tuple of start position and length of exon
	for NT in fasta_gene:
		# iterate through gene sequence
		exon_pos = ()
		for match in re.finditer(exon_format, fasta_gene[NT]): # "thing being searched", "thing you're searching"
			# search for exons
			start = match.start() # capture start position of exon
			end = match.end() # capture end position of exon
			diff = end-start # calculate length of the exon
			exon_pos = (start, diff)

		exons[NT] = exon_pos

	return(exons) # feed into DrawFig

def DrawFig(Motifs:dict, Exons:dict, Genes:dict, Legend_Motifs:list):
	'''
	Motifs: Dictionary with key as gene header and value as dictionary with
	key as motif and value as list of start positions.
	Exons: Dictionary with key as gene header and value as tuple of exon
	start position and length.
	Genes: Dictionary with key as gene header and value as the gene sequence.
	Legend_Motifs: List with motifs as originally listed in motif file.

	Draws one large .svg figure for all genes in the FASTA. Figure name will be
	name of the FASTA file with ".svg" appended to the end.
	'''
	color_list = []
	for m in range(len(Legend_Motifs)):
		r,g,b = ColorGenerator()
		color_list.append((r,g,b))

	# legend dimension start points
	legend_box_x = 50
	legend_box_y = 150
	motif_box_y = 175
	motif_text_y = 185
	legend_x_text = 80
	legend_y_text = 55

	longest_gene = 0
	for longest in Genes:
		# iterate through gene sequences to find longest in length
		current_gene = Genes[longest]
		if (len(current_gene)) >= longest_gene:
			longest_gene = len(current_gene)

	width = longest_gene * 1.5 # scale width to the longest gene in the FASTA file
	height = (len(Genes) + 500)*2 # scale height to the number of the sequences in FASTA file
	surface = cairo.SVGSurface(fasta_file + ".svg", width, height)
	context = cairo.Context(surface)

	longest_motif = 0
	for motif in Legend_Motifs:
		# iterate through motifs to find the longest in length
		current_motif = motif
		if len(current_motif) >= longest_motif:
			longest_motif = len(current_motif)

	# dimensions for legend
	context.set_line_width(2)
	context.move_to(50,legend_box_y)
	context.rectangle(50, legend_box_y, (longest_motif*2)+150, (len(Legend_Motifs)+100))
	context.set_source_rgb(0,0,0)
	context.stroke()

	# set up figure dimensions
	intron_y_coord = 300 + (len(Legend_Motifs)+30) # first intron y coord (horizontal line)
	exon_y_coord = 293 + (len(Legend_Motifs)+30) # first intron y coord (vertical black box)
	motif_y_coord = 293 + (len(Legend_Motifs)+30) # first motif y coord (vertical colored box)
	header_y_coord = 280 + (len(Legend_Motifs)+30) # first gene header y coord

	legend = False
	for sequence in Genes:
		# iterate through each gene sequence
		if legend == False:
			color_index = 0
			for m in Motifs[sequence]:
				#iterate through each motif

				# add legend title
				context.move_to(55, 165)
				context.set_source_rgb(0,0,0)
				context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
				context.show_text("Motif Legend")
				context.stroke()

				# add box with motif color
				context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
				context.set_line_width(2)
				context.move_to(55, motif_box_y)
				context.set_source_rgb(color_list[color_index][0], color_list[color_index][1], color_list[color_index][2])
				context.rectangle(55, motif_box_y, 8, 12)
				context.fill()
				context.stroke()

				# add original motif name next to box
				context.move_to(65, motif_text_y)
				context.show_text(Legend_Motifs[color_index])
				context.stroke()
				color_index += 1
				motif_box_y += 15
				motif_text_y += 15
			legend = True

		#add header above gene
		context.move_to(50, header_y_coord)
		context.set_source_rgb(0,0,0)
		context.show_text(sequence)
		header_y_coord += 80 # shift down for the next header

		# draw introns
		context.set_line_width(5)
		context.move_to(50,intron_y_coord) # start position for the intron, vertical coordinate
		context.line_to(len(Genes[sequence])+ 50, intron_y_coord) # length of the gene sequence, vertical coordinate
		context.set_source_rgb(0,0,0)
		context.stroke()
		intron_y_coord += 80 # shift down for next intron

		# draw exons
		context.rectangle(Exons[sequence][0]+ 50, exon_y_coord, Exons[sequence][1], 14)
		context.fill()
		exon_y_coord += 80 # shift down for next exon

		# draw motifs
		color_index = 0
		for mot in Motifs[sequence]:
			# iterate through each motif
			context.set_source_rgb(color_list[color_index][0], color_list[color_index][1], color_list[color_index][2])
			color_index += 1
			for start in Motifs[sequence][mot]:
				# iterate through each start position
				context.rectangle(start+50, motif_y_coord, len(Motifs[sequence]), 14)
				context.fill()
		motif_y_coord += 80

	# write out final figure
	surface.finish()

DrawFig(ProcessMotifs(motif_original, gene_dict), FindExons(gene_dict), gene_dict, motif_original)
