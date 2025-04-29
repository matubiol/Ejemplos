import os
import sys
import subprocess
import multiprocessing
import argparse
import re
import random
import datetime
import pandas as pd

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq


def main():

	# Primer_check script version
	primerCheck_version = 2.1

	print(f"Running PrimerCheck Version: {primerCheck_version}")

	### Input arguments
	options = parseArguments()
	
	# Get query sequence
	subprocess.call(f"qiime rescript get-ncbi-data --p-query '({options.search_term})' --output-dir {options.primer_name}", shell = True)

	# Export fasta file
	subprocess.call(f"qiime tools export --input-path {options.primer_name}/sequences.qza --output-path {options.primer_name}/", shell = True)
	subprocess.call(f"qiime tools export --input-path {options.primer_name}/taxonomy.qza --output-path {options.primer_name}/", shell = True)

	# Create dictionary of taxonomies
	taxon_dict = {}
	with open(f"{options.primer_name}/taxonomy.tsv") as taxonomy_in:
		next(taxonomy_in)
		for line in taxonomy_in:
			line = line.strip()
			linesplit = line.split("\t")
			species = linesplit[1].split("; ")[-1]
			taxon_dict.setdefault(species, [])
			taxon_dict[species].append(linesplit[0])

	print("From the below list, select the taxa required and input below as a comma delimited list, with no spaces after commas.")
	print("E.g s__Tomato brown rugose fruit virus,s__Tomato brown rugose fruit virus - Palestinian isolate,s__Tomato brown rugose fruit virus-israeli\n")
	for taxa in taxon_dict:
		print(taxa)

	taxa = input("\nPlease input comma delimited list, with no spaces after commas:\n")
	taxa = taxa.split(",")

	# Get accession numbers
	taxa_list = []
	for taxon in taxa:
		taxa_list.append(taxon_dict[taxon])

	# Flatten list
	taxa_list = [item for sublist in taxa_list for item in sublist]

	# Select accession numbers
	with open(f"{options.primer_name}/target_sequences.fasta", "w") as target_seqs:
		for seq_record in SeqIO.parse(open(f"{options.primer_name}/dna-sequences.fasta", mode = "r"), "fasta"):
			for taxon in taxa_list:
				if taxon == seq_record.id:
					SeqIO.write(seq_record, target_seqs, "fasta")

	# Align sequences
	subprocess.call(f"mafft --adjustdirection --thread -1 {options.primer_name}/target_sequences.fasta > {options.primer_name}/target_sequences.mafft.fasta", shell = True)

	# Extract primer and probe sequences
	primer_probe_dict = {}
	with open(options.assays) as assays_in:
		next(assays_in)
		for line in assays_in:
			line = line.strip()
			linesplit = line.split("\t")
			assay = linesplit[0]
			fprimer =linesplit[1]
			rprimer =linesplit[2]
			probe =linesplit[3]
			primer_probe_dict.setdefault(assay, [])
			primer_probe_dict[assay].append(linesplit[1])
			primer_probe_dict[assay].append(linesplit[2])
			primer_probe_dict[assay].append(linesplit[3])

	print("From the below list, select the assay required.")
	print("E.g Tomato brown rugose fruit virus\n")
	for assay in primer_probe_dict:
		print(assay)

	assay = input("\nPlease input assay:\n")
	
	# Create a primer and probe fasta file
	with open(f"{options.primer_name}/{options.primer_name}_primers.fasta", "w") as primer_file:
		primer_file.write(f">F_primer\n{primer_probe_dict[assay][0]}\n>R_primer\n{Seq(str(primer_probe_dict[assay][1])).reverse_complement()}\n>Probe\n{primer_probe_dict[assay][2]}\n")


	# # Add primers to alignment
	subprocess.call(f"mafft --multipair --addfragments {options.primer_name}/{options.primer_name}_primers.fasta --keeplength --thread -1 --mapout {options.primer_name}/target_sequences.mafft.fasta > {options.primer_name}/{options.primer_name}_primers.mafft.fasta", shell = True)
	
	# Redundant: Reorient any sequences that align better when in the reverse compliment
	# subprocess.call(f"mafft --adjustdirection {options.primer_name}/{options.primer_name}_primers.mafft.fasta > {options.primer_name}/{options.primer_name}_primers_adjusted.mafft.fasta", shell = True)

	# Create primer probe position dictionary
	map_dict = {}
	with open(f"{options.primer_name}/{options.primer_name}_primers.fasta.map") as map_in:
		primer_probe_id = ""
		for line in map_in:
			line = line.strip()
			if line.startswith(">"):
				primer_probe_id = line[1:]
				map_dict.setdefault(primer_probe_id, [])
			elif not line.startswith("#"):
				pos = line.split(" ")[2]
				map_dict[primer_probe_id].append(pos)

	# Read in alignment
	alignment = AlignIO.read(f"{options.primer_name}/{options.primer_name}_primers.mafft.fasta", "fasta")

	# Create slices of alignment over the primer/probe regions to export
	fprimer_alignment = alignment[:, int(map_dict["F_primer"][0]) - 1:int(map_dict["F_primer"][-1])]
	probe_alignment = alignment[:, int(map_dict["Probe"][0]) - 1:int(map_dict["Probe"][-1])]
	rprimer_alignment = alignment[:, int(map_dict["R_primer"][0]) - 1:int(map_dict["R_primer"][-1])]
	# final_alignment = fprimer_alignment + probe_alignment + rprimer_alignment
	# AlignIO.write(final_alignment, f"{options.primer_name}/{options.primer_name}_FINAL.msa.fasta", "fasta")

	output_frequency(fprimer_alignment, "Forward", options.primer_name)
	output_frequency(probe_alignment, "Probe", options.primer_name)
	output_frequency(rprimer_alignment, "Reverse", options.primer_name)

################################################################################################
def output_frequency(alignment, name, primer_name):
	frequency_dict = {}
	for position in range(len(alignment[0].seq) ):
		pos_x = alignment[: -3, position]
		nt_dict = {"A": 0, "T": 0, "C": 0, "G": 0, "-": 0, "Ambig": 0}
		for pos in pos_x:
			if pos.upper() in nt_dict:
				nt_dict[pos.upper()] += 1
			else:
				nt_dict["Ambig"] += 1
		frequency_dict[position] = nt_dict


	df = pd.DataFrame.from_dict(frequency_dict)
	df.to_csv(f"{primer_name}/{primer_name}_{name}_frequency.tsv", sep = "\t", index_label = name)

#################################### Get Arguments Function ####################################
def parseArguments():
	parser = argparse.ArgumentParser(description = "Runs the PrimerCheck pipeline.")

	# Main arguments
	parser.add_argument("--primer_name", help = "The name of the primer to check.", required = True)
	parser.add_argument("--search_term", help = "The term NCBI is searched for.", required = True)
	parser.add_argument("--assays", help = "The assay file.", required = True)

	return parser.parse_args()

################################################################################################

############################################################################################################################## Functions End ###############################################################################################################################

if __name__ == '__main__':
	main()
