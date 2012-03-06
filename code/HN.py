#Copyright 2010,2011 Brown University, Providence, RI.

                         #All Rights Reserved

#Permission to use, copy, modify, and distribute this software and its
#documentation for any purpose other than its incorporation into a
#commercial product is hereby granted without fee, provided that the
#above copyright notice appear in all copies and that both that
#copyright notice and this permission notice appear in supporting
#documentation, and that the name of Brown University not be used in
#advertising or publicity pertaining to distribution of the software
#without specific, written prior permission.

#BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
#INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
#PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
#ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
#WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
#ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
#OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#!/usr/bin/env python

import sys
import os
import random

def generate_num_alter_per_gene( mut_file, cna_file, out_file, tested_genes, gene_to_index ):
	#given name of files with mutations and cna, write on out_file the number of alterations per gene
	# mut_file contains name_sampleTABmutated_genes (mutated_genes are TAB separated). Even samples with no mutations should be reported
	# cna_file contains cna table: the gene name is in token 0, the samplesID/annotations starts at token 3
	# out_file: output, with geneTABnum_alterations
	# tested_genes: set of tested genes
	# gene_to_index: mapping between gene_name and gene_index in the network
	sequenced_samplesID = set()
	gene_altered_samples = dict()
	mut_f = open( mut_file , 'r' )
	for line in mut_f:
		v = line.split("\t")
		sample = v[0].strip("\n")
		sequenced_samplesID.add( sample )
		for i in range( len(v)-1 ):
			geneID = v[i+1].strip("\n")
			if geneID not in gene_altered_samples:
				gene_altered_samples[geneID] = set()
			gene_altered_samples[geneID].add( sample )
	mut_f.close()
	#rule for CNA: gene is altered IFF >= 90% of aberrated samples have the same sign
	file_cn = open( cna_file, 'r' )
	l_count = 0
	cna_samples = list()
	gene_alt_pos = dict()
	gene_alt_neg = dict()
	for line in file_cn:
		if l_count == 0:
			v = line.split('\t')
			for i in range( len(v) ):
				if i >2:
					r = v[i].strip('\",\n').split('-')
#					print v[i]
					sampleID = r[0]+'-'+r[1]+'-'+r[2]
					cna_samples.append( sampleID )
		else:
			v = line.split('\t')
			geneID = v[0]
			if geneID not in gene_alt_pos:
				gene_alt_pos[geneID] = set()
			if geneID not in gene_alt_neg:
				gene_alt_neg[geneID] = set()
			for i in range( len(v) ):
				if i>2 and v[i].strip() != 'NaN' and v[i].strip() != 'NA':
					if float(v[i].rstrip('\n')) > 1 :
						if cna_samples[i-3] in sequenced_samplesID:
							gene_alt_pos[ geneID ].add( cna_samples[i-3] )
					if float(v[i].rstrip('\n')) < -1 :
						if cna_samples[i-3] in sequenced_samplesID:
							gene_alt_neg[ geneID ].add( cna_samples[i-3] )
	#				if geneID not in alt_genes:
	#					alt_genes.append( geneID )
		l_count += 1
	
	file_cn.close()

	genes_sort = gene_alt_pos.keys()

	thresh = 0.9
	
	tot_cna_real = 0
	
	for i in range(len(genes_sort)):
		geneID = genes_sort[i]
		num_p = len(gene_alt_pos[geneID])
		num_n = len(gene_alt_neg[geneID])
		tot = num_p+num_n
		if tot > 0:
			if num_p >= num_n:
				ratio = float(num_p)/float(tot)
				if ratio >= thresh:
					if geneID not in gene_altered_samples:
						gene_altered_samples[ geneID ] = set()
					gene_altered_samples[ geneID ].update( gene_alt_pos[geneID] )
					if geneID in gene_to_index:
						if gene_to_index[ geneID ] in tested_genes:
							tot_cna_real += len(gene_alt_pos[geneID])
					
			else :
				ratio = float(num_n)/float(tot)
				if ratio >= thresh:
					if geneID not in gene_altered_samples:
						gene_altered_samples[ geneID ] = set()
					gene_altered_samples[ geneID ].update( gene_alt_neg[geneID] )
					if geneID in gene_to_index:
						if gene_to_index[ geneID ] in tested_genes:
							tot_cna_real += len(gene_alt_neg[geneID])
	
#	print "Total CNA real data: " +str(tot_cna_real)
	
	genes = gene_altered_samples.keys()
	
	out_f = open(out_file,'w')

	coutnn = 0
	num_alll = 0
	for i in range( len( genes ) ):
		if genes[i] in gene_to_index:
			coutnn += 1
			if gene_to_index[genes[i]] in tested_genes:
				out_f.write( str(gene_to_index[genes[i]]) + " 1 " + str(len( gene_altered_samples[genes[i]] ) )+"\n" )
		if len( gene_altered_samples[genes[i]] ) > 0:
			num_alll += 1
#		else:
#			print genes[i]+" is altered but not tested"
	out_f.close()
#	print "Number of altered genes: "+ str(num_alll)
	return sequenced_samplesID
	
def generateLengths( genes_length_file, tested_genes, index_to_gene, out_file ):

	#generate the lenght of the tested genes
	
	genes_length_f = open( genes_length_file , 'r' )
	
	out_f = open( out_file , 'w' )
	
	gene_to_length = dict()
	
	for line in genes_length_f:
		v = line.split()
		gene = v[0]
		length = int(v[1])
		gene_to_length[gene]=length

	for gene_id in tested_genes:
		gene = index_to_gene[ gene_id ]
		if gene not in gene_to_length:
			print "ERROR: " + gene + " is tested but has no length reported! Terminating."
			exit(0)
		gene_length = gene_to_length[ gene ]
		out_f.write(str(gene_id)+" 1 " + str(gene_length)+"\n")
	
	out_f.close()
	genes_length_f.close()
	return gene_to_length

def generateRandomMutations( tested_genes, index_to_gene, genes_length, BMR ):
	mut_genes = list()
	for gene in tested_genes:
		length = genes_length [ index_to_gene[gene] ]
		if length <= 0:
			print "Nonpositive length for " + index_to_gene[gene] + " - TERMINATING"
			exit(0)
		prob = 1 - pow( 1 - BMR, length)
		val = random.random()
		if val <= prob:
			mut_genes.append( gene )
#	print "Num. mutated genes: " + str( len( mut_genes ) )
	return mut_genes

def buildGeneOrder( ):
	config_f = open( "code/HN.config" , 'r')
	line_c = 0
	for line in config_f:
		line_c += 1
		if line_c == 5:
			order_file = line.strip("\n")
	config_f.close()
	order_f = open( order_file, 'r' )
	unit_counter = 0
	order_genes = dict()
	for line in order_f:
		v = line.strip("\n")
		order_genes[ unit_counter ] = list(v.split("\t"))
		unit_counter += 1
#	print "total units: " + str(unit_counter)
	return order_genes

def generateBlockRepresentation( order_genes, cna_file, sequenced_samples ):
	
	#build gene_cna_pos and gene_cna_neg tables
	
	file_cn = open( cna_file, 'r' )
	l_count = 0
	cna_samples = list()
	gene_alt_pos = dict()
	gene_alt_neg = dict()
	for line in file_cn:
		if l_count == 0:
			v = line.split('\t')
			for i in range( len(v) ):
				if i >2:
					r = v[i].strip('\",\n').split('-')
					sampleID = r[0]+'-'+r[1]+'-'+r[2]
					cna_samples.append( sampleID )
		else:
			v = line.split('\t')
			geneID = v[0]
			if geneID not in gene_alt_pos:
				gene_alt_pos[geneID] = set()
			if geneID not in gene_alt_neg:
				gene_alt_neg[geneID] = set()
			for i in range( len(v) ):
				if i>2 and v[i].strip() != 'NaN' and v[i].strip() != 'NA':
					if float(v[i].rstrip('\n')) > 1 :
						if cna_samples[i-3] in sequenced_samples:
							gene_alt_pos[ geneID ].add( cna_samples[i-3] )
					if float(v[i].rstrip('\n')) < -1 :
						if cna_samples[i-3] in sequenced_samples:
							gene_alt_neg[ geneID ].add( cna_samples[i-3] )
		l_count += 1
	
	file_cn.close()
	
	cna_blocks = dict()	
	#now build the table with the block representation of cna input data
	
	#loop on units (chrm or arms)
	for i in order_genes:
		cna_blocks[i] = dict()
		curr_order = order_genes[i]
		#loop on samples ID
		for sample in sequenced_samples:
			# build sequence of annotations
			curr_lists = list()
			#NB: if gene not tested for cna, it has a 0
			for j in range(len(curr_order)):
				curr_gene = curr_order[j]
				if curr_gene in gene_alt_pos:
					if sample in gene_alt_pos[ curr_gene ]:
						curr_annotation = 2
					else:
						if sample in gene_alt_neg[ curr_gene ]:
							curr_annotation = -2
						else:
							curr_annotation = 0
				else:
					curr_annotation = 0
				curr_lists.append( curr_annotation )
				#if curr_annotation != 0:
					#print curr_annotation
			#buil sequence of blocks
			curr_blocks = list()
			if len(curr_lists) > 0 :
				curr_annotation = curr_lists[0]
				curr_length = 1
				for c in range(len(curr_lists)-1):
					if curr_lists[c] == curr_lists[c+1]:
						curr_length += 1
					else:
						to_append = list([curr_annotation , curr_length])
						curr_blocks.append( to_append )
						curr_annotation = curr_lists[c+1]
						curr_length = 1
				to_append = list([curr_annotation , curr_length])
				curr_blocks.append( to_append )
			#now append the sequence of blocks
			(cna_blocks[i])[sample] = curr_blocks
	return cna_blocks

def generateRandomCNA( cna_blocks, gene_to_index, gene_order, sequenced_samples, tested_genes ):
	random_cna_blocks = dict()
	for sample in sequenced_samples:
		for j in cna_blocks:
			if j not in random_cna_blocks:
				random_cna_blocks[j] = dict()
			curr_block = (cna_blocks[j])[sample]
			zeros_list = list()
			alter_list = list()
			for i in range( len( curr_block ) ):
				if (curr_block[i])[0] == 0:
					zeros_list.append( curr_block[i] )
				else:
					alter_list.append( curr_block[i] )
			random_block = list()
			#discriminate on difference of lengths of alter_list and zero_list
			if len( zeros_list ) == len( alter_list ):
				#decide if start with zeros or alter:
				coin = random.randint(0 , 1)
				if coin == 0:
					#start with zero:
					to_choose = len( zeros_list)
					num_chosen = 0
					while num_chosen < to_choose:
						#print "len zeros: " + str(len(zeros_list))
						#print "len alter: " + str(len(alter_list))
						#print "chosen: " + str(num_chosen)
						#print "to_choose: " + str(to_choose)
						chosen = random.randint( 0, len(zeros_list) - 1)
						random_block.append( zeros_list.pop( chosen ) )
						chosen = random.randint( 0, len(alter_list) - 1)
						random_block.append( alter_list.pop( chosen ) )
						num_chosen += 1
				else:
					#start with alteration:
					to_choose = len( alter_list)
					num_chosen = 0
					while num_chosen < to_choose:
						chosen = random.randint( 0, len(alter_list) - 1)
						random_block.append( alter_list.pop( chosen ) )
						chosen = random.randint( 0, len(zeros_list) - 1)
						random_block.append( zeros_list.pop( chosen ) )
						num_chosen += 1
			else:
				if len( zeros_list ) > len( alter_list ):
					#start with zero:
					chosen = random.randint( 0, len(zeros_list) - 1)
					random_block.append( zeros_list.pop( chosen ) )
					to_choose = len( alter_list )
					num_chosen = 0
					while num_chosen < to_choose:
						chosen = random.randint( 0, len(alter_list) - 1)
						random_block.append( alter_list.pop( chosen ) )
						chosen = random.randint( 0, len(zeros_list) - 1)
						random_block.append( zeros_list.pop( chosen ) )
						num_chosen += 1
				else:
					#start with alteration:
					chosen = random.randint( 0, len(alter_list) - 1)
					random_block.append( alter_list.pop( chosen ) )
					to_choose = len(zeros_list)
					num_chosen = 0
					while num_chosen < to_choose:
						chosen = random.randint( 0, len(zeros_list) - 1)
						random_block.append( zeros_list.pop( chosen ) )
						chosen = random.randint( 0, len(alter_list) - 1)
						random_block.append( alter_list.pop( chosen ) )
						num_chosen += 1
			(random_cna_blocks[j])[sample] = random_block
	#for sample in sequenced_samples:
		#print sample
		#for i in random_cna_blocks:
			#print "---"+str(i)
			#print (random_cna_blocks[i])[sample]
	#now build the map of genes and samples with cna
	gene_cna_random_pos = dict()
	gene_cna_random_neg = dict()
	gene_cna_altered = list()
	for i in random_cna_blocks:
#		print i
		for sample in sequenced_samples:
			curr_block = (random_cna_blocks[i])[sample]
			first_elem = 0
			for j in range(len(curr_block)):
				curr_annotation = (curr_block[j])[0]
				length = (curr_block[j])[1]
				if curr_annotation == 2:
					for l in range( length ):
						gene = (gene_order[i])[first_elem + l]
						if gene in gene_to_index:
							index = gene_to_index[ gene ]
						else:
							index = -1
						if index in tested_genes:
							if index not in gene_cna_altered:
								gene_cna_altered.append(index)
							if index not in gene_cna_random_pos:
								gene_cna_random_pos[index] =set()
							gene_cna_random_pos[index].add( sample )
				else:
					if curr_annotation == -2:
						for l in range( length ):
							gene = (gene_order[i])[first_elem + l]
							if gene in gene_to_index:
								index = gene_to_index[ gene ]
							else:
								index = -1
							if index in tested_genes:
								if index not in gene_cna_altered:
									gene_cna_altered.append(index)
								if index not in gene_cna_random_neg:
									gene_cna_random_neg[index] =set()
								gene_cna_random_neg[index].add( sample )
				first_elem = first_elem + length
	#sort indices of altered genes
	gene_cna_random = dict()
	thresh = 0.9
	for index in gene_cna_altered:
		num_pos = 0.0
		num_neg = 0.0
		if index in gene_cna_random_pos:
			num_pos = float(len(gene_cna_random_pos[index]))
		if index in gene_cna_random_neg:
			num_neg = float(len(gene_cna_random_neg[index]))
		if num_pos != 0.0 and num_neg != 0.0:
			tot = num_pos + num_neg
			if num_pos > num_neg:
				ratio = num_pos / tot
				if ratio >= thresh:
					gene_cna_random[index]=set()
					gene_cna_random[index].update(gene_cna_random_pos[index])
			else:
				ratio = num_neg / tot
				if ratio >= thresh:
					gene_cna_random[index]=set()
					gene_cna_random[index].update(gene_cna_random_neg[index])
		else:
			if num_pos > 0:
				gene_cna_random[index]=set()
				gene_cna_random[index].update(gene_cna_random_pos[index])
			if num_neg > 0:
				gene_cna_random[index]=set()
				gene_cna_random[index].update(gene_cna_random_neg[index])
	return gene_cna_random

def generateRandomDatasets( cna_file, tested_genes, gene_to_index, index_to_gene, genes_length, BMR, num_rep, sequenced_samples ):
	
	random.seed()
	
	#build the order of genes in chr/arms
	order_genes = buildGeneOrder()
	
	#build the block representation of cna input data
	cna_blocks = generateBlockRepresentation( order_genes, cna_file, sequenced_samples )
	
	num_samples = len( sequenced_samples )
	#iterates on the number or random datasets
	for i in range(num_rep):
		
		#total random alterations; genes represented as indices
		gene_alter_random = dict()
		
		#generate random dataset of CNA; genes represented as indices
		gene_cna_random = generateRandomCNA( cna_blocks, gene_to_index, order_genes, sequenced_samples, tested_genes )
		
		for gene in gene_cna_random:
			if gene not in gene_alter_random:
				gene_alter_random[gene] = set()
			gene_alter_random[gene].update( gene_cna_random[gene] )
		
		#count num cna in dataset:
		tot_cna = 0
		for gene in gene_cna_random:
			tot_cna += len( gene_cna_random[gene] )
#		print "Total CNA in random: " + str(tot_cna)
		
		#iterates on the tested samples
		for sample in sequenced_samples:

			#generate random dataset of mutations; gene is represented as index
			mut_genes = generateRandomMutations( tested_genes, index_to_gene, genes_length, BMR )
			
			#print sample
			#print "Num. mut genes: "+str( len(mut_genes) )
			
			for gene in mut_genes:
				if gene not in gene_alter_random:
					gene_alter_random[gene] = set()
				gene_alter_random[gene].add(sample)
				
		alter_f = open( "tmp/edgesLists/permutation_test/"+str(i+1)+".sparse" , 'w' )
		genes = gene_alter_random.keys()
		genes.sort()
		for i in range(len(genes)):
			gene=genes[i]
			alter_f.write( str(gene)+" 1 "+str(len(gene_alter_random[gene]))+"\n" )
#			print gene_alter_random[gene]
		alter_f.close()
		
		
