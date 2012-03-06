#!/usr/bin/env python

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


import sys
import os, glob
sys.path.append("./code")
import HN

if len(sys.argv) < 9:
	print "Usage: python HotNet.py tested_file mutations_file CNA_file lengths_file BMR thresholds_file num_rep num_thresholds"
	print "tested_file: file with tested genes, one per line"
	print "mutations_file: file with mutations table. Format: sampleTABmutated_genes. mutated_genes are TAB separated"
	print "CNA_file: file with copy number aberrations. Format: see TCGA OV format"
	print "lengths_file: lenghts of tested genes"
	print "BMR: background mutation rate"
	print "thresholds_file: file containing thresholds to be tested, one per line"
	print "num_rep: number of repetions for permutation test"
	print "num_thresholds: number of thresholds to be tested"
	exit(0)

tested_file = sys.argv[1]
mutations_file = sys.argv[2]
CNA_file = sys.argv[3]
lengths_file = sys.argv[4]
BMR = float( sys.argv[5] )
thresh_file = sys.argv[6]
num_rep = int( sys.argv[7] )
num_thresholds = int( sys.argv[8] )

#load configuration info

config_f = open( "code/HN.config" , 'r')
line_c = 0
for line in config_f:
	line_c += 1
	if line_c == 1:
		last_id = int(line)
	if line_c == 2:
		next_id = int(line)
	if line_c == 3:
		gene_index_file = line.strip("\n")
	if line_c == 4:
		infl_matrix_f = line.strip("\n")
config_f.close()

gene_to_index = dict()
index_to_gene = dict()

gene_index_f = open( gene_index_file , 'r' )

for line in gene_index_f:
	v = line.split('\t')
	gene_to_index[v[1].strip('\n')] = int(v[0])
	index_to_gene[int(v[0])] = v[1].strip('\n')

gene_index_f.close()

tested_genes = list()
tested_f = open( tested_file , 'r' )
for line in tested_f:
	gene = line.strip() 
	if gene in gene_to_index:
		to_append = gene_to_index[gene]
	else:
		gene_to_index[gene] = next_id
		index_to_gene[next_id] = gene
		to_append = next_id
		next_id += 1
	if gene not in tested_genes:
		tested_genes.append( to_append )
tested_f.close()

tested_genes.sort()

#this is the index of the last analyzed gene in the network, where the indices starts from 1
last_tested_gene = 0

index_tested_f = open('tmp/index_tested.txt' , 'w')

for i in range( len( tested_genes ) ):
	index_tested_f.write(str(tested_genes[i])+"\n")
	if tested_genes[i] <= last_id:
		last_tested_gene = i+1

index_tested_f.close()

#generate file for matlab for real data; return the number of analyzed samples (= number of sequenced samples)
sequenced_samples = HN.generate_num_alter_per_gene( mutations_file, CNA_file, "tmp/num_alt_per_gene.sparse", tested_genes, gene_to_index )

#print "Number of analyzed samples: " + str( len(sequenced_samples) )

#now generate length of genes and obtain mapping gene_length for tested genes
gene_to_length = HN.generateLengths( lengths_file, tested_genes, index_to_gene, "tmp/index_length.sparse" )

#print parameters to generate cc on file input_val.dat

matlab_val_file = open( 'tmp/input_val.dat' , 'w' )
matlab_val_file.write( str(last_tested_gene)+"\n" )
matlab_val_file.write( str(num_rep) +"\n" )
matlab_val_file.write( str(BMR)+"\n" )
matlab_val_file.close()

matlab_files = open( 'tmp/input_files.dat' , 'w' )
matlab_files.write(infl_matrix_f)
matlab_files.close()

matlab_thresh_file = open( 'tmp/threshs.dat' , 'w' )
thresh_f = open( thresh_file, 'r' )
for line in thresh_f:
	matlab_thresh_file.write(line)
matlab_thresh_file.close()
thresh_f.close()

#generate all the random datasets in the permutation test for matlab
HN.generateRandomDatasets( CNA_file, tested_genes, gene_to_index, index_to_gene, gene_to_length, BMR, num_rep, sequenced_samples )

#now call the code for computing p-values
import connectedCompRandom_Distr

for infile in glob.glob( os.path.join('./tmp/', '*') ):
	if not os.path.isdir( infile) :
		os.remove(infile)
for infile in glob.glob( os.path.join('./tmp/edgesLists/', '*') ):
	if not os.path.isdir( infile) :
		os.remove(infile)
for infile in glob.glob( os.path.join('./tmp/edgesLists/permutation_test/', '*') ):
	if not os.path.isdir( infile) :
		os.remove(infile)
