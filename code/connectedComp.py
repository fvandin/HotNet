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
sys.path.append('..')
from pygraph.classes.graph import graph
from pygraph.algorithms.searching import breadth_first_search
from pygraph.algorithms.accessibility import connected_components

f_in = open("tmp/index_tested.txt" , 'r' )

nodes_list = list()

for line in f_in:
	nodes_list.append( int( line ) )
f_in.close()

num_lists = 1

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

f_genes = open( gene_index_file , 'r' )

index_to_gene = dict()
for line in f_genes:
	v = line.split()
	index_to_gene[ int(v[0]) ] = (v[1].strip("\n"))

indices_to_consider = [ 1 ]

f_mut = open("tmp/num_alt_per_gene.sparse" , 'r')

num_mut = dict() 
for line in f_mut:
	v = line.split()
	if int(v[0]) in index_to_gene:
		if index_to_gene[int(v[0])] in num_mut:
			print "ERROR: " + index_to_gene[int(v[0])] + " already in mutations"
		else:
			num_mut[ index_to_gene[int(v[0])] ] = v[2]

f_mut.close()

for z in range( num_lists ):

	# Graph creation
	gr = graph()

	gr.add_nodes( nodes_list )
	f_edges = open( 'tmp/edgesLists/edgesList.txt'+str( indices_to_consider[z] ) , 'r' )
	for line in f_edges:
		v = line.split()
		gr.add_edge( int(v[0]) , int(v[1]) )
#		gr.add_edge( (int(v[0]) , int(v[1])) )

	f_edges.close()

	connect_comp_dict = connected_components( gr )
	connect_comp = connect_comp_dict.keys()
	cc = dict()
	for i in range( len( connect_comp ) ):
		id_cc = connect_comp_dict[ connect_comp[ i ] ]
		if id_cc not in cc:
			cc[ id_cc ] = 0
		cc[ id_cc ] += 1
	cc_ids = cc.keys()
	cc_sizes = dict()
	for i in range( len( cc_ids ) ):
		curr_cc_size = cc[ cc_ids[i] ]
		if curr_cc_size not in cc_sizes:
			cc_sizes[ curr_cc_size ] = 0
		cc_sizes[ curr_cc_size ] += 1
	if z == 0:
		f_cc = open('OUTPUT/genes_cc.txt','w')
		cc_genes = dict()
		for i in range( len( connect_comp ) ):
			curr_cc_id = connect_comp_dict[ connect_comp[i] ]
			if cc[ curr_cc_id ] > 1:
				if curr_cc_id not in cc_genes:
					cc_genes[ curr_cc_id ] = set()
				cc_genes[ curr_cc_id ].add( index_to_gene [ connect_comp[i] ] )
				if index_to_gene [ connect_comp[i] ] == "-":
					print "size: " + str( cc[ curr_cc_id ] ) + ", index: " + str( connect_comp[i] )
		cc_g = cc_genes.keys()
		for i in range( len( cc_g ) ):
			curr_cc = cc_genes[ cc_g[i] ]
			l_curr_cc = list( curr_cc )
			for l in range( len( l_curr_cc ) -1):
				f_cc.write( str( l_curr_cc[l] )+"("+num_mut[l_curr_cc[l]]+")\t" )
			f_cc.write( str(l_curr_cc[ len( l_curr_cc) - 1] )+"("+num_mut[l_curr_cc[len(l_curr_cc)-1]]+")\n" )
		f_cc.close()
	sizes = cc_sizes.keys()
	sizes.sort()
	f_sizes = open( 'tmp/sizes/t' + str( indices_to_consider[z] ) , 'w' )

	to_test_f = open( 'code/to_test.txt','r' )
	to_test_list = list()
	for line in to_test_f:
		to_test_list.append( int(line.strip()) )
	to_test_f.close()

	for j in range( len( to_test_list ) ):
		to_test_s = to_test_list[j]
		to_p = 0
		for i in range(len( sizes )):
			if to_test_s <= sizes[i]:
				to_p += cc_sizes[sizes[i]]
		f_sizes.write( str(to_test_s) + "\t" + str( to_p ) +"\n" )
	f_sizes.close()
