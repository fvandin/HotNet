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
num_rep = int(sys.argv[7])

print "num_rep: " + str(num_rep)

size_file = open('tmp/sizes/t1' , 'r' )

size_index = list()
num_cc_index = list()

for line in size_file:
	v=line.split()
	if int(v[0]) > 1:
		size_index.append( int(v[0]) )
		num_cc_index.append( int(v[1]) )

size_file.close()

if len( size_index ) == 0:
	print "THERE ARE NO Connected Components. TERMINATING"
	exit(0)

exp_cc_index = list()
prob_cc_index = list()
for i in range( len( size_index ) ):
        exp_cc_index.append( 0 )
        prob_cc_index.append( 0 )

os.system( "matlab -nodesktop -nosplash -nodisplay -nojvm < code/matlab_for_generateEdges_random.m")

for rep in range( num_rep ):
	print rep

	curr_num_cc_index = list()

	for j in range( len( size_index ) ):
                curr_num_cc_index.append( 0 )

	for b in range( num_lists ):
		
		# Graph creation
		gr = graph()

		gr.add_nodes( nodes_list )
		f_edges = open( 'tmp/edgesLists/permutation_test/edgesList.txt'+str(rep+1) , 'r' )

		for line in f_edges:
			v = line.split()
			gr.add_edge( int(v[0]) , int(v[1]) )
#			gr.add_edge( (int(v[0]) , int(v[1])) )

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
		for i in range( len( cc_ids ) ):
			curr_cc_size = cc[ cc_ids[i] ]

			for j in range( len( size_index ) ):
                                if curr_cc_size >= size_index[j]:
                                        exp_cc_index[j] += 1
                                        curr_num_cc_index[j] += 1

                for j in range( len( size_index ) ):

			if curr_num_cc_index[j] >= num_cc_index[j]:
                                prob_cc_index[j] += 1

file_res = open( 'OUTPUT/randomTest.txt' ,'w')
for j in range( len( size_index ) ):
	prob = float( prob_cc_index[j] )/float( num_rep )
	file_res.write( "probability of size " + str( size_index[j] ) +": " + str(prob) + "\n")
	exp = float( exp_cc_index[j] )/float( num_rep )
        file_res.write( "expected: " + str(exp) + "\n" )
	file_res.write( "found: " + str(num_cc_index[j])+"\n\n" )
file_res.close()
