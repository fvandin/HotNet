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

num_rep = int(sys.argv[7])

num_lists = int( sys.argv[8] )

#percentage_alpha = float(sys.argv[9]);
#percentage_alpha = 0.01;

print "num_rep: " + str(num_rep)

os.system( "matlab -nodesktop -nosplash -nodisplay -nojvm < code/matlab_for_generateEdges_Distr.m")

max_sizes_cc = dict()
max_sizes_thre1cc = dict()
max_sizes_thre2cc = dict()
max_sizes_thre3cc = dict()
max_sizes_thre4cc = dict()
max_sizes_thre5cc = dict()
max_sizes_cc_stat = dict();
max_sizes_thre5cc_stat = dict();
max_sizes_thre4cc_stat = dict();
max_sizes_thre3cc_stat = dict();
max_sizes_thre2cc_stat = dict();
max_sizes_thre1cc_stat = dict();

# different thresholds for catching returning cc with size larger than these setting values
thre1 = 3
thre2 = 4
thre3 = 5
thre4 = 6
thre5 = 7

for rep in range( num_rep ):
	print rep

	if rep == 0:
		for j in range( num_lists ):
			max_sizes_cc[ j+1 ] = 0
			max_sizes_thre1cc[ j+1 ] = 0
			max_sizes_thre2cc[ j+1 ] = 0
			max_sizes_thre3cc[ j+1 ] = 0
			max_sizes_thre4cc[ j+1 ] = 0
			max_sizes_thre5cc[ j+1 ] = 0
			max_sizes_cc_stat[ j+1 ] = list();
			max_sizes_thre5cc_stat[ j+1 ] = list()
			max_sizes_thre4cc_stat[ j+1 ] = list()
			max_sizes_thre3cc_stat[ j+1 ] = list()
			max_sizes_thre2cc_stat[ j+1 ] = list()
			max_sizes_thre1cc_stat[ j+1 ] = list()

	for z in range( num_lists ):
		
		# Graph creation
		gr = graph()

		gr.add_nodes( nodes_list )
		f_edges = open( 'tmp/edgesLists/permutation_test/edgesList.txt'+str(rep+1)+'_'+str(z+1) , 'r' )
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
		max_cc_size = 0
		count_1=0
		count_2=0
		count_3=0
		count_4=0
		count_5=0
		for i in range( len( cc_ids ) ):
			curr_cc_size = cc[ cc_ids[i] ]
			# count the number of cc which is larger than our setting threshold
			if curr_cc_size > thre1:
				count_1 += 1
			if curr_cc_size > thre2:
				count_2 += 1
			if curr_cc_size > thre3:
				count_3 += 1
			if curr_cc_size > thre4:
				count_4 += 1
			if curr_cc_size > thre5:
				count_5 += 1
			if curr_cc_size > max_cc_size:
				max_cc_size = curr_cc_size

		if max_cc_size > max_sizes_cc[ z+1 ]: 
			max_sizes_cc[ z+1 ] = max_cc_size
		if count_1 > max_sizes_thre1cc[ z+1 ]:
			max_sizes_thre1cc[ z+1 ] = count_1
		if count_2 > max_sizes_thre2cc[ z+1 ]:
			max_sizes_thre2cc[ z+1 ] = count_2
		if count_3 > max_sizes_thre3cc[ z+1 ]:
			max_sizes_thre3cc[ z+1 ] = count_3
		if count_4 > max_sizes_thre4cc[ z+1 ]:
			max_sizes_thre4cc[ z+1 ] = count_4
		if count_5 > max_sizes_thre5cc[ z+1 ]:
			max_sizes_thre5cc[ z+1 ] = count_5

		max_sizes_cc_stat[ z+1 ].append(max_cc_size);
		max_sizes_thre5cc_stat[ z+1 ].append(count_5);
		max_sizes_thre4cc_stat[ z+1 ].append(count_4);
		max_sizes_thre3cc_stat[ z+1 ].append(count_3);
		max_sizes_thre2cc_stat[ z+1 ].append(count_2);
		max_sizes_thre1cc_stat[ z+1 ].append(count_1);

f_thresh = open( 'tmp/edgesLists/thresh_file.txt' , 'r')
file_res = open( 'OUTPUT/maxCC_VS_t.txt' ,'w')
mean_percentag = 0.5
for i in range( num_lists  ):
	line = f_thresh.readline()
	#file_res.write( line.rstrip() + "\t" + str( max_sizes_cc[i+1] ) + "\n")
	#nal = int(num_rep*percentage_alpha)
	mean = int(num_rep*mean_percentag)
	max_sizes_cc_stat[i+1].sort();

	inteNumsM = [float(x) for x in max_sizes_cc_stat[i+1]];
	meanM = sum(inteNumsM) / len(max_sizes_cc_stat[i+1]);

	inteNums1 = [float(x) for x in max_sizes_thre1cc_stat[i+1]];
	mean1 = sum(inteNums1) / len(max_sizes_thre1cc_stat[i+1]);
	
	inteNums2 = [float(x) for x in max_sizes_thre2cc_stat[i+1]];
	mean2 = sum(inteNums2) / len(max_sizes_thre2cc_stat[i+1]);

	inteNums3 = [float(x) for x in max_sizes_thre3cc_stat[i+1]];
        mean3 = sum(inteNums3) / len(max_sizes_thre3cc_stat[i+1]);

	inteNums4 = [float(x) for x in max_sizes_thre4cc_stat[i+1]];
        mean4 = sum(inteNums4) / len(max_sizes_thre4cc_stat[i+1]);

	inteNums5 = [float(x) for x in max_sizes_thre5cc_stat[i+1]];
        mean5 = sum(inteNums5) / len(max_sizes_thre5cc_stat[i+1]);

	#print max_sizes_cc_stat[i+1]
	#print max_sizes_cc_stat[i+1][num_rep - nal]
	#file_res.write( line.rstrip() + "\t" + str( max_sizes_cc_stat[i+1][num_rep - nal]) + "\n");
	#file_res.write( line.rstrip() + "\t" + str( max_sizes_cc_stat[i+1][num_rep - nal]) + "\t" + str(max_sizes_cc_stat[i+1][num_rep - mean]) + "\t" + str( max_sizes_cc[i+1] )+ "\t" + str(max_sizes_thre1cc[ i+1 ])+"\t"+str(max_sizes_thre2cc[ i+1 ])+"\t"+str(max_sizes_thre3cc[ i+1 ])+"\t"+str(max_sizes_thre4cc[ i+1 ])+"\t"+str(max_sizes_thre5cc[ i+1 ])+"\n")
	file_res.write( line.rstrip() + "\t" + str(meanM) + "\t" + str( max_sizes_cc[i+1] )+ "\t" + str(mean1)+"\t"+str(mean2)+"\t"+str(mean3)+"\t"+str(mean4)+"\t"+str(mean5)+"\n")
file_res.close()

