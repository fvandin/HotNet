%  Copyright 2010,2011 Brown University, Providence, RI.
%  
%                           All Rights Reserved
%  
%  Permission to use, copy, modify, and distribute this software and its
%  documentation for any purpose other than its incorporation into a
%  commercial product is hereby granted without fee, provided that the
%  above copyright notice appear in all copies and that both that
%  copyright notice and this permission notice appear in supporting
%  documentation, and that the name of Brown University not be used in
%  advertising or publicity pertaining to distribution of the software
%  without specific, written prior permission.
%  
%  BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
%  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%  PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
%  ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
%  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
%  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
%  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

function generateEdges_Distr( G, analyzed, last_real_an, num_rep)
close all

for z=1:1:num_rep

	name_file = strcat( '../tmp/edgesLists/permutation_test/' , int2str( z ) );
	name_file = strcat( name_file, '.sparse' );
	num_mut_per_gene = load(name_file);

	mut = spconvert( num_mut_per_gene );
	if size( mut ) < size( G )
		mut( length( G ) ) = 0;
	end
	if size( mut ) < max( analyzed )
		mut( max( analyzed ) ) = 0;
	end
	mut = mut( analyzed );
	
	NG = G( analyzed(1:last_real_an), analyzed(1:last_real_an) );
	maxt = max( max( NG ) );
	for i=1:length(NG)
		NG(i,i) = NG(i,i) - (maxt+1);
	end
	maxt = max( max( NG  ) );
	
	
	gene_lengths=load('../tmp/index_length.sparse');
	
	genes_l = spconvert( gene_lengths);
	if length(genes_l) < size(G)
		genes_l(length(G))=0;
	end
	if length( genes_l ) < max(analyzed)
		genes_l(max(analyzed))=0;
	end
	genes_l = genes_l(analyzed);
	total_l = sum( genes_l );
	length_analyzed=sum(genes_l);
	num_non_real_analyzed = length(find(genes_l == 0));
	%avg_length_non_real_analyzed = (total_length -length_analyzed)/num_non_real_analyzed;
	zer= find(genes_l == 0);
	for i=1:1:length(zer)
		genes_l(zer(i))=avg_length_non_real_analyzed;
	end
	mr=zeros(1,length(genes_l));
	for i=1:1:length(mr)
		mr(i)=mut(i);
	end

	threshs = load( '../tmp/threshs.dat');
	num_threshs = length( threshs );

	f_thresh = fopen( '../tmp/edgesLists/thresh_file.txt' , 'w' );

	num_edges = 0;
	indices_to_consider = [ 1 ];
	for th=1:1:length( threshs  )	
		name_file = strcat( '../tmp/edgesLists/permutation_test/edgesList.txt' , int2str( z ) );
		name_file = strcat( name_file , '_' );
		name_file = strcat( name_file , int2str(th) );
		f_edges = fopen( name_file ,'w');
		fprintf( f_thresh , '%g\n' , threshs(th) );
		for i=1:1:size(NG)
			for j=1:1:i
				if i ~= j
					if mr( i ) == 0 | mr( j ) == 0
						ANG = -1;
					else
						ANG = min( NG(i,j), NG(j,i) ).*max( mr(i),mr(j) );
					end
					if ANG >= threshs(th) 
						fprintf( f_edges , '%d\t%d\n', analyzed(i), analyzed(j) );
						num_edges = num_edges + 1;
					end
				end
			end
		end
		fclose( f_edges );
	end
	fclose( f_thresh );
	clear NG;
end
