Copyright 2010,2011 Brown University, Providence, RI.

                         All Rights Reserved

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose other than its incorporation into a
commercial product is hereby granted without fee, provided that the
above copyright notice appear in all copies and that both that
copyright notice and this permission notice appear in supporting
documentation, and that the name of Brown University not be used in
advertising or publicity pertaining to distribution of the software
without specific, written prior permission.

BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
http://cs.brown.edu/people/braphael/software.html

README file for HotNet:

Software for identifying mutated pathway by heat kernel.

Version: 0.9.1
Version Date: Apr 1 2011

contact: hotnet@cs.brown.edu

-------------------------------------------------------------------------------
Requirements to run Generalized HotNet:
1. Python
2. pygraph package for Python
3. MATLAB
4. Perl

MATLAB must be executable with the command "matlab".  

Generalized HotNet is known to be compatible with the following versions:
-Python 2.5, MATLAB R2010a, pygraph 1.6.2, perl v5.10.0

-------------------------------------------------------------------------------

INSTALL ====================================================

$ ./install.sh

USAGE ======================================================

0. Create the influence matrix for your network (e.g. protein protein interaction network):

To compute the influence matrix: let A be the adiacency matrix of the
undirected network (i.e, the entry of row i and column j is = 1 if there is
an interaction between i and j). The graph Laplacian L  is then given by L
= D-A, where D is a diagonal matrix with D(i,i)=degree of i in the network
(and D(i,j)=0 if i different from j). To compute the heat kernel at time t,
you need to find the exponential of the matrix L*t. For HotNet to work, you
need to assign the matrix the name 'Li', and then save on file.
For example in Matlab you can use:

Li=expm(-L*t);
save name_file.mat Li;

You will also need to create a gene name to index in this matrix mapping
file, associating the row/column index to the gene name. This is discussed
in point 1 below.

1. Pre-Setting
	* Edit your own configuration file: code/HN.config
	* Format:
		line (1) total number of genes in the network
		     (2) some number larger than the number of genes in the interaction network (used as starting index for tested genes not in the network
		     (3) location of index to gene mapping file (should be in format index \t gene name).		     
		     (4) Location of the matlab matrix containing the influence graph for the network (created in step 0).
		     (5) the file describes the order of genes in the genome 
			 - Each line represents a unit (e.g., a chromosome or one arm of a chromosome), and for each unit 
			 there is the list of genes in the order they are found in the unit. Ideally, all the genes in 
			 the genome should be represented in the order; all the genes that are tested must be present 
			 in the order.
	* Example: see code/HN.config (to run examples below, you can use that file and move corresponding files in INPUT/)

	
2. Find the threshold DELTA: we now describe the procedure we use to identify a first value for DELTA. We usually also test some other larger values close to that first value. Other procedures can be used as well.
	Step 1: Get the number of different size connected component found in random datasets using different delta
	
	* python findThreshold.py tested_file mutations_file CNA_file lengths_file BMR thresholds_file num_rep num_thresholds
	
		INPUT files:
		* tested_file: file with genes name, one gene per line (example: example/tested_genes.txt)
		* mutations_file: file with mutations table. Format: sampleTABmutated genes, where mutated genes are TAB 
				  separated(example: example/mutations_table.txt)
		* CNA_file: file with copy number aberrations. Format: i) first line has sample (IDs) starting from 4th token; ii) from second line: first token is gene name, from 4th token +2 if gene in focal amplification, -2 if gene in focal deletion 	(example example/CNA.txt)
		* lengths_file: number of bases analyzed for each gene, one line per gene, format: geneTABnumber 
			        (example INPUT/gene_lengths.txt)
		* BMR: background mutation rate (for random datasets)
		* thresholds_file: file with threshold to use, one threshold per line (example: example/deltas.txt)
		* num_rep: number of random datasets to generate for the permutation test
		* num_thresholds: number of thresholds to be tested in thresholds file

		OUTPUT files:
		* The output can be found in OUTPUT/maxCC_VS_t.txt. For each delta in thresholds file, there is one line 
		  reporting the delta, the maximum size of a pathway found in num_rep random datasets, and numbers of 
		  connected components found in random datasets which is larger than connected component size 3, 4, 5, 6, and 7.  

		Example: python findThreshold.py example/tested_genes.txt example/mutations_table.txt example/CNA.txt example/gene_lengths.txt 0.000001 example/deltas.txt 10 20

	Step 2: Determine a lower bound DELTA
	* perl DeltaSelection.pl OUTPUT/maxCC_VS_t.txt
		We select the DELTA which has the largest change among the distribution. The distribution is the list of 
		deltas and numbers of connected component which has size larger than 4.

3. Find mutated subnetworks

	Pre-step: Choose the subnetwork sizes to be tested, and put them in the file 'code/to_test.txt' (in increasing order). There is already one such file in the 'code' directory, for testing subnetworks of sizes from 3 to 10


	* python HotNet.py tested_file mutations_file CNA_file lengths_file BMR DELTA num_rep

	INPUT files:
	* tested_file: file with genes name, one gene per line (example: example/tested_genes.txt)
	* mutations_file: file with mutations table. Format: sampleTABmutated genes, where mutated genes are TAB
			  separated(example: example/mutations_table.txt)
	* CNA_file: file with copy number aberrations. For format, see 2. above  (example: example/CNA.txt)
	* lengths_file: number of bases analyzed for each gene, one line per gene, format: geneTABnumber (example INPUT/lengths.txt)
	* BMR: background mutation rate (for random datasets)
	* DELTA: the theshold use for analysis (see Section 1)
	* num_rep: number of random datasets to generate for the permutation test

	OUTPUT files:
	* OUTPUT/gene_cc.txt
		all (i.e., not only significant) pathways found in the input dataset using threshold DELTA
	* OUTPUT/randomTest.txt
		for each size s in code/to_test, three values are reported:
		(a) probability of size: the p-value for the number of pathways of size at least s in the input dataset 
		    (the p-value is for pathways of size at least s, not exactly s)
		(b) expected: the expected number of pathways of size at least s in a random dataset
		(c) found: the number of pathways of size at least s found in the input dataset

For example: python HotNet.py example/tested_genes.txt example/mutations_table.txt example/CNA.txt example/gene_lengths.txt 0.000001 0.08 10 

randomTest.txt (in OUTPUT/) would look like:

probability of size 3: 0.1
expected: 5.9
found: 9

probability of size 4: 0.1
expected: 3.5
found: 6

probability of size 5: 0.0
expected: 2.5
found: 5

probability of size 6: 0.0
expected: 1.8
found: 4

probability of size 7: 0.0
expected: 1.7
found: 3

probability of size 8: 0.5
expected: 1.4
found: 2

probability of size 9: 0.4
expected: 1.1
found: 2

probability of size 10: 0.0
expected: 0.5
found: 2

This tells you based on the random permutation trials, the number of
expected number of subnetworks of particular size AND GREATER, and the
ACTUAL number of the size OR GREATER that were found, along with the
calcuted probability of the event based on the random data sets generated.
Using these values you can perform the test described in the reference
below to identify significant subnetworks.

Visualization in Cytoscape: see files in visualization.zip

REFERENCES:

If you use HotNet in your research, please cite:

F. Vandin, E. Upfal, and B.J. Raphael. (2010) Algorithms for Detecting Significantly Mutated Pathways in Cancer. Proceedings of the 14th Annual International Conference on Research in Computational Molecular Biology (RECOMB 2010).
F. Vandin, E. Upfal, and B.J. Raphael. (2011) Algorithms for Detecting Significantly Mutated Pathways in Cancer. Journal of Computational Biology. 18(3):507-22.

WEBSITE:
http://cs.brown.edu/people/braphael/software.html
