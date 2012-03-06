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

fid = fopen('tmp/input_files.dat');
C = textscan(fid, '%s');
s = mat2str(cell2mat(C{1}));


load(s(2:end-1));

genes_tested = load('tmp/index_tested.txt');

load('tmp/input_val.dat');

delta = input_val(1);
last_tested = input_val(2);
num_rep =input_val(3);

cd code
generateEdges_random( Li, genes_tested, delta, last_tested, num_rep)
exit
