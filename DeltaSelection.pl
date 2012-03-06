#!/usr/bin/perl
# Usage: perl delta.pl [max file]

# Copyright 2010,2011 Brown University, Providence, RI.
# 
#                          All Rights Reserved
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose other than its incorporation into a
# commercial product is hereby granted without fee, provided that the
# above copyright notice appear in all copies and that both that
# copyright notice and this permission notice appear in supporting
# documentation, and that the name of Brown University not be used in
# advertising or publicity pertaining to distribution of the software
# without specific, written prior permission.
# 
# BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
# PARTICULAR PURPOSE.  IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR
# ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


my @valset;
my @delta_index;
open(IN, $ARGV[0]);
while(<IN>){
	chomp;
	my @tmp = split(/\t|\s+/, $_);
	push(@valset, $tmp[3]);
	push(@delta_index, $tmp[0]);
}
close(IN);

# Max point finding
my $max=0;
my $max_index;
for(my $i=0; $i<scalar(@valset);$i++){ 
	my $t = $valset[$i];
	if($t > $max){
		$max = $t;
		$max_index = $i;
	}
}


#print $max."	".$max_index*0.0002."\n";
#print $max_diff."	".$max_diff_index*0.0002."\n";
print "Delta starting at: ".$delta_index[$max_index]."\n";
