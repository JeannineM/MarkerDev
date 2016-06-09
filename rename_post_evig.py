##Usage python2.7 path/to/folder/rename_evig2.py <infile> <outfile>
#by Hannes Becher June 2014
#script to overwrite fasta file headers into '_int' and output second file
# .names containing former header information

import sys
import time
import re
inhandle=open(sys.argv[1],'r')
t0=time.time()
count=0
ht=0
out=open(sys.argv[2],'w')
out2=open('%s.names' % sys.argv[2], 'w')

for line in inhandle:
	if line[0]=='>':
		a=line.split()		
		out.write('%s\n' % (a[0]))
		out2.write('%s\t%s' % (a[0], line[1:]))
	else:
		out.write(line)