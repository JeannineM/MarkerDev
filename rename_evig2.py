##Usage python2.7 path/to/folder/rename_evig2.py <infile> <outfile> <prefix>
#by Hannes Becher June 2014
#script to rename fasta file headers into > prefix_int and output second file
#containing former header information

import sys
import time
import re
inhandle=open(sys.argv[1],'r')
t0=time.time()
count=0
ht=0
c=0
out=open(sys.argv[2],'w')
out2=open('%s.names' % sys.argv[2], 'w')

for line in inhandle:
	if line[0]=='>':
		c += 1
		out.write('>%s_%d\n' % (sys.argv[3],c))
		out2.write('%s_%d\t%s' % (sys.argv[3], c, line))
	else:
		out.write(line)