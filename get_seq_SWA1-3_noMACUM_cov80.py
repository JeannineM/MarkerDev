#!/usr/local/bin/python2.6
# -*- coding: utf-8 -*-
#USAGE python2.6 scriptname.py

import os as os

seq_name=open('SWA123_MACUM_cov80.tab','r')

list_of_nameProt=[]
list_of_nameSpec1=[]
list_of_nameSpec2=[]
list_of_nameSpec3=[]
#list_of_nameSpec4=[]

for line in seq_name:
  i=line.rstrip('\n')
  i=i.rsplit('\t')
  list_of_nameProt.append(i[0])
  list_of_nameSpec1.append(i[1])
  list_of_nameSpec2.append(i[2])
  list_of_nameSpec3.append(i[3])
 # list_of_nameSpec4.append(i[4])
  

  
def write_Fasta(input_file,list_of_name):
  fasta_file=open(input_file,'r')
  
  u=fasta_file.readline()

  
  while u:

    if u[0]=='>':
      if u.rstrip('\n')[1:] in list_of_name:
	f=open("gene_seq/"+list_of_nameProt[list_of_name.index(u.rstrip('\n')[1:])].replace('|','_')+".fas",'a')
	f.write(u)

	u=fasta_file.readline()
	while u[0]!='>':

	  f.write(u.rstrip('\n'))
	  u=fasta_file.readline()
	f.write('\n')
	f.close()
      else:
	u=fasta_file.readline()
    else:
      u=fasta_file.readline()
	
  fasta_file.close()
  
#write_Fasta('MACU_cDNA_codgenes.fas',list_of_nameProt)
write_Fasta('SWA1_ren.okay_ren.cds',list_of_nameSpec1)
write_Fasta('SWA2_ren.okay_ren.cds',list_of_nameSpec2)
write_Fasta('SWA3_ren.okay_ren.cds',list_of_nameSpec3)
#write_Fasta('SWA4_ren2.okay_ren.cds',list_of_nameSpec4)