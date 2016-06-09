#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# USAGE sript genelist.txt annotab.csv
# v 05/08/2014
# This script pulls out the entries from the banana gene annotation table based on a list of according banana gene names. 

import os as os
import sys

#genelist = open('02_get_common_gene_seq/sorted_genelist.txt','r')
#MACUanno = open('00D_DB/banana_genefunct_biomart.txt','r')

 


genelist = open(sys.argv[1],'r')
MACUanno = open(sys.argv[2],'r')
geneanno = open("parselength_mitoch_28oct.csv",'a')

#element=MACUanno.readlines().rsplit()
#gene=genelist.readline().rstrip("\n")
#element=element.split(',')
rows=[]
list_name=[]

for line in MACUanno:
    rows.append(line)
    list_name.append(line.rsplit(',')[0])


for i in genelist:
    j=i.rstrip('\n')
    marker=0
    while (list_name[marker:].count(j)>0): 
        index_g=list_name[marker:].index(j)
        geneanno.write(rows[marker+index_g])
        marker+=index_g+1
        