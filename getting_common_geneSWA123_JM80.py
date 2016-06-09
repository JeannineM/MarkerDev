#!/usr/local/bin/python2.6
# -*- coding: utf-8 -*-

import os as os
import numpy as nu


def Making_list(file1, file2):
  fileBanTo2 = open(file1,'r')

  fileBanTo2.next()
  outputlist1=[]
  seqname1=[]


  for line in fileBanTo2:
    element=line.rsplit('\t')
    if float(element[(len(element)-2)]) > 80:
      outputlist1.append(element[0])
      seqname1.append(element[1])
      
  fileBanTo2.close()
  
  #print len(outputlist)
  
  file2ToBan = open(file2,'r')
  file2ToBan.next()
  outputlist2=[]
  seqname2=[]

  for line in file2ToBan:
    element=line.rsplit('\t')
    if float(element[(len(element)-2)]) > 80:
      outputlist2.append(element[1])
      seqname2.append(element[0])

  
  file2ToBan.close()
  
  outputlist=[]
  seqname=[]
  
  for i in outputlist1:
    if i in outputlist2:
      if seqname1[outputlist1.index(i)]==seqname2[outputlist2.index(i)]:
	outputlist.append(i)
	seqname.append(seqname1[outputlist1.index(i)])
	
  
  return [outputlist,seqname]

List_Ban_Species1=Making_list("MACUMcDNA_to_DBSWA1aa.oft6.blx.w_pct_hit_length","SWA1aa_to_DBMACUMcDNA.oft6.tbln.w_pct_hit_length")
List_Ban_Species2=Making_list("MACUMcDNA_to_DBSWA2aa.oft6.blx.w_pct_hit_length","SWA2aa_to_DBMACUMcDNA.oft6.tbln.w_pct_hit_length")
List_Ban_Species3=Making_list("MACUMcDNA_to_DBSWA3aa.oft6.blx.w_pct_hit_length","SWA3aa_to_DBMACUMcDNA.oft6.tbln.w_pct_hit_length")
#List_Ban_Species4=Making_list("MACUMcDNA_to_DBSWA4aa.oft6.blx.w_pct_hit_length","SWA4aa_to_DBMACUMcDNA.oft6.tbln.w_pct_hit_length")

#List_Species1_Species2=Making_list("SWA1tr_to_DBtrSWA2.oft6.bln","SWA2tr_to_DBtrSWA1.oft6.bln")
#List_Species1_Species3=Making_list("SWA1tr_to_DBtrSWA3.oft6.bln","SWA3tr_to_DBtrSWA1.oft6.bln")
#List_Species2_Species3=Making_list("SWA2tr_to_DBtrSWA3.oft6.bln","SWA3tr_to_DBtrSWA2.oft6.bln")


List_Species1_Species2=Making_list("SWA1aa_to_DBaaSWA2.oft6.blp","SWA2aa_to_DBaaSWA1.oft6.blp")
List_Species1_Species3=Making_list("SWA1aa_to_DBaaSWA3.oft6.blp","SWA3aa_to_DBaaSWA1.oft6.blp")
#List_Species1_Species4=Making_list("SWA1aa_to_DBaaSWA4.oft6.blp","SWA4aa_to_DBaaSWA1.oft6.blp")
List_Species2_Species3=Making_list("SWA2aa_to_DBaaSWA3.oft6.blp","SWA3aa_to_DBaaSWA2.oft6.blp")
#List_Species2_Species4=Making_list("SWA2aa_to_DBaaSWA4.oft6.blp","SWA4aa_to_DBaaSWA2.oft6.blp")
#List_Species3_Species4=Making_list("SWA3aa_to_DBaaSWA4.oft6.blp","SWA4aa_to_DBaaSWA3.oft6.blp")

print len(List_Ban_Species1[0]), len(List_Ban_Species2[0]), len(List_Ban_Species3[0])

#len(List_Ban_Species4[0])

output=open("SWA123_MACUM_cov80.tab",'w')
output2=open("SWA123_MACUM_cov80_GSIDM_list.txt",'w')

# and i in List_Ban_Species4[0] 
#"\t"+List_Ban_Species4[1][List_Ban_Species4[0].index(i)]+

for i in List_Ban_Species1[0]:
  if (i in List_Ban_Species2[0] and i in List_Ban_Species3[0]):
    output.write(i+"\t"+List_Ban_Species1[1][List_Ban_Species1[0].index(i)]+"\t"+List_Ban_Species2[1][List_Ban_Species2[0].index(i)]+"\t"+List_Ban_Species3[1][List_Ban_Species3[0].index(i)]+"\n")
    output2.write(i.split('|')[0]+"\n")

print len(output2[0])

