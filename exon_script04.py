#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# USAGE ./sript.py 
# v 27/01/2015
# authors JM and AB
# This script extracts SNPs which fulfill certain criteria within certain boundaries given a sequence

# Input 
# blastresult exons vs ref outfmt 6 + manedt strandedness
# vcf.file containing all SNPs within loci
# vcf.file containing targetd SNPs for resequencing

###### input ######

#1 blast MUSAexons vs BB_refseq1042
file_blast2=open('MAevsrefseq1046.rmErr.sel1042.oft6.tblx.w_pct_hit_length','r')

#2 vcf file containing targeted SNPs 
#file_tarvcf=open('BBSWA12vs1046MusA.finalSnps_ef_AF.codanno_fixedSNPs_244.vcf','r')
file_tarvcf=open('BBSWA12vs1046MusA.finalSnps_ef_AF.codanno_fixedSNPs_sub730loci.wh.vcf','r')

#4 vcf file containing all SNPs within targeted loci

file_allvcf=open('BBSWA12vs1046MusA.finalSnps_ef_AF.codanno_sub730loci.wh.vcf','r')
#file_allvcf=open('BBSWA12vs1046MusA.finalSnps_ef_AF.codanno_sub244loci.vcf','r')

marker_neutral=open('neutralmarker139.greplist','r')
marker_GO=open('fix.GO.allNS153.greplist','r')

###### output ######

# output table
#out_tab=open('targetedSNPs_244loci.list','w') 
out_tab=open('targetedSNPs_730loci_80flank_09-02.list','w') 

#out_tab=open('targetedSNPs_test.list','w') 

out_tab.write('# LOCUS NAME \t SEQUENCE_ID \t SEQUENCE_INCLUDED_REGION \t SEQUENCE_EXCLUDED_REGION \t SEQUENCE_TARGET \t STRAND \t BORDER \t COD \n')

out_acc_snps=open('acceptedSNPs_730loci_80flank_09-02.vcf','w')
#out_acc_snps=open('acceptedSNPs_test.vcf','w')

###### script ######

#data_file=file_blast2.readlines()
#data_file=data_file[0].rsplit('\r')

## length of the flanking region around the targeted SNP

flank=80

## reading the file_blast2 and store in memory

dict_max_exon={}
data_exon=file_blast2.readlines()
for i in range(1,len(data_exon)):
  data_exon[i]=data_exon[i].rsplit('\t')
  if data_exon[i][1] not in dict_max_exon.keys():
    dict_max_exon[data_exon[i][1]]=1
#  exon_end
  dict_max_exon[data_exon[i][1]]=max(dict_max_exon[data_exon[i][1]],int(data_exon[i][9]))

dict_snp={}
for line in file_allvcf:
  if line[0]!='#':
    info=line.rsplit('\t')
    SNPpos=int(info[1])
    name_gene=info[0].rsplit('|')[0]
    if name_gene not in dict_snp.keys():
      dict_snp[name_gene]=[]
    dict_snp[name_gene].append(SNPpos)   

#dict_N={}
#for line in marker_neutral:
#  marker_name=line.rsplit('\n')[0]
#  char1='neutral'
#  if marker_name not in dict_N.keys():
#    dict_N[marker_name]=[]
#  dict_N[marker_name].append(char1)

#dict_GO={}
#for line in marker_GO:
#  marker_name=line.rsplit('\n')[0]
#  char2='flower'
#  if marker_name not in dict_GO.keys():
#    dict_GO[marker_name]=[]
#  dict_GO[marker_name].append(char2)

## 

for line in file_tarvcf:
  if line[0]!='#':
    info=line.rsplit('\t')
    vcfline=line.split('\t')
    SNPpos=int(info[1])
    COD=info[7].rsplit(';')[-1].rsplit('=')[1]    
    name_gene=info[0].rsplit('|')[0]
    bool_t=0
    #vwk=input()
    for j in range(1,len(data_exon)):
      if data_exon[j][1].count(name_gene):
# check strand positive
	if data_exon[j][10]=='+': 
	  #check SNPs have 100bp flanking region
	  if SNPpos<int(data_exon[j][9])-flank and SNPpos>int(data_exon[j][8])+flank:
	    #check if SNP in first or last exon and if true add 10bp to flanking region
	    if ((int(data_exon[j][8])==1 and SNPpos<int(data_exon[j][8])+flank+10)) or ((int(data_exon[j][9])==dict_max_exon[data_exon[j][1]] and SNPpos>int(data_exon[j][9])-10-flank)):
	      if name_gene in dict_snp.keys():
		list_excluded=''
		for i in dict_snp[name_gene]:
		  if i!= SNPpos and i <int(data_exon[j][9]) and i > int(data_exon[j][8]):
		    list_excluded=list_excluded+str(i)+', 1; '
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][8]))+', '+str(int(data_exon[j][13])-1)+'\t'+list_excluded+'\t'+str(SNPpos)+',1'+'\t'+data_exon[j][10]+'\t'+'border\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))
	      else:
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][8]))+', '+str(int(data_exon[j][13])-1)+'\t'+'\t'+str(SNPpos)+',1'+'\t'+data_exon[j][10]+'\t'+'border\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))
	    else:
	      if name_gene in dict_snp.keys():
		list_excluded=''
		for i in dict_snp[name_gene]:
		  if i!= SNPpos and i <int(data_exon[j][9]) and i > int(data_exon[j][8]):
		    list_excluded=list_excluded+str(i)+', 1; '
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][8]))+', '+str(int(data_exon[j][13])-1)+'\t'+list_excluded+'\t'+str(SNPpos)+',1\t'+data_exon[j][10]+'\t\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))
	      else:
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][8]))+', '+str(int(data_exon[j][13])-1)+'\t'+'\t'+str(SNPpos)+',1\t'+data_exon[j][10]+'\t\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))
# negative stranded		
	else: 
	  if SNPpos<int(data_exon[j][8])-flank and SNPpos>int(data_exon[j][9])+flank:	    
	    if ((int(data_exon[j][9])==1 and SNPpos<int(data_exon[j][9])+flank+10)) or ((int(data_exon[j][8])==dict_max_exon[data_exon[j][1]] and SNPpos>int(data_exon[j][8])-flank-10)):
	      if name_gene in dict_snp.keys():
		list_excluded=''
		for i in dict_snp[name_gene]:
		  if i!= SNPpos and i <int(data_exon[j][8]) and i > int(data_exon[j][9]):
		    list_excluded=list_excluded+str(i)+', 1; '
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][9]))+', '+str(int(data_exon[j][13])-1)+'\t'+list_excluded+'\t'+str(SNPpos)+',1'+'\t'+data_exon[j][10]+'\t'+'border\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))		  
	      else:
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][9]))+', '+str(int(data_exon[j][13])-1)+'\t'+'\t'+str(SNPpos)+',1'+'\t'+data_exon[j][10]+'\t'+'border\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))		  
	    else:
	      if name_gene in dict_snp.keys():
		list_excluded=''
		for i in dict_snp[name_gene]:
		  if i!= SNPpos and i <int(data_exon[j][8]) and i > int(data_exon[j][9]):
		    list_excluded=list_excluded+str(i)+', 1; '
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][9]))+', '+str(int(data_exon[j][13])-1)+'\t'+list_excluded+'\t'+str(SNPpos)+',1\t'+data_exon[j][10]+'\t\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))		  
	      else:
		out_tab.write(info[0]+'\t'+info[0]+'|exon.'+data_exon[j][0].rsplit('.')[1].rsplit('|')[0]+'\t'+str(int(data_exon[j][9]))+', '+str(int(data_exon[j][13])-1)+'\t'+'\t'+str(SNPpos)+',1\t'+data_exon[j][10]+'\t\t'+COD+'\n')
		#
		out_acc_snps.write(str(vcfline[0])+'\t'+str(vcfline[1])+'\t'+str(vcfline[2])+'\t'+str(vcfline[3])+'\t'+str(vcfline[4])+'\t'+str(vcfline[5])+'\t'+str(vcfline[6])+'\t'+str(vcfline[7])+'\t'+str(vcfline[8])+'\t'+str(vcfline[9])+'\t'+str(vcfline[10]))		  
	    
out_tab.close()
out_acc_snps.close()