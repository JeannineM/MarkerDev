#!/usr/bin/env Rscript
#USAGE Rscript scriptname depthtable outtable

args<-commandArgs(TRUE)
print(args)
#a=read.table('SWA2realign_vsens.coordSorted_samdepth.txt')
a=read.table(args[1])


#output=file('list_coverage_per_gene_per_quantile.txt',open='w')
output=file(args[2],open='w')

list_to_save=c('Gene name','Mean All','Mean Q1','Mean Q2','Mean Q3','Mean Q4','Min All','Min Q1','Min Q2','Min Q3','Min Q4','Max All','Max Q1','Max Q2','Max Q3','Max Q4','\n')
writeLines(list_to_save,con=output,sep='\t')
list_name=levels(a$V1)
# for (i in seq(length(list_name)))

for (i in seq(length(list_name)))
{
b=a[a$V1==list_name[i],]
c=floor(length(b$V1)/4)
d=length(b$V1)%%4

if (d==0) 
{b1=b$V3[(1):(c)];b2=b$V3[(c+1):(2*c)];b3=b$V3[(2*c+1):(3*c)];b4=b$V3[(3*c+1):(4*c)]} 
if (d==1) 
{b1=b$V3[(1):(c)];b2=b$V3[(c+1):(2*c)];b3=b$V3[(2*c+1):(3*c)];b4=b$V3[(3*c+1):(4*c+1)]} 
if (d==2) 
{b1=b$V3[(1):(c)];b2=b$V3[(c+1):(2*c)];b3=b$V3[(2*c+1):(3*c+1)];b4=b$V3[(3*c+2):(4*c+2)]} 
if (d==3) 
{b1=b$V3[(1):(c)];b2=b$V3[(c+1):(2*c+1)];b3=b$V3[(2*c+2):(3*c+2)];b4=b$V3[(3*c+3):(4*c+3)]} 

list_to_save=c(list_name[i],round(mean(b$V3),2),round(mean(b1),2),round(mean(b2),2),round(mean(b3),2),round(mean(b4),2),min(b$V3),min(b1),min(b2),min(b3),min(b4),max(b$V3),max(b1),max(b2),max(b3),max(b4),'\n')

writeLines(list_to_save,con=output,sep="\t")

if(i%%50==0) print(i)

}

close(output)