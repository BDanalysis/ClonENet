import pandas as pd
import numpy as np

def readfile_allect_states(path):
    dataframe=pd.read_csv(path,sep='\t',encoding='utf-8')
    dataframe = dataframe[['Chromosome','Start','End','CopyNr','A','B']]
    return dataframe

def readfile_snp(path):
    dataframe=pd.read_csv(path,sep='\t',header=None,comment='#',encoding='utf-8')
    dataframe.columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']
    dataframe = dataframe[['CHROM','POS','INFO']]
    dataframe['DP'] = dataframe['INFO'].str.split(';').str[0]
    dataframe['AF'] = dataframe['INFO'].str.split(';').str[1]
    dataframe['DP'] = dataframe['DP'].str.split('=').str[1]
    dataframe['AF'] = dataframe['AF'].str.split('=').str[1]
    dataframe['DP'] = dataframe['DP'].astype(int)
    dataframe['AF'] = dataframe['AF'].astype(float)
    dataframe['AD'] = (dataframe['DP']*dataframe['AF']).astype(int)
    dataframe = dataframe[['CHROM','POS','AD','DP']]
    return dataframe

def write_cnv(dataframe,filename):
    dataframe['Chromosome'] = dataframe['Chromosome'].str.split('chr').str[1]
    write = dataframe[['Chromosome','Start','End','A','B']]
    write.to_csv(filename,mode='w',index=0,header=0,sep='\t')

def write_snv(dataframe,filename):
    dataframe['CHROM'] = dataframe['CHROM'].str.split('chr').str[1]
    write = dataframe
    write.to_csv(filename,mode='w',index=0,header=0,sep='\t')

def write_purity(purity,filename):
    output = open(filename,'w')
    output.write(str(purity))
    output.close()
title = "EGAF00000098597"
path_vcf = 'example/'+title+'.sclust.vcf'
path_cnv = 'example/'+title+'_allelic_states.txt'

cnv = readfile_allect_states(path_cnv)
snp = readfile_snp(path_vcf)

write_purity(0.26,'./Sample_data/'+title+'_purity.txt')
write_cnv(cnv,'./Sample_data/'+title+'_cnv.txt')
write_snv(snp, 'Sample_data/'+title+'_snp.txt')