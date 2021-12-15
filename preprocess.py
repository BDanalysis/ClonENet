import numpy as np
import ElasticNet_GMM_process as Ep
import getopt
import sys
import gc

def theta(w,minor,major,cn):
  return(np.exp(w)*minor)/(cn+np.exp(w)*major)

def readfile_CNV_SNP(path):
    with open(path, 'r') as f1:
        CNV_data = f1.readlines()
    CNV_table = []
    for i in range(len(CNV_data)):
        CNV_data[i] = CNV_data[i].rstrip('\n')
        tmp=CNV_data[i].split("\t")
        CNV_table.append(tmp)
    table = np.array(CNV_table)
    return table

def main():
    CNV_path = ""
    SNV_path = ""
    purity = ""
    OutPath = ""
    opts, args = getopt.getopt(sys.argv[1:], "hp:c:s:o:x:", ["Help", "Purity=","CNV=", "SNP=","OutPath=","Prefix="])
    for opts, arg in opts:
        if opts == "-h" or opts == "--Help" or opts == "--help":
            print("Required parameters:")
            print("-c or --CNV:",end="\t")
            print("Path to record CNV information file")
            print("-p or --Purity:", end="\t")
            print("Tumor purity")
            print("-s or --SNP:", end="\t")
            print("Path to record SNP information file")
            print("-o or --OutPath:", end="\t")
            print("Path to save the output file")
            print("-x or --Prefix:", end="\t")
            print("The prefix of the output file")
            exit()
        elif opts == "-c" or opts == "--CNV" or opts == "--cnv":
            CNV_path = arg
        elif opts == "-s" or opts == "--SNP" or opts == "--snp":
            SNV_path = arg
        elif opts == "-p" or opts == "--Purity" or opts == "--purity":
            purity = float(arg)
        elif opts == "-o" or opts == "--OutPath" or opts == "--outpath":
            OutPath = arg
        elif opts == "-x" or opts == "--Prefix" or opts == "--prefix":
            Prefix = arg
    return CNV_path,SNV_path,purity,OutPath,Prefix

def check_path(CNV_path,SNV_path,purity,OutPath,Prefix):
    flag = 0
    if CNV_path=='':
        print("Error:\tMissing CNV path")
        flag = 1
    if SNV_path=='':
        print("Error:\tMissing SNP path")
        flag = 1
    if purity=='':
        print("Error:\tMissing purity")
        flag = 1
    if OutPath=='':
        print("Warning:\tMissing output path, use default path './output'")
        OutPath = './output'
    if Prefix=='':
        print("Warning:\tMissing prefix of the output-file, use default prefix 'sample0'")
        Prefix = 'sample0'
    if flag != 0:
        exit(flag)
    return OutPath,Prefix

def check_SNP(SNV_table):
    # check the chrom
    # check the read depth
    # check the minor and total > 0
    arm = np.array(
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22'])
    answer = np.array([-1,-1,-1,-1]).reshape([1,4])
    for char in arm:
        tmp =SNV_table[:,0]
        index_arm = np.argwhere(tmp == char).reshape(-1)
        answer= np.append(answer, SNV_table[index_arm,:], axis=0)
    answer = np.delete(answer, 0, axis=0)
    del tmp,index_arm
    gc.collect()
    tmp_value = answer[:,3].astype(int)#DP
    index_out = np.argwhere(tmp_value >0).reshape(-1)
    answer = answer[index_out,:]
    tmp_value = answer[:, 2].astype(int)#AD
    index_out = np.argwhere(tmp_value <= 0).reshape(-1)
    answer[index_out, 2]=1
    del index_out,tmp_value
    gc.collect()
    if answer.shape[0]<=0:
        print("Error:\tThe number of SNP detection results is insufficient")
        exit(1)
    return answer

def check_SNP_has_AB(SNV_table,CNV_table,col):
    arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X'])
    answer = np.array(['chrom','pos','AD','DP','major','minor']).reshape([1,6])
    for chrom in arm:
        CNV = np.argwhere((CNV_table[:,0]==chrom)).reshape(-1)
        SNV = np.argwhere(SNV_table[:,0]==chrom).reshape(-1)
        if np.size(CNV) <=0 or np.size(SNV) <=0:
            print("Information:\tThere is not any CNV or SNP in chrom "+str(chrom)+".")
            continue
        CNV = CNV_table[CNV,:]
        SNV = SNV_table[SNV,:]
        begins = CNV[:, 1]
        ends = CNV[:, 2]
        if col == 5:
            As = CNV[:, 3]
            Bs = CNV[:, 4]
        else:
            Ts = CNV[:, 3].astype(int)
        position = SNV[:,1].astype(int)
        AD = SNV[:, 2].astype(int)
        DP = SNV[:, 3].astype(int)
        a_chrom = np.array([-1,-1,-1,-1,-1,-1]).reshape([1,6])
        for i in range(np.size(begins)):
            begin = begins[i].astype(int)
            end = ends[i].astype(int)
            index_a = np.argwhere((position - begin)>0).reshape(-1)
            index_b = np.argwhere((position - end)<0).reshape(-1)
            index_merge = [val for val in index_a if val in index_b]
            if col == 5 and np.size(index_merge)>0:
                major = As[i].astype(int)
                minor = Bs[i].astype(int)
            elif np.size(index_merge)>0:
                vaf = (np.sum(DP[index_merge])-np.sum(AD[index_merge]))/np.sum(DP[index_merge])
                major_t = int(Ts[i].astype(int) * vaf)
                minor_t = Ts[i].astype(int) - major_t
                major = np.max([major_t,minor_t])
                minor = np.min([major_t,minor_t])
                #print(str(major)+","+str(minor))
            else:
                major = 2
                minor = 0
            major = (np.ones([np.size(index_merge), 1]) * major).astype(int)
            minor = (np.ones([np.size(index_merge), 1]) * minor).astype(int)
            tmp = np.c_[SNV[index_merge, :],major,minor]
            a_chrom = np.append(a_chrom,tmp, axis=0)
        a_chrom = np.delete(a_chrom,0,axis = 0)
        answer = np.append(answer,a_chrom, axis=0)
    answer = np.delete(answer, 0, axis=0)
    return answer

def calculate_major_minor(CNV_table):
    if CNV_table.shape[1]==5:
        major = CNV_table[:,3].astype(int)
        minor = CNV_table[:,4].astype(int)
        total = major + minor
        CNV_table = np.c_[CNV_table,total]
        return 5,CNV_table
    elif CNV_table.shape[1]==4:
        return 4,CNV_table
    else:
        print("Error:\tUnable to parse CNV file.")
        exit(1)

def calculate_multiplicity(SNV_table,purity):
    AD = SNV_table[:, 2].astype(int)
    DP = SNV_table[:, 3].astype(int)
    total = SNV_table[:, 4].astype(int)+SNV_table[:,5].astype(int)
    multiplicity = np.round(AD / DP / purity * (total * purity + (1 - purity) * 2)).astype(int)
    tmp = np.c_[SNV_table[:, 5].astype(int),multiplicity]
    tmp = np.max(tmp,axis=1)
    tmp = np.where(tmp<=1,1,tmp)
    SNV_table[:, 4] = SNV_table[:, 4].astype(int)+tmp#SNV_table[:,5].astype(int)
    SNV_table[:, 5] = tmp
    return SNV_table

def check_outlier(SNV_table):
    AD = SNV_table[:,2].astype(int)
    DP = SNV_table[:,3].astype(int)
    total = SNV_table[:,4].astype(int)
    multiplicity = SNV_table[:,5].astype(int)
    sorce = 2/(multiplicity / AD * DP - total + 2)
    pass_index = np.argwhere((sorce>0) & (sorce<2)).reshape(-1)
    #out_index = np.argwhere(sorce >1.5).reshape(-1)
    return pass_index#,out_index

CNV_path,SNV_path,purity,OutPath,Prefix = main()
OutPath,Prefix = check_path(CNV_path,SNV_path,purity,OutPath,Prefix)
SNV_table = readfile_CNV_SNP(SNV_path)
#SNV_table: chrom position AD DP
SNV_table = check_SNP(SNV_table)
CNV_table = readfile_CNV_SNP(CNV_path)
col,CNV_table = calculate_major_minor(CNV_table)

merge_table = check_SNP_has_AB(SNV_table,CNV_table,col)
#merge_table: chrom position AD DP major minor
merge_table_m = calculate_multiplicity(merge_table,purity)
#merge_table: chrom position AD DP total multiplicity
pass_index= check_outlier(merge_table_m)
merge_table_m = merge_table_m[pass_index,:]
merge_out = merge_table[pass_index,:]
res = Ep.Elastic_GMM_subclone(merge_table_m,purity,Prefix)
Ep.clust_stdout(merge_out,res,OutPath,Prefix)
