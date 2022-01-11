import numpy as np
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

def check_path(CNV_path,SNV_path,purity,OutPath,Prefix):
    flag = 0
    if CNV_path=='':
        print("Error:\tMissing CNV path")
        flag = 1
    if SNV_path=='':
        print("Error:\tMissing SNP path")
        flag = 1
    if purity=='':
        print("Warning:\tMissing purity")
        purity = -1
    if OutPath=='':
        print("Warning:\tMissing output path, use default path './output'")
        OutPath = './output'
    if Prefix=='':
        print("Warning:\tMissing prefix of the output-file, use default prefix 'sample0'")
        Prefix = 'sample0'
    if flag != 0:
        exit(flag)
    return OutPath,Prefix,purity

def check_SNP(SNV_table):
    # check the chrom
    # check the read depth
    # check the minor and total > 0
    arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X','Y'])
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
    arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X','Y'])
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
    if purity == -1:
        vaf = AD / DP
        index=np.argwhere(vaf>0)
        vaf = vaf[index]
        purity = np.mean(np.abs(2*vaf-1))
    total = SNV_table[:, 4].astype(int)+SNV_table[:,5].astype(int)
    multiplicity = np.round(AD / DP / purity * (total * purity + (1 - purity) * 2)).astype(int)
    tmp = np.c_[SNV_table[:, 5].astype(int),multiplicity]
    tmp = np.max(tmp,axis=1)
    tmp = np.where(tmp<=1,1,tmp)
    SNV_table[:, 4] = SNV_table[:, 4].astype(int)+tmp
    SNV_table[:, 5] = tmp
    return SNV_table,purity

def check_outlier(SNV_table):
    AD = SNV_table[:,2].astype(int)
    DP = SNV_table[:,3].astype(int)
    total = SNV_table[:,4].astype(int)
    multiplicity = SNV_table[:,5].astype(int)
    sorce = 2/(multiplicity / AD * DP - total + 2)
    pass_index = np.argwhere((sorce>0) & (sorce<2)).reshape(-1)
    out_index =  np.argwhere((sorce<0) | (sorce>2)).reshape(-1)
    return pass_index,out_index

def merge_two_sampe(merge_table_list):
    arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X','Y'])
    answer_0 = np.array(['chrom', 'pos', 'AD', 'DP', 'major', 'minor']).reshape([1, 6])
    answer_1 = np.array(['chrom', 'pos', 'AD', 'DP', 'major', 'minor']).reshape([1, 6])
    for chr in arm:
        index_0 = np.argwhere((merge_table_list[0][:, 0] == chr)).reshape(-1)
        index_1 = np.argwhere((merge_table_list[1][:, 0] == chr)).reshape(-1)
        pos_0 = merge_table_list[0][index_0, 1]
        pos_1 = merge_table_list[1][index_1, 1]
        pos_merge = [val for val in pos_0 if val in pos_1]
        if len(pos_merge) <= 0:
            continue
        for pos in pos_merge:
            index = np.argwhere((merge_table_list[0][:, 0] == chr) & (merge_table_list[0][:, 1] == pos)).reshape(-1)
            answer_0=np.append(answer_0,merge_table_list[0][index,:],axis=0)
            index = np.argwhere((merge_table_list[1][:, 0] == chr) & (merge_table_list[1][:, 1] == pos)).reshape(-1)
            answer_1=np.append(answer_1,merge_table_list[1][index,:],axis=0)
    answer_0 = np.delete(answer_0, 0, axis=0)
    answer_1 = np.delete(answer_1, 0, axis=0)
    return answer_0,answer_1



