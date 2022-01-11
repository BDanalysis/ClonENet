import numpy as np
import getopt
import sys

def main():
    OutputPath = ""
    InputPath = ""
    Prefix = ""
    Purity = -1
    CNV = ""
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:c:x:p:", ["Help", "InputPath=","Output=", "CNV=","prefix=","purity="])
    for opts, arg in opts:
        if opts == "-h" or opts == "--Help" or opts == "--help":
            print("Required parameters:")
            print("-i or --input:", end="\t")
            print("Path to read the input file")
            print("-o or --output:", end="\t")
            print("Path to save the output file")
            print("-c or --cnv:", end="\t")
            print("Path to the cnv file")
            print("-x or --prefix:", end="\t")
            print("The prefix of the output file")
            print("-p or --purity:", end="\t")
            print("The prefix of the output file")
            exit()
        elif opts == "-i" or opts == "--Input" or opts == "--input":
            InputPath = arg
        elif opts == "-o" or opts == "--Output" or opts == "--output":
            OutputPath = arg
        elif opts == "-c" or opts == "--CNV" or opts == "--cnv":
            CNV = arg
        elif opts == "-x" or opts == "--Prefix" or opts == "--prefix":
            Prefix = arg
        elif opts == "-p" or opts == "--Purity" or opts == "--purity":
            Purity = float(arg)
    if OutputPath == "" or InputPath == "" or Prefix == "":
        print("Error:\tloss the parameters")
        exit(1)
    return InputPath,OutputPath,CNV,Prefix,Purity

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

def find_max(ary,clone_clust):
    clust_sort =  clone_clust[:,0].astype(int)
    min_index = 10
    for i in ary:
        index = np.argwhere(clust_sort==i).reshape([-1])
        index =index[0]
        if min_index >= index:
            min_index = index
    return min_index

def make_subclone_CN(InputPath,clone_clust,CNV,Prefix,purity):
    snp_clust = readfile_CNV_SNP(InputPath+'/'+Prefix+".Clust_pos.txt")
    if purity == -1:
        AD = snp_clust[:,2].astype(int)
        DP = snp_clust[:, 3].astype(int)
        vaf = AD/DP
        index=np.argwhere(vaf>0)
        vaf = vaf[index]
        purity = np.mean(np.abs(2*vaf-1))
    cnv = readfile_CNV_SNP(CNV)
    arm = np.array(['1','2','3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X'])#
    answer = np.array(['#chrom','start','end','clust','major','minor']).reshape([1,6])
    for chrom in arm:
        print("segment and merge chrom: "+str(chrom)+" by "+str(1000)+" base...")
        CNV = np.argwhere(cnv[:,0]==chrom).reshape(-1)
        SNV = np.argwhere(snp_clust[:,0]==chrom).reshape(-1)
        if np.size(CNV) <= 0 or np.size(SNV) <= 0:
            print("Information:\tThere is not any CNV or SNP in chrom " + str(chrom) + ".")
            continue
        CNV = cnv[CNV, 0:5].astype(int)
        SNV = snp_clust[SNV, 0:5].astype(int)
        CNV = CNV[np.argsort(CNV[:, 1]), :]
        SNV = SNV[np.argsort(SNV[:, 1]), :]
        position = SNV[:, 1]/1000
        clusts = SNV[:, 4]
        begins = CNV[:, 1]/1000
        ends = CNV[:, 2]/1000
        base_end = int(np.max(ends))+1
        reads_pair = SNV[:,2:4]
        clones = []
        i = 0
        while i < base_end:
            index_snp = np.argwhere(((position - i) > 0) & ((position - i) < 1)).reshape(-1)
            if index_snp.shape[0]<=0: # none subclone snp
                if i > 0:
                    list_chr = clones[-1]
                    if chrom == list_chr[0] and -1 == list_chr[3] and 2 == list_chr[4] and 0 == list_chr[5]:
                        clones[-1][2] = (i + 1) * 1000 - 1
                    else:
                        clones.append([chrom, i * 1000, (i + 1) * 1000 - 1, -1, 2, 0])
                else:
                    clones.append([chrom, i * 1000, (i + 1) * 1000 - 1, -1, 2, 0])
                i += 1
                continue
            index_cnv = np.argwhere(((begins - i) < 0) & ((ends - i) < 1)).reshape(-1)
            if index_cnv.shape[0]>0: # have subclone snp, none CNV
                if CNV.shape[1]>4:
                    major = np.mean(CNV[index_cnv, 3]).astype(int)
                    minor = np.mean(CNV[index_cnv, 4]).astype(int)
                else:
                    major = np.mean(CNV[index_cnv, 3]).astype(int)
                    minor = 0
            else:
                major = 2
                minor = 0
            pair = reads_pair[index_snp].astype(int)
            AD = np.sum(pair[:,0])
            DP = np.sum(pair[:,1])
            clust = find_max(clusts[index_snp],clone_clust)
            theta = AD/DP
            tmp = theta*(purity*(major+minor)+2*(1-purity))+purity*(major+minor)
            alpha = int(np.floor(tmp/(2*purity)))
            beta = int((major+minor)-alpha)
            if beta < 0 :
                beta = minor
            if i>0 :
                list_chr = clones[-1]
                if chrom == int(list_chr[0]) and clust == int(list_chr[3]) and alpha == int(list_chr[4]) and beta == int(list_chr[5]):
                    clones[-1][2] = (i+1)*1000-1
                else:
                    clones.append([chrom, i * 1000, (i + 1) * 1000 - 1, clust, alpha, beta])
            i+=1
        clones = np.array(clones)
        clones = clones.astype(int)
        answer = np.append(answer, clones, axis=0)
    return answer

def makeout(clones,clone_clust):
    clust_sort =  clone_clust[:,0].astype(int)
    clust_No = clones[:,3]
    answer = np.ones([len(clust_No)]).astype(str)
    i = len(clust_sort)
    sub_index = np.argwhere(clust_No == '-1').reshape(-1)
    answer[sub_index] = "clone"
    for i_c in clust_sort:
        sub_index = np.argwhere(clust_No==str(i_c)).reshape(-1)
        answer[sub_index] = "subclone"+str(i)
        clones[sub_index, 3] = i
        i-=1
    answer = np.c_[clones[:,0:3],answer]
    return clones,answer

InputPath,OutputPath,CNV,Prefix,purity = main()
clone_clust = readfile_CNV_SNP(InputPath + '/' + Prefix + ".summary_table.txt")
clones = make_subclone_CN(InputPath,clone_clust,CNV,Prefix,purity)
clones,answer = makeout(clones,clone_clust)

np.savetxt(OutputPath+'/'+Prefix+".subclone.txt",answer,fmt='%s\t%s\t%s\t%s')
np.savetxt(OutputPath+'/'+Prefix+".CNi.txt",clones,fmt='%s\t%s\t%s\t%s\t%s\t%s')





"""
    while i < range(np.size(begins)):
        begin = begins[i]
        end = ends[i]
        if end-begin<=1:
            end = end+1
        index_a = np.argwhere((position - begin) > 0).reshape(-1)
        index_b = np.argwhere((position - end) < 0).reshape(-1)
        index_merge = [val for val in index_a if val in index_b]
        if np.size(index_merge)<=0:
            clones.append(np.array([chrom, (begin).astype(int), (end).astype(int), -1, 2, 0]))
            continue
        pos = position[index_merge]
        a_b = A_B[index_merge].astype(int)
        a_total = int(np.mean(a_b[:,0]))
        b_minor = int(np.mean(a_b[:,1]))
        clust1 = clusts[index_merge]
        clust = np.unique(clust1)
        if np.size(clust)>1:
            for ii in range(len(clust1)):
                thiss = clust1[ii]
                if ii + 1 == len(clust1):
                    end = ends[i].astype(int)
                    clones.append(np.array(
                        [chrom, (begin * 10).astype(int), (end * 10).astype(int), thiss, a_b[ii,0], a_b[ii,1]]))
                    break
                next = clust1[1 + ii]
                if thiss != next:
                    end = pos[ii] - 1
                    clones.append(np.array(
                        [chrom, (begin * 10).astype(int), (end * 10).astype(int), thiss, a_total, b_minor]))
                    begin = pos[ii]
        else:
            clones.append(np.array([chrom, (begin*10).astype(int), (end*10).astype(int), int(clust),a_total,b_minor]))
        i+=1
"""


