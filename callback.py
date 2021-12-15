import numpy as np
import getopt
import sys
#import matplotlib.pyplot as plt

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

def makeout(clones,clone_clust,purity):
    tmp = -1*clone_clust[:,2].reshape(-1).astype(float)
    index_clust = np.argsort(tmp)
    clust_No = clones[:,3]
    answer = np.ones([len(clust_No)]).astype(str)
    i = 1
    for i_c in range(len(index_clust)):
        clust = int(clone_clust[index_clust[i_c],0])
        sub_index = np.argwhere(clust_No==clust).reshape(-1)
        if float(clone_clust[index_clust[i_c],2])>purity:
            answer[sub_index] = "clone"
            clones[sub_index,3] = -1
        else:
            answer[sub_index] = "subclone"+str(i)
            clones[sub_index, 3] = i
            i+=1
    answer = np.c_[clones[:,0:3],answer]
    return clones,answer

InputPath,OutputPath,CNV,Prefix,purity = main()
snp_clust = readfile_CNV_SNP(InputPath+'/'+Prefix+".Clust_pos.txt")
clone_clust = readfile_CNV_SNP(InputPath+'/'+Prefix+".summary_table.txt")
cnv = readfile_CNV_SNP(CNV)

arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X'])#
#answer = np.array(['chrom','pos','AD','DP','major','minor']).reshape([1,6])
base_x = 0
answer = np.array(['#chrom','start','end','clust','major','minor']).reshape([1,6])
output = np.array(['#chrom','start','end','cloneType']).reshape([1,4])
#plt.figure(figsize=(84,6))
for chrom in arm:
    CNV = np.argwhere(cnv[:,0]==chrom).reshape(-1)
    SNV = np.argwhere(snp_clust[:,0]==chrom).reshape(-1)
    if np.size(CNV) <= 0 or np.size(SNV) <= 0:
        print("Information:\tThere is not any CNV or SNP in chrom " + str(chrom) + ".")
        continue
    CNV = cnv[CNV, 0:5].astype(int)
    SNV = snp_clust[SNV, 0:5].astype(int)
    CNV = CNV[np.argsort(CNV[:, 1]), :]
    SNV = SNV[np.argsort(SNV[:, 1]), :]
    position = SNV[:, 1]/10
    clusts = SNV[:, 4]
    begins = CNV[:, 1]/10
    ends = CNV[:, 2]/10
    A_B = SNV[:,2:4]
    base_end = np.max(ends)
    clones = []
    for i in range(np.size(begins)):
        begin = begins[i]
        end = ends[i]
        index_a = np.argwhere((position - begin) > 0).reshape(-1)
        index_b = np.argwhere((position - end) < 0).reshape(-1)
        index_merge = [val for val in index_a if val in index_b]
        if np.size(index_merge)<=0:
            clones.append(np.array([chrom, (begin*10).astype(int), (end*10).astype(int), -1, 2, 0]))
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
                    #segment.append(np.array([begin, end, thiss]))
                    break
                next = clust1[1 + ii]
                if thiss != next:
                    end = pos[ii] - 1
                    clones.append(np.array(
                        [chrom, (begin * 10).astype(int), (end * 10).astype(int), thiss, a_total, b_minor]))
                    #segment.append(np.array([begin, end, thiss]))
                    begin = pos[ii]
            #segment = np.array(segment)
            #for seg in segment:
            #    pass
        else:
            clones.append(np.array([chrom, (begin*10).astype(int), (end*10).astype(int), int(clust),a_total,b_minor]))
    clones = np.array(clones)
    clones = clones.astype(int)
    clones,output_tmp = makeout(clones,clone_clust,purity)
    answer = np.append(answer, clones, axis=0)
    output=np.append(output,output_tmp, axis=0)
output = np.array(output)
np.savetxt(OutputPath+'/'+Prefix+".subclone.txt",output,fmt='%s\t%s\t%s\t%s')
np.savetxt(OutputPath+'/'+Prefix+".CNi.txt",answer,fmt='%s\t%s\t%s\t%s\t%s\t%s')


