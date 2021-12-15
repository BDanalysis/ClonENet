import numpy as np
import matplotlib.pyplot as plt
def readfile_CNV_SNP(path):
    with open(path, 'r') as f1:
        CNV_data = f1.readlines()
    CNV_table = []
    for i in range(len(CNV_data)):
        if i == 0:
            continue
        CNV_data[i] = CNV_data[i].rstrip('\n')
        tmp=CNV_data[i].split("\t")
        CNV_table.append(tmp)
    table = np.array(CNV_table)
    return table

def split_clone_sub(clone_clust):
    clust_No = clone_clust[:,0]
    balance = clone_clust[:,2]
    index_max = np.argmax(balance)
    return clust_No[index_max]

snp_clust = readfile_CNV_SNP("./CLiP/hx_019_CLiP_position.txt")
clone_clust = readfile_CNV_SNP("./CLiP/hx_019_clust_summary.txt")
cnv = readfile_CNV_SNP("./WGS_input/hx_005_cnv.txt")
clone_main = split_clone_sub(clone_clust)
arm = np.array(['1','2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22','X'])#
answer = np.array(['chrom','pos','AD','DP','major','minor']).reshape([1,6])
base_x = 0
plt.figure(figsize=(84,6))
for chrom in arm:
    CNV = np.argwhere(cnv[:,0]==chrom).reshape(-1)
    SNV = np.argwhere(snp_clust[:,0]==chrom).reshape(-1)
    if np.size(CNV) <= 0 or np.size(SNV) <= 0:
        print("Information:\tThere is not any CNV or SNP in chrom " + str(chrom) + ".")
        continue
    CNV = cnv[CNV, 0:5].astype(int)
    SNV = snp_clust[SNV, 0:3].astype(int)
    CNV = CNV[np.argsort(CNV[:, 1]), :]
    SNV = SNV[np.argsort(SNV[:, 1]), :]
    position = SNV[:, 1]/10
    A_B = CNV[:, 3]
    clusts = SNV[:, 2]
    begins = CNV[:, 1]/10
    ends = CNV[:, 2]/10
    base_end = np.max(ends)
    clones = []

    for i in range(np.size(begins)):
        begin = begins[i].astype(int)
        end = ends[i].astype(int)
        index_a = np.argwhere((position - begin) > 0).reshape(-1)
        index_b = np.argwhere((position - end) < 0).reshape(-1)
        index_merge = [val for val in index_a if val in index_b]
        if np.size(index_merge)<=0:
            continue
        pos = position[index_merge]
        a_b = A_B[i].astype(int)
        clust1 = clusts[index_merge]
        clust = np.unique(clust1)
        segment = []
        if np.size(clust)>1:
            for ii in range(len(clust1)):
                thiss = clust1[ii]
                if ii+1 == len(clust1):
                    end = ends[i].astype(int)
                    segment.append(np.array([begin,end,thiss]))
                    break
                next = clust1[1 + ii]
                if thiss != next:
                    end = pos[ii]-1
                    segment.append(np.array([begin,end,thiss]))
                    begin = pos[ii]
            max_length = 0
            max_clust = -1
            segment = np.array(segment)
            for lis in clust:
                indexs = np.argwhere(segment[:, 2] == lis).reshape(-1)
                length = segment[indexs,0:2]
                length = np.sum(length[:,1]-length[:,0])
                if length > max_length:
                    max_length = length
                    max_clust = lis
            #    print(str(chrom)+'\t'+str(lis[0])+'\t'+str(lis[1])+'\t'+str(lis[2]))
            clones.append(np.array([chrom,begins[i].astype(int),ends[i].astype(int),max_clust,a_b,1]))
        else:
            #print(str(chrom) + '\t' + str(begin) + '\t' + str(end) + '\t' + str(clust))
            clones.append(np.array([chrom, begin, end, int(clust),a_b,1]))
    clones = np.array(clones)
    clones = clones.astype(int)
    #print(clones)
    begins = clones[:, 1]
    ends = clones[:, 2]-1
    clusts = clones[:,3]
    A = clones[:,4]
    B = clones[:,5]
    point1_A = np.c_[begins,A,B,clusts].reshape((-1,4))
    point2_A = np.c_[ends,A,B,clusts].reshape((-1,4))
    points_A = np.c_[point1_A,point2_A].reshape((-1,4))
    flag = False
    colors_1 = ["darkgreen", "sienna"]
    colors_2 = ["green", "goldenrod"]
    for i in range(np.shape(points_A)[0]-1):
        if flag:
            flag = False
            continue
        line = points_A[i:i+2,:]
        print(line)
        x = line[:,0]+base_x
        y_1 = line[:,1]
        y_2 = line[:,2]
        if int(clone_main) == line[0, 3]:
            c = 0
        else:
            c = 1
        plt.plot(x,y_1, linestyle='-', marker='.', color=colors_1[c], alpha=1,linewidth=5)
        plt.plot(x,y_2, linestyle='-', marker='.', color=colors_2[c], alpha=1, linewidth=5)
        flag = True
    plt.text(base_x,0.8, "chr" + str(chrom), rotation=45, verticalalignment='center')
    base_x = base_x + base_end
    print(base_x)
    plt.axvline(x=base_x, ls=":", c="black")
plt.xticks([])
plt.savefig("./output/hx_005.pdf", format='pdf')
plt.show()
plt.close()

