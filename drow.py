import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc
import numpy as np
import GMM_clust as GMM
from scipy.stats import norm
import getopt
import sys

def main():
    OutputPath = ""
    InputPath = ""
    Prefix = ""
    Purity = -1
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:x:p:", ["Help", "Path=","CNV=", "SNP=","OutPath=","Prefix="])
    for opts, arg in opts:
        if opts == "-h" or opts == "--Help" or opts == "--help":
            print("Required parameters:")
            print("-i or --input:", end="\t")
            print("Path to read the input file")
            print("-o or --output:", end="\t")
            print("Path to save the output file")
            print("-x or --prefix:", end="\t")
            print("The prefix of the output file")
            print("-p or --purity:", end="\t")
            print("The prefix of the output file")
            exit()
        elif opts == "-i" or opts == "--Input" or opts == "--input":
            InputPath = arg
        elif opts == "-o" or opts == "--Output" or opts == "--output":
            OutputPath = arg
        elif opts == "-x" or opts == "--Prefix" or opts == "--prefix":
            Prefix = arg
        elif opts == "-p" or opts == "--Purity" or opts == "--purity":
            Purity = arg
    if OutputPath == "" or InputPath == "" or Prefix == "":
        print("Error:\tloss the parameters")
        exit(1)
    return InputPath,OutputPath,Prefix,Purity

def read_profile(path1):
    with open(path1, 'r') as f1:
        CNV_data = f1.readlines()
    data1 = []
    for i in range(len(CNV_data)):
        CNV_data[i] = CNV_data[i].rstrip('\n')
        if '#' in CNV_data[i] or 'c' in CNV_data[i]:
            continue
        tmp=CNV_data[i].split('\t')
        data1.append(tmp)
    table = np.array(data1)
    return table

def read_plug(path):
    with open(path, 'r') as f1:
        CNV_data = f1.readlines()
    data1 = []
    for i in range(len(CNV_data)):
        CNV_data[i] = CNV_data[i].rstrip('\n')
        tmp=float(CNV_data[i])
        data1.append(tmp)
    table = np.array(data1)
    return table

def drow_copynum(clones,path):
    base_x = 0
    clone_main = -1  # split_clone_sub(clone_clust)
    arm = np.array(
        ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22', 'X','Y'])
    subplots = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,23,24]
    clones_tmp = clones
    plt.figure(figsize=(24, 36))
    iis = 0
    for chrom in arm:
        print(chrom)
        idx = np.argwhere(clones_tmp[:, 0] == chrom).reshape(-1)
        if len(idx) <= 0 :
            continue
        clones = clones_tmp[idx,:].astype(int)
        begins = clones[:, 1] / 1000
        ends = (clones[:, 2] - 1) / 1000
        base_end = np.max(ends)
        clusts = clones[:, 3]
        A = clones[:, 4]+clones[:, 5]
        B = clones[:, 5]
        point1_A = np.c_[begins, A, B, clusts].reshape((-1, 4))
        point2_A = np.c_[ends, A, B, clusts].reshape((-1, 4))
        points_A = np.c_[point1_A, point2_A].reshape((-1, 4))
        flag = False
        colors_1 = ["darkgreen", "lawngreen"]
        colors_2 = ["sienna", "yellow"]
        for i in range(np.shape(points_A)[0] - 1):
            if flag:
                flag = False
                continue
            line = points_A[i:i + 2, :]
            x = line[:, 0] + base_x
            y_1 = line[:, 1]
            y_2 = line[:, 2]
            if int(clone_main) == line[0, 3]:
                c = 0
            else:
                continue
            plt.subplot(6,4,subplots[iis])
            plt.plot(x, y_1, linestyle='-',  color=colors_1[c], alpha=1, linewidth=3)#marker='.',
            plt.plot(x, y_2, linestyle='-',color=colors_2[c], alpha=1, linewidth=5) #marker='.',
            flag = True

        for i in range(np.shape(points_A)[0] - 1):
            if flag:
                flag = False
                continue
            line = points_A[i:i + 2, :]
            x = line[:, 0] + base_x
            y_1 = line[:, 1]
            y_2 = line[:, 2]
            if int(clone_main) == line[0, 3]:
                continue
            else:
                c = 1
            plt.subplot(6,4,subplots[iis])
            plt.plot(x, y_1, linestyle='-',  color=colors_1[c], alpha=1, linewidth=3)#marker='.',
            plt.plot(x, y_2, linestyle='-',color=colors_2[c], alpha=1, linewidth=5) #marker='.',
            flag = True
        plt.text(base_x, -0.2, "chr" + str(chrom), rotation=45, verticalalignment='center')
        plt.axvline(x=base_x, ls=":", c="black")
        iis += 1
    plt.xticks([])
    plt.savefig(path, format='pdf')
    plt.close()

def drow_clust(path,result,vaf,phi,clust=-1):
    gs = gsc.GridSpec(2, 5)
    lab = np.unique(result)
    plt.figure(1, figsize=(20,10))
    plt.subplots_adjust(wspace=0.3,hspace=0.0)
    as_raw = plt.subplot(gs[0,0:2])
    num_bins = int((np.max(vaf)-np.min(vaf))/0.01)
    #as_raw.plot(vaf.reshape([-1,1]),linestyle='',marker='.')
    as_raw.hist(vaf, num_bins, alpha=0.6,color="blue")
    #as_raw.set_xlim([-0.01,1.01])
    as_raw.yaxis.set_ticks_position('left')
    as_raw.set_title("CCF of raw data")
    #as_raw.set_xlabel("CCF")
    as_phi = plt.subplot(gs[1, 0:2])
    as_phi.hist(phi, num_bins,alpha=0.8, color="navy")#alpha=0.6,
    #as_phi.plot(phi.reshape([-1,1]),linestyle='',marker='.')
    #as_phi.set_xlim([-0.01, 1.01])
    as_phi.yaxis.set_ticks_position('left')
    #as_phi.set_title("Corrected CCF")
    as_phi.set_xlabel("Corrected CCF")
    colors = GMM.colors(np.max(lab)+1)
    as_clust = plt.subplot(gs[0:3, 2:5])
    as_clust.yaxis.set_ticks_position('right')
    #n, bins, patches = as_clust.hist(vaf, num_bins,color="dimgrey")
    for i in np.unique(result):
        c_list = vaf[result == i]
        if np.max(c_list)-np.min(c_list)==0:
            c_list = np.r_[c_list,np.min(c_list)-0.01]
        x = np.array(c_list)
        num_bins = int((np.max(c_list)-np.min(c_list))/0.01) # 直方图柱子的数量
        n, bins, patches = as_clust.hist(c_list, num_bins,color='k')#,color=colors[i],alpha=0.6
        position = np.linspace(bins[0],bins[-1],np.size(n))
        mu = float(np.mean(position[n==np.max(n)]))
        #mu = np.mean(x)
        if clust !=-1:
            mu = float(clust[i])
        sigma = np.std(x)*0.7#round(np.max([np.std(x),0.01]),2)
        x = np.linspace(-0.01,1.01, 102)
        y = norm(mu, sigma).pdf(x)
        print(mu, end='\t')
        print(sigma)
        y = y*(np.max(n)/np.max(y))
        as_clust.plot(x, y, color=colors[i],linestyle='-')  # 绘制y的曲线
        as_clust.axvline(x=mu,ls=":",c=colors[i])#colors[i]
    as_clust.set_xlim([-0.01,1.01])  # 绘制x轴
    as_clust.set_ylabel('Probability')  # 绘制y轴
    as_clust.set_xlabel("CCF")
    as_clust.set_title("Clustering results")
    plt.savefig(path)
    #plt.show()
    plt.close()

InputPath,OutputPath,Prefix,Purity = main()
table1 = read_profile(InputPath+'/'+Prefix+'.Clust_pos.txt')
clones = read_profile(InputPath+'/'+Prefix+'.CNi.txt')
summary = read_profile(InputPath+'/'+Prefix+'.summary_table.txt')
result = table1[:,4].reshape([-1]).astype(int)
phi = table1[:,5].reshape([-1]).astype(float)
major = table1[:,2].reshape([-1]).astype(int)
minor = table1[:,3].reshape([-1]).astype(int)
table2 = read_plug(InputPath+'/'+Prefix+'.vaf.txt')
vaf = table2.reshape([-1])
print("========================Elastic_net clust===========================")
drow_clust(OutputPath+'/clust_'+Prefix+".png",result,vaf,phi)
drow_copynum(clones,OutputPath+'/CN_'+Prefix+".pdf")
