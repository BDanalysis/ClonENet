import numpy as np
import getopt
import sys

def main():
    truth=""
    pres=""
    opts, args = getopt.getopt(sys.argv[1:], "g:p:", ["ground", "perpresice"])
    for opts, arg in opts:
        if opts == "-g":
            truth = arg
        elif opts == "-p":
            pres = arg
    if truth=="" or pres=="":
        exit(1)
    return truth,pres

def readfile_data(path):
    with open(path, 'r') as f1:
        data = f1.readlines()
    table = np.array([-1,-1,-1,-1]).reshape([1,4])
    points = []
    type = "clone0"
    for i in range(len(data)):
        if "->" in data[i]:
            data[i] = data[i].rstrip('\n')
            point = data[i].split("->")
            points = np.append(points,point, axis=0)
            point = int(point[-1])
            if point != 0:
                type = "subclone" + str(point)
        else:
            data[i] = data[i].rstrip('\n')
            position= data[i].split(" ")[0]
            position = position.lstrip('(')
            position = position.rstrip(')')
            position = position.split(",")
            chrom = data[i].split(" ")[-1]
            chrom = chrom.split("-")[-1].lstrip('chr')
            position.append(chrom)
            position.append(type)
            position = np.array(position).reshape([1,4])
            table = np.r_[table,position]
    table = np.array(table)
    return table

def readfile_prep(path):
    with open(path, 'r') as f1:
        data = f1.readlines()
    table = np.array([-1,-1,-1,-1]).reshape([1,4])
    for i in range(len(data)):
            data[i] = data[i].rstrip('\n')
            llise = data[i].split("\t")
            line = np.array([llise[1],llise[2],llise[0],llise[3]]).reshape(1,4)
            table = np.r_[table,line]
    table = np.array(table)
    return table

def mapping(truth,prep):
    arm = np.array(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
         '21', '22', 'X'])
    #prep_num = prep.shape[0]
    #truth_num = truth.shape[0]
    #num = 0
    truth_l = 0
    prep_l = 0
    mapping_l = 0
    for chr in arm:
        index_truth = np.argwhere(str(chr) == truth[:,2]).reshape(-1)
        index_prep = np.argwhere(str(chr) == prep[:,2]).reshape(-1)
        T = truth[index_truth,:]
        T = T[np.argsort(T[:,0]),:]
        P = prep[index_prep,:]
        for i in range(T.shape[0]):
            start_t = T[i,0].astype(int)
            end_t = T[i,1].astype(int)
            type_t= T[i,3]
            for ii in range(P.shape[0]):
                start_p = P[ii,0].astype(int)
                end_p = P[ii, 1].astype(int)
                type_p = P[ii,3]
                if start_p < end_t and start_t < end_p:
                    start = np.max([start_t,start_p])
                    end = np.min([end_p,end_t])
                    a = (end-start)/(end_t - start_t)
                    b = (end-start)/(end_p - start_p)
                    F1 = 2 * a * b / (a + b)
                    if b > 0.9:
                        truth_l = truth_l + (end_t - start_t)/10
                        prep_l = prep_l + (end_p - start_p)/10
                        if "sub" in type_p and "sub" in type_t:
                            mapping_l = mapping_l + (end-start)/10
                        if ("sub" not in type_t) and ("sub" not in type_p):
                            mapping_l = mapping_l + (end - start)/10
                        break
    a = mapping_l / truth_l
    b = mapping_l/prep_l
    print(a)
    print(b)
    print(2 * a * b / (a + b))
    print("================================================")


#truth,prep = main()
preps = ["sim2_np1c0p7c1p2_clone.txt","sim2_np2c1p6c2p2_clone.txt","sim4_np1c0p5c2p2c3p2_clone.txt","sim4_np1c1p4c2p3c3p2_clone.txt","sim4_np2c0p4c1p3c2p1_clone.txt","sim4_np3c0p7_clone.txt"]
truths=["sim2","sim2","sim4","sim4","sim4","sim4"]
for i in range(len(preps)):
    prep = "./output/"+preps[i]
    truth = "./example/fasta_"+truths[i]+"/tumor.dot"
    table = readfile_data(truth)
    prep = readfile_prep(prep)
    mapping(table,prep)
