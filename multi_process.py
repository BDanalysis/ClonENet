from preprocess import *
import ElasticNet_GMM_process as Ep
import getopt
import sys

def main():
    CNV_file_list = []
    SNV_file_list = []
    purity = []
    Prefix = ""
    OutPath = ""
    opts, args = getopt.getopt(sys.argv[1:], "hc:s:p:o:n:", ["help", "CNV=","SNV=","Name=","Purity=","OutPath=","Name="])
    if len(opts) == 0:
        print("Use analysis_Tree -h to get help")
        exit()
    for opts, arg in opts:
        if opts == "-h" or opts == "--help":
            print("Required parameters:")
            print("-c or --CNV:", end="\t")
            print("The path and file of input CNV, e.g. : /the/way/to/file/a.cnv.txt")
            print("-s or --SNV:", end="\t")
            print("The path and file of input SNV, e.g. : /the/way/to/file/a.snv.txt")
            print("-p or --Purity:", end="\t")
            print("Tumor purity")
            print("\t\t\t\tThe program allows you to enter multiple sample files")
            print("-o or --OutPath:", end="\t")
            print("Path to save the output file")
            print("-n or --Name:", end="\t")
            print("Prefix of the output file")
            exit()
        elif opts == "-c" or opts == "--CNV":
            CNV_file_list.append(str(arg))
        elif opts == "-p" or opts == "--Purity":
            purity.append(float(arg))
        elif opts == "-s" or opts == "--SNV":
            SNV_file_list.append(str(arg))
        elif opts == "-n" or opts == "--Name":
            Prefix=str(arg)
        elif opts == "-o" or opts == "--OutPath":
            OutPath = arg
    return CNV_file_list,SNV_file_list,purity,Prefix,OutPath

CNV_file_list,SNV_file_list,Puritys,Prefix,OutPath = main()
merge_table_list = []
i = 0
while i < len(CNV_file_list):
    CNV_file = CNV_file_list[i]
    SNV_file = SNV_file_list[i]
    if len(Puritys)<=0:
        Puritys.append(-1)
        Puritys.append(-1)
        purity = -1
    else:
        purity = Puritys[i]
    OutPath, Prefix, purity = check_path(CNV_file, SNV_file, purity, OutPath, Prefix)
    SNV_table = readfile_CNV_SNP(SNV_file)
    SNV_table = check_SNP(SNV_table)
    CNV_table = readfile_CNV_SNP(CNV_file)
    col, CNV_table = calculate_major_minor(CNV_table)
    merge_table = check_SNP_has_AB(SNV_table, CNV_table, col)
    merge_table_list.append(merge_table)
    i+=1
    # merge_table: chrom position AD DP major minor
merge_table_0,merge_table_1 = merge_two_sampe(merge_table_list)
merge_table0_m, purity0 = calculate_multiplicity(merge_table_0, Puritys[0])
merge_table1_m, purity1 = calculate_multiplicity(merge_table_1, Puritys[1])
distance0, phi_new0, mutation_num0,n0= Ep.Elastic_GMM_subclone(merge_table0_m, purity0,Prefix)
distance1, phi_new1, mutation_num1,n1= Ep.Elastic_GMM_subclone(merge_table1_m, purity1,Prefix)
phi_corrected0, class_label0 = Ep.corrected(distance0,phi_new0,n0,mutation_num0)
phi_corrected1, class_label1 = Ep.corrected(distance1,phi_new1,n1,mutation_num1)
phi_new_merge = np.c_[phi_new0.reshape([-1,1]),phi_new1.reshape([-1,1])]
phi_res_merge = np.c_[phi_corrected0.reshape([-1,1]),phi_corrected1.reshape([-1,1])]
res = Ep.clust_GMM_multisample(phi_new_merge,phi_res_merge)
Ep.clust_stdout_muliti(merge_table_0,res,OutPath,Prefix)
#res1 = Ep.Elastic_GMM_subclone(merge_table0_m,purity,Prefix)
# merge_table: chrom position AD DP total multiplicity
#pass_index, out_index = check_outlier(merge_table_m)
#merge_table_m = merge_table_m[pass_index, :]
#merge_out = merge_table[pass_index, :]
