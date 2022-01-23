from preprocess import *
import ElasticNet_GMM_process as Ep
import getopt
import sys

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

CNV_path,SNV_path,purity,OutPath,Prefix = main()
OutPath,Prefix,purity = check_path(CNV_path,SNV_path,purity,OutPath,Prefix)
SNV_table = readfile_CNV_SNP(SNV_path)
#SNV_table: chrom position AD DP
SNV_table = check_SNP(SNV_table)
CNV_table = readfile_CNV_SNP(CNV_path)
col,CNV_table = calculate_major_minor(CNV_table)
merge_table = check_SNP_has_AB(SNV_table,CNV_table,col)
#merge_table: chrom position AD DP major minor
merge_table_m,purity = calculate_multiplicity(merge_table,purity)
#merge_table: chrom position AD DP total multiplicity
pass_index,out_index= check_outlier(merge_table_m)
merge_table_m = merge_table_m[pass_index,:]
merge_out = merge_table[pass_index,:]
print("preprocess end")
print("ElasticNet_GMM clust start...")
distance, phi_new, mutation_num,n= Ep.Elastic_GMM_subclone(merge_table_m,purity,Prefix)
phi_corrected, class_label = Ep.corrected(distance,phi_new,n,mutation_num)
res = Ep.clust_GMM(phi_new,phi_corrected,mutation_num, has_sklearn=True)
Ep.clust_stdout(merge_out,phi_new,res,OutPath,Prefix)
print("ElasticNet_GMM clust end")
