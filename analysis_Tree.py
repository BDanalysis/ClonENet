import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc
import getopt
import sys
from scipy.stats import norm
from matplotlib import rcParams

config = {
    "font.family": 'serif', # 衬线字体
    "font.size": 15, # 相当于小四大小
    "font.serif": ['SimSun'], # 宋体
    "mathtext.fontset": 'stix', # matplotlib渲染数学字体时使用的字体，和Times New Roman差别不大
    'axes.unicode_minus': False # 处理负号，即-号
}
rcParams.update(config)
np.set_printoptions(suppress=True)
def main():
    Clust_pos_list=[]
    summary_table_list = []
    purity = ""
    OutPath = ""
    opts, args = getopt.getopt(sys.argv[1:], "hf:o:", ["help", "Purity=","=File","OutPath="])
    if len(opts) == 0:
        print("Use analysis_Tree -h to get help")
        exit()
    for opts, arg in opts:
        if opts == "-h" or opts == "--help":
            print("Required parameters:")
            print("-f or --File:", end="\t")
            print("The path and prefix of input file, e.g. : /path/to/file/sample0")
            print("\t\t\t\tThe program allows you to enter multiple sample files")
            print("-o or --OutPath:", end="\t")
            print("Path to save the output file")
            exit()
        elif opts == "-f" or opts == "--CNV":
            Clust_pos_list.append(str(arg)+".Clust_pos.txt")
            summary_table_list.append(str(arg)+".summary_table.txt")
        elif opts == "-p" or opts == "--Purity":
            purity = float(arg)
        elif opts == "-o" or opts == "--OutPath":
            OutPath = arg
    return Clust_pos_list,summary_table_list,purity,OutPath

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

def pos_Hashmap(path,pfx):
    with open(path, 'r') as f1:
        CNV_data = f1.readlines()
    position = {}
    for i in range(len(CNV_data)):
        if '#' in CNV_data[i]:
            continue
        CNV_data[i] = CNV_data[i].rstrip('\n')
        tmp = CNV_data[i].split("\t")
        key = str(tmp[0])+":"+str(tmp[1])
        value=[str(pfx)+":"+str(int(tmp[4])+1)+","+str(tmp[5])]
        dict={key:value}
        position.update(dict)
    return position

def calling_clone_type(path_list):
    clone_position={}
    for i in range(len(path_list)):
        path = path_list[i]
        sample = pos_Hashmap(path,i)
        tmp_set1 = set(clone_position)
        tmp_set2 = set(sample)
        key_list=tmp_set1.intersection(tmp_set2)
        if len(key_list)==0:
            clone_position.update(sample)
            continue
        for key in key_list:
            tmp_list = clone_position[key]+sample[key]
            sample.pop(key)
            clone_position[key] = tmp_list
        clone_position.update(sample)
    return len(path_list),clone_position

def sort_chr_position(clone_position,sample_num):
    keys = clone_position.keys();
    answer = []
    for key in keys:
        pos_list = key.split(':')
        pos_list = np.array(pos_list).astype(int)
        CCF = np.zeros([sample_num*2])
        values = clone_position[key]
        for i in range(len(values)):
            value = np.array(values[i].split(':'))
            tmp_message = value[-1].split(',')
            CCF[int(value[0])*2]= int(tmp_message[0])
            CCF[int(value[0])*2+1] = float(tmp_message[1])
        pos_list=np.r_[pos_list,CCF].reshape(-1)
        answer.append(pos_list)
    answer = np.array(answer)
    index  = np.lexsort(np.transpose(answer)[::-1])
    answer = answer[index]
    return answer

def get_clust_center(path_list):
    answer = []
    for i in range(len(path_list)):
        path = path_list[i]
        position={}
        with open(path, 'r') as f1:
            Clust_data = f1.readlines()
        for i in range(len(Clust_data)):
            if '#' in Clust_data[i]:
                continue
            Clust_data[i] = Clust_data[i].rstrip('\n')
            tmp = Clust_data[i].split("\t")
            key = str(tmp[0])
            value = [int(tmp[1]),float(tmp[2])]
            dict = {key: value}
            position.update(dict)
        answer.append(position)
    return answer

def merge_matrix(pos_value,table_center):
    a = pos_value[:,2:6]
    #plt.scatter(a,b,alpha=0.2)
    table = []
    points = []
    zeros = np.zeros([int(np.max(a[:,0]))+1,int(np.max(a[:,2]))+1])
    for i_x,value_x,i_y,value_y in a:
        i_x = int(i_x)
        i_y = int(i_y)
        zeros[i_x, i_y] += 1
    for i_x,value_x,i_y,value_y in a:
        i_x = int(i_x)-1
        i_y = int(i_y)-1
        type = i_x+i_y * 1j
        tmp_point = np.array([value_x,value_y,type])
        if i_x < 0 :
            tmp_clust = np.array([0, table_center[1][str(i_y)][1],zeros[i_x+1,i_y+1]])
        elif i_y < 0 :
            tmp_clust = np.array([table_center[0][str(i_x)][1], 0,zeros[i_x+1,i_y+1]])
        else:
            tmp_clust = np.array([table_center[0][str(i_x)][1], table_center[1][str(i_y)][1],zeros[i_x+1,i_y+1]])
        table.append(tmp_clust)
        points.append(tmp_point)
    table = np.array(table)
    points = np.array(points)
    return table,points

def merge_clust_by_matrix(matrix):
    x = matrix[:, 0] + matrix[:, 1] * 1j
    idx = np.unique(x,return_index=True)[1]
    c = matrix[idx]
    leave = []
    for lis in c:
        x = matrix[:, 0] + matrix[:, 1] * 1j
        map = lis[0]+lis[1]* 1j
        num = np.size(np.argwhere(x==map))
        if num >=5:
            leave.append([lis[0],lis[1],num])
    leave = np.array(leave)
    for i in range(leave.shape[0]):
        line = leave[i,:]
        if line[0]==0:
            col1 = leave[:,1]
            idx = np.argwhere(col1 == line[1])
            if np.max(leave[idx,0])>0:
                leave[i,:]=-1
        elif line[1]==0:
            col1 = leave[:, 0]
            idx = np.argwhere(col1 == line[0])
            if np.max(leave[idx, 1]) > 0:
                leave[i, :] = -1
    leave = leave[np.all(leave >= 0, axis=1), :]
    #leave = np.where(leave == 0,2,leave)
    return leave

def make_Tree(leave):
    x = np.sort(np.unique(leave[:,0]))
    y = np.sort(np.unique(leave[:,1]))
    matrix = np.zeros([np.size(x),np.size(y)])
    matrix_dis = np.zeros([np.size(x),np.size(y)])-1
    map = leave[:, 0] + leave[:, 1] * 1j
    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            idx = np.argwhere(map==(x[i]+y[j]*1j)).reshape(-1)
            if np.size(idx)==1:
                matrix[i,j] = np.mean(leave[idx,0:2])
                matrix_dis[i, j] = i+j
    Tree={}
    for child in np.unique(matrix_dis):
        if child >= 0 and child < np.max(matrix_dis):
            lists = np.argwhere(matrix_dis==child)
            for i_x,i_y in lists:
                matrix_tmp = np.where(matrix_dis<=child,100,matrix)
                dis_me = matrix_tmp-matrix[i_x,i_y]
                #dis = np.where(dis_me<=0,100,dis_me)
                father = np.argwhere(dis_me == np.min(dis_me))[0]
                key=str(i_x)+","+str(i_y)
                value = str(father[0])+","+str(father[1])
                Tree.update({key:value})
    for key in Tree.keys():
        ix = int(key.split(",")[0])
        iy = int(key.split(",")[1])
        value = Tree[key]
        fx = int(value.split(",")[0])
        fy = int(value.split(",")[1])
        print(str(x[ix])+","+str(y[iy]),end=" ")
        print("->", end=" ")
        print(str(x[fx]) + "," + str(y[fy]))

def fq(list):
    max = np.max(list)
    min = np.min(list)
    list = list - min
    space = np.arange(min,max,0.01)
    fq = []
    for i in range(np.size(space)):
        tmp_list= list-((i+1)*0.01)
        index = np.argwhere(tmp_list>0).reshape(-1)
        fq.append(np.size(list)-np.size(index))
        if np.size(index) > 0:
            list = list[index]
    return fq

def drow_bar(Clust_pos_list,colors,as_plot1,as_plot2):
    table1 = read_profile(Clust_pos_list[0])
    clusts = table1[:,4].astype(int)
    phi_x = table1[:, 5].astype(float)
    fqlist = []
    for c_id in np.unique(clusts):
        indx = np.argwhere(c_id == clusts).reshape(-1)
        x = phi_x[indx]
        sigma = np.std(x)
        mu = np.mean(x)
        x = np.linspace(-0.01, 1.51, 152)
        y = norm(mu, sigma).pdf(x)
        fqlist.append(np.max(y))
    for c_id in np.unique(clusts):
        indx = np.argwhere(c_id==clusts).reshape(-1)
        x = phi_x[indx]
        sigma = np.std(x)
        mu = np.mean(x)
        x = np.linspace(-0.01,1.51, 152)
        y = norm(mu, sigma).pdf(x)

        as_plot1.plot(x, y, color=colors[c_id], linestyle='-.')
    table2 = read_profile(Clust_pos_list[1])
    clusts = table2[:, 4].astype(int)
    phi_x = table2[:, 5].astype(float)
    fqlist = []
    for c_id in np.unique(clusts):
        indx = np.argwhere(c_id == clusts).reshape(-1)
        x = phi_x[indx]
        sigma = np.std(x)
        mu = np.mean(x)
        x = np.linspace(-0.01, 1.51, 152)
        y = norm(mu, sigma).pdf(x)
        fqlist.append(np.max(y))
    for c_id in np.unique(clusts):
        indx = np.argwhere(c_id == clusts).reshape(-1)
        x = phi_x[indx]
        sigma = np.std(x)
        mu = np.mean(x)
        x = np.linspace(-0.01, 1.51, 152)
        y = norm(mu, sigma).pdf(x)
        as_plot2.plot(y, x, color=colors[c_id], linestyle='-.')

plt.figure(1, figsize=(12,12))
gs = gsc.GridSpec(3, 3)
plt.subplots_adjust(wspace=0.0,hspace=0.0)
as_raw = plt.subplot(gs[1:3,0:2])
as_raw.set_xlim(-0.01, 1)
as_raw.set_ylim(-0.01, 1)
as_raw.set_xlabel("CCF of sample A",fontsize=18)
as_raw.set_ylabel("CCF of sample B",fontsize=18)
as_plot1 = plt.subplot(gs[0,0:2])
as_plot1.set_ylabel("Frequency",fontsize=18)
as_plot1.spines["right"].set_visible(False)
as_plot1.spines["top"].set_visible(False)
as_plot1.spines["left"].set_visible(False)
as_plot1.set_xticks([])
as_plot1.set_yticks([])
as_plot2 = plt.subplot(gs[1:3,2])
as_plot2.set_xlabel("Frequency")
as_plot2.spines["right"].set_visible(False)
as_plot2.spines["bottom"].set_visible(False)
as_plot2.spines["top"].set_visible(False)
as_plot2.set_xticks([])
as_plot2.set_yticks([])
colors = sns.color_palette(n_colors=30)
Clust_pos_list,summary_table_list,purity,OutPath=main()
drow_bar(Clust_pos_list,colors,as_plot1,as_plot2)
sample_num,clone_position = calling_clone_type(Clust_pos_list)
pos_value = sort_chr_position(clone_position,sample_num)
table_center = get_clust_center(summary_table_list)
table,points = merge_matrix(pos_value,table_center)
leave = merge_clust_by_matrix(table)
make_Tree(leave)
x2 = points[:, 0].astype(float)
y2 = points[:, 1].astype(float)
types = points[:,2]
typ = np.unique(types)
clas = np.zeros([len(types)])
for i in range(len(typ)):
    index = np.argwhere(types==typ[i]).reshape(-1)
    x2_1 = x2[index]
    y2_1 = y2[index]
    as_raw.scatter(x2_1,y2_1,c=colors[i],alpha=0.8,s=4)
    if len(index) > 5:
        cen_x = np.mean(x2_1)
        cen_y = np.mean(y2_1)
        as_raw.scatter(cen_x, cen_y, c="black", alpha=0.8, s=12)

plt.savefig("./output_real/merge6.pdf",format="pdf")
plt.close()