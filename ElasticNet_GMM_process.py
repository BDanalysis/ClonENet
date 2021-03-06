import numpy as np
#import GMM_clust as GMM
from sklearn.mixture import GaussianMixture as GMM
import scipy
from scipy.special import logit
from scipy.special import expit
import gc

def update(x,alpha,lam,tau,r=0.5):
    flag = np.where(x>0,-1,1)
    tau_line = tau.flatten()
    val = x*alpha + flag * lam * r + tau_line
    val = np.where(x==0,0,val)
    return val/(2*lam*(1-r)*alpha)
"""
r: tumor RD
n: total RD
minor: tumor CN
total: total CN
ploidy: ploidy -> normal:2
"""

def fq(list):
    max = np.max(list)
    min = np.min(list)
    if max == min:
        return list[0]
    list = list-min
    space = np.arange(min,max+0.01,0.01)
    fq = []
    for i in range(np.size(space)):
        tmp_list= list-((i+1)*0.01)
        index = np.argwhere(tmp_list>0).reshape(-1)
        fq.append(np.size(list)-np.size(index))
        if np.size(index) > 0:
            list = list[index]
    if len(fq)>0:
        top = int(np.mean(np.argmax(fq)))
        top = space[top]
    else:
        return round(np.mean(list), 2)
    return top

def ElasticNet_correct(r, n, multiplicity, total, ploidy, Lambda, alpha, rho, Run_limit, precision,control_large, post_th,purity,Prefix):
    mutation_num = len(r)
    print("ElasticNet:"+str(mutation_num))
    vaf = r/n  # 次等位基因数/总数 thate
    phi = vaf * ((ploidy - purity * ploidy + purity * total) / multiplicity) #  CCF 文章中用 phi 表示， CP = CCF*purity
    tmp = np.where(phi>=1,phi,1)
    phi_new = phi/tmp                                               # CCF 最大为1
    phi_new[phi_new > expit(control_large)] = expit(control_large)  # 极大控制--由置信区间控制？
    phi_new[phi_new < expit(-control_large)] = expit(-control_large)  # 负极大控制
    omiga_new = logit(phi_new)
    omiga_new[omiga_new > control_large] = control_large  # 极大控制
    omiga_new[omiga_new < -control_large] = -control_large  # 极大控制
    distance = np.subtract.outer(omiga_new, omiga_new)
    index_dis = np.triu_indices(distance.shape[1], 1)
    eta = distance[index_dis] #  ωi − ωj − ηij = 0

    tau = np.ones((int(mutation_num * (mutation_num - 1) / 2), 1))   # 拉格朗日乘子的参数？
    col_id = np.append(np.array(range(int(mutation_num * (mutation_num - 1) / 2))),
                       np.array(range(int(mutation_num * (mutation_num - 1) / 2))))
    row1 = np.zeros(int(mutation_num * (mutation_num - 1) / 2))
    row2 = np.zeros(int(mutation_num * (mutation_num - 1) / 2))
    starting = 0
    for i in range(mutation_num - 1):
        row1[starting:(starting + mutation_num - i - 1)] = i
        row2[starting:(starting + mutation_num - i - 1)] = np.array(range(mutation_num))[(i + 1):]
        starting = starting + mutation_num - i - 1
    row_id = np.append(row1, row2)
    vals = np.append(np.ones(int(mutation_num * (mutation_num - 1) / 2)),
                     -np.ones(int(mutation_num * (mutation_num - 1) / 2)))

    DELTA = scipy.sparse.coo_matrix((vals, (row_id, col_id)),
                                 shape=(mutation_num, int(mutation_num * (mutation_num - 1) / 2)))
    del row_id, col_id,row1,row2,vals,phi
    gc.collect()
    k = 0  # iterator
    residual = 100
    res_last = 100
    eta_new = eta
    while residual > precision and k < Run_limit :
        k = k+1
        omiga_old = omiga_new
        #eta_new = eta
        theta = np.exp(omiga_old) * multiplicity / (2 + np.exp(omiga_old) * total)
        A = np.sqrt(n) *(theta-vaf)/np.sqrt(theta * (1 - theta))
        B = np.sqrt(n) * theta/np.sqrt(theta * (1 - theta))
        omiga_tmp_2 = (DELTA.tocsr() * np.array((alpha * eta_new + tau.T).T)).flatten() - (B * A)
        Minv = 1 / (B ** 2 + mutation_num * alpha)
        trace_g = -alpha * np.sum(Minv)

        Minv_outer = np.outer(Minv, Minv)
        omiga_tmp_1 = np.diag(Minv) + (1 / (1 + trace_g) * alpha * Minv_outer)# add?
        del A,B,trace_g
        gc.collect()
        omiga_new = np.matmul(omiga_tmp_1, omiga_tmp_2.T)
        distance = np.subtract.outer(omiga_new, omiga_new)
        delt = (distance[index_dis] - 1 / alpha * tau.T).ravel()  #
        eta_new = update(delt, alpha, Lambda, tau)
        del delt,omiga_tmp_1,omiga_tmp_2
        gc.collect()
        omiga_new[omiga_new > control_large] = control_large
        omiga_new[omiga_new < -control_large] = -control_large
        tau = tau - np.array([alpha * (distance[index_dis] - eta_new)]).T
        residual = np.max(distance[index_dis] - eta_new)
        v = res_last-residual
        res_last = residual
        next_step = 1 + rho*(v/np.abs(v))
        alpha = alpha * next_step
        print(residual)
    eta_new = np.where(np.abs(eta_new) < post_th,0,eta_new)
    distance[index_dis] = eta_new
    return distance,phi_new,mutation_num

def corrected(distance,phi_new,n,mutation_num):
    class_label = -np.ones(mutation_num)
    class_label[0] = 0
    group_size = [1]
    labl = 1
    for i in range(0, mutation_num):
        for j in range(i+1):
            if distance[j, i] == 0:
                class_label[i] = class_label[j]
                group_size[int(class_label[j])] += 1
                break
        if class_label[i] == -1:
            class_label[i] = labl
            labl += 1
            group_size.append(1)
    labels = np.unique(class_label)
    phi_out = np.zeros([len(labels)])-1
    index = 0
    for i in range(len(labels)):
        ind = np.where(class_label == labels[i])[0]
        phi_tmp = np.sum(phi_new[ind] * n[ind]) / np.sum(n[ind])
        indexs= index
        if phi_tmp in phi_out:
            indexs = np.argwhere(phi_out==phi_tmp)
        else:
            index+=1
        class_label[ind] = indexs
        phi_out[indexs]=phi_tmp
    np.delete(phi_out, -1, axis=0)
    phi_res = np.zeros(mutation_num)
    for lab in range(len(phi_out)):
        phi_res[class_label == lab] = phi_out[lab]
    phi_res = np.round(phi_res,3)
    return phi_res,class_label

def clust_GMM(phi_new,phi_res,mutation_num, has_sklearn=True):
    limit_diff = 0.02
    limit_clust = 0.03
    Xmoon = phi_res.reshape([-1,1])
    if has_sklearn == True:
        n_components = np.arange(1, 5)
        models = [GMM(n, random_state=0).fit(Xmoon)
                  for n in n_components]
        bic_list = []
        for m in models:
            bic_list.append(m.aic(Xmoon))
        num_clust = np.argwhere(bic_list==np.min(bic_list))
        num_clust = np.max(num_clust)+1
        print(num_clust)
        gmm = GMM(n_components=num_clust,random_state=0).fit(phi_res.reshape([-1,1]))
        result = gmm.predict(phi_res.reshape([-1,1]))
    else:
        model1 = GMM.GMM(phi_res,5)
        result = model1.fit()

    result = np.array(result)
    lab = np.unique(result)
    balance = []
    for i in lab:
        c_list = phi_new[result == i]
        mu = fq(c_list)#round(np.mean(c_list), 2)
        sigma = np.std(c_list)
        print(mu,end="\t")
        print(sigma, end="\t")
        print(np.size(c_list))
        balance.append(mu)
    index_sort = np.argsort(balance)
    print("======================================================")
    if len(index_sort)==1:
        return {'phi': phi_res, 'balance': balance, 'label': result}
    for i in range(0, len(index_sort)):
        previous = i - 1
        the = i
        next = i + 1
        if i == 0:
            previous = the
            diff_after = balance[index_sort[next]] - balance[index_sort[the]]
            diff_befor = 1
        elif i == (len(index_sort)-1):
            next = the
            diff_befor = balance[index_sort[the]] - balance[index_sort[previous]]
            diff_after = 1
        else:
            diff_after = balance[index_sort[next]] - balance[index_sort[the]]
            diff_befor = balance[index_sort[the]] - balance[index_sort[previous]]
        if diff_befor == 0 :
            meger = 1
        elif diff_after == 0:
            meger = -1
        elif diff_befor > diff_after:
            meger = 1
        else:
            meger = -1
        c_1 = lab[index_sort[previous]]
        c1 = lab[index_sort[the]]
        c2 = lab[index_sort[next]]
        sigma = np.std(phi_new[result == c1])
        num = np.size(phi_new[result == c1])
        if meger == -1 and (diff_befor < limit_diff or num <= int(limit_clust * mutation_num) or sigma<0.002):
            print("merge:"+str(c1)+'\t'+str(c_1))
            clust1 = phi_new[result == c1]
            clust2 = phi_new[result == c_1]
            mu = np.mean(np.r_[clust1, clust2])
            balance[index_sort[the]] = mu
            balance[index_sort[previous]] = mu
            result = np.where(result == c_1, c1, result)
            lab[index_sort[previous]] = c1
        elif meger == 1 and (diff_after < limit_diff or num <= int(limit_clust * mutation_num) or sigma<0.002):
            print("merge:" + str(c1) + '\t' + str(c2))
            clust1 = phi_new[result == c1]
            clust2 = phi_new[result == c2]
            mu = np.mean(np.r_[clust1, clust2])
            balance[index_sort[the]] = mu
            balance[index_sort[next]] = mu
            result = np.where(result == c2, c1, result)
            lab[index_sort[next]] = c1
    return {'phi': phi_res,'balance':balance, 'label': result}

def clust_GMM_multisample(phi_new,phi_res, has_sklearn=True):
    Xmoon = phi_res.reshape([-1,2])
    if has_sklearn == True:
        n_components = np.arange(1, 6)
        models = [GMM(n, random_state=0).fit(Xmoon)
                  for n in n_components]
        bic_list = []
        for m in models:
            bic_list.append(m.aic(Xmoon))
        num_clust = np.argwhere(bic_list==np.min(bic_list))
        num_clust = np.max(num_clust)+1
        print(num_clust)
        gmm = GMM(n_components=4,random_state=0).fit(phi_res.reshape([-1,2]))
        result = gmm.predict(phi_res.reshape([-1,2]))
    else:
        model1 = GMM.GMM(phi_res,5)
        result = model1.fit()

    result = np.array(result)
    lab = np.unique(result)
    balance = []
    for i in lab:
        c_list = phi_res[result == i,:]
        mu = np.mean(c_list,axis=0)
        sigma = np.std(c_list,axis=0)
        print(mu,end="\t")
        print(sigma, end="\t")
        print(np.size(c_list))
        balance.append(mu)
    print("======================================================")
    return {'phi': phi_res,'balance':balance, 'label': result}

def Elastic_GMM_subclone(SNV_pass,purity,Prefix):
    r = SNV_pass[:,2].astype(int)
    n = SNV_pass[:,3].astype(int)
    total = SNV_pass[:,4].astype(int)
    multiplicity = SNV_pass[:,5].astype(int)
    #--------------------------------epoch canshu--------------------------------#
    Lambda = float(1)
    ploidy = 2
    alpha = 0.2
    rho = 0.02
    precision = 0.01
    Run_limit = 500
    control_large = 5
    #--------------------------------jingdu canshu-------------------------------#
    post_th = 0.01
    # -------------------------------run   --------------------------------------#
    distance, phi_new, mutation_num = ElasticNet_correct(r, n, multiplicity, total, ploidy, Lambda, alpha,
                                        rho, Run_limit, precision, control_large,post_th,purity,Prefix)
    return distance, phi_new, mutation_num,n

def clust_stdout(SNVtable,phi_new,res,path,name):
    labl = np.unique(res['label'])
    summary = np.zeros([len(labl), 3])
    for i in range(len(labl)):
        summary[i, 0] = labl[i]
        tmp = np.mean(res['phi'][np.where(res['label'] == labl[i])])
        #tmp = np.mean(res['balance'][labl[i]])
        print(tmp)
        summary[i, 2] = np.round(tmp, 3)
        summary[i, 1] = len(np.where(res['label'] == labl[i])[0])
    index = np.argsort(summary[:,2]).reshape(-1)
    summary = summary[index,:]
    answer = SNVtable[:,0:4].astype(int)
    #answer2 = SNVtable[:, 4:6].astype(int)
    answer = np.c_[answer,res['label'].astype(int),res['phi'].astype(float)]
    #with open(path+'/'+str(name)+'.Clust_pos.txt','w') as f:
    #    f.write('#chrom\tposition\tmajor\tminor\tclass\tvaf_phi \n')
    np.savetxt(path+'/'+str(name)+'.vaf.txt', phi_new, fmt='%.3f', delimiter=',')
    np.savetxt(path+'/'+str(name)+'.Clust_pos.txt', answer, fmt='%d\t%d\t%d\t%d\t%d\t%.3f', delimiter=',')
    np.savetxt('%s/%s.summary_table.txt' % (path, str(name)), summary, fmt='%d\t%d\t%.3f')

def clust_stdout_muliti(SNVtable,res,path,name):
    labl = np.unique(res['label'])
    summary = np.zeros([len(labl), 4])
    for i in range(len(labl)):
        summary[i, 0] = labl[i]
        #tmp = np.mean(res['phi'][np.where(res['label'] == labl[i]),:],axis=0)
        tmp = res['balance'][labl[i]]
        print(tmp)
        summary[i, 2:4] = np.round(tmp, 3)
        summary[i, 1] = len(np.where(res['label'] == labl[i])[0])
    index = np.argsort(summary[:,2]).reshape(-1)
    summary = summary[index,:]
    answer = SNVtable[:,0:4].astype(int)
    #answer2 = SNVtable[:, 4:6].astype(int)
    answer = np.c_[answer,res['label'].astype(int),res['phi'][:,0],res['phi'][:,1]]
    #with open(path+'/'+str(name)+'.Clust_pos.txt','w') as f:
    #    f.write('#chrom\tposition\tmajor\tminor\tclass\tvaf_phi \n')
    np.savetxt(path+'/'+str(name)+'.Clust_pos.txt', answer, fmt='%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f', delimiter=',')
    np.savetxt('%s/%s.summary_table.txt' % (path, str(name)), summary, fmt='%d\t%d\t%.3f\t%.3f')