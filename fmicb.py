import pandas as pd
import numpy as np
import math

# 设置常量
r_m=1
alpha=0.2
cutoff=1e-6

df = pd.read_csv("HMDAD.csv")                   # 读取csv
df.drop(labels=['Position','Evidence','PMID'],axis=1,inplace=True)      # 删除后三列
norepeat_df=df.drop_duplicates(ignore_index=True)            # 删除重复数据
print(df.head(5))                   # 查看前5行

disease=pd.unique(norepeat_df['Disease'])    # 获取所有疾病
nd=len(disease)
microbe=pd.unique(norepeat_df['Microbe'])    # 获取所有微生物
nm=len(microbe)
data=np.array(norepeat_df)                   # 转成np矩阵便于矩阵运算

#===========================================================计算A和H============================================================
A=np.zeros((nd,nm))         # 初始化adjacency matrix
for i in range(nd):
    disease_index=list(norepeat_df[norepeat_df.Disease==disease[i]].index)        # 获取疾病出现的位置
    for j in disease_index:
        microbe_index=int(np.where(microbe==data[j][1])[0])
        A[i][microbe_index]=1
H=A.transpose()             # incidence matrix是adjacency matrix的转置

#===========================================================计算GM=================================================================
GM=np.zeros((nm,nm))        # 初始化GM矩阵
for i in range(nm):
    for j in range(nm):
        GM[i][j]=math.exp(-r_m*(np.sum((A.transpose()[i]-A.transpose()[j])**2)))           # 根据公式1计算相似度矩阵（是对称矩阵）

#===========================================================计算D_v=================================================================
D_v=np.zeros((nm,nm))       # 初始化顶点的degree matrix
for i in range(nm):
    D_v[i][i]=(1/nd)*np.sum(H[i])          # 是一个对角阵，对角线元素是顶点的度，即这个微生物和多少种疾病关联

#===========================================================计算W=================================================================
W=np.zeros((nm,nd))       # 初始化顶点的权重矩阵
for i in range(nm):
    for j in range(nd):
        if H[i][j]==1:      # 如果第i种微生物和第j种疾病相关
            W[i][j]=np.sum(GM[i])

#===========================================================计算D_ve=================================================================
D_ve=np.zeros((nd,nd))      # 初始化超边的degree matrix（根据W得到）
for i in range(nd):
    D_ve[i][i]=np.sum(W.transpose()[i])       # 是一个对角阵，对角线元素是超边的度。因为和该病无关的微生物权重为0，所以可以转置后对行求和

#===========================================================计算W_e=================================================================
W_e=np.zeros((nd,nd))       # 初始化边的权重矩阵
for i in range(nd):
    W_e[i][i]=1/nd

#=====================================================计算P==========================================================================
P=(np.linalg.inv(D_v)).dot(H).dot(W_e).dot(np.linalg.inv(D_ve)).dot(W.transpose())

#=====================================================计算v，排序=======================================================================
pre_dieases=['Asthma','Type 2 diabetes',"Crohn's disease(CD)"]
for pre_diease in pre_dieases:
    v=np.zeros((nm,1))
    pre_diease_index=int(np.where(disease==pre_diease)[0])          # 获取该疾病的编号
    disease_index=list(np.where(A[pre_diease_index]==1)[0])      # 获取与该疾病有关的微生物列表
    for j in disease_index:
        v[j]=1/len(disease_index)       # 设置初始值
    # 更新
    v_old=v
    v_new=(1-alpha)*P.transpose().dot(v_old)+alpha*v
    while(np.linalg.norm(v_old-v_new,ord=1)>=cutoff):
        v_old=v_new
        v_new = (1 - alpha) * P.transpose().dot(v_old) + alpha * v
        print(np.linalg.norm(v-v_new,ord=1))
    rank=np.argsort(-v_new,axis=0)
    count=0
    for k in range(len(rank)):
        if v[rank[k]]!=0:               # 不等于0就是seed microbe
            continue
        else:                           # 等于0才是candidate microbe
            count += 1
            print("[{}] disease:{},  microbe:{}".format(count,pre_diease,microbe[rank[k]]))
        if count==10:
            break