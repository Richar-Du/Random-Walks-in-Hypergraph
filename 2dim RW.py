# 二维随机游走
import matplotlib.pyplot as plt
import numpy as np
# ========================================================生成游走矩阵======================================
a=[0,0]
A=[a,]
np.random.seed(1)
while len(A)<10000:
    p=np.array([[0,1],[1,0],[0,-1],[-1,0]])/10          # 初始化四个方向
    b=np.random.randint(len(p))                         # 随机选择一个
    a=np.array(a)+p[b]                                  # 获取位置
    c=list(a)
    A.append(c)                                         # 追加到位置矩阵中
# =================================================作图====================================================
B=np.array(A).T
plt.title('One dimension Random Walk',size=20)
plt.scatter(0,0,c='r')
plt.axis('off')
plt.plot(B[0],B[1])
plt.savefig('2dim RW.png')
plt.show()