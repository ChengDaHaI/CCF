#coding=utf-8 
from sage.all import *
import math
import numpy

#信道维度
L=2
M=2
#两跳信道矩阵
is_set_H=True
if is_set_H == True:
    #Ha = matrix(RR, M, L, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
    Ha=matrix(RR,M,L,[[-0.541978155712295 ,0.740073351426688],[-0.773073785628476 ,0.584325217080305]])
    #Ha=matrix(RR,M,L,[[ -0.129901655668989 ,0.803654841233371],[  0.180173049695095,-0.0675284022380773]])
    #Hb = matrix(RR, M, L, [[ 0.806026835557602,-0.267139360752616], [0.455755216914796, 0.590419325969173]])
    Hb = matrix(RR, 1, 2, [0.806026835557602,-0.267139360752616])
else:
    Ha = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    Hb = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))

#计算AmplifyForward系统速率
#即等效MIMO信道速率
#默认L=M,若不等，会出错
'''
def Compute_AF_rate(Ps,Pr,Ha,Hb,L=2,M=2):

    alpha=[0]*M
    for i in range(M):
        alpha[i]=math.sqrt(Pr/(Ha.row(i).norm()**2*Ps+1))
    #calculate diagonal matrix alpha
    A=matrix.diagonal([alpha[i] for i in range(M)])
    #calculate equivalent matrix H
    H=Hb*A*Ha
    #source power diagonal matrix
    Q=matrix.diagonal([Ps for i in range(L)])
    #the equivalent noise matrix
    N=(Hb*A)*(Hb*A).transpose()
    I=matrix.diagonal([1]*M)
    N=N+I
    #calculate MIMO channel sum rate
    sum_rate=max(0,0.5*log((I+H*Q*H.transpose()*N.inverse()).determinant(),2))
    return sum_rate
'''

#终端为单天线
#故Hb为一横向量
#此处L表示中继数量，M表示终端天线数量
def Compute_AF_rate(Ps,Pr,Ha,Hb,L=2,M=2):

    alpha=[0]*L
    for i in range(L):
        alpha[i]=math.sqrt(Pr/(Ha.row(i).norm()**2*Ps+1))
    #calculate diagonal matrix alpha
    A=matrix.diagonal([alpha[i] for i in range(L)])
    #calculate equivalent matrix H
    H=Hb*A*Ha
    #source power diagonal matrix
    Q=matrix.diagonal([Ps for i in range(L)])
    #the equivalent noise matrix
    N=(Hb*A)*(Hb*A).transpose()
    I=matrix.diagonal([1]*M)
    N=N+I
    #calculate MIMO channel sum rate
    sum_rate=max(0,0.5*log((1+H*Q*H.transpose()*N.inverse()).determinant(),2))
    return sum_rate

if __name__=="__main__" :
    Ps=1023
    Pr=0.25*Ps
    print Compute_AF_rate(Ps, Pr, Ha, Hb,L=2,M=1)