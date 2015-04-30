#coding=utf-8 
'''
This is for computing the second hop rate tuple with a MAC channel in CoF
'''
from sage.all import *
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
import itertools

if is_set_H == True:
    H_b = set_H_b
else:
    set_random_seed() # to avoid producing the same H_a in different threads
    H_b = (matrix.random(RR, L,M , distribution=RealDistribution('gaussian', 1)))

#用于计算第二跳多天线MAC 信道容量rate region
#第二跳为终端多天线的MAC信道，N*N
#输入信道矩阵H_b，relay转发功率（假设都为P_relay）
#输出第二跳信道achievable region constraint tuple
def ComputeSecRate(M,P_relay,H_b):
    rate_sec_hop=[0]*M
    P=[0]*M
    constraint=[]
    for i in range(M):
        P[i]=(H_b.column(i).norm()**2)*P_relay
        rate_sec_hop[i]=0.5*log(1+P[i],2)
    constraint.extend(rate_sec_hop)
    if M==2:
        MaxConstraint=0.5*log(1+P[0]+P[1],2)
        constraint.append(MaxConstraint)
        return constraint
    if M==3:
        rate_bound=[]
        for i in range(M):#实质为一个组合
            for j in range(i+1,M):
                rate_bound.append(0.5*log(1+P[i]+P[j],2))
        constraint.extend(rate_bound)
        MaxConstraint=0.5*log(1+P[0]+P[1]+P[2],2)
        constraint.append(MaxConstraint)
        return constraint
    

