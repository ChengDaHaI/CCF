#coding=utf-8 
'''
This is for computing the second hop rate tuple with a MAC channel in CoF
'''
from sage.all import *
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
from itertools import chain, combinations

#produce the subsets
def powerset(iterable):
  xs = list(iterable)
  # note we return an iterator rather than a list
  return chain.from_iterable( combinations(xs,n) for n in range(len(xs)+1) )

#用于计算第二跳多天线MAC 信道容量rate region
#第二跳为终端多天线的MAC信道，N*N
#输入信道矩阵H_b，relay转发功率（假设都为P_relay）
#输出第二跳信道achievable region constraint list
def ComputeSecRate(M,P_relay,H_b):
    rate_sec_hop=[0]*M
    P=[0]*M
     #calculate the second hop channel capacity
    constraint=[]
    for i in range(M):
        P[i]=(H_b.column(i).norm()**2)*P_relay
        rate_sec_hop[i]=0.5*log(1+P[i],2)
    constraint.extend(rate_sec_hop)
    list_M=range(1,M+1)
    subsets_list=list(powerset(set(list_M)))
    for i in range(M+1,pow(2,M)):
        pow_forward=0
        for j in subsets_list[i]:
            pow_forward+=P[j-1]
        constraint.append(0.5*log(1+pow_forward,2))
    return constraint
    
if  __name__=="__main__":
    M=3
    P_relay=25
    #set_random_seed(1)
    H_b= matrix.random(RR, M,M, distribution=RealDistribution('gaussian', 1))
    constraint_list=ComputeSecRate(M, P_relay, H_b)
    print constraint_list

    
