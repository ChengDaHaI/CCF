#coding=utf-8 
'''
This includes some functions related with the second hop in CoF.
That seems to be wrong!
'''
from sage.all import *
from scipy import optimize
from numpy import arange
from sage.parallel.all import *
import time
import copy
from CoF_basic import *
from itertools import chain, combinations
from scipy.optimize import linprog

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


#mod_scheme=asym_mod,quan_scheme=asym_quan
#assume all matrix invertible
#return the two corner point in slepian-wolf coding
def Relay_Forward_Rates(relay_fine_lattices,trans_coarse_lattices,A):
    (M, L) = (A.nrows(), A.ncols())
    if M != L:
        raise Exception("L and M should be the same in destination's perspective.")
    
    relay_compute_fine_lattices = list(relay_fine_lattices)
    # determine the fine lattice of the l-th transmitter according to computation constraints
    trans_compute_fine_lattices = [float(0)]*L
    for i_L in range(0, L):
        for i_M in range(0, M):
            if (A[i_M, i_L]!=0) and (relay_compute_fine_lattices[i_M]>trans_compute_fine_lattices[i_L]):
                trans_compute_fine_lattices[i_L] = relay_compute_fine_lattices[i_M]
    #sort the fine lattices in the transmitters
    trans_compute_fine_lattices.sort()
    sorted_trans_fine_lattices=trans_compute_fine_lattices
    #determine the compressing fine lattices at the relays
    relay_compress_fine_lattices=sorted_trans_fine_lattices
    
    #sort the coarse lattices in the transmitters
    trans_coarse_lattices.sort()
    sorted_trans_coarse_lattices=trans_coarse_lattices
    #determine the compressing coarse lattices at the relays
    sorted_trans_coarse_lattices.reverse()
    relay_compress_coarse_lattices =sorted_trans_coarse_lattices
    
    #calculate the forwarding rates
    r_f=[0]*M
    for i_M in range(0, M):
        r_f[i_M]=max(0,0.5*log(relay_compress_coarse_lattices[i_M]/relay_compress_fine_lattices[i_M],2))
    #return forwarding rate tuple/list
    return r_f


#input the forwarding rates list and the second hop channel capacity constraint list
#using linear programming to calculate the max rates the second hop can support
def SecHop_Support_Rates(ForwardRate,SecHopConstraint):
    #the sum of forwarding rates exceed the second hop capacity
    if sum(ForwardRate)>=SecHopConstraint[-1]:
        Support_rate=SecHopConstraint[-1]
        print 'Forwarding rates exceed the second hop capacity\n'
        return Support_rate
    else:#two rate region intersect
        #the left inequality coefficient in second hop capacity
        A1_un=[0]*(pow(2,M)-1-M)
        #the left inequality coefficient in forwarding rates
        A2_un=[0]*(pow(2,M)-1-M)
        list_M=range(1,M+1)
        subsets_list=list(powerset(set(list_M)))
        for i in range(M+1,pow(2,M)):
            a1_un=[0]*M
            a2_un=[0]*M
            for j in subsets_list[i]:
               a1_un[j-1]=1
               a2_un[j-1]=-1
            A1_un[i-M-1]=a1_un
            A2_un[i-M-1]=a2_un
        
        A1_un.extend(A2_un)
        A_un=A1_un
        #the right inequality coefficient in second hop capacity
        b1_un=SecHopConstraint[M:]
        #the left inequality coefficient in forwarding rates
        b2_un=[0]*(pow(2,M)-1-M)
        ForRate=copy(ForwardRate)
        ForRate.reverse()
        for i in range(M+1,pow(2,M)):
            j=len(subsets_list[i])
            b2_un[i-M-1]=sum(ForRate[0:j])
        b1_un.extend(b2_un)
        b_un=b1_un
        #the rate bounds
        bound=[0]*M
        for i in range(M):
            bound[i]=(ForwardRate[i],SecHopConstraint[i])
        bound=tuple(bound)
        
        #what's the object function???
        C=[1]*M
        result=linprog(C,A_ub=A_un,b_ub=b_un,bounds=bound,options={"disp":True})
        return result

if __name__=="__main__":
    print '-----------------------------------\n'+ \
    'testing CoF_second_hop\n'
    M=3
    P_relay=15
    set_random_seed(1)
    H_b= matrix.random(RR, M,M, distribution=RealDistribution('gaussian', 1))
    constraint_list=ComputeSecRate(M, P_relay, H_b)
    print 'calculated the constraint_list \n'
    #R = [1, 2]
    A = matrix(ZZ, 3, 3, [[1, 2,1], [2,1, 1],[2,1,3]])
    p = 3
    relay_fine_lattices = [0.5, 0.2,0.4]
    trans_coarse_lattices = [1.5, 2.0,1.8]
    beta = [1]*M
    Forward_rate=Relay_Forward_Rates(relay_fine_lattices, trans_coarse_lattices, A)
    print 'calculated the Forward rates \n'
    Res=SecHop_Support_Rates(Forward_rate,constraint_list)
    print Res
    print 'test ended! \n'
    