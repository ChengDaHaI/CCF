'''
Simulation of Relay Compressing in CCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
'''
from sage.all import *
from NEW_CCF_Modle import Relay_Forward_Rate, powerset
from NewSecondHopChannel import ComputeSecRate
from scipy import optimize
from CoF_LLL import Find_A_and_Rate
from NEW_basic import *
import math
import time
import itertools
import copy

#the varables to optimize are the all rate piece beta1~2L-1
#per_s is [1,,,L]'s permutation, per_c is [L+1,2*L]'s permutation
def Linear_Program(entropy_coefficient_list,secChannel_constiant,source_rate_upbound_list,per_s,per_c):
    #object Function in the linear programming problem 
    C=[0]*(2*L-1)
    for i in range(0,2*L-1):
        if i <=L-1:
            C[i]=-(i+1)#Attention, the real object function coefficient should be positive
        elif i>=L:
            C[i]=-(2*L-1-i)#Attention, the real object function coefficient should be positive
    # source rate is the coefficient list of rate pieces
    SourseRate=[]
    for i in range(L):
        SourseRate.extend([[0]*(2*L-1)])
    for i in range(0,L):
        #piece_mount=per_c[i]-per_s[i]
        for j in range(0,2*L-1):
            #SourseRate[i][j]=[0]*(2*L-1)
            if (j>=per_s[i])&(j<=per_c[i]+L-1):
                SourseRate[i][j]=1
    #construct the linear programming equation
    channel_mode="parallel"
    if channel_mode=="parallel":
        A_ConstriantMatrix=SourseRate+entropy_coefficient_list
        b_ConstriantVector=source_rate_upbound_list+secChannel_constiant[0:M]
        #change the parallel channel capacity constraints
        subset_list=list(powerset(range(0,L)))
        for i in range(L+1,len(subset_list)):
            bound_sum=0
            for j in subset_list[i]:
                bound_sum=bound_sum+secChannel_constiant[i-1]
            b_ConstriantVector.append(bound_sum)
    elif channel_mode=="MAC":
        b_ConstriantVector=source_rate_upbound_list+secChannel_constiant
        A_ConstriantMatrix=SourseRate+entropy_coefficient_list
    # the default bound is nonnegative, that is (0,None)
    bound=[0]*(2*L-1)
    for i in range(2*L-1):
        bound[i]=(0, None)
    bound=tuple(bound)
    result=optimize.linprog(C, A_ub=A_ConstriantMatrix, b_ub=b_ConstriantVector, bounds=bound, options={"disp": False})
    print result.x
    return result

#compute the source rate upbound and matrix A when given variable beta
# the output would be the input of linear programme function
def CCF_fix_pow_sourceRate_upbound(P_con,H_a,beta=[]):
    (L, L) = (H_a.nrows(), H_a.ncols())#Assuming the H_a matrix is L by L
    if beta == []:
        beta = vector(RR, [1]*L)
    for be in list(beta):
        if be <= 0:
            return 0
    B = diagonal_matrix(beta)
    P_t=P_con
    try:
        P_t[0]
    except:
        P_t = [P_t]
    for i_P in range(0, L):
        if math.isnan(P_t[i_P]):
            print 'P', str(i_P), ' should not be NaN!'
            raise Exception('Invalid power setting reached.')
        '''
        if P_t[i_P] <= 0 or P_t[i_P] > (P_con+0.1):
            # print 'P', str(i_P), ' should be positive and less than P_con'
            return 0
        '''
    P_vec = vector(RR, P_t)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, source_rate_list, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, True, beta)
    except:
        print 'error in seeking A and source rate upbound list'
        raise
    A_best_LLL_F = matrix(GF(p), A_best_LLL)
    if A_best_LLL_F.rank() != min(L, M):
        source_rate_list=0
    return (source_rate_list, A_best_LLL_F)

#compute the new CCF system sum rate
def CCF_new_sumrate_func(betaScale, H_a, H_b, P_con, P_relay,per_c):
    #compute the source rate upbound and matrix A
    source_rate_upbound_list, A  =CCF_fix_pow_sourceRate_upbound([P_con]*L,H_a,betaScale)
    #the second hop channel capacity constriant
    SecChannel_constiant=ComputeSecRate(L,P_relay,H_b)
    #test program running time cost
    Max_sumrate=0
    t1=time.time()
    #=======#
    #compute the proper nested shaping lattice order from the betaScale
    beta=copy.copy(betaScale)
    beta=list(beta)
    per_s=[]
    max_beta=0
    for i in range(L):
        max_beta=max(beta)
        max_beta_index=beta.index(max_beta)
        per_s.append(max_beta_index)
        beta[max_beta_index]=0
    # a larger beta corresponds to a coarse shping lattice
    '''
    per_s.reverse()
    '''
    #compute the coefficient of rate pieces of the coditional entropy 
    entropy_coefficient_list=Relay_Forward_Rate(per_s,per_c,A)
    Res=Linear_Program(entropy_coefficient_list,SecChannel_constiant,source_rate_upbound_list,per_s,per_c)
    #========#
    t2=time.time()
    t=t2-t1
    return Res.fun
    
def RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay,per_c=[]):
    '''
    CCF_beta_func=lambda x: CCF_sumrate_compute(vector(RR, [1,]+list(x[0:L-1])), H_a, H_b, P_con, P_relay, per_s, per_c)
    '''
    if per_c==[]:
        #compute the proper coding lattice order
        H_a_col_min=[]
        H_a_trans=H_a.transpose()
        H_a_trans=list(H_a_trans)
        for i in range(L):
            for j in range(L):
                temp=math.fabs(H_a_trans[i][j])
                H_a_trans[i][j]=copy.copy(temp)
            H_a_col_min.append(min(H_a_trans[i]))
        per_c=[]
        for i in range(L):
            H_a_colmin_max=max(H_a_col_min)
            H_a_colmin_max_index=H_a_col_min.index(H_a_colmin_max)
            per_c.append(H_a_colmin_max_index)
            H_a_col_min[H_a_colmin_max_index]=0
        per_c.reverse()#a larger channel coefficient corresponds to a finer coding lattice
        
    #perform differential evolution before computing two permutation 
    CCF_beta_func=lambda x: CCF_new_sumrate_func(vector(RR, [1,]+list(x[0:L-1])), H_a, H_b, P_con, P_relay,per_c)
    Pranges=((0.1,betaScale_max),)*(L-1)
    if P_Search_Alg=='differential_evolution':
        #test program running time cost
        t1=time.time()
        ResSearch=optimize.differential_evolution(CCF_beta_func,Pranges,strategy="best1bin",maxiter=10)
        t2=time.time()
        t=t2-t1
        beta_opt=ResSearch.x
        sum_rate_opt=-ResSearch.fun
    else:
        Exception("error: Not Such Search Algorithm!")
    return beta_opt, sum_rate_opt
    

if __name__=="__main__":
    '''
    L=3
    A=Matrix([[1,2,3],[1,2,2],[2,1,3]])#matrix A must be full-rank
    beta_s=[25,20,18]
    beta_c=[15,10,8]
    per_s=[0,2,1]
    per_c=[0,1,2]
    P_con=1000
    P_relay=0.25*P_con
    '''
    beta_opt, sum_rate_opt=RandomSearch('differential_evolution', H_a, H_b, P_con, P_relay, per_s, per_c)
    print "optimal beta list:", beta_opt
    print "optimal sum rate:", sum_rate_opt
    '''
    #just for test the source rate up bound
    source_rate_upbound_list=[5,4.5,3.5]
    #produce the conditional entropy and the coefficient of rate piece
    conditional_entropy_list,entropy_coefficient_list=Relay_Forward_Rate(beta_s,beta_c,per_s,per_c,A)
    #the second hop channel capacity constriant
    SecChannel_constiant=ComputeSecRate(L,P_relay,H_b)
    Res=Linear_Program(entropy_coefficient_list, SecChannel_constiant, source_rate_upbound_list, per_s, per_c)
    print Res
    '''
    