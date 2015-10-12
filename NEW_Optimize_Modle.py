'''
Simulation of Relay Compressing in CCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
'''
from sage.all import *
from NEW_CCF_Modle import Relay_Forward_Rate
from NewSecondHopChannel import ComputeSecRate
from scipy import optimize
from CoF_LLL import Find_A_and_Rate
import math
from NEW_basic import *

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
    SourseRate=[[0]*(2*L-1)]*L
    for i in range(0,L):
        #piece_mount=per_c[i]-per_s[i]
        for j in range(0,2*L-1):
            #SourseRate[i][j]=[0]*(2*L-1)
            if j>=per_s[i]-1&j<=per_c[i]-1:
                SourseRate[i][j]=1
    #construct the linear programming equation
    A_ConstriantMatrix=SourseRate+entropy_coefficient_list
    b_ConstriantVector=source_rate_upbound_list+secChannel_constiant
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

def CCF_sumrate_compute(betaScale, H_a, H_b, P_con, P_relay, per_s, per_c):
    #compute the source rate upbound and matrix A
    source_rate_upbound_list, A                                 =CCF_fix_pow_sourceRate_upbound([P_con]*L,H_a,betaScale)
    #compute the coditional entropy and the coefficient of rate pieces
    conditional_entropy_list,entropy_coefficient_list=Relay_Forward_Rate(beta_s,beta_c,per_s,per_c,A)
    #the second hop channel capacity constriant
    SecChannel_constiant=ComputeSecRate(L,P_relay,H_b)
    Res=Linear_Program(entropy_coefficient_list,SecChannel_constiant,source_rate_upbound_list,per_s,per_c)
    #x=Res.x
    #fun=Res.fun
    #return the optimizer and the object value(the real object value should be -Res.fun)
    return Res.fun
    
def RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c):
    CCF_beta_func=lambda x: CCF_sumrate_compute(vector(RR, [1,]+list(x[0:L-1])), H_a, H_b, P_con, P_relay, per_s, per_c)
    Pranges=((0.1,betaScale_max),)*(L-1)
    if P_Search_Alg=='differential_evolution':
        ResSearch=optimize.differential_evolution(CCF_beta_func,Pranges)
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
    