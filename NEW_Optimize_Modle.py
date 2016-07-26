'''
Simulation of Relay Compressing in CCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
'''
from sage.all import *
from NEW_CCF_Modle import Relay_Forward_Rate, Powerset
#from NewSecondHopChannel import ComputeSecRate
from scipy import optimize
from CoF_LLL import Find_A_and_Rate
from NEW_general_optimize_model import GCCF_new_sumrate_func
from NEW_basic import *
import math
import time
import itertools
import copy
import numpy as np

#the varables to optimize are the all rate piece beta1~2L-1
#per_s is [1,,,L]'s permutation, per_c is [L+1,2*L]'s permutation

#scaling factor beta decide the rate pieces in shaping lattice part
def Linear_Program(entropy_coefficient_list,secChannel_constriant,source_rate_upbound_list,per_s,per_c,beta):
    
    #object Function in the linear programming problem 
    C=[0]*(L)#the last L rate pieces
    for i in range(0,L):
        C[i]=-(L-i)#Attention, the real object function coefficient should be positive
    #the first (L-1) rate pieces
    Part_ratepiece=[0]*(L-1)
    for i in range(L-1):
        Part_ratepiece[i]=0.5*np.log2(beta[per_s[i]]/beta[per_s[i+1]])
    # source rate is the coefficient list of rate pieces
    SourseRate=[]
    for i in range(L):
        SourseRate.extend([[0]*(2*L-1)])
    for i in range(0,L):
        #piece_mount=per_c[i]-per_s[i]
        for j in range(0,2*L-1):
            #SourseRate[i][j]=[0]*(2*L-1)
            if (j>=per_s.index(i))&(j<=per_c.index(i)+L-1):
                SourseRate[i][j]=1
    #construct the linear programming equation
    channel_mode="parallel"
    if channel_mode=="parallel":
        A_ConstriantMatrix=[SourseRate[i][L-1:2*L-1] for i in range(len(SourseRate))]+[entropy_coefficient_list[i][L-1:2*L-1] for i in range(len(entropy_coefficient_list))]
        #SourseRate[:][L-1:2*L-1]+entropy_coefficient_list[:][L-1:2*L-1]
        b_ConstriantVector=source_rate_upbound_list+secChannel_constriant[0:M]
        #change the parallel channel capacity constraints
        subset_list=list(Powerset(set(range(0,L))))
        for i in range(L+1,len(subset_list)):
            bound_sum=0
            for j in subset_list[i]:
                bound_sum=bound_sum+secChannel_constriant[j]
            b_ConstriantVector.append(bound_sum)
        #substract the known part form the b_ConstriantVector
        Part_Aconstriant=np.array([SourseRate[i][0:L-1] for i in range(len(SourseRate))]+[entropy_coefficient_list[i][0:L-1] for i in range(len(entropy_coefficient_list))])
        Part_ratepiece=np.array(Part_ratepiece)
        To_Be_Sub=np.dot(Part_Aconstriant,Part_ratepiece)
        b_ConstriantVector=np.subtract(np.array(b_ConstriantVector),To_Be_Sub)
        b_ConstriantVector=b_ConstriantVector.tolist()
    elif channel_mode=="MAC":
        b_ConstriantVector=source_rate_upbound_list+secChannel_constriant
        A_ConstriantMatrix=SourseRate+entropy_coefficient_list
    # the default bound is nonnegative, that is (0,None)
    bound=[0]*(L)
    for i in range(L):
        bound[i]=(0, None)
    bound=tuple(bound)
    result=optimize.linprog(C, A_ub=A_ConstriantMatrix, b_ub=b_ConstriantVector, bounds=bound, options={"disp": False})
    Part_SourceRate=np.dot(np.array([-(i+1) for i in range(0,L-1)]),Part_ratepiece)
    if result.success == False:
        # print 'optimization failure'
        return 0
    else:
        #print 'source rate pieces:', Part_ratepiece, result.x
        return result.fun+Part_SourceRate
    #return the true max summation of source rates
        

#compute the source rate upbound and matrix A when given variable beta
# the output would be the input of linear programme function
def CCF_fix_pow_sourceRate_upbound(P_source,H_a,beta=[]):
    (L, L) = (H_a.nrows(), H_a.ncols())#Assuming the H_a matrix is L by L
    if beta == []:
        beta = vector(RR, [1]*L)
    for be in list(beta):
        if be <= 0:
            return 0
    B = diagonal_matrix(beta)
    P_t=P_source
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
    return (source_rate_list, A_best_LLL)



#compute the new CCF system sum rate
#-------------------
#    Let Source Power be L variables to be optimized!!!
#----------------------
def CCF_new_sumrate_func(betaScale, P_source, H_a, rate_sec_hop, per_c):
    #compute the source rate upbound and matrix A
    source_rate_upbound_list, A  =CCF_fix_pow_sourceRate_upbound(P_source,H_a,betaScale)
    #the second hop channel capacity constriant
    #SecChannel_constiant=ComputeSecRate(M,P_relay,H_b)
    SecChannel_constiant=rate_sec_hop

    beta2_power=list([betaScale[i]**2*P_source[i] for i in range(0,L)])
    
    beta=copy.copy(beta2_power)
    per_s=[]#nested shaping lattice order
    max_beta=0
    for i in range(L):
        max_beta=max(beta)
        max_beta_index=beta.index(max_beta)
        per_s.append(max_beta_index)
        beta[max_beta_index]=0
    
    #compute the coefficient of rate pieces of the conditional entropy

    if per_c == []:
        max_sum_rate = 0
        for code_order in itertools.permutations(list(range(0, L)), L):
            per_c = list(code_order)
            entropy_coefficient_list=Relay_Forward_Rate(per_s,per_c,A)
            t1=time.time()
            Res=Linear_Program(entropy_coefficient_list,SecChannel_constiant,source_rate_upbound_list,per_s,per_c,beta2_power)
            #========#
            t2=time.time()
            #print 'linear programming time cost:', t2 - t1
            if max_sum_rate > Res:
                max_sum_rate = copy.copy(Res)
        Res = copy.copy(max_sum_rate)
    else:
        entropy_coefficient_list = Relay_Forward_Rate(per_s, per_c, A)
        t1 = time.time()
        Res = Linear_Program(entropy_coefficient_list, SecChannel_constiant, source_rate_upbound_list, per_s, per_c,
                             beta2_power)
        # ========#
        t2 = time.time()
        print 'linear programming time cost:', t2 - t1
    return Res
    
    
    #-------------------------------------------
    #        The main optimization function
    #------------------------------------------
def RandomSearch(P_Search_Alg, H_a, rate_sec_hop, P_con, per_c=[]):
    '''
    CCF_beta_func=lambda x: CCF_sumrate_compute(vector(RR, [1,]+list(x[0:L-1])), H_a, H_b, P_con, P_relay, per_s, per_c)
    '''
    fix_pow = True
    GCCF = True
    if fix_pow:  # fixed source power
        if GCCF:
            CCF_beta_func = lambda x: GCCF_new_sumrate_func(vector(RR, [1, ] + list(x[0:L - 1])), [P_con] * L, H_a,
                                                       rate_sec_hop, per_c)
        else:
            CCF_beta_func = lambda x: CCF_new_sumrate_func(vector(RR, [1, ] + list(x[0:L - 1])), [P_con] * L, H_a,
                                                       rate_sec_hop, per_c)
        Pranges = ((0.01, betaScale_max),) * (L - 1)  # L beta and L source power
    else:  # optimize source power
        if GCCF:
            CCF_beta_func = lambda x: GCCF_new_sumrate_func(vector(RR, list(x[0:L])), x[L:2 * L], H_a, rate_sec_hop, per_c)
        else:
            CCF_beta_func = lambda x: CCF_new_sumrate_func(vector(RR, list(x[0:L])), x[L:2 * L], H_a, rate_sec_hop, per_c)
        Pranges = ((0.01, betaScale_max),) * (L) + ((0.01, P_con),) * L  # L beta and L source power

    if P_Search_Alg == 'differential_evolution':
        # test program running time cost
        t1 = time.time()
        try:
            # set_random_seed()
            # seed_int = np.random.randint(1,100)
            seed_int = randint(1, 100)
            #             print 'NCCF seed: ', seed_int
            # return [0]*(L),0
            ResSearch = optimize.differential_evolution(CCF_beta_func, Pranges, maxiter=50, seed=seed_int,
                                                        disp=False)
        except:
            print 'error in differential evolution algorithm'
            raise
        t2 = time.time()
        t = t2 - t1
        print 'New CCF Differential Evolution:', ResSearch.success
        beta_opt = ResSearch.x
        # print 'optimal beta and source power(new)', beta_opt
        sum_rate_opt = -ResSearch.fun
    return vector(RR, [1,] + list(beta_opt)), sum_rate_opt
    

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
    