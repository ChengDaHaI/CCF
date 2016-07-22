'''
Simulation of Relay Compressing in GCCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
V 1.0
'''
from sage.all import *
from NEW_CCF_Modle import Relay_Forward_Rate, Powerset
from scipy import optimize
from CoF_LLL import Find_A_and_Rate
from NEW_basic import *
import math
import time
import itertools
import copy
import numpy as np

'''
Function to compute the sum rate of generalized CCF 
input: source power P_source, scale factor betaScale, first hop channel matrix H_a, second hop channel capacity region rate_sec_hop
'''

#compute the source rate upbound and matrix A when given variable beta
# the output would be the input of linear programme function
def CCF_computation_SourceRate_upbound(P_source,H_a,beta=[]):
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
    return (A_best_LLL, source_rate_list, relay_fine_lattices)

def GCCF_coding_search_func(coding_lattice_voronoi, shaping_lattice_voronoi, per_s, A, rate_sec_hop):
    
    sum_rate = 0
    for i in range(len(coding_lattice_voronoi)):
        sum_rate = sum_rate + 0.5 * math.log(shaping_lattice_voronoi[i]/coding_lattice_voronoi[i], 2)
    # compute per_c
    per_c = [0]*L
    fine_voronoi = copy.copy(list(coding_lattice_voronoi))
    for i in range(L):
        per_c[i] = fine_voronoi.index(max(fine_voronoi))
        fine_voronoi[per_c[i]] = 0
    
    # rate component
    r_comp = [0]*(2*L - 1)
    
    # determine nested lattice chaine
    if shaping_lattice_voronoi[per_s[1]] > coding_lattice_voronoi[per_c[0]]:
        r_comp[0] = 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2)
        r_comp[1] = 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/coding_lattice_voronoi[per_c[0]], 2)
        r_comp[2] = 0.5 * math.log(coding_lattice_voronoi[per_c[0]]/coding_lattice_voronoi[per_c[1]], 2)
        comp_source_set = [[per_s[0]], [0,1], [per_c[1]]]
    else:
        r_comp[0] = 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/coding_lattice_voronoi[per_c[0]], 2)
        r_comp[1] = 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/coding_lattice_voronoi[per_c[1]], 2)
        r_comp.pop(2)
        comp_source_set = [[per_s[0]], [per_c[1]]]
        
    #produce subset list
    subset_list=list(Powerset(range(0,L)))
    set_L=set(range(0,L))
    
    #coefficients of compression rate region bounds
    entropy_coefficient=[0]*(pow(2,L)-1)
    
    #for every subset, we first calculate the sub-matrix 
    #then calculate the coefficients
    for i in range(1,len(subset_list)):
        piece_coefficient=[]
        '''
        conditional_entropy[i-1]=rate_total
        '''
        subset=set(subset_list[i])
        complement_set=set_L.difference(subset)
        row=[]#sub-matrix row index
        row.extend(list(complement_set))
        
        rank_value=0
        for j in range(len(comp_source_set)):
            colum = comp_source_set[j] #sub-matrix colum index
            sub_A=A[row,colum]
            rank_value=rank(A[list(set_L), colum]) - rank(sub_A)
            #record the coefficient of rate pieces
            #that is also the i row of transform matrix between  
            #conditional_entropy list and rate_piece list
            piece_coefficient.append(rank_value)
        entropy_coefficient[i-1]=piece_coefficient
    
    # compute the compression rate region bound
    compr_bound = [np.dot(entropy_coefficient[i], r_comp) for i in range(pow(2,L)-1)]
    for i in range(pow(2,L)-1):
        if compr_bound[i] > rate_sec_hop[i]:
            sum_rate = 0
            break
        
    return -sum_rate

def GCCF_new_sumrate_func(betaScale, P_source, H_a, rate_sec_hop):
    
    (A, source_rate_list, coding_lattice_lowerbound) = CCF_computation_SourceRate_upbound(P_source, H_a, betaScale)
    shaping_lattice_voronoi = [P_source[i]*betaScale[i]**2 for i in range(L)]
    # compute per_s
    per_s = [0]*L
    coarse_voronoi = copy.copy(shaping_lattice_voronoi)
    for i in range(L):
        per_s[i] = coarse_voronoi.index(max(coarse_voronoi))
        coarse_voronoi[per_s[i]] = 0
    
    # optimize the voronoi of coding lattice 
    GCCF_coding_func = lambda x: GCCF_coding_search_func(vector(RR, [1, ] + list(x[0:L - 1])), shaping_lattice_voronoi, per_s, A,
                                                       rate_sec_hop)
    Pranges = [(coding_lattice_lowerbound[i], shaping_lattice_voronoi[i]) for i in range(L)]
    t1 = time.time()
    try:
        seed_int = randint(1, 100)
        ResSearch = optimize.differential_evolution(GCCF_coding_func, Pranges, maxiter=50, seed=seed_int,
                                                    disp=False)
    except:
        print 'error in differential evolution algorithm'
        raise
    t2 = time.time()
    t = t2 - t1
    print 'New CCF Differential Evolution:', ResSearch.success
    coding_lattice_opt = ResSearch.x
    sum_rate_opt = ResSearch.fun
    
    return sum_rate_opt
    