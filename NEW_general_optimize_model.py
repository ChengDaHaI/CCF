'''
Simulation of Relay Compressing in GCCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
V 1.0
'''
from sage.all import *
from NEW_CCF_Modle import Relay_Forward_Rate, Powerset
from NewSecondHopChannel import ComputeSecRate
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

# 2*2 channel
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

# For a 3*3 channel, we compute the sum rate 
def GCCF_new_sumrate_func(betaScale, P_source, H_a, rate_sec_hop, per_c = []):

    (A, source_rate_list, coding_lattice_lowerbound) = CCF_computation_SourceRate_upbound(P_source, H_a, betaScale)
    shaping_lattice_voronoi = [P_source[i]*betaScale[i]**2 for i in range(L)]
    # compute per_s
    per_s = [0]*L
    coarse_voronoi = copy.copy(shaping_lattice_voronoi)
    for i in range(L):
        per_s[i] = coarse_voronoi.index(max(coarse_voronoi))
        coarse_voronoi[per_s[i]] = 0
    
    #produce subset list
    subset_list=list(Powerset(range(0,L)))
    set_L=set(range(0,L))
    
    # list for linear program
    # fixed rate component
    rate_comp_fix = []
    # -1 represent none
    rate_comp_fix_index = [-1]
    # objective function coefficient
    obj_func = []
    obj_func_fix = [0]
    # coefficients of source rate to be optimize
    source_coef = [0]*L
    # fixed part
    source_coef_fix = [0]*L 
    # coefficients of compression rate region bounds to be optimize
    entropy_coef = [0]*(pow(2,L)-1)
    # fixed part
    entropy_coef_fix = [0]*(pow(2,L)-1) 
    
    if per_c == [0,0,0]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append(0)
        obj_func = [1,1,1]
        source_coef = [[0]*3]*3
        for i in range(L):
            temp = copy.copy(source_coef[i])
            temp[per_s.index(i)] = 1
            source_coef[i] = copy.copy(temp)
        # sub matrix column index, i.e., the source index corresponding to each rate component
        comp_source_set = [[per_s[0]], [per_s[1]], [per_s[2]]]
        
    elif per_c == [1,1,1]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [1]
        obj_func = [1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[1]] = 1
        source_coef = [[0]*3]*3
        for i in range(L):
            temp = copy.copy(source_coef[i])
            if per_s.index(i) != 2:
                temp[per_s.index(i)] = 1
                source_coef[i] = copy.copy(temp)
            else:
                temp[1:3] = [1,1]
                source_coef[i] = copy.copy(temp)
        comp_source_set = [[per_s[0]], [per_s[1]], [per_s[1], per_s[2]], [per_s[2]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             
#             if rate_comp_fix_index != []:
#                 temp = []
#                 for k in range(len(rate_comp_fix_index)):
#                     temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#                 entropy_coef_fix[i-1] = temp
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [2,2,2]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [1]
        obj_func = [1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[1]] = 1
        source_coef = [[0]*3]*3
        source_coef[per_s[0]] = copy.copy([1,0,0])
        source_coef[per_s[1]] = copy.copy([0,1,1])
        source_coef[per_s[2]] = copy.copy([0,1,0])
        comp_source_set = [[per_s[0]], [per_s[1]], [per_s[1], per_s[2]], [per_s[1]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [3,3,3]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = 1
        source_coef = [[0]*3]*3
        source_coef[per_s[0]] = copy.copy([1,0,0])
        source_coef[per_s[1]] = copy.copy([1,1,0])
        source_coef[per_s[2]] = copy.copy([0,0,1])
        comp_source_set = [[per_s[0]], [per_s[0], per_s[1]], [per_s[1]], [per_s[2]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [4,4,4]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = 1
        source_coef = [[0]*3]*3
        source_coef[per_s[0]] = copy.copy([1,1,0])
        source_coef[per_s[1]] = copy.copy([1,0,0])
        source_coef[per_s[2]] = copy.copy([0,0,1])
        comp_source_set = [[per_s[0]], [per_s[0], per_s[1]], [per_s[0]], [per_s[2]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [5,5,5]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        #rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = [1]
        source_coef_fix[per_s[1]] = [0]
        source_coef_fix[per_s[2]] = [0]
        source_coef = [[0]*4]*3
        source_coef[per_s[0]] = copy.copy([1,0,0,0])
        source_coef[per_s[1]] = copy.copy([1,1,1,0])
        source_coef[per_s[2]] = copy.copy([0,0,1,1])
        comp_source_set = [[per_s[0]], [per_s[0],per_s[1]], [per_s[1]], [per_s[1],per_s[2]], [per_s[2]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [6,6,6]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        #rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = [1]
        source_coef_fix[per_s[1]] = [0]
        source_coef_fix[per_s[2]] = [0]
        source_coef = [[0]*4]*3
        source_coef[per_s[0]] = copy.copy([1,0,0,0])
        source_coef[per_s[1]] = copy.copy([1,1,1,1])
        source_coef[per_s[2]] = copy.copy([0,0,1,0])
        comp_source_set = [[per_s[0]], [per_s[0],per_s[1]], [per_s[1]], [per_s[1],per_s[2]], [per_s[1]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [7,7,7]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        #rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = [1]
        source_coef_fix[per_s[1]] = [0]
        source_coef_fix[per_s[2]] = [0]
        source_coef = [[0]*4]*3
        source_coef[per_s[0]] = copy.copy([1,1,1,0])
        source_coef[per_s[1]] = copy.copy([1,0,0,0])
        source_coef[per_s[2]] = copy.copy([0,0,1,1])
        comp_source_set = [[per_s[0]], [per_s[0],per_s[1]], [per_s[1]], [per_s[0],per_s[2]], [per_s[2]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    elif per_c == [8,8,8]:
        print 'when per_c =', per_c, ':\n'
        rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) )
        #rate_comp_fix.append( 0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) )
        rate_comp_fix_index = [0]
        obj_func = [2,1,2,1]
        obj_func_fix = [1]
        source_coef_fix[per_s[0]] = [1]
        source_coef_fix[per_s[1]] = [0]
        source_coef_fix[per_s[2]] = [0]
        source_coef = [[0]*4]*3
        source_coef[per_s[0]] = copy.copy([1,1,1,1])
        source_coef[per_s[1]] = copy.copy([1,0,0,0])
        source_coef[per_s[2]] = copy.copy([0,0,1,0])
        comp_source_set = [[per_s[0]], [per_s[0],per_s[1]], [per_s[1]], [per_s[0],per_s[2]], [per_s[0]]]
#         #for every subset, we first calculate the sub-matrix 
#         #then calculate the coefficients
#         for i in range(1,len(subset_list)):
#             piece_coefficient=[]
#             subset=set(subset_list[i])
#             complement_set=set_L.difference(subset)
#             row=[]#sub-matrix row index
#             row.extend(list(complement_set))
#             
#             rank_value=0
#             for j in range(len(comp_source_set)):
#                 colum = comp_source_set[j] #sub-matrix colum index
#                 sub_A = A[row,colum]
#                 rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
#                 #record the coefficient of rate pieces
#                 #that is also the i row of transform matrix between  
#                 #conditional_entropy list and rate_piece list
#                 piece_coefficient.append(rank_value)
#             temp = []
#             for k in range(len(rate_comp_fix_index)):
#                 temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
#             entropy_coef_fix[i-1] = copy.copy(temp)
#             entropy_coef[i-1] = piece_coefficient
    
    #for every subset, we first calculate the sub-matrix 
    #then calculate the coefficients
    for i in range(1,len(subset_list)):
        piece_coefficient=[]
        subset=set(subset_list[i])
        complement_set=set_L.difference(subset)
        row=[]#sub-matrix row index
        row.extend(list(complement_set))
        rank_value=0
        for j in range(len(comp_source_set)):
            colum = comp_source_set[j] #sub-matrix colum index
            sub_A = A[row,colum]
            rank_value = rank(A[list(set_L), colum]) - rank(sub_A)
            #record the coefficient of rate pieces
            #that is also the i row of transform matrix between  
            #conditional_entropy list and rate_piece list
            piece_coefficient.append(rank_value)
        if rate_comp_fix_index != [-1]:
            temp = []
            for k in range(len(rate_comp_fix_index)):
                temp.append(piece_coefficient.pop(rate_comp_fix_index[0]))
            entropy_coef_fix[i-1] = temp
        entropy_coef[i-1] = piece_coefficient
    
    opt_sum_rate = GCCF_linear_prog(obj_func, obj_func_fix, rate_comp_fix, rate_comp_fix_index, source_coef, source_coef_fix, \
                     entropy_coef, entropy_coef_fix, source_rate_list, shaping_lattice_voronoi, per_c, per_s, rate_sec_hop)
    
    print 'Done!\n'
    return opt_sum_rate

# applying linear programming to optimize the sum rate
def GCCF_linear_prog(obj_func, obj_func_fix, rate_comp_fix, rate_comp_fix_index, source_coef, source_coef_fix, entropy_coef, entropy_coef_fix, source_rate_bound, shaping_lattice_voronoi, per_c, per_s, rate_sec_hop):

    obj_C = np.subtract([0]*len(obj_func), obj_func)
    # computation region constraints
    per_c_constr_A = []
    per_c_constr_b = []
    if per_c == [0,0,0]:
        per_c_constr_A = source_coef[0:2]
        per_c_constr_b = [0.5 * math.log(shaping_lattice_voronoi[per_s[i]]/shaping_lattice_voronoi[per_s[i+1]], 2) for i in [0,1]]
    elif (per_c == [1,1,1])|(per_c == [2,2,2]) :
        per_c_constr_A = [source_coef[0]]
        per_c_constr_b = [0.5 * math.log(shaping_lattice_voronoi[per_s[0]]/shaping_lattice_voronoi[per_s[1]], 2) ]
    elif (per_c == [3,3,3])|(per_c == [4,4,4]) :
        per_c_constr_A = [[1,1,0]]
        per_c_constr_b = [0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2) ]
    else:
        pass
    
    # coefficient matrix of inequality constranit
    constr_A_un = source_coef + entropy_coef + per_c_constr_A

    Part_SourceRate = np.dot(np.transpose(np.array(source_coef_fix)[np.newaxis]), np.array(rate_comp_fix))
    Part_entropy = np.dot(np.transpose(np.array(entropy_coef_fix)[np.newaxis]), np.array(rate_comp_fix))
    if len(Part_SourceRate.shape) == 2:
        Part_SourceRate = Part_SourceRate[0]
    if len(Part_entropy.shape) == 2:
        Part_entropy = Part_entropy[0]
    # vector of inequality constranit
    constr_b_un = (np.subtract(np.array(source_rate_bound), Part_SourceRate.tolist())).tolist() + (np.subtract(np.array(rate_sec_hop), Part_entropy.tolist())).tolist() + per_c_constr_b
    
    # coefficient matrix and vector of equality constranit
    # it is part of fixed source rate, result from per_c
    if per_c[0] >= 5:
        A_eq = [[1,1,0,0]]
        b_eq = [0.5 * math.log(shaping_lattice_voronoi[per_s[1]]/shaping_lattice_voronoi[per_s[2]], 2)]
    else:
        A_eq = None
        b_eq = None
    
    bound = [(0,None)]*len(obj_C.tolist())

    opt_result = optimize.linprog(c = obj_C.tolist(), A_ub = constr_A_un, b_ub = constr_b_un, A_eq = A_eq, b_eq = b_eq, bounds = bound, options = {"disp": False})
    
    return opt_result.fun + ( -sum(Part_SourceRate.tolist()))

if __name__ == '__main__':
    betaScale = vector(RR, [1,1,1])
    P_source = [10,30,18]
    P_relay = 8
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution = RealDistribution('gaussian', 1))
    H_b = matrix.random(RR, 1, M, distribution = RealDistribution('gaussian', 1))
#     H_a = matrix(RR, 3, 3, [[  0.699276348994144,  0.979966803608800,   0.731095879215959],
#                             [-0.0467540769081729,  0.253952649489866,  -0.775007249377646],
#                             [  0.669860797898610,  0.183608409840893,  -0.343789835395397]])
    rate_sec_hop=ComputeSecRate(M,P_relay,H_b)
    #rate_sec_hop = [2]*7
    for i in range(9):
        per_c = [i]*3
        print  GCCF_new_sumrate_func(betaScale, P_source, H_a, rate_sec_hop, per_c)
        
    
    