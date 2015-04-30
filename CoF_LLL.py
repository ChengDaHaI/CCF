'''
This is something related with lattice reduction algorithms.
'''
from sage.all import *
import copy
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
from CoF_second_hop import *


# If the rates at the transmitters are not supported by the second hop, set the rate
# as the largest rates that the second hop can support.
'''
def CoF_compute_fixed_pow_flex_fine_lattice(P_t, H_a, rate_sec_hop, beta=[]):
    (M, L) = (H_a.nrows(), H_a.ncols())
    if beta == []:
        beta = vector(RR, [1]*L)
    try:
        P_t[0]
    except:
        P_t = [P_t]
    
    for i_P in range(0, len(P_t)):
        if math.isnan(P_t[i_P]):
            print 'P', str(i_P), ' should not be NaN!'
            raise Exception('Invalid power setting reached.')
        if P_t[i_P] <= 0:
            # print 'P', str(i_P), ' should be positive'
            return 0
        
    P_vec = vector(RR, P_t)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, rate_list_A_LLL, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, True, beta)
    except:
        print 'error in seeking A and rate'
        raise
    rate_list = list(rate_list_A_LLL)
    for i_l in range(0, L):
        for i_m in range(0, M):
            if (A_best_LLL[i_m, i_l]!=0) and (rate_sec_hop[i_m]<rate_list[i_l]):
                rate_list[i_l] = rate_sec_hop[i_m]
    sum_rate = sum(rate_list)
    return sum_rate
'''

def CoF_compute_fixed_pow_flex(P_t, P_con, is_return_A, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sym_mod', quan_scheme='sym_quan', beta=[]):
    
    # print 'P: '; print P_t
    # print 'beta = ', list(beta)
    
    (M, L) = (H_a.nrows(), H_a.ncols())
    if beta == []:
        beta = vector(RR, [1]*L)
    for be in list(beta):
        if be <= 0:
            # print 'beta = ', list(beta), 'should be positive'
            if is_return_A == True:
                return (0, zero_matrix(M, L))
            else:
                return 0
    B = diagonal_matrix(beta)
    try:
        P_t[0]
    except:
        P_t = [P_t]
        
    for i_P in range(0, L):
        if math.isnan(P_t[i_P]):
            print 'P', str(i_P), ' should not be NaN!'
            raise Exception('Invalid power setting reached.')
        if P_t[i_P] <= 0 or P_t[i_P] > (P_con+0.1):
            # print 'P', str(i_P), ' should be positive and less than P_con'
            if is_return_A == True:
                return (0, zero_matrix(M, L))
            else:
                return 0
    
    P_vec = vector(RR, P_t)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, sum_rate_A_LLL, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, False, beta)
    except:
        print 'error in seeking A and rate'
        raise
    
    A_best_LLL_F = matrix(GF(p), A_best_LLL)
    if A_best_LLL_F.rank() == min(L, M):
        try:
            if is_dual_hop == True:
                '''constraints of the second hop'''
                # relay_fine_lattices is already obtained
                # compute the coarse lattice of the l-th transmitter
                
                # The true coarse lattices have scaling factor beta.
                trans_coarse_lattices = list(P_vec.pairwise_product(vector([b**2 for b in beta]))) # copy
                
                # check whether the second-hop constraint rate_sec_hop can support the first-hop rate r
                try:
                    support_rates = RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A_best_LLL, rate_sec_hop, mod_scheme, quan_scheme))
                except:
                    print 'error in second hop'
                    raise
            else:
                support_rates = sum_rate_A_LLL
        except:
            print 'error in checking second hop'
            raise
    else:
        support_rates = 0
    if is_return_A == True:
        return (support_rates, A_best_LLL)
    else:
        return support_rates
        

def Find_A_m_list(P, H, beta):
    # P is a LxL matrix, H is a MxL matrix
    (M, L) = (H.nrows(), H.ncols())
    A_m_list = []
    B = diagonal_matrix(beta)
    for i_a in range(0, M):
        h = H.row(i_a)
        D_m = B*P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm())**2)))*P*B
        # D_m = P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm())**2)))*P
        D_m_rdf = D_m.change_ring(RDF)
        if D_m_rdf.is_positive_definite() == False:
            # raise Exception('G should be positive definite!')
            # FIXME
            print 'WARNNING! G should be positive definite!'
            return []
        F = D_m_rdf.cholesky()
        if F.rank() != L:
            raise Exception('F should be full rank. Check it!')
        amp = 10**5
        F_a = amp*F
        F_a_rdf = F_a.change_ring(RDF)
        F_a_z = matrix(ZZ, L, L, [round(x) for x in F_a_rdf.list()])
        F_a_reduced = F_a_z.LLL(use_givens=True)
        F_reduced = F_a_reduced/amp
        F_reduced_list = F_reduced.rows()
        F_reduced_list = sorted(F_reduced_list, key=lambda x:x.norm(), reverse=False)
        F_reduced = matrix(F_reduced_list)
        FF = F_reduced*F.inverse()
        FF_z = matrix(ZZ, L, L, [round(x) for x in FF.list()])
        A_m_list += [FF_z]
    return A_m_list

def Find_A_and_Rate(P_mat, P_vec, H, is_return_rate_list=False, beta=[]):
    # Use LLL to find a good A matrix
    (M, L) = (H.nrows(), H.ncols())
    if beta == []:
        beta = vector(RR, [1]*L)
    A_list = Find_A_m_list(P_mat, H, beta)
    if not A_list:
        print 'WARNING! empty A_list. Will return zero rate and matrix.'
        if is_return_rate_list == True:
            return (zero_matrix(RR, M, L), [0]*L, float('inf')*M)
        else:
            return (zero_matrix(RR, M, L), 0, float('inf')*M)
    A = zero_matrix(ZZ, M, L)
    for i_a in range(0, M):
        A.set_row(i_a, A_list[i_a].row(0))
    A_F = matrix(GF(p), A)
    rank_first_row = A_F.rank()
#     alpha_opt = zero_vector(RR, M)
    relay_fine_lattices = [0]*M
    sum_rate = 0
    rate_list = [0]*L
    A_best = zero_matrix(ZZ, M, L)
    if rank_first_row == min(L, M):
        # full rank
        try:
            (rate, relay_fine_lattices) = rate_computation_MMSE_alpha(L, M, P_vec, H, A, beta)
        except:
            print 'error in rate_computation_MMSE_alpha: pos 0'
            raise
        sum_rate = sum(rate)
        rate_list = list(rate)
        A_best = A
    else:
        # rank deficiency: search for full-rank matrix A
        D = A.nrows() # search in [0, D-1] in each matrix in the list
        for i_s in range(0, D**M):
            A = zero_matrix(RR, M, L)
            for i_a in range(0, M):
                idx_row_i_a = (i_s%(D**(i_a+1)))/(D**i_a)
                A.set_row(i_a, A_list[i_a].row(idx_row_i_a))
            A_F = matrix(GF(p), A)
            if A_F.rank() == min(L, M):
                try:
                    (rate, relay_fine_lattices_i_row_search) = rate_computation_MMSE_alpha(L, M, P_vec, H, A, beta)
                except:
                    print 'error in rate_computation_MMSE_alpha: pos 1'
                    raise
                if sum_rate < sum(rate):
                    sum_rate = sum(rate)
                    rate_list = list(rate)
#                     alpha_opt = alpha_opt_i_row_search
                    relay_fine_lattices = relay_fine_lattices_i_row_search
                    A_best = A
            else:
                pass
#     A_best_F = matrix(GF(p), A_best)
#     if A_best_F.rank() < min(L, M):
#         pass # for test
#     if A_best == matrix(GF(p), 2, 2, [[1, 2], [2, 1]]):
#         pass
    A_best_F = matrix(GF(p), A_best)
    if A_best_F.rank() < min(L, M):
        #raise Exception('Even the best coefficient matrix is not full-rank in finite field')
        rate_list = [0]*L
        sum_rate = 0
    
    if is_return_rate_list == True:
        return (A_best, rate_list, relay_fine_lattices)
    else:
        return (A_best, sum_rate, relay_fine_lattices)

    
@parallel(ncpus=Cores)
def CoF_rank_deficiency(P_con):
    iter_H = 2000
    rank_deficiency = 0
    M = 2
    L = 2
    
    iter_d = 10000
    for i_d in range(0, iter_d):
        H = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        P = diagonal_matrix((vector(RR, [1, 1])+random_vector(RR, 2))*100000)
        h = H.row(1)
        beta = vector(RR, [1, 1])
        B = diagonal_matrix(beta)
        D_m = B*P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm())**2)))*P*B
        # D_m = P*(identity_matrix(RR, L)-(P.T*h.column()*h.row()*P/(1+((h.row()*P).norm())**2)))*P
        D_m_rdf = D_m.change_ring(RDF)
        if D_m_rdf.is_positive_definite() == False:
            # raise Exception('G should be positive definite!')
            # FIXME
            print 'WARNNING! G should be positive definite!'
            print 'h = ', h, 'D_m = ', D_m
            return []
    print 'finish: no error found'
    return []
    
#     for i_H in range(0, iter_H):
#         set_random_seed() # to avoid producing the same H in different threads
#         H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
#         # H_a = matrix(RR, M, L, [[-0.869243498832414, 0.762538785616652], [-0.914996809982352, 0.801032403084523]])
#         # Fixed power
#         P_mat = sqrt(P_con)*identity_matrix(L)
#         P_vec = P_con*vector([1]*L)
#         
#         # Use LLL to find a good A matrix
#         (A, dummy, dummy) = Find_A_and_Rate(P_mat, P_vec, H_a, is_return_rate_list=False)
#         
#         A_F = matrix(GF(p), A)
#         alpha_opt = zero_vector(RR, M)
#         if A_F.rank() == min(L, M):
#             pass
#         else:
#             # rank deficiency: outage
#             sum_rate_A = 0
#             rank_deficiency += 1
#     
#     return float(rank_deficiency)/iter_H


if __name__ == "__main__":
    P_eq_dB_Min = float(20)
    P_eq_dB_Max = float(40)
    P_delta = 5
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    t1 = time.ctime()
    result = list(CoF_rank_deficiency(Pl_con))
    t2 = time.ctime()
    
    print 'Simulation started at ', t1
    print 'Simulation ended at ', t2
    rank_deficiency = [result[i][1] for i in range(0, len(Pl_con))]
    plot_rank_deficiency = list_plot(zip(P_eq_dB, rank_deficiency), plotjoined=True, marker='o', \
                              rgbcolor=Color('red'), linestyle="--", \
                              legend_label= 'probability', \
                              title = 'Rank deficiency probability')
    print 'rank_deficiency = ', rank_deficiency
    show(plot_rank_deficiency)
    
    raw_input() # stop Sage from shutting down
    
    
