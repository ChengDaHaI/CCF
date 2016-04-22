'''
This includes some functions related with the second hop in CoF.
'''
from sage.all import *
from scipy import optimize
from numpy import arange
from sage.parallel.all import *
import time
from CoF_basic import *
import itertools
import copy

# A: matrix, mod_order: permutation
def check_feasible_permutation(A, decode_order, mod_order):
    mod_order_list = list(mod_order) # note: not equal to pi_e()
    is_decodable = True
    A_F = matrix(GF(p), A)
    A_i_L = A_F
    x = list(decode_order) # copy
    for i_L in range(0, L):
        if A_i_L.rank() < L-i_L:
            is_decodable = False
            break
        elif A_i_L.rank() == L-i_L:
            x_min_idx, x_min_val = min(enumerate(x), key=lambda x:x[1])
            x.pop(x_min_idx)
            if A_i_L.rank() > 1:
                A_i_L = A_i_L.delete_columns([x_min_idx])
                row_thresh = mod_order_list[x_min_idx]
                A_i_L = A_i_L.delete_rows([mod_order_list.pop(x_min_idx)])
                for i_row in range(0, len(mod_order_list)):
                    if mod_order_list[i_row] > row_thresh:
                        mod_order_list[i_row] -= 1
                    pass
            else:
                # don't need to reduce since no further operations need it
                pass
        else:
            raise Exception('error: the rank should not exceed L-i_L')
    return is_decodable


def second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, rate_sec_hop, mod_scheme, quan_scheme):
    (M, L) = (A.nrows(), A.ncols())
    if M != L:
        raise Exception("L and M should be the same in destination's perspective.")
    
    if mod_scheme == 'asym_mod':
        # iterate all permutations 
        sum_rate_max = 0
        r_opt=[0]*L# optimal source rate
        trans_shaping_opt = [0]*L# optimal shaping lattice
        trans_coding_opt  = [0]*L# optimal coding lattice
        has_feasible_mod = False
        for mod_order in itertools.permutations(list(range(0, M)), L):
            # each element in mod_order: the relay that the l-th lattice should be assigned to
            is_mod_decodable = check_feasible_permutation(A, trans_coarse_lattices, mod_order)
            if is_mod_decodable == True:
                has_feasible_mod = True
                # if decodable, then calculate how much sum rate it can support
                mod_order_list = list(mod_order)
                # determine the coarse lattices at the relays
                relay_coarse_lattices = [0]*M
                for i_m in range(0, M):
                    relay_coarse_lattices[i_m] = trans_coarse_lattices[mod_order_list.index(i_m)]

                if quan_scheme == 'sym_quan':
                    relay_actual_fine_lattices = list(relay_fine_lattices)
                    # determine the achievable fine lattices at the relays
                    for i_m in range(0, M):
                        relay_actual_fine_lattices[i_m] = max(relay_actual_fine_lattices[i_m], relay_coarse_lattices[i_m]/(2**(2*rate_sec_hop[i_m])))
                    # determine the fine lattice of the l-th transmitter
                    trans_fine_lattices = [float(0)]*L
                    for i_L in range(0, L):
                        for i_M in range(0, M):
                            if (A[i_M, i_L]!=0) and (relay_actual_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                                trans_fine_lattices[i_L] = relay_actual_fine_lattices[i_M]
                                
                    r = [0]*L
                    for i_l in range(0, L):
                        r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
                    sum_rate = sum(r)
                    if sum_rate_max < sum_rate:
                        sum_rate_max = sum_rate
                elif quan_scheme == 'asym_quan':
                    has_feasible_quan = False
                    
                    if False: # comment Tan's method
                        # Tan's method to find source coding lattices and quan_order seems to be wrong!
                        for quan_order in itertools.permutations(list(range(0, M)), L):
                            relay_compute_fine_lattices = list(relay_fine_lattices)
                            is_quan_decodable = check_feasible_permutation(A, [-relay_compute_fine_lattices[i] for i in range(0, L)], quan_order)
                            if is_quan_decodable == True:
                                has_feasible_quan = True
                                quan_order_list = list(quan_order)
                                # determine the fine lattice of the l-th transmitter according to computation constraints
                                trans_compute_fine_lattices = [float(0)]*L
                                for i_L in range(0, L):
                                    for i_M in range(0, M):
                                        if (A[i_M, i_L]!=0) and (relay_compute_fine_lattices[i_M]>trans_compute_fine_lattices[i_L]):
                                            trans_compute_fine_lattices[i_L] = relay_compute_fine_lattices[i_M]
                                # map the trans_compute_fine_lattices to the relays
                                relay_map_fine_lattices = [0]*M
                                for i_m in range(0, M):
                                    relay_map_fine_lattices[i_m] = trans_compute_fine_lattices[quan_order_list.index(i_m)]
                                
                                # determine the quantization lattice at the m-th relay
                                relay_quan_fine_lattices = [0]*M
                                for i_m in range(0, M):
                                    relay_quan_fine_lattices[i_m] = max(relay_map_fine_lattices[i_m], relay_coarse_lattices[i_m]/(2**(2*rate_sec_hop[i_m])))
                                
                                # map the relay_quan_fine_lattices back to the transmitters
                                trans_fine_lattices = [0]*L
                                for i_l in range(0, L):
                                    trans_fine_lattices[i_l] = relay_quan_fine_lattices[quan_order_list[i_l]]
                                
                                # calculate transmission sum rates
                                r = [0]*L
                                for i_l in range(0, L):
                                    r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
                                sum_rate = sum(r)
                                if sum_rate_max < sum_rate:
                                    sum_rate_max = sum_rate
                                    r_opt = copy.copy(r)
                                    trans_shaping_opt = copy.copy(trans_coarse_lattices)
                                    trans_coding_opt  = copy.copy(trans_fine_lattices)
                                    mod_order_opt = copy.copy(mod_order_list)
                                    quan_order_opt = copy.copy(quan_order_list)
                    else:
                        # My revised method to find source coding lattices and quan_order.
                        # I just change to check the feasibility of quan_order at the last step.
                        for quan_order in itertools.permutations(list(range(0, M)), L):
                            relay_compute_fine_lattices = list(relay_fine_lattices)
                            '''
                            is_quan_decodable = check_feasible_permutation(A, [-relay_compute_fine_lattices[i] for i in range(0, L)], quan_order)
                            '''
                            quan_order_list = list(quan_order)
                            # determine the fine lattice of the l-th transmitter according to computation constraints
                            trans_compute_fine_lattices = [float(0)]*L
                            for i_L in range(0, L):
                                for i_M in range(0, M):
                                    if (A[i_M, i_L]!=0) and (relay_compute_fine_lattices[i_M]>trans_compute_fine_lattices[i_L]):
                                        trans_compute_fine_lattices[i_L] = relay_compute_fine_lattices[i_M]
                            # map the trans_compute_fine_lattices to the relays
                            relay_map_fine_lattices = [0]*M
                            for i_m in range(0, M):
                                relay_map_fine_lattices[i_m] = trans_compute_fine_lattices[quan_order_list.index(i_m)]
                            
                            # determine the quantization lattice at the m-th relay
                            relay_quan_fine_lattices = [0]*M
                            for i_m in range(0, M):
                                relay_quan_fine_lattices[i_m] = max(relay_map_fine_lattices[i_m], relay_coarse_lattices[i_m]/(2**(2*rate_sec_hop[i_m])))
                            
                            # map the relay_quan_fine_lattices back to the transmitters
                            trans_fine_lattices = [0]*L
                            for i_l in range(0, L):
                                trans_fine_lattices[i_l] = relay_quan_fine_lattices[quan_order_list[i_l]]
                            
                            # check the feasibility of quan_order with given trans_fine_lattices
                            is_quan_decodable = check_feasible_permutation(A, [-trans_fine_lattices[i] for i in range(0, L)], quan_order)
                            
                            if is_quan_decodable:
                                has_feasible_quan = True
                                # calculate transmission sum rates
                                r = [0]*L
                                for i_l in range(0, L):
                                    r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
                                sum_rate = sum(r)
                                if sum_rate_max < sum_rate:
                                    sum_rate_max = sum_rate
                                    r_opt = copy.copy(r)
                                    trans_shaping_opt = copy.copy(trans_coarse_lattices)
                                    trans_coding_opt  = copy.copy(trans_fine_lattices)
                                    mod_order_opt = copy.copy(mod_order_list)
                                    quan_order_opt = copy.copy(quan_order_list)
                                
                    if has_feasible_quan == False:
                        raise Exception('If Q is full rank, then there must be at leat one feasible quantizatioin way. But no one found.')
                else:
                    raise Exception('quan_scheme should take value as "asym_quan" or "sym_quan"')
                
        # if no successful scheme found,
        if has_feasible_mod == False: 
            raise Exception('If Q is full rank, then there must be at leat one feasible modulo way. But no one found.')
        # return the optimal sum rate , rate tuple, shaping lattice , and coding lattice
        return sum_rate_max, r_opt, trans_shaping_opt, trans_coding_opt, mod_order_opt, quan_order_opt
    
    elif mod_scheme == 'sym_mod' and quan_scheme == 'asym_quan':
        relay_coarse_lattice = max(trans_coarse_lattices)
        # iterate all permutations 
        sum_rate_max = 0
        has_feasible_quan = False
        for quan_order in itertools.permutations(list(range(0, M)), L):
            relay_compute_fine_lattices = list(relay_fine_lattices)
            is_quan_decodable = check_feasible_permutation(A, [-relay_compute_fine_lattices[i] for i in range(0, L)], quan_order)
            if is_quan_decodable == True:
                has_feasible_quan = True
                quan_order_list = list(quan_order)
                # determine the fine lattice of the l-th transmitter according to computation constraints
                trans_compute_fine_lattices = [float(0)]*L
                for i_L in range(0, L):
                    for i_M in range(0, M):
                        if (A[i_M, i_L]!=0) and (relay_compute_fine_lattices[i_M]>trans_compute_fine_lattices[i_L]):
                            trans_compute_fine_lattices[i_L] = relay_compute_fine_lattices[i_M]
                # map the trans_compute_fine_lattices to the relays
                relay_map_fine_lattices = [0]*M
                for i_m in range(0, M):
                    relay_map_fine_lattices[i_m] = trans_compute_fine_lattices[quan_order_list.index(i_m)]
                
                # determine the quantization lattice at the m-th relay
                relay_quan_fine_lattices = [0]*M
                for i_m in range(0, M):
                    relay_quan_fine_lattices[i_m] = max(relay_map_fine_lattices[i_m], RR(relay_coarse_lattice/(2**(2*rate_sec_hop[i_m]))))
                
                # map the relay_quan_fine_lattices back to the transmitters
                trans_fine_lattices = [0]*L
                for i_l in range(0, L):
                    trans_fine_lattices[i_l] = relay_quan_fine_lattices[quan_order_list[i_l]]
                
                # calculate transmission sum rates
                r = [0]*L
                for i_l in range(0, L):
                    r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
                sum_rate = sum(r)
                if sum_rate_max < sum_rate:
                    sum_rate_max = sum_rate
        if has_feasible_quan == False:
            raise Exception('If Q is full rank, then there must be at leat one feasible quantizatioin way. But no one found.')
        return sum_rate_max
        
    elif mod_scheme == 'sym_mod' and quan_scheme == 'sym_quan':
        relay_actual_fine_lattices = list(relay_fine_lattices)
        relay_coarse_lattice = max(trans_coarse_lattices)
        # determine the achievable fine lattices at the relays
        for i_m in range(0, M):
            R_m = max(0, 0.5*log(relay_coarse_lattice/relay_actual_fine_lattices[i_m], 2))
            if rate_sec_hop[i_m] < R_m:
                relay_actual_fine_lattices[i_m] = relay_coarse_lattice/(2**(2*rate_sec_hop[i_m]))
        # determine the fine lattice of the l-th transmitter
        trans_fine_lattices = [float(0)]*L
        for i_L in range(0, L):
            for i_M in range(0, M):
                if (A[i_M, i_L]!=0) and (relay_actual_fine_lattices[i_M]>trans_fine_lattices[i_L]):
                    trans_fine_lattices[i_L] = relay_actual_fine_lattices[i_M]
        # calculate transmission sum rates
        r = [0]*L
        for i_l in range(0, L):
            r[i_l] = max(0, 0.5*log(trans_coarse_lattices[i_l]/trans_fine_lattices[i_l], 2))
        sum_rate = sum(r)
        return sum_rate
        #return r #return rate tuple
    else:
        raise Exception("mod_scheme should take value as 'asym_mod' or 'sym_mod'! \
            quan_scheme should take value as 'asym_quan' or 'sym_quan'!")



if __name__ == "__main__":
    print '-----------------------------------\n'+ \
        'testing CoF_second_hop\n'
#     fine_lattices = [0.5, 0.2]
#     coarse_lattices = [1.5, 2.1]
#     R = [1.5, 1.7]
#     A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
#     p = 3 
#     print second_hop(fine_lattices, coarse_lattices, A, R, 'asym_mod')
#     
#     fine_lattices = [0.5, 0.2, 0.4]
#     coarse_lattices = [1.5, 2.1, 1]
#     R = [1.5, 1.7, 3]
#     A = matrix(ZZ, 3, 3, [[0, 0, 2], [2, 2, 0], [0, 2, 2]])
#     print second_hop(fine_lattices, coarse_lattices, A, R, 'asym_mod')
    
#     R = [1, 2]
#     A = matrix(ZZ, 2, 2, [[1, 2], [1, 1]])
#     p = 3
#     relay_fine_lattices = [0.5, 0.2]
#     trans_coarse_lattices = [1.5, 2.0]
#     beta = [1]*2
#     print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'asym_mod', 'asym_quan'))
    # should be 1.792
    
#     R = [1, 10]
#     A = matrix(ZZ, 2, 2, [[1, 2], [0, 1]])
#     p = 3
#     relay_fine_lattices = [0.5, 0.2]
#     trans_coarse_lattices = [1.5, 2.0]
#     print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'asym_mod', 'asym_quan'))
#     # should be 1.792
# 
#     R = [1, 1]
#     A = matrix(ZZ, 2, 2, [[0, 2], [1, 1]])
#     p = 3
#     relay_fine_lattices = [0.5, 0.2]
#     trans_coarse_lattices = [1.5, 2.0]
#     print RR(second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A, R, 'asym_mod', 'asym_quan'))
    # should be 1.792
    
    #test 
    
    A = matrix(ZZ, 3, 3, [[1,0,0], [-1,0,1], [1,-1,1]])
    print 'Integer matrix A: ', A
    trans_fine_lattices = [17,70, 10.45, 17.74]
    print 'coding lattice at transmitter:', trans_fine_lattices
    quan_order = [2, 0, 1]
    print 'quan_order at relays: ', quan_order
    is_quan_decodable = check_feasible_permutation(A, [-trans_fine_lattices[i] for i in range(0, L)], quan_order)
    print is_quan_decodable
    