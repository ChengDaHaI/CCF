import copy
from sage.all import *
import math


L = 2 # L transmitters
M =2# M relays
p = 17 # The prime number

Cores = 8 # The number of CPU cores used in parallel computing
DEBUG_H = False # When this value is True, the channel matrix H is set as certain matrices
P_MIN = 25
P_MAX = 40
P_DEL = 5
P_Search_Alg='brute'
# 'brute', 'TNC', 'anneal', 'brute_fmin', 'brute_fmin_cobyla', 'brute_fmin_beta', 'brute_fmin_cobyla_beta'
# 'brute_brute'
#genetic
#differential_evolution
brute_number = 100
brute_fmin_number = 20
brute_fmin_maxiter = 70
brute_brute_first_number = 25
brute_brute_second_number = 25

beta_max = 4
if mod(brute_number, beta_max) != 0:
    raise Exception('Error: beta=1 will not be included!')

is_alternate = False # True or False
is_set_H = False # given channel

if is_set_H == True:
    set_H_a = matrix(RR, 2, 2, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
    set_H_b = vector(RR, [0.5, 0.5])
is_set_beta = False # set beta for given channel
if is_set_beta == True:
    set_beta = vector(RR, [1, 1.8])
iter_H = 8
batch_H = 1

#This is the key function in terms of rate!!!
# P is a LxL matrix P_mat
def rate_computation_MMSE_alpha(L, M, P_t, H, A, beta):
    for i_P in range(0, len(P_t)):
        if math.isnan(P_t[i_P]):
            print 'P', str(i_P), ' should not be NaN!'
            return ([0]*L, [float('inf')]*M)
        if P_t[i_P] <= 0:
            print 'P', str(i_P), ' should be positive'
            return ([0]*L, [float('inf')]*M)
    P = matrix.diagonal([sqrt(x) for x in P_t])
    
    r = zero_vector(RR, L)
    phi = [float('inf')]*M
    try:
        for i_m in range(0, M):
            h_m = H.row(i_m)
            a_m = A.row(i_m)
            a_tilde_m = a_m.pairwise_product(beta)
            phi[i_m] = (a_tilde_m.row()*P).norm()**2- \
                ((h_m.row()*P*P.T*a_tilde_m.column()).norm()**2)/(1+(h_m*P).norm()**2)
        for i_l in range(0, L):
            phi_max = 0
            for i_m in range(0, M):
                if A[i_m, i_l] != 0:
                    phi_max = max(phi[i_m], phi_max)
            r[i_l] = 0.5*log(max(1, beta[i_l]**2*(P[i_l,i_l]**2/phi_max)), 2)
    except:
        print 'error in rate_computation_MMSE_alpha'
        raise

    return (r, phi)
    # phi is exactly the fine lattices at the relays


def sum_rate_computation_MMSE_alpha(L, M, P_t, H, A, beta):
    r, relay_fine_lattices = rate_computation_MMSE_alpha(L, M, P_t, H, A, beta)
    return sum(r)


# Simple structure, for containing simulation result
class CoF_Sim_Result:
    def __init__(self, sum_rate, sum_rate_var):
        self.sum_rate = sum_rate
        self.sum_rate_var = sum_rate_var
        
class CoF_Dual_Hops_Sim_Result:
    def __init__(self, sum_rate_fixed_pow_sym_mod, sum_rate_sym_mod, sum_rate_asym_mod):
        self.sum_rate_fixed_pow_sym_mod = sum_rate_fixed_pow_sym_mod
        self.sum_rate_sym_mod = sum_rate_sym_mod
        self.sum_rate_asym_mod = sum_rate_asym_mod

class CoF_Sim_Result_for_Fixed_H_Fixed_P:
    def __init__(self, H_a, sum_rate_i_H, A, alpha_opt, P_vec):
        self.H_a = H_a
        self.sum_rate_i_H = sum_rate_i_H
        self.A = A
        self.alpha_opt = alpha_opt
        self.P_vec = P_vec

    def __repr__(self):
        return '---------------Fixed_H_Fixed_P---------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H = '+self.sum_rate_i_H.__str__()+'\n' \
            +'A = \n'+self.A.__str__()+'\n' \
            +'alpha_opt = '+self.alpha_opt.__str__()+'\n' \
            +'P_vec = '+self.P_vec.__str__()+'\n' \
            +'--------------------------------------------\n'
        
class CoF_Sim_Result_for_Fixed_H_Variable_P:
    def __init__(self, H_a, sum_rate_i_H_var, A, \
                 alpha_opt, P_vec_best_var):
        self.H_a = H_a
        self.sum_rate_i_H_var = sum_rate_i_H_var
        self.A = A
        self.alpha_opt = alpha_opt
        self.P_vec_best_var = P_vec_best_var

    def __repr__(self):
        return '-------------Fixed_H_Variable_P--------------\n' \
            +'H_a = \n'+self.H_a.__str__()+'\n' \
            +'sum_rate_i_H_var = '+self.sum_rate_i_H_var.__str__()+'\n' \
            +'A = \n'+ self.A.__str__()+'\n' \
            +'alpha_opt = '+self.alpha_opt.__str__()+'\n' \
            +'P_vec_best_var = '+self.P_vec_best_var.__str__()+'\n' \
            +'--------------------------------------------\n'