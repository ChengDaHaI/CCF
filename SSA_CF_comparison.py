'''
This file aims to generate simulation result of CCF in MIMO C-RAN system.
We will compare CCF with SSA PNC scheme.
Date :2017-5-4
'''

from sage.all import *
from NEW_basic import *
from CoF_LLL import Find_A_and_Rate
from scipy import optimize
from math import log, log10, fabs
import time
import copy
import matplotlib.pyplot as plt
import numpy as np
from SSA_fading_channel_model import chanMatrix


# Compute the sum rate of CCF in a MIMO C-RAN system
# A user with N antenna is viewed as N independent users
# input: P_tran -- transmiter power list, H_a -- first hop channel,
#        capacity_bh -- capcaity of backhual, N -- number of antenna in transmiter
# output:sum rate of a C-RAN system
def Rate_CCF_CRAN( P_tran,beta, H_a, capacity_bh):
    # transmitter power
    P_vec = vector(RR, P_tran)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, source_rate_upbound_list, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, True,beta)
    except:
        print 'error in seeking A and rate'
        raise
    # assuming SW coding is used and thus the forwarding rate is the sum of source message rates
    return max(min(sum(source_rate_upbound_list), capacity_bh), 0)


# optimized sum rate of a C-RAN system
# input: P_con -- power budget of transmiter, C_BH -- coefficient of capacity_bh, alg -- used algorithm
# output:optimized sum rate
@parallel(ncpus=Cores)
def rate_opt(P_con,  C_BH,  alg = 'None'):

    set_random_seed()
    # generate channel matrix with M*N rows and L*N columns
    # H = Matrix(RR, M, L,  lambda i, j: normalvariate(0, 1))
    # large scale fading channel
    H = chanMatrix(3, 3, 2, 2) # 2 * 2 system with 3 antenna
    H = Matrix(H )
    # BH capcaity
    # capacity_bh = C_BH * log(P_con,2)
    capacity_bh = C_BH * K
    sum_rate = 0
    if alg == 'DE':
        opt_fun = lambda x: - Rate_CCF_CRAN([P_con]*L, vector(RR, [1, ] + list(x[0:L - 1])), H, capacity_bh)
        bounds = ((0.1, betaScale_max), ) * (L-1)
        opt_res = optimize.differential_evolution(opt_fun, bounds, maxiter = 20, disp = False)
        sum_rate = -opt_res.fun
        print 'CCF Differential Evolution:', opt_res.success
    elif alg == 'None':
        sum_rate = Rate_CCF_CRAN([P_con]*L, H, capacity_bh)
        print 'CCF None Optimization!'

    return sum_rate

if __name__ == '__main__':

    # C_BH = 3
    alg = 'DE'
    num_batch = 500
    sum_rate_list = []
    New_sum_rate = []
    # PI_con = [10 ** 1, 10 ** 1.25, 10 ** 1.5, 10 ** 1.75, 10 ** 2, 10 ** 2.25, 10 ** 2.5, 10 ** 2.75, 10 ** 3, 10 ** 3.25,
    #             10 ** 3.5]
    # PI_con = [ 10**3, 10**3.5,  1e4, 10**4.5, 10**5, 10**5.5, 10**6, 10**6.5,  1e7, 10**7.5, 10**8, 10**8.5, 10**9]
    PI_con = [  1e7]
    # Pow = 10 ** 2
    # BH_list = [2, 4, 6, 8, 10, 12, 14, 16] # BH capacity per user
    print 'Simulation Starts!\n'
    t1 = time.time()
    for Pow in PI_con:
        C_BH = 100 # set to be very large
        print 'Transmitter Power:', Pow
        # the power of each antenna user is Pi/N
        # test_rate = rate_opt(Pi/N, C_BH, alg)
        # print 'Test'
        # Rate_list = []
        # for i in range(num_batch):
        #     Rate_list.append(rate_opt(Pi, C_BH, alg))
        result_list = list(rate_opt([(Pow/N, C_BH,  alg )]* num_batch))
        Rate_list = [result_list[i][1] for i in range(0, num_batch)]
        sum_rate = sum(Rate_list)/num_batch
        sum_rate_list.append(sum_rate)
    t2 = time.time()
    print 'Total Time Cost: ', (t2 - t1)
    PI_dB = [10 * log10(P_con) for P_con in PI_con]
    Full_Result = np.column_stack((PI_dB, sum_rate_list))
    print 'sum_rate_list:', sum_rate_list
    np.savetxt('/home/haizi/PycharmProjects/SSA/Simu_result/' + 'CCF_OPT' +' iter =' + num_batch.__str__() + ' K=' + K.__str__() + ' N=' + N.__str__() + time.ctime()
               + 'Simu_Data.txt', Full_Result, fmt='%1.5e')
