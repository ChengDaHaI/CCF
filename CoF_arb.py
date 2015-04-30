'''
Simulation of Compute and Forward. In alpha test. 
Author: Yihua Tan
Email: adherentx.tyh@gmail.com
The Chinese University of Hong Kong
'''

import copy
import sys
from sage.all import *
from numpy import arange
import time
from CoF_LLL import *
from CoF_basic import *
from CoF_second_hop import *
from alternate_opt import *
from CoF_GA import *
from scipy import optimize
import math
import pickle

#from ComputeRate import *
from random import *
#import gc
#sys.path.insert(1, '/usr/local/lib/python2.7/dist-packages/')
#from guppy import hpy
#my_hpy = hpy()

print 'Hello, this is the simulation of CoF.'

# This is for test.


def CoF_compute_search_pow_flex(P_con, H_a, is_dual_hop, rate_sec_hop=[], mod_scheme='sym_mod', quan_scheme='sym_quan', beta=[]):
    (M, L) = (H_a.nrows(), H_a.ncols())
    global P_Search_Alg
    if beta == []:
        beta = vector(RR, [1]*L)
    cof_pow = lambda x: -CoF_compute_fixed_pow_flex(x, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
    cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex(x[0:L], P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, x[L:L+M]))
    #Pranges = ((P_con/brute_number, P_con), )*L # (slice(0, P_con+0.1, P_con/brute_number), )*L
    Pranges = ((0, P_con), )*L
    initial_guess = [0.5*P_con]*L
    try:
        if P_Search_Alg == 'brute':
            res_cof = optimize.brute(cof_pow, Pranges, Ns=brute_number, full_output=True, finish=None)
            P_opt = res_cof[0]
            sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow
        elif P_Search_Alg == 'TNC':
            #res_cof = optimize.minimize(cof_pow, initial_guess, method='TNC', bounds=Pranges, options={'maxiter': 400, 'approx_grad': True})
            #P_opt = list(res_cof.x)
            #sum_rate_opt = -res_cof.fun # negative! see minus sign in cof_pow
            res_cof = optimize.fmin_tnc(cof_pow, initial_guess, bounds=list(Pranges), approx_grad=True, epsilon=1, stepmx=10)
            P_opt = res_cof[0]
            sum_rate_opt = CoF_compute_fixed_pow_flex(P_opt, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
        elif P_Search_Alg == 'anneal':
            res_cof = optimize.anneal(cof_pow, initial_guess, schedule='cauchy', T0=1, Tf=1e-6, \
                      full_output=True, maxiter=30, lower=[1, 1], upper=[P_con, P_con], dwell=30, disp=True)
            P_opt = list(res_cof[0])
            sum_rate_opt = -res_cof[1]
        elif P_Search_Alg == 'brute_fmin':
            res_brute = optimize.brute(cof_pow, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow
            res_fmin = optimize.fmin(cof_pow, P_brute_opt, xtol=1, ftol=0.01, maxiter=brute_fmin_maxiter, full_output=True)
            #P_fmin_opt = res_fmin[0]
            P_opt = res_fmin[0]
            sum_rate_opt = -res_fmin[1]
        elif P_Search_Alg == 'brute_brute':
            res_brute1 = optimize.brute(cof_pow, Pranges, Ns=brute_brute_first_number, full_output=True, finish=None)
            P_brute_opt1 = res_brute1[0]
            sum_rate_brute1 = -res_brute1[1] # negative! see minus sign in cof_pow
            Pranges_brute_2 = tuple([(max(0,P_i-P_con/brute_brute_first_number), min(P_con,P_i+P_con/brute_brute_first_number)) for P_i in P_brute_opt1])
            res_brute2 = optimize.brute(cof_pow, Pranges_brute_2, Ns=brute_brute_second_number, full_output=True, finish=None)
            P_brute_opt2 = res_brute2[0]
            sum_rate_brute2 = -res_brute2[1] # negative! see minus sign in cof_pow
            sum_rate_opt = sum_rate_brute2
        elif P_Search_Alg == 'brute_fmin_beta':
            res_brute = optimize.brute(cof_pow, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow
            res_fmin_beta = optimize.fmin(cof_pow_beta, list(P_brute_opt)+[1]*M, xtol=0.01, ftol=0.01, maxiter=brute_fmin_maxiter*50, full_output=True)
            P_fmin_opt = res_fmin_beta[0]
            sum_rate_opt = -res_fmin_beta[1]
        elif P_Search_Alg == 'brute_fmin_cobyla':
            res_brute = optimize.brute(cof_pow, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            def pow_constraint(x):
                return x
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow
            p_cobyla = optimize.fmin_cobyla(cof_pow, P_brute_opt, pow_constraint, maxfun=100)
            sum_rate_fmin_cobyla = CoF_compute_fixed_pow_flex(p_cobyla, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
            sum_rate_opt = sum_rate_fmin_cobyla
        elif P_Search_Alg == 'brute_fmin_cobyla_beta':
            res_brute = optimize.brute(cof_pow, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            def pow_beta_constraint(x):
                return x
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow
            p_cobyla = optimize.fmin_cobyla(cof_pow_beta, list(P_brute_opt)+[1]*M, pow_beta_constraint, maxfun=200)
            sum_rate_fmin_cobyla = CoF_compute_fixed_pow_flex(p_cobyla, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
            sum_rate_opt = sum_rate_fmin_cobyla
        #Add differential evolution
        elif P_Search_Alg =="differential_evolution":
            bounds=((0.1,P_con),)*L
            res_brute=optimize.differential_evolution(cof_pow,bounds)
            P_opt=res_brute.x
            sum_rate_opt=-res_brute.fun
        #Add Genetic Algorithm
        elif P_Search_Alg=="genetic":
            res_cof=GeneticAlgorithm(P_con,H_a)
            P_opt=res_cof[0]
            sum_rate_opt=res_cof[1]
        #The Genetic Algorithm End
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    return sum_rate_opt
    #return (sum_rate_opt,P_opt)#return tuple,for comparing

def CoF_compute_search_pow_flex_beta(P_con, H_a, is_fixed_power, is_dual_hop, rate_sec_hop=[], mod_scheme='sym_mod', quan_scheme='sym_quan'):
    (M, L) = (H_a.nrows(), H_a.ncols())
    '''
    def cof_pow_beta(x):
        power = x[0:L]
        beta = vector(RR, [1,]+list(x[L:2*L-1]))
        -CoF_compute_fixed_pow_flex(power, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
    '''
    if is_fixed_power == False:
        cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex(x[0:L], P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, [1,]+list(x[L:2*L-1])))
        Pranges = ((P_con/brute_number, P_con), )*L + ((float(beta_max)/brute_number, beta_max), )
    else:
        cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex((P_con,P_con), P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, [1,x]))
        Pranges = ((float(beta_max)/brute_number, beta_max), )
        
    try:
        if P_Search_Alg == 'brute':
            res_cof = optimize.brute(cof_pow_beta, Pranges, Ns=brute_number, full_output=True, finish=None)
            P_opt = res_cof[0]
            sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow_beta
        elif P_Search_Alg == 'brute_fmin':
            res_brute = optimize.brute(cof_pow_beta, Pranges, Ns=brute_fmin_number, full_output=True, finish=None)
            P_brute_opt = res_brute[0]
            sum_rate_brute = -res_brute[1] # negative! see minus sign in cof_pow_beta
            res_fmin = optimize.fmin(cof_pow_beta, P_brute_opt, xtol=1, ftol=0.01, maxiter=brute_fmin_maxiter, full_output=True)
            P_fmin_opt = res_fmin[0]
            sum_rate_opt = -res_fmin[1]
        elif P_Search_Alg == 'brute_brute':
            res_brute1 = optimize.brute(cof_pow_beta, Pranges, Ns=brute_brute_first_number, full_output=True, finish=None)
            P_brute_opt1 = res_brute1[0]
            sum_rate_brute1 = -res_brute1[1] # negative! see minus sign in cof_pow_beta
            Pranges_brute_2 = tuple([(max(0,P_i-P_con/brute_brute_first_number), min(P_con,P_i+P_con/brute_brute_first_number)) for P_i in P_brute_opt1])
            res_brute2 = optimize.brute(cof_pow_beta, Pranges_brute_2, Ns=brute_brute_second_number, full_output=True, finish=None)
            P_brute_opt2 = res_brute2[0]
            sum_rate_brute2 = -res_brute2[1] # negative! see minus sign in cof_pow_beta
            sum_rate_opt = sum_rate_brute2
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    return sum_rate_opt
    

@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_first_hop(P_con, M, L):
    sum_rate = 0
    sum_rate_var = 0
    sum_rate_var_beta = 0
    sum_rate_beta = 0
    for i_H in range(0, batch_H):
        if is_set_H == True:
            H_a = set_H_a
        else:
            set_random_seed() # to avoid producing the same H_a in different threads
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
            
        sum_rate_i_H = 0
        sum_rate_i_H_var = 0
        sum_rate_i_H_var_beta = 0
        
        # Fixed power
        sum_rate_i_H = CoF_compute_fixed_pow_flex((P_con, )*L, P_con, False, H_a, False)
        
        # Variable power
        if is_alternate == True:
            sum_rate_i_H_var = alternate_optimize(P_con, H_a, is_dual_hop=False)
        else:
            sum_rate_i_H_var = CoF_compute_search_pow_flex(P_con, H_a, is_dual_hop=False)
            
        sum_rate_i_H_beta = CoF_compute_search_pow_flex_beta(P_con, H_a, is_fixed_power=True, is_dual_hop=False)
        sum_rate_i_H_var_beta = CoF_compute_search_pow_flex_beta(P_con, H_a, is_fixed_power=False, is_dual_hop=False)
        
        if sum_rate_i_H_var_beta > 1.01*sum_rate_i_H_beta:
            print 'var_beta is better than sole beta when channel H = ', H_a
            print 'sum_rate_i_H_beta = ', sum_rate_i_H_beta
            print 'sum_rate_i_H_var_beta', sum_rate_i_H_var_beta
        
        sum_rate += sum_rate_i_H
        sum_rate_var += sum_rate_i_H_var
        sum_rate_beta += sum_rate_i_H_beta
        sum_rate_var_beta += sum_rate_i_H_var_beta
    # for i_H
    sum_rate /= batch_H
    sum_rate_var /= batch_H
    sum_rate_beta /= batch_H
    sum_rate_var_beta /= batch_H
    
#     print P_con, sum_rate, sum_rate_var
    return {'sum_rate': sum_rate, 'sum_rate_var': sum_rate_var, 'sum_rate_var_beta': sum_rate_var_beta, 'sum_rate_beta': sum_rate_beta}

'''
The Following Are Dual Hop Functions
'''
@parallel(ncpus=Cores)
def CoF_compute_eq_pow_con_dual_hops(P_con, M, L):
    sum_rate_fixed_pow_sym_mod = 0
    sum_rate_fixed_pow_asym_quan = 0
    sum_rate_sym_mod = 0
    sum_rate_asym_mod = 0
    sum_rate_asym_quan = 0
    sum_rate_asym_mod_asym_quan = 0 
    '''How to determine the power of the relays?'''
    P_relay = 0.25*P_con
    #P_relay = P_con
    #P_relay = sqrt(P_con) 
    for i_H in range(0, batch_H):
        if is_set_H == True:
            H_a = set_H_a
            H_b = set_H_b
        else:
            set_random_seed() # to avoid producing the same H_a in different threads
            H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
            set_random_seed()
            H_b = (matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))).column(0)
        if is_set_beta == True:
            beta = set_beta
        else:
            beta = vector(RR, [1]*L)
        rate_sec_hop = [0]*M # ? bit/s for each parallel channel in the second hop
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
        
        try:
            (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_fixed_pow_asym_quan, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_quan, sum_rate_i_H_asym_mod_asym_quan) = \
                CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, True, rate_sec_hop, beta)
        except:
            print 'error in searching for good power'
            raise
        
        sum_rate_fixed_pow_sym_mod += sum_rate_i_H_fixed_pow_sym_mod
        sum_rate_fixed_pow_asym_quan += sum_rate_i_H_fixed_pow_asym_quan
        sum_rate_sym_mod += sum_rate_i_H_sym_mod
        sum_rate_asym_mod += sum_rate_i_H_asym_mod
        sum_rate_asym_quan += sum_rate_i_H_asym_quan
        sum_rate_asym_mod_asym_quan += sum_rate_i_H_asym_mod_asym_quan
    # for i_H
    sum_rate_fixed_pow_sym_mod /= batch_H
    sum_rate_fixed_pow_asym_quan /= batch_H
    sum_rate_sym_mod /= batch_H
    sum_rate_asym_mod /= batch_H
    sum_rate_asym_quan /= batch_H
    sum_rate_asym_mod_asym_quan /= batch_H
    
    return {'sum_rate_fixed_pow_sym_mod': sum_rate_fixed_pow_sym_mod,
            'sum_rate_fixed_pow_asym_quan': sum_rate_fixed_pow_asym_quan, 
            'sum_rate_sym_mod': sum_rate_sym_mod, 
            'sum_rate_asym_mod': sum_rate_asym_mod, 
            'sum_rate_asym_quan': sum_rate_asym_quan, 
            'sum_rate_asym_mod_asym_quan': sum_rate_asym_mod_asym_quan}


@parallel(ncpus=Cores)
def CoF_beta_search_dual_hops(P_con, M, L, beta):
    H_a = set_H_a
    H_b = set_H_b
    rate_sec_hop = [0]*M # ? bit/s for each parallel channel in the second hop
    P_relay = 0.25*P_con
    for i_h_b in range(0, M):
        rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
    (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_fixed_pow_asym_quan, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_quan, sum_rate_i_H_asym_mod_asym_quan) = \
        CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, True, rate_sec_hop, beta)
    
    return {'sum_rate_i_H_fixed_pow_sym_mod': sum_rate_i_H_fixed_pow_sym_mod, 
            'sum_rate_i_H_fixed_pow_asym_quan': sum_rate_i_H_fixed_pow_asym_quan, 
            'sum_rate_i_H_sym_mod': sum_rate_i_H_sym_mod, 
            'sum_rate_i_H_asym_mod': sum_rate_i_H_asym_mod, 
            'sum_rate_i_H_asym_quan': sum_rate_i_H_asym_quan, 
            'sum_rate_i_H_asym_mod_asym_quan': sum_rate_i_H_asym_mod_asym_quan}
    

def CoF_compute_eq_pow_con_dual_hops_fixed_H(is_alternate, P_con, H_a, is_dual_hop, rate_sec_hop, beta):
    # Fixed power
    sum_rate_i_H_fixed_pow_sym_mod = CoF_compute_fixed_pow_flex([P_con]*L, P_con, False, H_a, is_dual_hop, rate_sec_hop, 'sym_mod', 'sym_quan', beta)
    sum_rate_i_H_fixed_pow_asym_quan = CoF_compute_fixed_pow_flex([P_con]*L, P_con, False, H_a, is_dual_hop, rate_sec_hop, 'sym_mod', 'asym_quan', beta)
    
    # Variable power
    sum_rate_i_H_sym_mod = 0
    sum_rate_i_H_asym_mod = 0
    sum_rate_i_H_asym_quan = 0
    sum_rate_i_H_asym_mod_asym_quan = 0
    
    if is_alternate == True:
        sum_rate_i_H_sym_mod = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_quan = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'asym_quan', beta)
        sum_rate_i_H_asym_mod_asym_quan = alternate_optimize(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'asym_quan', beta)
    else:
        sum_rate_i_H_sym_mod = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_mod = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'sym_quan', beta)
        sum_rate_i_H_asym_quan = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'sym_mod', 'asym_quan', beta)
        sum_rate_i_H_asym_mod_asym_quan = CoF_compute_search_pow_flex(P_con, H_a, True, rate_sec_hop, 'asym_mod', 'asym_quan', beta)
    
    
    if (sum_rate_i_H_fixed_pow_sym_mod<0) or (sum_rate_i_H_fixed_pow_asym_quan<0) or (sum_rate_i_H_sym_mod<0) or (sum_rate_i_H_asym_quan<0) or (sum_rate_i_H_asym_mod<0) or (sum_rate_i_H_asym_mod_asym_quan<0):
        raise Exception('Function CoF_compute_eq_pow_con_dual_hops_fixed_H() gets negative result!')
        # return (0, 0, 0, 0)
    
    return (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_fixed_pow_asym_quan, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_quan, sum_rate_i_H_asym_mod_asym_quan)



if __name__ == "__main__": 
    '''Equal Power Constraint'''
    P_eq_dB_Min = float(P_MIN)
    P_eq_dB_Max = float(P_MAX)
    P_delta = P_DEL
    P_eq_dB = arange(P_eq_dB_Min, P_eq_dB_Max, P_delta)
    P_eq = [10**(P_eq_dB_i/10) for P_eq_dB_i in P_eq_dB]
    Pl_con = P_eq
    
    num_batch = iter_H/batch_H # number of batch
    # assert mod(num_batch, Cores)==0, 'For the sake of speed, the number of batches should be a multiple of Cores.'
    
    t1 = time.ctime()
    tf1 = time.time()
    if is_alternate == True:
        dir_name = '/home/chenghai/Research/CoF_Sim/'+str(tf1)+'-'+'alternate'+'-M=L='+str(M)+'-iter_H='+str(iter_H)
    else:
        dir_name = '/home/chenghai/Research/CoF_Sim/'+str(tf1)+'-'+P_Search_Alg+'-M=L='+str(M)+'-iter_H='+str(iter_H)
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    original_dir = os.getcwd()
    os.chdir(dir_name) # change to the directory where simulation results should be placed
    
    '''First Hop'''
    if False:
        sum_rate = [0]*len(Pl_con)
        sum_rate_var = [0]*len(Pl_con)
        sum_rate_beta = [0]*len(Pl_con)
        sum_rate_var_beta = [0]*len(Pl_con)
        for i_P in range(0, len(Pl_con)):
            print 'First Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_first_hop([(Pl_con[i_P], M, L)]*num_batch))
            sum_rate_bat = [result_bat[i][1]['sum_rate'] for i in range(0, num_batch)]
            sum_rate_var_bat = [result_bat[i][1]['sum_rate_var'] for i in range(0, num_batch)]
            sum_rate_beta_bat = [result_bat[i][1]['sum_rate_beta'] for i in range(0, num_batch)]
            sum_rate_var_beta_bat = [result_bat[i][1]['sum_rate_var_beta'] for i in range(0, num_batch)]
            sum_rate[i_P] = sum(sum_rate_bat)/num_batch
            sum_rate_var[i_P] = sum(sum_rate_var_bat)/num_batch
            sum_rate_beta[i_P] = sum(sum_rate_beta_bat)/num_batch
            sum_rate_var_beta[i_P] = sum(sum_rate_var_beta_bat)/num_batch
            
        sum_rate = [RR(sum_rate[i]) for i in range(0, len(Pl_con))]
        sum_rate_var = [RR(sum_rate_var[i]) for i in range(0, len(Pl_con))]
        sum_rate_beta = [RR(sum_rate_beta[i]) for i in range(0, len(Pl_con))]
        sum_rate_var_beta = [RR(sum_rate_var_beta[i]) for i in range(0, len(Pl_con))]
        
        t2 = time.ctime()
        
        print 'Simulation of the first hop started at ', t1
        print 'Simulation of the first hop ended at ', t2
                
        plot_sum_rate = list_plot(zip(P_eq_dB, sum_rate), plotjoined=True, marker='o', \
                                  rgbcolor=Color('black'), linestyle="--", \
                                  legend_label= 'Sum rate with fixed power', \
                                  title = 'Comparison in the First Hop')
        plot_sum_rate_var = list_plot(zip(P_eq_dB, sum_rate_var), plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'Sum rate with variable power')
        plot_sum_rate_beta = list_plot(zip(P_eq_dB, sum_rate_beta), plotjoined=True, marker='D', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'Sum rate with variable beta')
        plot_sum_rate_var_beta = list_plot(zip(P_eq_dB, sum_rate_var_beta), plotjoined=True, marker='x', \
                                      rgbcolor=Color('red'), linestyle='-', \
                                      legend_label = 'Sum rate with variable power and beta')
        plot_compare = plot_sum_rate+plot_sum_rate_var+plot_sum_rate_beta+plot_sum_rate_var_beta
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.set_legend_options(loc='upper left')
        plot_compare.save('Comparison_in_the_First_Hop-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'-M=L='+str(M)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, {'sum_rate': sum_rate, 
                               'sum_rate_var': sum_rate_var,
                               'sum_rate_beta': sum_rate_beta,
                               'sum_rate_var_beta': sum_rate_var_beta}), 
                    open('First_Hop.pkl', 'w'))
        print 'sum_rate: '; print sum_rate
        print 'sum_rate_var: '; print sum_rate_var
        print 'sum_rate_beta: '; print sum_rate_beta
        print 'sum_rate_var_beta: '; print sum_rate_var_beta
    
    '''Dual Hops'''
    if True:
        sum_rate_fixed_pow_sym_mod = [0]*len(Pl_con)
        sum_rate_fixed_pow_asym_quan = [0]*len(Pl_con)
        sum_rate_sym_mod = [0]*len(Pl_con)
        sum_rate_asym_mod = [0]*len(Pl_con)
        sum_rate_asym_quan = [0]*len(Pl_con)
        sum_rate_asym_mod_asym_quan = [0]*len(Pl_con)
        t3 = time.ctime()
        # debug:
        #CoF_compute_eq_pow_con_dual_hops(100)
        
        for i_P in range(0, len(Pl_con)):
            print 'Dual Hop Simulation: P_con=', Pl_con[i_P]
            result_bat = list(CoF_compute_eq_pow_con_dual_hops([(Pl_con[i_P], M, L)]*num_batch))
            
            sum_rate_fixed_pow_sym_mod_bat = [result_bat[i][1]['sum_rate_fixed_pow_sym_mod'] for i in range(0, num_batch)]
            sum_rate_fixed_pow_asym_quan_bat = [result_bat[i][1]['sum_rate_fixed_pow_asym_quan'] for i in range(0, num_batch)]
            sum_rate_sym_mod_bat = [result_bat[i][1]['sum_rate_sym_mod'] for i in range(0, num_batch)]
            sum_rate_asym_mod_bat = [result_bat[i][1]['sum_rate_asym_mod'] for i in range(0, num_batch)]
            sum_rate_asym_quan_bat = [result_bat[i][1]['sum_rate_asym_quan'] for i in range(0, num_batch)]
            sum_rate_asym_mod_asym_quan_bat = [result_bat[i][1]['sum_rate_asym_mod_asym_quan'] for i in range(0, num_batch)]
            
            sum_rate_fixed_pow_sym_mod[i_P] = sum(sum_rate_fixed_pow_sym_mod_bat)/num_batch
            sum_rate_fixed_pow_asym_quan[i_P] = sum(sum_rate_fixed_pow_asym_quan_bat)/num_batch
            sum_rate_sym_mod[i_P] = sum(sum_rate_sym_mod_bat)/num_batch
            sum_rate_asym_mod[i_P] = sum(sum_rate_asym_mod_bat)/num_batch
            sum_rate_asym_quan[i_P] = sum(sum_rate_asym_quan_bat)/num_batch
            sum_rate_asym_mod_asym_quan[i_P] = sum(sum_rate_asym_mod_asym_quan_bat)/num_batch
            
        sum_rate_fixed_pow_sym_mod = [RR(sum_rate_fixed_pow_sym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_fixed_pow_asym_quan = [RR(sum_rate_fixed_pow_asym_quan[i]) for i in range(0, len(Pl_con))]
        sum_rate_sym_mod = [RR(sum_rate_sym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_asym_mod = [RR(sum_rate_asym_mod[i]) for i in range(0, len(Pl_con))]
        sum_rate_asym_quan = [RR(sum_rate_asym_quan[i]) for i in range(0, len(Pl_con))]
        sum_rate_asym_mod_asym_quan = [RR(sum_rate_asym_mod_asym_quan[i]) for i in range(0, len(Pl_con))]
        t4 = time.ctime()
        
        print 'Simulation of dual hops started at ', t3
        print 'Simulation of dual hops ended at ', t4
        
        plot_sum_rate_fixed_pow_sym_mod = list_plot(zip(P_eq_dB, sum_rate_fixed_pow_sym_mod), plotjoined=True, marker='D', \
                                  rgbcolor=Color('green'), linestyle="-.", \
                                  legend_label= 'Fixed power(Original CoF)', \
                                  title = 'Comparison in Two-Hop System')
        plot_sum_rate_fixed_pow_asym_quan = list_plot(zip(P_eq_dB, sum_rate_fixed_pow_asym_quan), plotjoined=True, marker='d', \
                                  rgbcolor=Color('tomato'), linestyle="-.", \
                                  legend_label= 'Fixed power with asymmetric quantization')
        plot_sum_rate_sym_mod = list_plot(zip(P_eq_dB, sum_rate_sym_mod), plotjoined=True, marker='H', \
                                  rgbcolor=Color('black'), linestyle=":", \
                                  legend_label= 'Symmetric modulo approach with variable power')
        plot_sum_rate_asym_mod = list_plot(zip(P_eq_dB, sum_rate_asym_mod), plotjoined=True, marker='<', \
                                      rgbcolor=Color('blue'), linestyle='-', \
                                      legend_label = 'Asymmetric modulo approach with variable power')
        plot_sum_rate_asym_quan = list_plot(zip(P_eq_dB, sum_rate_asym_quan), plotjoined=True, marker='>', \
                                      rgbcolor=Color('seagreen'), linestyle='-', \
                                      legend_label = 'Asymmetric quantization approach with variable power')
        plot_sum_rate_asym_mod_asym_quan = list_plot(zip(P_eq_dB, sum_rate_asym_mod_asym_quan), plotjoined=True, marker='o', \
                                        rgbcolor=Color('red'), linestyle='--', \
                                        legend_label = 'Asymmetric modulo and asymmetric quantization approach with variable power')
        
        plot_compare = plot_sum_rate_fixed_pow_sym_mod+plot_sum_rate_fixed_pow_asym_quan+plot_sum_rate_sym_mod+\
            plot_sum_rate_asym_mod+plot_sum_rate_asym_quan+plot_sum_rate_asym_mod_asym_quan
        plot_compare.set_axes_range(ymax=max(sum_rate_fixed_pow_sym_mod+sum_rate_fixed_pow_asym_quan+sum_rate_sym_mod+\
                                             sum_rate_asym_mod+sum_rate_asym_quan+sum_rate_asym_mod_asym_quan)*1.4)
        plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
        plot_compare.set_legend_options(loc='upper left')
        plot_compare.save('Comparison_in_Dual_Hops_System-' \
                          +P_Search_Alg+'-is_alternate='+str(is_alternate)+'-M=L='+str(M)+'.eps')
        show(plot_compare)
        pickle.dump((P_eq_dB, {'sum_rate_fixed_pow_sym_mod': sum_rate_fixed_pow_sym_mod, 
                               'sum_rate_fixed_pow_asym_quan': sum_rate_fixed_pow_asym_quan, 
                               'sum_rate_sym_mod': sum_rate_sym_mod, 
                               'sum_rate_asym_mod': sum_rate_asym_mod, 
                               'sum_rate_asym_quan': sum_rate_asym_quan, 
                               'sum_rate_asym_mod_asym_quan': sum_rate_asym_mod_asym_quan}), 
                    open('Dual_Hops.pkl', 'w'))
        
        print 'fixed power: '; print sum_rate_fixed_pow_sym_mod
        print 'fixed power with asymmetric quantization: '; print sum_rate_fixed_pow_asym_quan
        print 'symmetric mod with variable power: '; print sum_rate_sym_mod
        print 'asymmetric mod with variable power: '; print sum_rate_asym_mod
        print 'asymmetric quan with variable power: '; print sum_rate_asym_quan
        print 'asym mod and asym quan with variable power: '; print sum_rate_asym_mod_asym_quan
        


    '''Beta optimization'''
    if False:
        beta_min = 1
        beta_max = 2
        beta_delta = 1
        beta_number = int(round((beta_max-beta_min)/beta_delta)+1)
        #assert beta_number%Cores==0, 'wrong beta_number setting!'
        P_con = 10**(20/10) # 40 dB
        result = []
        if L==M==2:
            result = list(CoF_beta_search_dual_hops([(P_con, M, L, vector(RR, [beta_min+j0*beta_delta, beta_min+j1*beta_delta])) \
                                                     for j0 in range(0,beta_number) for j1 in range(0,beta_number)]))
        else:
            pass
        # (sum_rate_i_H_fixed_pow_sym_mod, sum_rate_i_H_sym_mod, sum_rate_i_H_asym_mod, sum_rate_i_H_asym_mod_asym_quan)
        for re in result:
            print re
        pickle.dump(result, open('beta_search.pkl', 'w'))
        best_result_0 = max(result, key=lambda x: x[1]['sum_rate_i_H_asym_mod_asym_quan'])
        print 'The best result for asym_mod_asym_quan is: \n', best_result_0
        best_result_1 = max(result, key=lambda x: x[1]['sum_rate_i_H_asym_mod'])
        print 'The best result for asym_mod is: \n', best_result_1
        best_result_2 = max(result, key=lambda x: x[1]['sum_rate_i_H_sym_mod'])
        print 'The best result for sym_mod is: \n', best_result_2
        best_result_3 = max(result, key=lambda x: x[1]['sum_rate_i_H_fixed_pow_sym_mod'])
        print 'The best result for fixed_pow_sym_mod is: \n', best_result_3
        pickle.dump({'sum_rate_i_H_asym_mod_asym_quan': best_result_0, 
                     'sum_rate_i_H_asym_mod': best_result_1, 
                     'sum_rate_i_H_sym_mod': best_result_2, 
                     'sum_rate_i_H_fixed_pow_sym_mod': best_result_3},
                    open('best_beta.pkl', 'w'))
        result_all_asym_mod_asym_quan = [tuple(x[0][0][3])+(x[1]['sum_rate_i_H_asym_mod_asym_quan'],) for x in result]
        print result_all_asym_mod_asym_quan
        result_beta_plot = list_plot3d(result_all_asym_mod_asym_quan)
        sum_rate_max = max([x[1]['sum_rate_i_H_asym_mod_asym_quan'] for x in result])
        sum_rate_min = min([x[1]['sum_rate_i_H_asym_mod_asym_quan'] for x in result])
        beta_label_plot = text3d('beta1', (beta_max,beta_min,sum_rate_min), color='red') + \
            text3d('beta2', (beta_min,beta_max,sum_rate_min), color='red') + \
            text3d('sum_rate', (beta_min,beta_min,sum_rate_max), color='red')
        
        show(result_beta_plot+beta_label_plot)
    # P_con, M, L, beta
    
    os.chdir(original_dir) # recover directory
    raw_input() # stop Sage from shutting down
