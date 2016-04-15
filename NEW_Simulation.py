'''
This file is to simulate our new CCF modle 
and get a convinced result through 
Monte-Carlo simulation.
'''
from sage.all import *
from NEW_Optimize_Modle import RandomSearch, CCF_new_sumrate_func
from NEW_basic import *
from NEW_CCF_Modle import Relay_Forward_Rate, powerset
from NewSecondHopChannel import ComputeSecRate
from ComputeRate import CoF_compute_search_pow_flex_beta
from CoF_LLL import Find_A_and_Rate
from CoF_second_hop import second_hop_support_rates
from math import log10, fabs
import time
import copy
import itertools
import matplotlib.pyplot as plt
import numpy as np

@parallel(ncpus=Cores)
def CCF_Model_Comparison(P_Search_Alg,P_con,P_relay):
    set_random_seed()
    set_HaHb=False
    if set_HaHb==True:
        H_a=matrix(RR, M, L, [[-0.612642983158090, 0.293476637932474, -0.335563118861373],\
                    [ 0.377058724156508, -0.206528716605357, -0.615919802428871],\
                    [ 0.264570807268208, -0.592605902143905,  0.699237149229798]])
        '''
        H_b=matrix(RR, M, L, [[ -0.844849781483391 , -0.678659125685948  ,-0.484271670880304],\
                [-0.0729932845848398   ,0.609420751701606   ,0.846395865838560],\
                [ 0.0645367208093419  ,-0.205375774175623 , -0.480734935684002]])
        '''
        H_b=matrix(RR, 1, L, [-0.995387738148530, 0.746656886108230, -0.926365502487133])
        print 'H_a:', H_a
        print 'H_b:', H_b
        print 'Transmitter Power:', P_con
    else:   
        H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        #second hop channel is parallel 
        H_b = (matrix.random(RR, 1, M, distribution=RealDistribution('gaussian', 1)))
        print 'H_a:', H_a
        print 'H_b:', H_b
        print 'Transmitter Power:', P_con
    #second hop channel capacity, 2**L-1 inequalities
    rate_sec_hop=ComputeSecRate(M,P_relay,H_b)
    Max_New_sum_rate=0
    t1=time.time()
    
    sum_rate_opt, beta_pow_opt = CoF_compute_search_pow_flex_beta(P_con,H_a,True, True, P_Search_Alg,rate_sec_hop[0:M],'asym_mod','asym_quan')
    
    if True:
        # check the feasibility of beta_pow_opt
        try:
            is_feasible = Opt_feasible_check(beta_pow_opt, sum_rate_opt, P_con,  H_a, rate_sec_hop)
        except:
            print 'Error in Opt_feasible_check function !'
            raise
        
        if is_feasible:
            
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
            
            # put the optimal solution beta_pow_opt into our NEW CCF system
            try:
                true_beta = vector(RR, [1,] + list(beta_pow_opt))
                LP_res = CCF_new_sumrate_func(true_beta, [P_con]*L, H_a, rate_sec_hop, P_relay, per_c)
                
            except:
                print 'Error In Check The Optimal SOlution to NEW CCF system!'
                raise
            print 'When put beta_pow_opt into OUR NEW CCF system, the sum rate is:', -LP_res
            
        else:
            # if the optimal solution is infeasible in NEW CCF system, we set the sum rate to zero.
            return 0, 0 ,0, 0 #New_sum_rate_opt, sum_rate_opt, (t3-t2) ,(t2-t1)
        
    t2=time.time()
    per_c_search=False
    if per_c_search==True:
        for code_order in itertools.permutations(list(range(0, L)), L):
            per_c=list(code_order)
            (beta_opt, New_sum_rate_opt)=RandomSearch(P_Search_Alg, H_a, rate_sec_hop, P_con, P_relay, per_c)
            if Max_New_sum_rate<New_sum_rate_opt:
                Max_New_sum_rate=New_sum_rate_opt
                New_sum_rate_opt=Max_New_sum_rate
    elif per_c_search==False: 
        '''
        #global per_s, per_c
        (beta_opt, New_sum_rate_opt)=RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c)
        '''
        #compute two permutation after differential evolution operation
        (beta_opt, New_sum_rate_opt)=RandomSearch(P_Search_Alg, H_a, rate_sec_hop, P_con, P_relay)
        
    
    t3=time.time()
    return New_sum_rate_opt, sum_rate_opt, (t3-t2) ,(t2-t1)
    
    
# chech the feasibility of optimal solution of originl CCF in the NEW CCF system
# Input: beta_opt shoulde be L-1 beta factor when with fixed transimitter power
#        rate_sec_hop should be the all constriants
# Output: the feasibility logical value feasible_flag
def Opt_feasible_check(beta_opt, sum_rate_opt, P_con, H_a, rate_sec_hop):
    
    # transmitter power
    P_t = [P_con]*L
    # all beta factor
    beta_opt = vector(RR, [1] + list(beta_opt))
    #compute the shaping lattice and coding lattice of original CCF
    
    P_vec = vector(RR, P_t)
    P_mat = matrix.diagonal([sqrt(x) for x in P_vec])
    # Use LLL to find a good A matrix
    # determine the fine lattice of m-th relay at the same time
    try:
        (A_best_LLL, source_rate_upbound_list, relay_fine_lattices) = Find_A_and_Rate(P_mat, P_vec, H_a, True, beta_opt)
    except:
        print 'error in seeking A and rate'
        raise
    
    A_best_LLL_F = matrix(GF(p), A_best_LLL)
    if A_best_LLL_F.rank() == min(L, M):
        '''constraints of the second hop'''
        # relay_fine_lattices is already obtained
        # compute the coarse lattice of the l-th transmitter
        # The true coarse lattices have scaling factor beta.
        trans_coarse_lattices = list(P_vec.pairwise_product(vector([b**2 for b in beta_opt]))) # copy
        # check whether the second-hop constraint rate_sec_hop can support the first-hop rate r
        try:
            support_result = second_hop_support_rates(relay_fine_lattices, trans_coarse_lattices, A_best_LLL, rate_sec_hop[0:M], 'asym_mod','asym_quan')
        except:
            print 'error in second hop'
            raise
        
        support_rates = support_result[0]
        source_rate   = support_result[1]
        shaping_lattice= support_result[2]
        coding_lattice = support_result[3]
        
        if np.abs(support_rates - sum_rate_opt) > 1* 10**(-3): # prevent the numerical error:
            print 'Something Wrong When recovery the Original CCF sum rate!'
            #raise
    else:
        print 'rank problem in check function'
        #raise
    
    # coding latice nested order
    per_c=[]
    H_a_col_min = copy.copy(coding_lattice)
    for i in range(L):
        H_a_colmin_max=max(H_a_col_min)
        H_a_colmin_max_index=H_a_col_min.index(H_a_colmin_max)
        per_c.append(H_a_colmin_max_index)
        H_a_col_min[H_a_colmin_max_index]=0
    #a larger voronoi corresponds to a coarser coding lattice
    
    #compute the proper nested shaping lattice order from the betaScale
    #------------------
    #nested shaping lattice order is decided by betaScale and source power P_source
    #-------------------
    per_s=[]#nested shaping lattice order
    beta=copy.copy(beta_opt)
    # a larger (beta^2*power) corresponds to a coarse shping lattice
    beta2_power=list([beta[i]**2*P_t[i] for i in range(0,L)])
    beta=copy.copy(beta2_power)
    max_beta=0
    for i in range(L):
        max_beta=max(beta)
        max_beta_index=beta.index(max_beta)
        per_s.append(max_beta_index)
        beta[max_beta_index]=0
    
    
    if max(coding_lattice) >= min(shaping_lattice):
        print 'The lattice nested order in Original CCF is mixed! Not the Same as New CCF'
    
    # compute the all 2*L-1 rate pieces
    rate_piece = [0]*(2*L-1)
    for i in range(2*L-1):
        if i <= L-2:# rate piece from 1 to L-1
            rate_piece[i] = max(0.5*np.log2(shaping_lattice[per_s[i]]/shaping_lattice[per_s[i+1]]), 0)
        elif i == L-1:
            rate_piece[i] = max(0.5*np.log2(shaping_lattice[per_s[i]]/coding_lattice[per_c[i-(L-1)]]), 0)
        elif i >= L:
            rate_piece[i] = max(0.5*np.log2(coding_lattice[per_c[i-L]]/coding_lattice[per_c[i-(L-1)]]), 0)
    
    
    #compute the coefficient of rate pieces of the conditional entropy 
    entropy_coefficient_list = Relay_Forward_Rate(per_s,per_c,A_best_LLL)
    
    # Main Part: check those result whethe satisfy My NEW_CCF system constriants
    
    #---------------------------
    #     Construct the coefficients of constriants
    #----------------------------
    
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
        A_ConstriantMatrix=SourseRate+entropy_coefficient_list
        b_ConstriantVector=source_rate_upbound_list+rate_sec_hop[0:M]
        #change the parallel channel capacity constraints
        subset_list=list(powerset(set(range(0,L))))
        for i in range(L+1,len(subset_list)):
            bound_sum=0
            for j in subset_list[i]:
                bound_sum=bound_sum+rate_sec_hop[j]
            b_ConstriantVector.append(bound_sum)
    elif channel_mode=="MAC":
        b_ConstriantVector=source_rate_upbound_list+rate_sec_hop
        A_ConstriantMatrix=SourseRate+entropy_coefficient_list
    
    feasible_flag = True
    for i in range(len(A_ConstriantMatrix)):
        temp = np.dot(np.array(A_ConstriantMatrix[i]), np.array(rate_piece))
        if temp <= b_ConstriantVector[i] + temp * 10**(-5): # prevent the numerical error
            continue
        else:
            feasible_flag = False
            print 'constriant', i+1, 'is not satisfied!'
            #raise
            #break
    
    if feasible_flag:
        print 'All constriants are satisfied!'
    
    return feasible_flag# return the logical value
    
        
if __name__=="__main__":
    num_batch = 240
    sum_rate=[]
    New_sum_rate=[]
    New_sum_time=[]
    sum_time=[]
    #ratelist
    #result_list=[]
    #PI_con=[10**1,10**1.5,10**2,10**2.5,10**3,10**3.5]
    #PI_con=[10**2, 10**2.2, 10**2.4, 10**2.6, 10**2.8, 10**3.0]
    PI_con=[10**1.8, 10**2.0, 10**2.2, 10**2.4, 10**2.6, 10**2.8, 10**3.0]
    #PI_con=[10**2.6, 10**2.8, 10**3.0, 10**3.2, 10**3.4]
    #PI_con=[10**1.5]
    print 'Simulation Starts!\n'
    t1=time.time()
    for Pi in PI_con:
        result_list=list(CCF_Model_Comparison([(SearchAlgorithm,Pi,k_P_ratio*Pi)]*num_batch))
        New_Rate_list=[result_list[i][1][0] for i in range(0,num_batch)]
        Rate_list=[result_list[i][1][1] for i in range(0,num_batch)]
        New_time_list=[result_list[i][1][2] for i in range(0,num_batch)]
        time_list=[result_list[i][1][3] for i in range(0,num_batch)]
        ##
        New_ratelist=[New_Rate_list[i] for i in range(0,num_batch)]
        New_sum_rate.append(sum(New_ratelist)/num_batch)
        ratelist=[Rate_list[i] for i in range(0,num_batch)]
        sum_rate.append(sum(ratelist)/num_batch)
        New_timelist=[New_time_list[i] for i in range(0,num_batch)]
        New_sum_time.append(sum(New_timelist)/num_batch)
        timelist=[time_list[i] for i in range(0,num_batch)]
        sum_time.append(sum(timelist)/num_batch)
    t2=time.time()
    print 'Total Time Cost: ' ,(t2-t1)
    print 'New CCF Model Time Cost:' , New_sum_time
    print 'CCF Model Time Cost:' , sum_time
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    plot_rate=list_plot(zip(PI_dB,sum_rate),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'CCF_Modle',gridlines=True)
#     plot_rate.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
#     plot_rate.set_legend_options(loc='upper left')
    plot_new_rate=list_plot(zip(PI_dB,New_sum_rate),plotjoined=True, marker='o', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'New_CCF_Modle',gridlines=True)
#     plot_new_rate.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
#     plot_new_rate.set_legend_options(loc='upper left')
    plot_compare=plot_new_rate+plot_rate
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    #plot_compare.title('Comparision of Two CCF')
    plot_compare.set_legend_options(loc='upper left')
    
    plot_compare.save('/home/haizi/Pictures/Results/' + 'L=M=' + L.__str__() + time.ctime() +'.png')
    plot_compare.show()
    raw_input()