'''
This file is to simulate our new CCF modle 
and get a convinced result through 
Monte-Carlo simulation.
'''
from sage.all import *
from NEW_Optimize_Modle import RandomSearch
from NEW_basic import *
from NewSecondHopChannel import ComputeSecRate
from ComputeRate import CoF_compute_search_pow_flex_beta
from math import log10
import time
import itertools


@parallel(ncpus=Cores)
def CCF_Model_Comparison(P_Search_Alg,P_con,P_relay):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    H_b = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    rate_sec_hop=ComputeSecRate(M,P_relay,H_b)
    Max_New_sum_rate=0
    t1=time.time()
    per_search=False
    if per_search==True:
        for shape_order in itertools.permutations(list(range(0, L)), L):
            per_s=list(shape_order)
            for code_order in itertools.permutations(list(range(0, L)), L):
                per_c=list(code_order)
                (beta_opt, New_sum_rate_opt)=RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c)
                if Max_New_sum_rate<New_sum_rate_opt:
                    Max_New_sum_rate=New_sum_rate_opt
        New_sum_rate_opt=Max_New_sum_rate
    elif per_search==False: 
        global per_s, per_c
        (beta_opt, New_sum_rate_opt)=RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c)
    t2=time.time()
    sum_rate_opt=CoF_compute_search_pow_flex_beta(P_con,H_a,True,True,P_Search_Alg,rate_sec_hop,'asym_mod','asym_quan')
    t3=time.time()
    return New_sum_rate_opt, sum_rate_opt,(t2-t1),(t3-t2)

if __name__=="__main__":
    num_batch=8
    sum_rate=[]
    New_sum_rate=[]
    New_sum_time=[]
    sum_time=[]
    #ratelist
    #result_list=[]
    PI_con=[10,100,1000,1000]
    print 'Simulation Starts!\n'
    for Pi in PI_con:
        t1=time.time()
        result_list=list(CCF_Model_Comparison([(SearchAlgorithm,Pi,k_P_ratio*Pi)]*num_batch))
        t2=time.time()
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
    print 'Total Time Cost: ' ,(t2-t1)
    print 'New CCF Model Time Cost:' , New_sum_time
    print 'CCF Model Time Cost:' , sum_time
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    plot_rate=list_plot(zip(PI_dB,sum_rate),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'CCF_Modle',gridlines=True)
    plot_rate.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_rate.set_legend_options(loc='upper left')
    plot_new_rate=list_plot(zip(PI_dB,New_sum_rate),plotjoined=True, marker='o', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'New_CCF_Modle',gridlines=True)
    plot_new_rate.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_new_rate.set_legend_options(loc='upper left')
    plot_compare=plot_new_rate+plot_rate
    #plot_compare.save("/home/chenghai/pictures/foo1.png")
    show(plot_compare)
    raw_input()