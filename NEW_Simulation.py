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


@parallel(ncpus=Cores)
def NONNameFunc(P_Search_Alg,P_con,P_relay,per_s,per_c):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    H_b = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    rate_sec_hop=ComputeSecRate(M,P_relay,H_b)
    beta_opt, New_sum_rate_opt=RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c)
    sum_rate_opt=CoF_compute_search_pow_flex_beta(P_con,H_a,True,True,P_Search_Alg,rate_sec_hop,'asym_mod','asym_quan')
    return New_sum_rate_opt, sum_rate_opt

if __name__=="__main__":
    num_batch=4
    sum_rate=[]
    New_sum_rate=[]
    #ratelist
    #result_list=[]
    PI_con=[10,100,1000]
    for Pi_c in PI_con:
        t1=time.time()
        result_list=list(NONNameFunc([(SearchAlgorithm,Pi_c,k_P_ratio*Pi_c,per_s,per_c)]*num_batch))
        t2=time.time()
        New_Rate_list=[result_list[i][1][0] for i in range(0,num_batch)]
        Rate_list=[result_list[i][1][1] for i in range(0,num_batch)]
        ratelist=[Rate_list[i] for i in range(0,num_batch)]
        sum_rate.append(sum(ratelist)/num_batch)
        New_ratelist=[New_Rate_list[i] for i in range(0,num_batch)]
        New_sum_rate.append(sum(New_ratelist)/num_batch)
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