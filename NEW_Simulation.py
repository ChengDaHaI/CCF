'''
This file is to simulate our new CCF modle 
and get a convinced result through 
Monte-Carlo simulation.
'''
from sage.all import *
from NEW_Optimize_Modle import RandomSearch
from NEW_basic import *
from math import log10
import time


@parallel(ncpus=Cores)
def NONNameFunc(P_Search_Alg,P_con,P_relay,per_s,per_c):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    H_b = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    beta_opt, sum_rate_opt=RandomSearch(P_Search_Alg, H_a, H_b, P_con, P_relay, per_s, per_c)
    return beta_opt, sum_rate_opt

if __name__=="__main__":
    num_batch=120
    sum_rate=[]
    #ratelist
    #result_list=[]
    PI_con=[10,100,1000]
    for Pi_c in PI_con:
        t1=time.time()
        result_list=list(NONNameFunc([(SearchAlgorithm,Pi_c,k_P_ratio*Pi_c,per_s,per_c)]*num_batch))
        t2=time.time()
        Rate_list=[result_list[i][1][1] for i in range(0,num_batch)]
        ratelist=[Rate_list[i] for i in range(0,num_batch)]
        sum_rate.append(sum(ratelist)/num_batch)
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    plot_compare=list_plot(zip(PI_dB,sum_rate),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'New_CCF_Modle')
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    #show(plot_compare)
    plot_compare.show(gridlines=True)
    #plot(sum_rate,PI_dB)
    raw_input()