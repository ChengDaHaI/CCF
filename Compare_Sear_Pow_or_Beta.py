from sage.all import *
from CoF_basic import *
from ComputeRate import CoF_compute_search_pow_flex
from ComputeRate import CoF_compute_search_pow_flex_beta
from AF_basic import Compute_AF_rate
from NewSecondHoP import ComputeSecRate
from math import log10
import time

@parallel(ncpus=Cores)
def Compare_Pow_Beta(P_con,is_dual_hop):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    #print 'First Hop channel matrix H_a:\n',H_a
    H_b = matrix.random(RR, 2, M, distribution=RealDistribution('gaussian', 1)).column(0)
    #print 'Second Hop channel matrix H_b:\n',H_b
    P_relay=0.5*P_con#set the P_relay
    sum_rate=[0]*2
    rate_sec_hop=[]
    if is_dual_hop==True:
        #rate_sec_hop.extend(ComputeSecRate(M,P_relay,H_b))
        rate_sec_hop = [0]*M # ? bit/s for each parallel channel in the second hop
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
    P_search_Alg='differential_evolution'
    sum_rate[0],P_opt=CoF_compute_search_pow_flex(P_con,H_a,is_dual_hop,P_search_Alg,rate_sec_hop)
    is_fixed_power=True
    sum_rate[1]=CoF_compute_search_pow_flex_beta(P_con,H_a,is_fixed_power,is_dual_hop,P_search_Alg,rate_sec_hop)
    return sum_rate

if __name__=="__main__":
    print 'Simulation Start!\n'
    t1=time.time()
    num_batch=120
    is_dual_hop=False
    P_I=[1,10,100,1000,10000,100000]
    sum_rate_list=[]#save the sum_rate list of two scheme
    Power_rate_list=[]
    Beta_rate_list=[]
    for P_con in P_I:
        result=list(Compare_Pow_Beta([(P_con, is_dual_hop)]*num_batch))
        Power_list=[result[i][1][0] for i in range(0,num_batch)]
        Beta_list=[result[i][1][1] for i in range(0,num_batch)]
        Power_rate=sum(Power_list)/num_batch
        Beta_rate=sum(Beta_list)/num_batch
        Power_rate_list.append(Power_rate)
        Beta_rate_list.append(Beta_rate)
    t2=time.time()
    print 'Simulation cost %i s'%(t2-t1)
    #Plot the figure
    PI_dB=[10*log10(P_con) for P_con in P_I]
    plot_CF=list_plot(zip(PI_dB,Power_rate_list),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'Power Search')
    plot_AF=list_plot(zip(PI_dB,Beta_rate_list),plotjoined=True, marker='D', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'Beta Search')
    plot_compare=plot_CF+plot_AF
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    plot_compare.show(gridlines=True)
    
    raw_input()