from sage.all import *
from CoF_basic import *
from ComputeRate import CoF_compute_search_pow_flex
from AF_basic import Compute_AF_rate
from NewSecondHoP import ComputeSecRate
from math import log10
import time
#compare AF and CCF
#return sum_rate list
def Compare_fix_H(H_a,H_b,P_con,is_dual_hop):
    P_relay=0.25*P_con#set the P_relay
    sum_rate=[0]*L
    rate_sec_hop=[]
    rate_sec_hop.extend(ComputeSecRate(M,P_relay,H_b))
    P_search_Alg='differential_evolution'
    sum_rate[0],P_opt=CoF_compute_search_pow_flex(P_con,H_a,is_dual_hop,P_search_Alg,rate_sec_hop)
    #P_search_Alg='genetic'
    sum_rate[1]=Compute_AF_rate(P_con,P_relay,H_a,H_b)
    return sum_rate

@parallel(ncpus=Cores)
def Compare_CF_AF(P_con,is_dual_hop):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    #print 'First Hop channel matrix H_a:\n',H_a
    H_b = matrix.random(RR, 2, 2, distribution=RealDistribution('gaussian', 1))
    #print 'Second Hop channel matrix H_b:\n',H_b
    P_relay=0.15*P_con#set the P_relay
    sum_rate=[0]*2
    rate_sec_hop=[]
    rate_sec_hop.extend(ComputeSecRate(M,P_relay,H_b))
    P_search_Alg='differential_evolution'
    sum_rate[0],P_opt=CoF_compute_search_pow_flex(P_con,H_a,is_dual_hop,P_search_Alg,rate_sec_hop)
    #P_search_Alg='genetic'
    sum_rate[1]=Compute_AF_rate(P_con,P_relay,H_a,H_b,L=2,M=2)
    return sum_rate

if __name__=="__main__":
    print 'Simulation Start!\n'
    t1=time.time()
    num_batch=120
    is_dual_hop=True
    P_I=[1,10,100,1000,10000,100000]
    sum_rate_list=[]#save the sum_rate list of two scheme
    CF_rate_list=[]
    AF_rate_list=[]
    for P_con in P_I:
        result=list(Compare_CF_AF([(P_con, is_dual_hop)]*num_batch))
        CF_list=[result[i][1][0] for i in range(0,num_batch)]
        AF_list=[result[i][1][1] for i in range(0,num_batch)]
        CF_rate=sum(CF_list)/num_batch
        AF_rate=sum(AF_list)/num_batch
        CF_rate_list.append(CF_rate)
        AF_rate_list.append(AF_rate)
    t2=time.time()
    print 'Simulation cost %i s'%(t2-t1)
    #Plot the figure
    PI_dB=[10*log10(P_con) for P_con in P_I]
    plot_CF=list_plot(zip(PI_dB,CF_rate_list),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'CF')
    plot_AF=list_plot(zip(PI_dB,AF_rate_list),plotjoined=True, marker='D', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'AF')
    plot_compare=plot_CF+plot_AF
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    plot_compare.show(gridlines=True)
    '''
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    print 'First Hop channel matrix H_a:\n',H_a
    H_b = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    print 'Second Hop channel matrix H_b:\n',H_b
    for P_con in P_I:
        sum_rate_list.append(Compare_fix_H(H_a, H_b, P_con, is_dual_hop))
    AF_list=[sum_rate_list[i][1] for i in range(0,len(P_I))]
    CF_list=[sum_rate_list[i][0] for i in range(0,len(P_I))]
    
    #Plot the figure
    PI_dB=[10*log10(P_con) for P_con in P_I]
    plot_genetic=list_plot(zip(PI_dB,CF_list),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'CF')
    plot_AF=list_plot(zip(PI_dB,AF_list),plotjoined=True, marker='D', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'AF')
    plot_compare=plot_genetic+plot_AF
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    show(plot_compare)
    '''
    raw_input()