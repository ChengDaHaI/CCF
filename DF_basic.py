from sage.all import *
import numpy as np
from NEW_basic import *
from AF_basic import Compute_AF_rate
import time
from math import log10, fabs
# channel matrix
# H_a = matrix(RR, 3, 3, [ [-0.604774174080910, -0.516611703927027, 0.0251878692137226],\
#                          [-0.350171195717287,  0.814517492278491, -0.236238019733556],\
#                          [ 0.232228528459459, -0.518860603180491,  0.973647111105997]])
# H_b = matrix(RR, 1, 3, [ 0.227086968428515,  0.682635663808828, -0.814728906414353])


# sum rate of decode and forward Scheme
# input : Ha, Hb, Ps, Pr
# output: DF sum rate
def Compute_DF_rate(Ha, Hb, Ps, Pr):
    
    # first hop  channel rates 
    RelayDecodeRate = [0] * M
    for i in range(0,M):
        RelayDecodeRate[i] = 0.5 * log(1 + Ha.row(i).norm()**2*Ps, 2)
    
    # second hop channel rates
    DestinationRate = [0] * M
    for i in range(0,M):
        DestinationRate[i] = 0.5 * log(1 + Hb.column(i).norm()**2*Pr, 2)
    
    # the DF sum rate is the minimal of { first hop sum rate(MISO MAC channel), maximal of second hop rate }
    # this is the original DF, the poorest performance.
    # DF_sumrate = min(RelayDecodeRate + [max(DestinationRate)])
    DF_rate_perlink = [0]*M
    for i in range(M):
        DF_rate_perlink[i] = min(DestinationRate[i],RelayDecodeRate[i])

    # DF scheme only need one relay decoding successfully or need all relay decoding successfully.
    DF_sumrate = max(DF_rate_perlink) # one relay
    # DF_sumrate = min(DF_rate_perlink) # all relay
    return DF_sumrate

@parallel(ncpus=Cores)
def Compare_DF_AF(P_s, P_r):
    set_random_seed()
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    #print 'First Hop channel matrix H_a:\n',H_a
    H_b = matrix.random(RR, 1, M, distribution=RealDistribution('gaussian', 1))
    #print 'Second Hop channel matrix H_b:\n',H_b

    sum_rate=[0]*2
    
    sum_rate[0] = Compute_DF_rate(H_a, H_b, P_s, P_r)
    sum_rate[1] = Compute_AF_rate(P_s, P_r,H_a, diagonal_matrix(list(H_b.row(0))), M, M)

    return sum_rate

if __name__ == '__main__':
#     P_s = 10**(28.0/10)
#     P_r = 0.25 * P_s
#     DF_rate = Compute_DF_rate(H_a, H_b, P_s, P_r)
#     AF_rate = Compute_AF_rate(P_s, P_r,H_a, diagonal_matrix(list(H_b.row(0))), 3,3)
#     print 'sum rate of decode and forward Scheme is \n', DF_rate, AF_rate
    print 'Simulation Start!\n'
    t1=time.time()
    num_batch = 800
    #PI_con=[10**2.0, 10**2.2, 10**2.4, 10**2.6, 10**2.8, 10**3.0, 10**3.2, 10**3.4, 10**3.6, 10**3.8, 10**4.0]
    #PI_con=[10**2.0, 10**2.25, 10**2.5, 10**2.75, 10**3.0, 10**3.25, 10**3.5]
    PI_con=[10**1.0, 10**1.2, 10**1.4, 10**1.6, 10**1.8, 10**2.0, 10**2.2, 10**2.4, 10**2.6, 10**2.8, 10**3.0]
    sum_rate_list=[]#save the sum_rate list of two scheme
    DF_rate_list=[]
    AF_rate_list=[]
    for P_con in PI_con:
        result=list(Compare_DF_AF([(P_con, 0.25 * P_con)]*num_batch))
        DF_list=[result[i][1][0] for i in range(0,num_batch)]
        AF_list=[result[i][1][1] for i in range(0,num_batch)]
        DF_rate=sum(DF_list)/num_batch
        AF_rate=sum(AF_list)/num_batch
        DF_rate_list.append(DF_rate)
        AF_rate_list.append(AF_rate)
    t2=time.time()
    print 'Simulation cost %i s'%(t2-t1)
    # save result
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    Full_Result = np.column_stack((PI_dB, DF_rate_list, AF_rate_list))
    np.savetxt('/home/haizi/Pictures/Results/TxtFile/' + 'DF_AF' + time.ctime() + 'L=' + L.__str__() + 'iter = ' + num_batch.__str__() + 'Full_Result.txt', Full_Result ,fmt = '%1.5e')
    
    #Plot the figure
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    plot_DF=list_plot(zip(PI_dB,DF_rate_list),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'DF')
    plot_AF=list_plot(zip(PI_dB,AF_rate_list),plotjoined=True, marker='D', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'AF')
    plot_compare=plot_DF+plot_AF
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    plot_compare.save('/home/haizi/Pictures/Results/' + 'DF_AF' + time.ctime() + 'L=' + L.__str__() + 'iter = ' + num_batch.__str__() + '.eps')
    plot_compare.show(gridlines=True)