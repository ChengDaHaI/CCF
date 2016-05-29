#coding=utf-8
from CoF_basic import *
#from CoF_arb import CoF_compute_search_pow_flex 
from CoF_LLL import CoF_compute_fixed_pow_flex
from NewSecondHopChannel import ComputeSecRate
from sage.all import *
from scipy import optimize
import time
from CoF_GA import *

#H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
#H_a=matrix(RR,M,L,[[-0.811303682004443 ,0.580911387597112],[ 0.807900039113608,-0.451669881492410]])
#H_a=matrix(RR, M, L, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
#H_b=(matrix(RR,M,L,[[ 0.325586153061723, -0.194076144797190],[-0.568119458140173, -0.312043789566336]])).column(0)
#H_a=matrix(RR,M,L,[[ 0.325586153061723, -0.194076144797190],[-0.568119458140173, -0.312043789566336]])
#print 'Channel Matrix:\n',H_a


@parallel(ncpus=Cores)
def CoF_compute_search_pow_flex(P_con, H_a, is_dual_hop, P_Search_Alg, rate_sec_hop=[], mod_scheme='asym_mod', quan_scheme='asym_quan', beta=[]):
    (M, L) = (H_a.nrows(), H_a.ncols())
    set_random_seed()#avoid producing same random number from other thread
    if beta == []:
        beta = vector(RR, [1]*L)
    cof_pow = lambda x: -CoF_compute_fixed_pow_flex(x, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
    cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex(x[0:L], P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, x[L:L+M]))
    Pranges = ((P_con/brute_number, P_con), )*L # (slice(0, P_con+0.1, P_con/brute_number), )*L
    initial_guess = [0.5*P_con]*L
    try:
        if P_Search_Alg == 'brute':
            res_cof = optimize.brute(cof_pow, Pranges, Ns=brute_number, full_output=True, finish=None)
            P_opt = res_cof[0]
            sum_rate_opt = -res_cof[1] # negative! see minus sign in cof_pow
        elif P_Search_Alg == 'fmin_cg':
            res_cof=optimize.fmin_cg(cof_pow,initial_guess,gtol=1e-4,epsilon=1e-4)
            P_opt = res_cof[0]
            sum_rate_opt = -res_cof[1]
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
        #Add differential evolution
        elif P_Search_Alg =="differential_evolution":
            bounds=((0.1,P_con),)*L
            res_brute=optimize.differential_evolution(cof_pow,bounds)
            P_opt=res_brute.x
            sum_rate_opt=-res_brute.fun
        #Add Genetic Algorithm
        elif P_Search_Alg=="genetic":
            res_cof=GeneticAlgorithm(P_con,H_a,is_dual_hop,rate_sec_hop,mod_scheme,quan_scheme)
            P_opt=res_cof[0]
            sum_rate_opt=res_cof[1]
        #The Genetic Algorithm End
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    #return sum_rate_opt
    return (sum_rate_opt,P_opt)#return tuple,for comparing

@parallel(ncpus=Cores)
def CoF_compute_search_pow_flex_beta(P_con, H_a, is_fixed_power, is_dual_hop, P_Search_Alg, rate_sec_hop=[], mod_scheme='asym_mod', quan_scheme='asym_quan'):
    (M, L) = (H_a.nrows(), H_a.ncols())
    '''
    def cof_pow_beta(x):
        power = x[0:L]
        beta = vector(RR, [1,]+list(x[L:2*L-1]))
        -CoF_compute_fixed_pow_flex(power, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, beta)
    '''
    if is_fixed_power == False:
        cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex(x[0:L], P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, [1,]+list(x[L:2*L-1])))
        Pranges = ((P_con/brute_number, P_con), )*L + ((float(beta_max)/brute_number, beta_max), )*(L-1)
    else:
        cof_pow_beta = lambda x: -CoF_compute_fixed_pow_flex((P_con,)*L, P_con, False, H_a, is_dual_hop, rate_sec_hop, mod_scheme, quan_scheme, vector(RR, [1,]+list(x[0:L-1])))
        Pranges = ((0.01, beta_max), )*(L-1)
        
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
        #add differential evolution algorithm
        elif P_Search_Alg=='differential_evolution':
            #Pranges=[(float(beta_max)/brute_number, beta_max)]
            #test program running time cost
            t1=time.time()
            res_brute=optimize.differential_evolution(cof_pow_beta,Pranges, maxiter = 50, disp = True)
            t2=time.time()
            t=t2-t1
            P_opt=res_brute.x
            # print 'optimal beta and source power(original)', P_opt
            sum_rate_opt=-res_brute.fun
            # print 'maximum sum source rate(original)', sum_rate_opt
            print 'CCF Differential Evolution:', res_brute.success
        #end
        else:
            raise Exception('error: algorithm not supported')
    except:
        print 'error in search algorithms'
        raise
    # retrun beta to check the feasiblility
    is_return_beta = True
    if is_return_beta:
        return sum_rate_opt, P_opt
    else:
        return sum_rate_opt


#传入P_con,H_a,is_dual_hop,rate_sec_hop等参数
#计算相同参数下三种算法的运行结果
#其中遗传算法用并行计算以减少运算时间
#返回sum_rate_list, time_list
#返回顺序为：brute,differential_evolution,genetic
@parallel(ncpus=Cores)
def CoF_compare_algorithm(P_con,is_dual_hop,rate_sec_hop=[],mod_scheme='asym_mod', quan_scheme='asym_quan'):
    set_random_seed()#避免产生相同信道矩阵
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
    #H_a=matrix(RR,M,L,[[-0.541978155712295 ,0.740073351426688],[-0.773073785628476 ,0.584325217080305]])
    #H_a=matrix(RR,M,L,[[ 0.0582433959114153 ,-0.0331203286619994  ,-0.298110577723396],[  0.204499188754200  ,0.516092269060844 ,-0.388103212165450],[ -0.235001703409716 ,0.0189319096867924  ,0.935347565704346]])

    print 'First Hop channel matrix H_a:\n',H_a
    if is_dual_hop==True:
        #parallel channel
        H_b = (matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))).column(0)
        print 'Second Hop channel matrix H_b:\n',H_b
        rate_sec_hop=[0]*M
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
        '''
        #产生MIMO信道forwarding rate bounds
        H_b = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
        print 'Second Hop channel matrix H_b:\n',H_b
        rate_sec_hop=[]
        rate_sec_hop.extend(ComputeSecRate(M,P_relay,H_b))
        '''
    #Algorithm=['brute','differential_evolution','genetic']
    #Algorithm=['TNC','differential_evolution','genetic']
    #Algorithm=['differential_evolution']
    sum_rate_list=[]
    time_list=[]
    #for P_Search_Alg in Algorithm:
    P_Search_Alg='differential_evolution'
    print P_Search_Alg,':'
    t1=time.time()
    sum_rate,P_opt=CoF_compute_search_pow_flex(P_con,H_a,is_dual_hop,P_Search_Alg,rate_sec_hop,mod_scheme, quan_scheme)
    t2=time.time()
    sum_rate_list.append(sum_rate)
    time_list.append(t2-t1)
    #减小种群规模，并行遗传算法取各次结果最优值，缩短运行时间
    '''
    P_Search_Alg='genetic'
    print P_Search_Alg,':'
    t1=time.time()
    result=list(CoF_compute_search_pow_flex([(P_con,H_a,is_dual_hop,P_Search_Alg,rate_sec_hop,mod_scheme,quan_scheme)]*8))
    t2=time.time()
    rate_list=[result[i][1][0] for i in range(0,8)]
    #P_opt_list=[result[i][1][1] for i in range(0,num_batch)]
    sum_rate=max(rate_list)
    sum_rate_list.append(sum_rate)
    time_list.append(t2-t1)
    '''
    #返回顺序为：brute,differential_evolution,genetic
    return sum_rate_list,time_list

if __name__=='__main__':
    num_batch=8
    #PI_con=[1023,2047,3071,4095]
    PI_con=[1000]
    is_dual_hop=True
    sum_rate_brute=[]
    sum_rate_genetic=[]
    sum_rate_differential=[]
    time_brute=[]
    time_genetic=[]
    time_differential=[]
    for P_con in PI_con:
        P_relay=0.5*P_con
        t1=time.time()
        result_list=list(CoF_compare_algorithm([(P_con, is_dual_hop)]*num_batch))
        t2=time.time()
        print 'When power constraint is %i,time spend: %i s'%(P_con,t2-t1)
        Rate_list=[result_list[i][1][0] for i in range(0,num_batch)]
        Time_list=[result_list[i][1][1] for i in range(0,num_batch)]
        '''
        #brute
        ratelist=[Rate_list[i][0] for i in range(0,num_batch)]
        sum_rate_brute.append(sum(ratelist)/num_batch)
        timelist=[Time_list[i][0] for i in range(0,num_batch)]
        time_brute.append(sum(timelist)/num_batch)
        #differential
        ratelist=[Rate_list[i][1] for i in range(0,num_batch)]
        sum_rate_differential.append(sum(ratelist)/num_batch)
        timelist=[Time_list[i][1] for i in range(0,num_batch)]
        time_differential.append(sum(timelist)/num_batch)
        #genetic
        ratelist=[Rate_list[i][2] for i in range(0,num_batch)]
        sum_rate_genetic.append(sum(ratelist)/num_batch)
        timelist=[Time_list[i][2] for i in range(0,num_batch)]
        time_genetic.append(sum(timelist)/num_batch)
        '''
        ratelist=[Rate_list[i][0] for i in range(0,num_batch)]
        sum_rate_differential.append(sum(ratelist)/num_batch)
        timelist=[Time_list[i][0] for i in range(0,num_batch)]
        time_differential.append(sum(timelist)/num_batch)
    #画图
    #表现各算法性能之间的关系
    PI_dB=[10*log10(P_con) for P_con in PI_con]
    '''
    #Rate Comparison
    plot_brute=list_plot(zip(PI_dB,sum_rate_brute), plotjoined=True, marker='o', \
                                  rgbcolor=Color('black'), linestyle="--", \
                                  legend_label= 'TNC', \
                                  title = 'Rate Comparison')
    plot_differential=list_plot(zip(PI_dB,sum_rate_differential),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'differential_evolution')
    plot_genetic=list_plot(zip(PI_dB,sum_rate_genetic),plotjoined=True, marker='<', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'genetic')
    plot_compare=plot_brute+plot_genetic+plot_differential
    '''
    plot_compare=list_plot(zip(PI_dB,sum_rate_differential),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'differential_evolution')
    plot_compare.axes_labels(['SNR(dB)', 'Sum rate(bps)'])
    plot_compare.set_legend_options(loc='upper left')
    #show(plot_compare)
    plot_compare.show(gridlines=True)
    '''
    #Time Comparison
    plot_brute=list_plot(zip(PI_dB,time_brute), plotjoined=True, marker='o', \
                                  rgbcolor=Color('black'), linestyle="--", \
                                  legend_label= 'TNC', \
                                  title = 'Time Comparison')
    plot_differential=list_plot(zip(PI_dB,time_differential),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'differential_evolution')
    plot_genetic=list_plot(zip(PI_dB,time_genetic),plotjoined=True, marker='<', \
                                      rgbcolor=Color('green'), linestyle='-.', \
                                      legend_label = 'genetic')
    plot_compare=plot_brute+plot_genetic+plot_differential
    '''
    plot_compare=list_plot(zip(PI_dB,time_differential),plotjoined=True, marker='d', \
                                      rgbcolor=Color('blue'), linestyle='-.', \
                                      legend_label = 'differential_evolution')
    plot_compare.axes_labels(['SNR(dB)', 'time(s)'])
    plot_compare.set_legend_options(loc='upper left')
    #show(plot_compare)
    plot_compare.show(gridlines=True)
    raw_input()
    
