#coding=utf-8 
'''This is a genetic algorithm for Compute-and-Forward scheme'''
'''We ever try to utilize the DEAP software package'''
'''
1、初始化，产生初始个体和种群，确定交叉、变异和进化代数等数值
2、计算每一个体适应度
3、选择操作
4、交叉与变异操作
5、满足循环停止条件则退出（？）
'''
#import sys
#sys.path.append("/usr/local/lib/python2.7/dist-packages")
#print sys.path
#import random
import math
import time
from copy import copy
from copy import deepcopy
from operator import itemgetter

from sage.all import *
from random import *
from numpy import *
from CoF_basic import *
from CoF_LLL import CoF_compute_fixed_pow_flex#评估函数

#指定迭代代数NGEN，判断终止迭代参数BREAK
#变异概率MUPB，能量上限P_con
#每个维度染色体数目mount，
NGEN=30
BREAK=NGEN/3
MUPB=0.05
mount=8
#产生随机信道矩阵
'''
P_con=100
chromolength=int(math.ceil(math.log(P_con+1,2)))
if is_set_H == True:
    H_a = set_H_a
else:
    set_random_seed() # to avoid producing the same H_a in different threads
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
#H_a = matrix(RR, 2, 2, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
'''
#由功率P_con计算染色体长度
def ComputeChromolength(P_con):
    return int(math.ceil(math.log(P_con+1,2)))

#十进制转二进制数序列
#输入为十进制数，输出为转换后二进制0,1序列
'''
def int2bin(pop):
    popbin=[]
    str=bin(pop)[2:]
    for i in range(0,chromolength-len(str)):
        popbin.append(0)
    for i in range(0,len(str)):
        popbin.append(int(str[i]))
    return popbin
'''
#实数转化为精度为0.125的二进制序列
#序列长度加三
def int2bin(pop):
    popbin=[]
    pop_int=int(pop)#pop整数部分
    pop_dec=pop-pop_int#pop小数部分
    str=bin(pop_int)[2:]
    for i in range(0,chromolength-len(str)):
        popbin.append(0)
    for i in range(0,len(str)):
        popbin.append(int(str[i]))
    temp=int(pop_dec/0.125)
    str=bin(temp)[2:]
    for i in range(0,len(str)):
        popbin.append(int(str[i]))
    for i in range(0,3-len(str)):
        popbin.append(0)
    return popbin
    
#二进制数转化十进制数
#输入为二进制序列，输出十进制数
def bin2int(popbin):
    pop=0
    for i in range(0,len(popbin)):
        pop+=popbin[i]*pow(2,chromolength-i-1)
    return pop

#种群产生函数，种群个体采用均匀分布
#P_con为种群个体染色体数值上限（能量上限）
#chromolength为染色体长度（二进制位数）
#mount为每个维度染色体数目
def initialpop(P_con,mount):
    Population=[]#存储初始种群
    #chromolength=math.ceil(math.log(P_con+1,2))#向上取整
    for i in range(int(pow(mount,L))):
        pop=[]
        for j in range(L):
            #l=int(uniform(0,P_con))#产生整数个体
            l=uniform(0,P_con)#产生实数个体
            pop.append(int2bin(l))
            #print pop
            #print bin2int(pop[0])
        Population.append(pop)    
    return Population

#print initialpop(100, 2, 10, 20)

#计算个体适应度函数
#仅仅返回最大支持速率
#输入Power
def evaluate(P_t,P_con,H_a):
     return CoF_compute_fixed_pow_flex(P_t, P_con, False, H_a,is_dual_hop,rate_sec_hop,mod_scheme, quan_scheme)
'''
print evaluate((1023,150))
H_a=set_H_a
print evaluate((bin2int(Population[0][0]),bin2int(Population[0][1])))
'''
 
#计算种群个体适应度函数
#返回tuple，即种群中各个个体适应度
def fitvalue(pop,P_con,H_a):
    fitness=[]
    for i in range(int(pow(mount,L))):
        popvalue=[]
        for j in range(L):
            popvalue.append(bin2int(pop[i][j]))
        fitness.append(evaluate(popvalue,P_con,H_a))
    return fitness

#定义交叉函数
'''
#两点交叉
def crossover(ind1,ind2):
    #for i in range(2):
    #chromolength=math.ceil(math.log(P_con+1,2))
    p1=randint(0,chromolength)
    p2=randint(0,chromolength-1)
    if p2>=p1:
        p2+=1
    else:
        p1,p2=p2,p1
    p3=randint(0,chromolength)
    p4=randint(0,chromolength-1)
    if p4>=p3:
        p4+=1
    else:
        p3,p4=p4,p3        
    #ind1[i][p1:p2],ind2[i][p1:p2]=ind2[i][p1:p2],ind1[i][p1:p2]
    if L==2:
        return (ind1[0][:p1]+ind2[0][p1:p2]+ind1[0][p2:],ind1[1][:p3]+ind2[1][p3:p4]+ind1[1][p4:])\
            ,(ind2[0][:p1]+ind1[0][p1:p2]+ind2[0][p2:],ind2[1][:p3]+ind1[1][p3:p4]+ind2[1][p4:])
    elif L==3:
        p5=randint(0,chromolength)
        p6=randint(0,chromolength-1)
        if p6>=p5:
            p6+=1
        else:
            p5,p6=p6,p5
        return (ind1[0][:p1]+ind2[0][p1:p2]+ind1[0][p2:],ind1[1][:p3]+ind2[1][p3:p4]+ind1[1][p4:],ind1[2][:p5]+ind2[2][p5:p6]+ind1[2][p6:])\
            ,(ind2[0][:p1]+ind1[0][p1:p2]+ind2[0][p2:],ind2[1][:p3]+ind1[1][p3:p4]+ind2[1][p4:],ind2[2][:p5]+ind1[2][p5:p6]+ind2[2][p6:])
    else:
        raise Exception("error: cannot support  a channel more than 3 input &output")
'''

#一点交叉函数
def crossover(ind1,ind2):
    p1=randint(0,chromolength+3)
    p2=randint(0,chromolength+3)
    if L==2:
        return (ind1[0][:p1]+ind2[0][p1:],ind1[1][:p2]+ind2[1][p2:])\
            ,(ind2[0][:p1]+ind1[0][p1:],ind2[1][:p2]+ind1[1][p2:])
    elif L==3:
        p3=randint(0,chromolength+3)
        return (ind1[0][:p1]+ind2[0][p1:],ind1[1][:p2]+ind2[1][p2:],ind1[2][:p3]+ind2[2][p3:])\
            ,(ind2[0][:p1]+ind1[0][p1:],ind2[1][:p2]+ind1[1][p2:],ind2[2][:p3]+ind1[2][p3:])
    else:
        raise Exception("error: cannot support  a channel more than 3 input &output")

#定义变异函数，ind为个体，indpb为变异概率
def mutate(individual,indpb):
    for i in range(0,L):
        for j in range(0,chromolength+3):
            if random() < indpb:
                individual[i][j] = type(individual[i][j])(not individual[i][j])
    #return individual

#排序函数，按照适应度对种群进行排序
#返回排序后的种群序列(zip,sorted)
#函数中对fitscore乘以100,是为了提高各个体之间适应度差异
def rankPop(pop,P_con,H_a):
    fitness=fitvalue(pop,P_con,H_a)
    sum_fit=sum(fitness)
    fitscore=[]
    for i in range(len(fitness)):
        fitscore.append(100*fitness[i]/sum_fit)
    pairedpop=zip(pop,fitscore)
    rankedpop=sorted(pairedpop,key=itemgetter(-1),reverse=True)
    return rankedpop

#使用轮盘赌选择法
#输入归一化的适应度
#输出相应选择项的下标
def Roulette(fitscore):
    caculatefitscore=0
    for i in range(len(fitscore)):
        caculatefitscore+=fitscore[i]
        if caculatefitscore>100*random():
            return i
        
#产生新种群函数
#采用精英选择和轮盘赌的方法
#输入当代种群
#输出新一代种群
def NewPop(pop,P_con,H_a):
    lenpop=len(pop)
    rankedpop=rankPop(pop,P_con,H_a)
    fitscore=[item[-1] for item in rankedpop]
    rankedchromes=[item[0] for item in rankedpop]    
    newpop=[]
    #精英选择，选择当代适应度最高的前百分之十
    newpop.extend(rankedchromes[0:lenpop/10])
    #轮盘赌选择
    #选择的个体参与遗传操作：交叉与变异
    while len(newpop)<lenpop:
        ind1,ind2=[],[]
        index1=Roulette(fitscore)
        index2=Roulette(fitscore)
        while index2==index1:
            index2=Roulette(fitscore)
        ind1=rankedchromes[index1]
        ind2=rankedchromes[index2]
        ind1,ind2=crossover(ind1, ind2)
        #crossover(ind1, ind2, CXPB)
        #ind1=mutate(ind1,MUPB)
        #ind2=mutate(ind2,MUPB)
        mutate(ind1, MUPB)
        mutate(ind2, MUPB)
        newpop.extend([ind1,ind2])
        ##newpop.append(ind1)
    if len(newpop)!=lenpop:
        newpop=newpop[0:lenpop]
    return newpop
    
#遗传算法主程序
#输入功率上限和信道矩阵
#输出最大适应度个体数值及其适应度tuple
def GeneticAlgorithm(P_con,H_a,is_double_hop,rate_sechop=[],modscheme='sym_mod', quanscheme='sym_quan'):
    global chromolength#产生全局变量
    chromolength=ComputeChromolength(P_con)
    global is_dual_hop#便于参数传递，产生全局变量
    is_dual_hop=is_double_hop
    global rate_sec_hop#便于参数传递，产生全局变量
    rate_sec_hop=rate_sechop
    global mod_scheme
    mod_scheme=modscheme
    global quan_scheme
    quan_scheme=quanscheme
    #产生初始种群
    pop=initialpop(P_con,mount)
    #记录每代最大适应度
    #当适应度在多代均保持不变后
    #终止迭代
    maxfit1,maxfit2,max_equal=0,0,0
    for i in range(0,NGEN):
        fitness=fitvalue(pop,P_con,H_a)
        '''
        print "Generation: ", i
        print "该代最大适应度：",max(fitness)
        print "最大适应度对应个体：",pop[fitness.index(max(fitness))]
        '''
        maxfit2=max(fitness)
        if maxfit2>maxfit1:
            maxfit1=maxfit2
        elif maxfit2==maxfit1:
            max_equal+=1
        #判断是否终止迭代
        if max_equal>BREAK:
            break
        #选择并交叉变异产生下一代种群
        newpop=[]
        newpop=NewPop(pop,P_con,H_a)
        pop=copy(newpop)
    fitness=fitvalue(pop,P_con,H_a)#最终代个体适应度tuple
    MaxValue=max(fitness)#最大适应度
    BestPop=pop[fitness.index(max(fitness))]#最大适应度对应个体
    Bestpopvalue=[]#最大适应度个体十进制数值
    for j in range(L):
        Bestpopvalue.append(bin2int(BestPop[j]))
    print "最大适应度：",MaxValue
    print"最大适应度对应个体：",Bestpopvalue
    return (Bestpopvalue,MaxValue)

if __name__=="__main__":
    #H_a = matrix(RR, M, L, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
    H_a=matrix(RR,M,L,[[-0.541978155712295 ,0.740073351426688],[-0.773073785628476 ,0.584325217080305]])
    #H_a= matrix(RR,M,L,[[ 0.624353542234016 ,-0.650505767940636],[ 0.916985921698907  ,0.539574915164442]])
    H_b=(matrix(RR,M,L,[[ 0.325586153061723, -0.194076144797190],[-0.568119458140173, -0.312043789566336]])).column(0)
    P_con=1023
    is_dual_hop=True
    rate_sec_hop = [0]*M
    if is_dual_hop==True:
        P_relay=0.25*P_con
        for i_h_b in range(0, M):
            rate_sec_hop[i_h_b] = 0.5*log(1+H_b[i_h_b]**2*P_relay, 2)
    t1=time.time()
    print "程序开始运行："
    GeneticAlgorithm(P_con,H_a,is_dual_hop,rate_sec_hop)
    t2=time.time()
    print "程序运行完成！"
    print "程序耗时：%i 秒"%(t2-t1)