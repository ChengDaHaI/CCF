from sage.all import *

N = 3 # number of antenna, used in SSA simulation
K = 2    # number of user
L = N * K
M = N * K
p=17
#A=Matrix([[1,2,3],[1,2,2],[2,1,3]])#matrix A must be full-rank
#beta_s=[2.5,2.0,1.8]
#beta_c=[1.5,1.0,0.8]
#per_s=[0,1,2]
#per_c=[2,1,0]
#P_con=1000
k_P_ratio=0.25
betaScale_max=4
Cores = 20

set_HaHb = False


'''
this channel corresponds to a mixed nested order, but there isn't any redundant.
set_H_a = matrix(RR, 3, 3, [[   0.820261848804884,    0.140458547512844, -0.00857116868370333],
							[  0.0252808891663803,    0.745716024139580,   -0.540385424646098],
							[  0.0888656378230059,    0.548100597640509,   -0.163573195375633]])
set_H_b = matrix(RR, 1, 3,  [-0.936314065730628,  0.618026134156727, -0.299471215058300])

#this channel can prove that Tan's compression is wrong!
set_H_a = matrix(RR, 3, 3, [[ 0.730043250944212, -0.161588668372594, -0.420512555770676],
							[-0.731239625346022, -0.667746791944543, -0.156364457993504],
							[ 0.783181077170343,  0.227800143283188,  0.821617367387947]])
set_H_b = matrix(RR, 1, 3,  [ 0.122985097913771,  0.832646255255913, -0.105196243509679])



#new_sum_rate = -0.0
set_H_a = matrix(RR, 3, 3, [[  0.307225944671425,  -0.585803434882211,  -0.840564151862982],
                            [ -0.478108939306364,   0.926766353326669, -0.0792394740282352],
                            [ -0.542984347795854,  0.0384952002592667,   0.347943766191883]])
set_H_b = matrix(RR, 1, 3, [-0.666210873369349,  0.994714108828944, -0.302557252996096])

set_H_a = matrix(RR, 3, 3, [[-0.348375475413612,  0.204705322673501, -0.172233857483020],
                            [-0.947128695185850,  0.802209214990737, -0.124109607989854],
                            [ 0.870413887174075, -0.136209452578445,  0.141487070140212]])
set_H_b = matrix(RR, 1, 3,  [-0.0460360313180472,  0.0274367396953652,   0.251385721396213])
'''
SearchAlgorithm="differential_evolution"
#H_a=matrix.random(RR, L, M, distribution=RealDistribution('gaussian', 1))
#produce the second channel matrix
#H_b= matrix.random(RR, 1, L, distribution=RealDistribution('gaussian', 1))