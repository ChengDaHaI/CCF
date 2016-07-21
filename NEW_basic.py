from sage.all import *
L = 3
M = 3
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

set_H_a = matrix(RR, 3, 3, [[  0.699276348994144,   0.979966803608800,   0.731095879215959],
	                        [-0.0467540769081729,  0.253952649489866,  -0.775007249377646],
	                        [  0.669860797898610,   0.183608409840893,  -0.343789835395397]])

set_H_b = matrix(RR, 1, 3, [-0.622194984940395, -0.563158590610974, -0.777483868327708])

'''
set_H_a = matrix(RR, 3, 3, [[-0.420339821302539, -0.256727346251365, -0.720655914254134],
                            [-0.402509674288030, -0.682940704430745, -0.203551145932588],
                            [ 0.665845601783994, -0.211669229448582, -0.180877130764904]])

set_H_b = matrix(RR, 1, 3, [-0.936758795240084, 0.880013370030291, -0.369571658106388])
'''
SearchAlgorithm="differential_evolution"
#H_a=matrix.random(RR, L, M, distribution=RealDistribution('gaussian', 1))
#produce the second channel matrix
#H_b= matrix.random(RR, 1, L, distribution=RealDistribution('gaussian', 1))