from sage.all import *
L=3
M=3
p=17
#A=Matrix([[1,2,3],[1,2,2],[2,1,3]])#matrix A must be full-rank
#beta_s=[2.5,2.0,1.8]
#beta_c=[1.5,1.0,0.8]
#per_s=[0,1,2]
#per_c=[2,1,0]
#P_con=1000
k_P_ratio=0.25
betaScale_max=4
Cores=20
SearchAlgorithm="differential_evolution"
#H_a=matrix.random(RR, L, M, distribution=RealDistribution('gaussian', 1))
#produce the second channel matrix
#H_b= matrix.random(RR, 1, L, distribution=RealDistribution('gaussian', 1))