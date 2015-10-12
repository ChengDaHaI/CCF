from sage.all import *
L=2
M=2
p=17
#A=Matrix([[1,2,3],[1,2,2],[2,1,3]])#matrix A must be full-rank
beta_s=[2.5,2.0]
beta_c=[1.5,1.0]
per_s=[0,1]
per_c=[0,1]
P_con=1000
P_relay=0.5*P_con
betaScale_max=4
H_a=matrix.random(RR, L, M, distribution=RealDistribution('gaussian', 1))
#produce the second channel matrix
H_b= matrix.random(RR, 1, L, distribution=RealDistribution('gaussian', 1))