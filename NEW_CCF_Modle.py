'''
Simulation of Relay Compressing in CCF Scheme.
Author: ChengHai
Email: chenghai@shanghaitech.edu.cn
The ShangHaiTech University
'''
from sage.all import *
#from CoF_basic import *
import math
from NEW_basic import *
from itertools import chain, combinations, permutations

#generate the subsets
def powerset(iterable):
  xs = list(iterable)
  # note we return an iterator rather than a list
  return chain.from_iterable( combinations(xs,n) for n in range(len(xs)+1) )

#Given matrix H and A, coefficients shaping lattice beta_s and coding lattice beta_c,
#the nested order of shaping lattice per_s and coding lattice per_c.
#beta_s&beta_c are degressive, that is to say the first element 
#in beta is the coarsest and the last one is the finest .
#Compute the relay forwarding rates constraints.
#input:the first four lists and the last one matrix
#output: the conditional entropy list
def Relay_Forward_Rate(beta_s,beta_c,per_s,per_c,A):
    if len(beta_c)!=len(beta_s):
        raise Exception('beta_c and beta_s should have the same length')
    elif len(beta_c)!=L:
        raise Exception('beta_c and beta_s should have %i elements'%L)
    #get  beta_s and beta_c after the permutation operation
    per_beta_s=[beta_s[per_s[i]] for i in per_s]
    per_beta_c=[beta_c[per_c[i]] for i in per_c]
    #incorrect code about list
    #per_beta=[per_beta_s,per_beta_c]
    #beta=[beta_s,beta_c]
    per_beta=per_beta_s+per_beta_c
    beta=beta_s+beta_c
    #compute the total rate in sources
    rate_total=0#the total rate in sources
    for i in range(0,L):
        rate_total+=log(per_beta_s[i]/per_beta_c[i],2)
    #compute the rate pieces
    #suppose there isn't two same lattice.
    rate_piece=[0]*(2*L-1)#the rate pieces divided by nested lattices
    for i in range(0,2*L-1):
        rate_piece[i]=0.5*log(beta[i]/beta[i+1],2)
    #produce subset list
    subset_list=list(powerset(range(0,L)))
    set_L=set(range(0,L))
    #conditional entropy list
    conditional_entropy=[0]*(pow(2,L)-1)
    #coefficient of rate pieces
    entropy_coefficient=[0]*(pow(2,L)-1)
    
    #for every subset, we first calculate the sub-matrix 
    #then calculate the conditional entropy
    for i in range(1,len(subset_list)):
        piece_coefficient=[]
        conditional_entropy[i-1]=rate_total
        subset=set(subset_list[i])
        complement_set=set_L.difference(subset)
        row=[]#sub-matrix row index
        row.extend(list(complement_set))
        colum=[]#sub-matrix colum index
        for j in range(0,2*L-1):
            if j<=L-1:
                #lattice_s=per_beta[j]
                #lattice_c=per_beta[j+1]
                colum.append(per_s[j])
            elif j>L-1:
                index=colum.index(per_c[j-L])
                colum.pop(index)
            sub_A=A[row,colum]
            rank_value=0
            rank_value=rank(sub_A)
            if rank_value>rank(A):
                raise Exception('sub-matrix rank cannot great than the original matrix!')
            conditional_entropy[i-1]-=rank_value*rate_piece[j]
            #record the coefficient of rate pieces
            #that is also the i row of transform matrix between  
            #conditional_entropy list and rate_piece list
            if j<=L-1:
                piece_coefficient.append(j+1-rank_value)
            elif j>L-1:
                piece_coefficient.append(2*L-1-j-rank_value)
        entropy_coefficient[i-1]=piece_coefficient
        
    return conditional_entropy,entropy_coefficient


#given the conditional entropy piece rate coefficient list
#find all valid  vertex for the polytope     
#input : conditional entropy rate pieces coefficient list
#output: every vertex rate's piece rate coefficient  list 
def Vertex_RatePiece_Coefficient(coefficient_list):
    subset_list=list(powerset(range(1,L+1)))
    vertex_list=[]#store all vertex
    #generate object function weight order
    for weight_order in permutations(list(range(1,L+1))):
        weight_list=list(weight_order)
        vertex_rate_list=[0]*L#a vertex rate coordinates list
        for i in range(1,L+1):
            subset_Te=range(1,L+1)[0:i]
            index1=subset_list.index(tuple(subset_Te))-1#subset_list contain none subset
            if i==1:
               # vertex_rate_list[i]=entropy_list[index1]
               vertex_rate_list[i-1]=coefficient_list[index1]
            else:
                #index2=subset_list.index(tuple(subset_Te.pop(-1)))-1#subset_list contain none subset
                subset_Te.pop(-1)
                index2=subset_list.index(tuple(subset_Te))-1
                #subtraction of element in list
                vertex_rate_list[i-1]=list(map(lambda x: x[0]-x[1],zip(coefficient_list[index1],coefficient_list[index2])))
                #v = list(map(lambda x: x[0]-x[1], zip(v2, v1)))
        ordered_rate_list=[0]*L
        for i in range(L):
            #there is a permutation between vertex_rate_list and the true vertex rate list
            ordered_rate_list[weight_list[i]-1]=vertex_rate_list[i]
        temp=[]
        temp.extend(ordered_rate_list)
        vertex_list.append(temp)
    
    return vertex_list
    


#given the coefficient of  pieces rate  which form the vertex rate
#with this coefficient we can calculate which lattice tuple the relay should adopt
#input: shaping lattice beta_s, coding lattice beta_c, vertex rate coefficient list
#output: relay's lattice tuple list 
def Relay_Compress_Lattice_Tuple(beta_s,beta_c,vertex_rate_coefficient_list):
    beta=beta_s+beta_c
    vertex_mount=len(vertex_rate_coefficient_list)
    if vertex_mount!=factorial(L):
        raise Exception('there should be L s factorial vertex')
    relay_lattice_pair_list=[0]*vertex_mount
    
    for i in range(vertex_mount):
        relay_lattice_pair=[0]*L
        vertex_coefficient_list=vertex_rate_coefficient_list[i]
        for j in range(L):
            index1=vertex_coefficient_list[j].index(1)#the first nonzero element in the coefficient list
            #index2=copy(index1)
            temp=copy(vertex_coefficient_list[j])
            temp.reverse()
            index2=temp.index(1)
            index2=2*L-2-index2
            if index1==index2:
                relay_lattice_pair[j]=(beta[index1],beta[1+index1])
            else:
                relay_lattice_pair[j]=(beta[index1],beta[1+index2])
        relay_lattice_pair_list[i]=relay_lattice_pair
    
    return relay_lattice_pair_list
    
if __name__=='__main__':
    L=3
    A=Matrix([[1,2,3],[1,2,2],[2,1,3]])#matrix A must be full-rank
    beta_s=[2.5,2.0,1.8]
    beta_c=[1.5,1.0,0.8]
    per_s=[0,2,1]
    per_c=[0,1,2]
    conditional_entropy_list,entropy_coefficient_list=Relay_Forward_Rate(beta_s,beta_c,per_s,per_c,A)
    vertex_rate_coefficient_list                                    =Vertex_RatePiece_Coefficient(entropy_coefficient_list)
    relay_lattice_pair_list                                              =Relay_Compress_Lattice_Tuple(beta_s,beta_c,vertex_rate_coefficient_list)
    print relay_lattice_pair_list
    