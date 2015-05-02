__author__ = 'andreas'

import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import *
import scipy.special
from math import *

#vec1 = coo_matrix([[1], [2], [3]])
#vec2 = coo_matrix([[4], [5], [6]])
#mat1 = hstack([vec1, vec2])
# print mat1
#print mat1.toarray()
#mat1 = np.transpose(mat1)
#print mat1.toarray()

def findLowestNZIndex(v):
    i=len(v)-1
    while i>=0 and v[i]==0:
        i-=1
    return i

def H_int_ofState(v):
    result=0
    for n in v:
        result+=n*(n-1)
    return result

class BHState:
    rep=np.array([])
    def __init__(self, dim, nB):
        self.rep=np.array([nB]+[0]*(dim-1)) # initial state
    def next(self):
        i=findLowestNZIndex(self.rep[0:len(self.rep)-1])
        if i<0: # cannot increment any further
            self.rep=np.array([np.sum(self.rep)]+[0]*(len(self.rep)-1)) # initial state
            return False
        self.rep[i]-=1
        self.rep[i+1:]=np.array([np.sum(self.rep[i+1:])+1]+[0]*(len(self.rep[i+1:])-1))
        return True

def hashState(v):
    result=0.
    for i in range(len(v)):
        result+=sqrt(100.*i+3)*v[i]
    return result

def H_kin_row(stateRep, dictionary):
    result=np.array([0.]*len(dictionary))
    for i in range(len(stateRep)):
        if stateRep[i]!=0:
            stateRep[i]-=1
            stateRep[(i+1)%len(stateRep)]+=1
            # THE FOLLOWING IS THE WRONG FORMULA!!!
            if hashState(stateRep) in dictionary:
                result[dictionary[hashState(stateRep)]]=sqrt(float(stateRep[i]*stateRep[(i+1)%len(stateRep)]))
            stateRep[i]+=1
            stateRep[(i+1)%len(stateRep)]-=1
        if stateRep[(i+1)%len(stateRep)]!=0:
            stateRep[i]+=1
            stateRep[(i+1)%len(stateRep)]-=1
            # THE FOLLOWING IS THE WRONG FORMULA!!!
            if hashState(stateRep) in dictionary:
                result[dictionary[hashState(stateRep)]]=sqrt(float(stateRep[i]*stateRep[(i+1)%len(stateRep)]))
            stateRep[i]-=1
            stateRep[(i+1)%len(stateRep)]+=1
    return csc_matrix(result)

print "Setting up the vectors of the problem..."

def dimension(N, M): # problem size for N bosons on M lattice sites
    return scipy.special.binom(N+M-1,M-1)

def getIndex(v):
    M=len(v)
    if M==0:
        return 0
    N=sum(v)
    n0=N
    result=0
    while n0>v[0]:
        result+=dimension(N-n0,M-1)
        n0-=1
    return int(result+getIndex(v[1:]))

state=BHState(9,9)
i=0 # gives the indices and ultimately the size of the problem
allStates={hashState(state.rep):i}

while state.next():
    # add to allStates; set condition here
    if True:
        i+=1
        #allStates.update({hashState(state.rep):i})
        #getIndex(state.rep)

print "...Done. Setting up matrices..."

H_kin=csc_matrix(H_kin_row(state.rep, allStates))
H_int=np.array([H_int_ofState(state.rep)]+[0]*len(allStates))
j=0
while state.next():
    j+=1
    if False:#True: # same as above
        H_kin=vstack((H_kin, H_kin_row(state.rep, allStates))) # Very inefficient
        H_int[j]=H_int_ofState(state.rep)

H_i=dia_matrix(([H_int],[0]),shape=(len(H_int)-1,len(H_int)-1))
#print H_i.toarray()
#print H_kin.toarray()
print "...Done. Getting eigenvalues..."
vals, vecs = eigsh(1.*H_i+H_kin)
print "...Done! Eigenvalues are:"
print vals
print dimension(3, 3)
print dimension(4,4)