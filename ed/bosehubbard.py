__author__ = 'andreas'

import numpy as np
from scipy.sparse import *
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
            # THE FOLLOWING IS THE WRONG FORUMULA!!!
            if hashState(stateRep) in dictionary:
                result[dictionary[hashState(stateRep)]]=sqrt(float(stateRep[i]*stateRep[(i+1)%len(stateRep)]))
            stateRep[i]+=1
            stateRep[(i+1)%len(stateRep)]-=1
        if stateRep[(i+1)%len(stateRep)]!=0:
            stateRep[i]+=1
            stateRep[(i+1)%len(stateRep)]-=1
            # THE FOLLOWING IS THE WRONG FORUMULA!!!
            if hashState(stateRep) in dictionary:
                result[dictionary[hashState(stateRep)]]=sqrt(float(stateRep[i]*stateRep[(i+1)%len(stateRep)]))
            stateRep[i]-=1
            stateRep[(i+1)%len(stateRep)]+=1
    return csc_matrix(result)

print "Setting up the vectors of the problem..."

state=BHState(3,3)
i=0 # gives the indices and ultimately the size of the problem
allStates={hashState(state.rep):i}
while state.next():
    # add to allStates; set condition here
    if True:
        i+=1
        allStates.update({hashState(state.rep):i})

print "... Done. Setting up matrices..."

H_kin=csc_matrix(H_kin_row(state.rep, allStates))
H_int=np.array([H_int_ofState(state.rep)])
j=0
while state.next():
    print j
    j+=1
    if True: # same as above
        H_kin=vstack((H_kin, H_kin_row(state.rep, allStates)))
        H_int=np.append(H_int, H_int_ofState(state.rep))

H_i=dia_matrix(([H_int],[0]),shape=(len(H_int),len(H_int)))
print H_i.toarray()
print H_kin.toarray()
print "...Done."
#H_i=dia_matrix(([H_int],[0]),shape=(len(H_int),len(H_int)))
#print H_i.toarray()