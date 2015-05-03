__author__ = 'andreas'

import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import *
import scipy.special
from math import *

def findLowestNZIndex(v):
    i = len(v) - 1
    while i >= 0 and v[i] == 0:
        i -= 1
    return i

def H_int_ofState(v):
    result = 0
    for n in v:
        result += n * (n - 1)
    return result

class BHState:
    rep = np.array([])

    def __init__(self, dim, nB):
        self.rep = np.array([nB] + [0] * (dim - 1))  # initial state

    def next(self):
        i = findLowestNZIndex(self.rep[0:len(self.rep) - 1])
        if i < 0:  # cannot increment any further
            self.rep = np.array([np.sum(self.rep)] + [0] * (len(self.rep) - 1))  # initial state
            return False
        self.rep[i] -= 1
        self.rep[i + 1:] = np.array([np.sum(self.rep[i + 1:]) + 1] + [0] * (len(self.rep[i + 1:]) - 1))
        return True

def hashState(v):
    result = 0.
    for i in range(len(v)):
        result += sqrt(100. * i + 3) * v[i]
    return result

def aiaj(i, j, stateRep, dictionary, currentIndex, x, y, val): # currentIndex is for computational convenience only
    if stateRep[i] != 0:
        # apply H_kin to state
        stateRep[i] -= 1
        stateRep[j] += 1
        # find state in dictionary and compute matrix element
        if hashState(stateRep) in dictionary:
            x += [currentIndex]
            y += [dictionary[hashState(stateRep)]]
            val += [sqrt(float((stateRep[i] + 1) * stateRep[j]))]
        # return to former state
        stateRep[i] += 1
        stateRep[j] -= 1

def H_kin_row(stateRep, dictionary, currentIndex, x, y, val):
    for i in range(len(stateRep)):
        aiaj(i, (i+1)%len(stateRep), stateRep, dictionary, currentIndex, x, y, val)
        aiaj((i+1)%len(stateRep), i, stateRep, dictionary, currentIndex, x, y, val)
        # if stateRep[i] != 0:
        #     # apply H_kin to state
        #     stateRep[i] -= 1
        #     stateRep[(i + 1) % len(stateRep)] += 1
        #     # find state in dictionary and compute matrix element
        #     if hashState(stateRep) in dictionary:
        #         x += [currentIndex]
        #         y += [dictionary[hashState(stateRep)]]
        #         val += [sqrt(float((stateRep[i] + 1) * stateRep[(i + 1) % len(stateRep)]))]
        #     # return to former state
        #     stateRep[i] += 1
        #     stateRep[(i + 1) % len(stateRep)] -= 1
        # # repeat this for the h.c.
        # if stateRep[(i + 1)%len(stateRep)] != 0:
        #     stateRep[i] += 1
        #     stateRep[(i + 1) % len(stateRep)] -= 1
        #     if hashState(stateRep) in dictionary:
        #         x += [currentIndex]
        #         y += [dictionary[hashState(stateRep)]]
        #         val += [sqrt(float(stateRep[i] * (stateRep[(i + 1) % len(stateRep)] + 1)))]
        #     stateRep[i] -= 1
        #     stateRep[(i + 1) % len(stateRep)] += 1

def dimension(N, M):  # problem size for N bosons on M lattice sites
    return scipy.special.binom(N + M - 1, M - 1)


def getBasisAndExpectationValues(state, operators=[], conditional=True):
    i=0
    allStates={hashState(state.rep):i}
    expectationValues=[[]*len(operators)]
    while state.next():
        if conditional:
            i+=1
            allStates.update({hashState(state.rep):i})
            for j in range(len(operators)):
                expectationValues[i]+=[operators[i](state.rep)]
    if expectationValues[0]:
        return allStates, expectationValues
    else:
        return allStates

def getCorrelationMatrix(i, stateRep, dictionary):
    x=[]
    y=[]
    val=[]
    currentIndex=0
    aiaj(0, i, stateRep, dictionary, currentIndex, x, y, val)
    while state.next():
        currentIndex+=1
        aiaj(0, i, stateRep, dictionary, currentIndex, x, y, val)
    return coo_matrix((val,(x,y)), shape=(len(dictionary),len(dictionary)))

# def getIndex(v):
#     M = len(v)
#     if M == 0:
#         return 0
#     N = sum(v)
#     n0 = N
#     result = 0
#     while n0 > v[0]:
#         result += dimension(N - n0, M - 1)
#         n0 -= 1
#     return int(result + getIndex(v[1:]))

print "Setting up the vectors of the problem..."

nrSites=10
nrBosons=nrSites

state = BHState(nrSites, nrBosons)
allStates = getBasisAndExpectationValues(state)

print "...Done. Setting up matrices..."

H_int = [H_int_ofState(state.rep)] + [0] * (len(allStates)-1)
H_int_x = []
H_kin_x = []
H_kin_y = []
H_kin_val = []
matsForSPDM = [[]*(nrSites/2+1)]

j = 0
H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
while state.next():
    j += 1
    if True:  # same as above
        H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
        H_int[j] = H_int_ofState(state.rep)

H_i = dia_matrix(([H_int], [0]), shape=(len(H_int), len(H_int)))
H_kin=coo_matrix((H_kin_val,(H_kin_x, H_kin_y)))
print "...Done. Getting eigenvalues..."
# Add Parameters U, J here
U=1.
J=1.
vals, vecs = eigsh(U/2. * H_i -J*H_kin)
print "...Done! Eigenvalues are:"
print vals
gs=np.transpose(vecs)[0]
print "Occupation number variance (groundstate):"
O_n=[]
i=nrBosons
while i>=0:
    O_n=O_n+[i]*dimension(nrBosons-i, nrSites-1)
    i-=1
N_0=dia_matrix(([O_n],[0]), shape=(len(O_n), len(O_n)))
N2_0=dia_matrix(([np.power(np.array(O_n),2)],[0]), shape=(len(O_n), len(O_n)))
occupationNumberVariance=sqrt(np.transpose(gs).dot(N2_0.dot(gs))-pow(np.transpose(gs).dot(N_0.dot(gs)),2))
print occupationNumberVariance
print "Computing SPDM..."
SPDM=np.transpose(gs).dot(N_0.dot(gs))*np.identity(nrSites)
for j in range(1, nrSites/2+1): # a_0^*a_0 is <N_0>, see above
    # get matrix a_0^*a_j
    mat=getCorrelationMatrix(j, state.rep, allStates)
    temp = np.transpose(gs).dot(mat.dot(gs))
    for a in range(nrSites):
        SPDM[a, (a+j)%nrSites]=temp
        SPDM[a, (a-j)%nrSites]=temp
print "...Done. Eigenvalues of SPDM:"
SPDMeig = np.linalg.eigvalsh(SPDM)
print SPDMeig
#print np.amax(SPDMeig)