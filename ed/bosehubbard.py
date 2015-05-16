__author__ = 'andreas'

import numpy as np
from scipy.sparse import *
from scipy.sparse.linalg import *
import scipy.special
from math import *
import matplotlib.pylab as plt

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

def getCorrelationMatrix(i, state, dictionary):
    x=[]
    y=[]
    val=[]
    currentIndex=0
    aiaj(0, i, state.rep, dictionary, currentIndex, x, y, val)
    while state.next():
        currentIndex+=1
        aiaj(0, i, state.rep, dictionary, currentIndex, x, y, val)
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

def plotOccupationNumberVarianceAndCondensateFraction(nrSites, nrBosons):
    plotX=[]
    plotOccupationNumberVariance=[]
    plotCondensateFraction=[]
    print "Setting up the vectors of the problem..."

    # nrSites=5
    # nrBosons=nrSites

    state = BHState(nrSites, nrBosons)
    allStates = getBasisAndExpectationValues(state)

    print "...Done. Setting up matrices..."

    H_int = [H_int_ofState(state.rep)] + [0] * (len(allStates)-1)
    H_kin_x = []
    H_kin_y = []
    H_kin_val = []

    j = 0
    H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
    while state.next():
        j += 1
        if True:  # impose a condition, such as a maximum occupation number
            H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
            H_int[j] = H_int_ofState(state.rep)

    H_i = dia_matrix(([H_int], [0]), shape=(len(H_int), len(H_int)))
    H_kin=coo_matrix((H_kin_val,(H_kin_x, H_kin_y)))
    plt.spy(H_kin)
    plt.show()
    # print "...Done. Getting eigenvalues..."
    # Add Parameters U, J here
    J=1.
    UoverJ=0.
    print "Calculating for U/J="
    while UoverJ<20:
        print "...", UoverJ
        plotX+=[UoverJ]
        U=UoverJ*J # J=1, but this is for correctness...
        vals, vecs = eigsh(U/2. * H_i -J*H_kin, which='SA', k=3)
    # print "...Done! Eigenvalues are:"
    # print vals
        gs=np.transpose(vecs)[0]
    # print gs
    # print "Occupation number variance (groundstate):"
        O_n=[]
        i=nrBosons
        while i>=0:
            O_n=O_n+[i]*dimension(nrBosons-i, nrSites-1)
            i-=1
        N_0=dia_matrix(([O_n],[0]), shape=(len(O_n), len(O_n)))
        N2_0=dia_matrix(([np.power(np.array(O_n),2)],[0]), shape=(len(O_n), len(O_n)))
        occupationNumberVariance=sqrt(np.transpose(gs).dot(N2_0.dot(gs))-pow(np.transpose(gs).dot(N_0.dot(gs)),2))
        plotOccupationNumberVariance+=[occupationNumberVariance]
    # print occupationNumberVariance
    # print "Computing SPDM..."
        SPDM=np.transpose(gs).dot(N_0.dot(gs))*np.identity(nrSites)
        for j in range(1, nrSites/2+1): # a_0^*a_0 is <N_0>, see above
            # get matrix a_0^*a_j
            mat=getCorrelationMatrix(j, state, allStates)
            temp = np.transpose(gs).dot(mat.dot(gs))
            for a in range(nrSites):
                SPDM[a, (a+j)%nrSites]=temp
                SPDM[a, (a-j)%nrSites]=temp
    # print "...Done. Eigenvalues of SPDM:"
        SPDMeig = np.linalg.eigvalsh(SPDM)
        plotCondensateFraction+=[np.nanmax(SPDMeig)/nrBosons]
        UoverJ+=1
        # print SPDM
        # print SPDMeig
    plt.plot(plotX, plotOccupationNumberVariance, "b-", label="Occupation Number Variance")
    plt.plot(plotX, plotCondensateFraction, "r-", label="Condensate Fraction")
    plt.legend()
    plt.show()

def plotOccupationNumber(nrSites):
    nrBosons=nrSites
    # plot probabilities of occupation number of site 0
    print "Setting up the vectors of the problem..."
    state = BHState(nrSites, nrBosons)
    allStates = getBasisAndExpectationValues(state)

    print "...Done. Setting up matrices..."

    H_int = [H_int_ofState(state.rep)] + [0] * (len(allStates)-1)
    H_kin_x = []
    H_kin_y = []
    H_kin_val = []
    H_NN=[[]]*(nrBosons+1) # need a list of operators for occupation number

    j = 0
    H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
    while state.next():
        j += 1
        if True:  # impose a condition, such as a maximum occupation number
            H_kin_row(state.rep, allStates, j, H_kin_x, H_kin_y, H_kin_val)
            H_int[j] = H_int_ofState(state.rep)
            for i in range(len(H_NN)):
                if i==state.rep[0]:
                    H_NN[i]=np.append(H_NN[i], [1])
                else:
                    H_NN[i]=np.append(H_NN[i], [0])
    H_i = dia_matrix(([H_int], [0]), shape=(len(H_int), len(H_int)))
    H_kin=coo_matrix((H_kin_val,(H_kin_x, H_kin_y)))
    # plt.spy(H_kin)
    # plt.show()
    print "U/J<<0..."
    U=.1
    J=1.
    vals, vecs = eigsh(U/2. * H_i -J*H_kin, which='SA', k=3)
    gs=np.transpose(vecs)[0]
    plotX=range(len(H_NN))
    plotY=[]
    for i in plotX:
        tempMat=dia_matrix((H_NN[i],[0]), shape=(len(gs),len(gs)))
        plotY+=[np.transpose(gs).dot(tempMat.dot(gs))]
    plt.plot(plotX, plotY, label="U<<J")
    print "U/J>>0..."
    U=1.
    J=.1
    vals, vecs = eigsh(U/2. * H_i -J*H_kin, which='SA', k=3)
    gs=np.transpose(vecs)[0]
    plotY=[]
    for i in plotX:
        tempMat=dia_matrix((H_NN[i],[0]), shape=(len(gs),len(gs)))
        plotY+=[np.transpose(gs).dot(tempMat.dot(gs))]
    plt.plot(plotX, plotY, label="U>>J")
    plt.legend()
    plt.show()

def plotMomentumDistribution(nrSites, UoverJ):
    nrBosons=nrSites  # unit filling
    J=1.
    U=UoverJ*J  # J=1., but for clearity

# plotOccupationNumberVarianceAndCondensateFraction(7, 7)
plotOccupationNumber(7)