/*
 * BoseHubbardState.cpp
 *
 *  Created on: 26.04.2015
 *      Author: andreas
 */

#include "BoseHubbardState.h"

#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

namespace {
	unsigned int addH_int_atomic(unsigned int i, unsigned int j) {return i+j*(j-1);}
}

BoseHubbardState::BoseHubbardState(unsigned int nrSites, unsigned int nrBosons, double U, double J):
nrSites(nrSites),
nrBosons(nrBosons),
rep(nrSites, 0),
kCurr(--rep.rend())
{
	setInitialState();
}

BoseHubbardState::~BoseHubbardState() {
	// TODO Auto-generated destructor stub
}

bool BoseHubbardState::next() {
	if(nextIndex()) {
		--*kCurr; // decrement n_k
		(*--kCurr)=getNrBosonsFromKCurr()+1; //Attention! Left-hand side should be evaluated first!
		setTailToZero();
		return true;
	}
	else {return false;}
}

//void BoseHubbardState::addToH_int(unsigned int i) {
//	H_int(i,i)=Uover2*accumulate(rep.begin(), rep.end(), 0, [](unsigned int i, unsigned int j){return i+j*(j-1);});
//}

//void BoseHubbardState::addToH_kin(unsigned int i) {
//	vector<unsigned int> tempState(rep);
//	map<vector<unsigned int>, long int>::const_iterator it;
//	for (unsigned int j(0); j!=tempState.size(); ++j) {
//		if (tempState[j]!=0) {
//			--tempState[j]; ++tempState[(j+1)%tempState.size()];
//			it = allStates.find(tempState);
//			if (it!=allStates.end()) {
//				H_kin(i, it->second)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
//				H_kin(it->second, i)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
//			}
//		}
//		if (tempState[(j+1)%tempState.size()]!=1) {
//			tempState[j]+=2; tempState[(j+1)%tempState.size()]-=2;
//			it = allStates.find(tempState);
//			if (it!=allStates.end()) {
//				H_kin(i, it->second)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
//				H_kin(it->second, i)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
//			}
//		}
//	}
//}


void BoseHubbardState::setInitialState() {
	rep[0]=nrBosons;
}

bool BoseHubbardState::nextIndex() { // returns false if invoked on the last state
	kCurr=find_if(++rep.rbegin(), rep.rend(), nonZero); // find first non-zero from back
	return kCurr!=rep.rend();
}

unsigned int BoseHubbardState::getNrBosonsFromKCurr() const {
	// I'm sorry, this line looks horrible:
	return accumulate<vector<unsigned int>::const_reverse_iterator, unsigned int>(rep.rbegin(), kCurr, 0)+*kCurr;
}

void BoseHubbardState::setTailToZero(){
	replace_if(rep.rbegin(), kCurr, nonZero, 0);
}
