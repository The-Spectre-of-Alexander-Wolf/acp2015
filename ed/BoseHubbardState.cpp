/*
 * BoseHubbardState.cpp
 *
 *  Created on: 26.04.2015
 *      Author: andreas
 */

#include "BoseHubbardState.h"

using namespace std;

BoseHubbardState::BoseHubbardState(unsigned int nrSites, unsigned int nrBosons):
nrSites(nrSites),
nrBosons(nrBosons),
rep(nrSites, 0),
kCurr(--rep.rend())
{
	rep[0]=nrBosons;
}

BoseHubbardState::~BoseHubbardState() {
	// TODO Auto-generated destructor stub
}

bool BoseHubbardState::nextIndex() { // returns false if invoked on the last state
	kCurr=find_if(++rep.rbegin(), rep.rend(), nonZero); // find first non-zero from back
	return kCurr!=rep.rend();
}

unsigned int BoseHubbardState::getNrBosonsFromKCurr() {
	return accumulate(rep.rbegin(), kCurr, 0)+*kCurr;
}

void BoseHubbardState::setTailToZero() {
	replace_if(rep.rbegin(), kCurr, nonZero, 0);
}

