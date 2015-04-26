/*
 * BoseHubbardState.h
 *
 *  Created on: 26.04.2015
 *      Author: andreas
 */
#ifndef BOSEHUBBARDSTATE_H_
#define BOSEHUBBARDSTATE_H_

#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

namespace {
	inline bool nonZero(unsigned int n) {return n!=0;}
}

class BoseHubbardState {
public:
	BoseHubbardState(unsigned int nrSites, unsigned int nrBosons);
	~BoseHubbardState();
	bool next() {
		if(nextIndex()) {
			--*kCurr; // decrement n_k
			(*--kCurr)=getNrBosonsFromKCurr()+1; //Attention! Left-hand side should be evaluated first!
			setTailToZero();
			return true;
		}
		else {return false;}
	}
private:
	const unsigned int nrSites, nrBosons;
public:
	std::vector<unsigned int> rep; // representation of state
private:
	std::vector<unsigned int>::reverse_iterator kCurr;

	bool nextIndex();
	unsigned int getNrBosonsFromKCurr();
	void setTailToZero();
};

#endif /* BOSEHUBBARDSTATE_H_ */
