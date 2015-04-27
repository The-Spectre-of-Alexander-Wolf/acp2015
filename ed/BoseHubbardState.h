/*
 * BoseHubbardState.h
 *
 *  Created on: 26.04.2015
 *      Author: andreas
 */
#ifndef BOSEHUBBARDSTATE_H_
#define BOSEHUBBARDSTATE_H_

#include <vector>
#include <map>
#include <armadillo>

namespace {
	inline bool nonZero(unsigned int n) {return n!=0;}
}

class BoseHubbardState {
public:
	BoseHubbardState(unsigned int nrSites, unsigned int nrBosons, double U=1., double J=1.);
	~BoseHubbardState();
	bool next();
	void addToH_int(unsigned long i);
	void addToH_kin(unsigned long i);
private:
	const unsigned int nrSites, nrBosons;
	const double Uover2, minusJ;
public:
	std::map<std::vector<unsigned int>, long int> allStates;
	std::vector<unsigned int> rep; // representation of state; only public for testing
	arma::sp_mat H_int, H_kin;
private:
	std::vector<unsigned int>::reverse_iterator kCurr;

	void setInitialState();
	bool nextIndex();
	unsigned int getNrBosonsFromKCurr() const;
	void setTailToZero();
};

#endif /* BOSEHUBBARDSTATE_H_ */
