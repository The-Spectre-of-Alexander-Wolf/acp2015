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
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace {
	inline bool nonZero(unsigned int n) {return n!=0;}
}

// To use the unordered map (constant-time lookup) with non-standard type vector<unsigned int>,
// we need to declare our own hash class
template<typename T>
struct MyHash;
template<> struct MyHash<std::vector<unsigned int> >
{
	size_t operator()(const std::vector<unsigned int>& v) const {
		size_t t(0);
		// NB: does NOT return double, but is rounded to a size_t!
		return accumulate(v.begin(), v.end(), t, [](unsigned int i, unsigned int j){static double a(0.); return i+=10000000*j*sqrt(100.*(a+=1)+3.);});
	}
};

class BoseHubbardState {
public:
	BoseHubbardState(unsigned int nrSites, unsigned int nrBosons, double U=1., double J=1.);
	~BoseHubbardState();
	bool next();
private:
	const unsigned int nrSites, nrBosons;
public:
	std::unordered_map<std::vector<unsigned int>, int, MyHash<std::vector<unsigned int> > > allStates;
	std::vector<unsigned int> rep; // representation of state; only public for testing
	std::vector<int> H_int_irow, H_int_pcol, H_kin_irow, H_kin_pcol;
	std::vector<double> H_int_A, H_kin_A; // data in matrices
private:
	std::vector<unsigned int>::reverse_iterator kCurr;

	void setInitialState();
	bool nextIndex();
	unsigned int getNrBosonsFromKCurr() const;
	void setTailToZero();
};

#endif /* BOSEHUBBARDSTATE_H_ */
