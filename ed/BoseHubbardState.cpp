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
#include <armadillo>

#ifndef ARMA_64BIT_WORD
#define ARMA_64BIT_WORD
#endif

using namespace std;
using namespace arma;

namespace {
	unsigned int addH_int_atomic(unsigned int i, unsigned int j) {return i+j*(j-1);}
}

BoseHubbardState::BoseHubbardState(unsigned int nrSites, unsigned int nrBosons, double U, double J):
nrSites(nrSites),
nrBosons(nrBosons),
Uover2(U/2.),
minusJ(-J),
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

void BoseHubbardState::addToH_int(unsigned long i) {
	H_int(i,i)=Uover2*accumulate(rep.begin(), rep.end(), 0, [](unsigned int i, unsigned int j){return i+j*(j-1);});
}

void BoseHubbardState::addToH_kin(unsigned long i) {
	vector<unsigned int> tempState(rep);
	map<vector<unsigned int>, long int>::const_iterator it;
	for (unsigned int j(0); j!=tempState.size(); ++j) {
		if (tempState[j]!=0) {
			--tempState[j]; ++tempState[(j+1)%tempState.size()];
			it = allStates.find(tempState);
			if (it!=allStates.end()) {
				H_kin(i, it->second)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
				H_kin(it->second, i)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
			}
		}
		if (tempState[(j+1)%tempState.size()]!=1) {
			tempState[j]+=2; tempState[(j+1)%tempState.size()]-=2;
			it = allStates.find(tempState);
			if (it!=allStates.end()) {
				H_kin(i, it->second)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
				H_kin(it->second, i)=minusJ*sqrt((tempState[j]+1)*tempState[(j+1)%tempState.size()]);
			}
		}
	}
}


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

template<unsigned int N>
unsigned long int myFac() {
	return myFac<N-1>()*N;
}
template<>
unsigned long int myFac<2>() {return 2;}
template<>
unsigned long int myFac<1>() {return 1;}
template<>
unsigned long int myFac<0>() {return 1;}

int main(int argc, char const *argv[]){
    //runAllTests(argc,argv);


    cout << "Indexing all basis vectors and computing the Hamiltonian...\n";

    const unsigned int nrSites(7), nrBosons(7);
    const unsigned long int dim(myFac<nrBosons+nrSites-1>()/myFac<nrBosons>()/myFac<nrSites-1>());

//    vector<double> H_int(dim, 0.); // adjust for trimmed Hilbert space
//    vector<pair<vector<unsigned int>, vector<double> > > H_kin(dim-1); // adjust for trimmed Hilbert space

    BoseHubbardState s(nrSites, nrBosons);

    s.H_int=sp_mat(dim, dim);
    s.H_kin=sp_mat(dim, dim);
    long int i(0);
    do { // for all states
    	if (true) {s.addToH_int(i); s.addToH_kin(i);}
//    	H_int[i]=s.H_int();
//    	H_kin[i]=s.H_kin();
    	s.allStates.insert(pair<vector<unsigned int>, long int>(s.rep, i));
    	++i;
    } while (s.next());
    cout << "...done after " << i << " vectors. First eigenvalues:\n";
//    vec eigval=eigs_sym(s.H_int+s.H_kin, 5);
//    cout << eigval << endl;
//    cout << s.H_int << endl << s.H_kin << endl;
}

