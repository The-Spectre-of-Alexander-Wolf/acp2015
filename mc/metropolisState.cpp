//
// Created by andreas on 14.05.15.
//

#include "metropolisState.h"

#include <random>
#include <iostream>

using namespace std;

MetropolisState::MetropolisState(double K, xSizeType sizeX, ySizeType sizeY)
        : K(K),
          rep(vector<vector<bool> >(sizeX, vector<bool>(sizeY>0?sizeY:sizeX, true))) {}



void MetropolisState::step(unsigned int nSteps)
{
    for (unsigned int i(0); i!=nSteps; ++i) {step();}
}


void MetropolisState::step(void)
{
    // pick random site
    static auto seed(1);
    static default_random_engine re(seed);
    static uniform_int_distribution<xSizeType> uniXDist(0, rep.size()-1);
    static uniform_int_distribution<ySizeType> uniYDist(0, rep[0].size()-1);
    auto xIndex(uniXDist(re));
    auto yIndex(uniYDist(re));
    // check if we take the step
    static uniform_real_distribution<double> rng(0.,1.);
    double prob(getProbablity(xIndex, yIndex));
    if (prob>=1 or prob>rng(re)) {rep[xIndex][yIndex]=not rep[xIndex][yIndex];}
}


double MetropolisState::getProbablity(xSizeType x, xSizeType y)
{
    // CONSIDER MEMOIZATION!
    // there are only three cases: 1., exp(-4*K), exp(-8*K)
    int sum(sumAdjacentSpins(x, y));
    return exp(rep[x][y]?-K*sum:K*sum); // INEFFICIENT! THE FACTORS ARE ALWAYS THE SAME!
}


int MetropolisState::sumAdjacentSpins(xSizeType x, xSizeType y)
{
    int neighboursCount(0);
    for (int i(-1); i!=1; i+=2) {
        for (int j(-1); j!=1; j+=2) {
            // attention: periodic boundary conditions
            if (rep[(x+i+rep.size())%rep.size()][(y+j+rep[0].size())%rep[0].size()]) {
                neighboursCount+=1;
            }
            else {
                neighboursCount-=1;
            }
        }
    }
    return neighboursCount;
}


void MetropolisState::show(void)
{
    for (auto it(rep.cbegin()); it!=rep.cend(); ++it) {
        for (auto jt(it->cbegin()); jt!=it->cend(); ++jt) {
            *jt?cout << "O":cout << "X";
        }
        cout << endl;
    }
}