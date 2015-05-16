//
// Created by andreas on 14.05.15.
//

#ifndef ACP_METROPOLISSTATE_H
#define ACP_METROPOLISSTATE_H

#include <vector>

/*!
 * class MetropolisState represents a spin configuration of the Ising model
 */
class MetropolisState
{
private:
    typedef std::vector<std::vector<bool> >::size_type xSizeType;
    typedef std::vector<bool>::size_type ySizeType;
public:
    //! initialize (rectangular) spin configuration with all spins up/1
    explicit MetropolisState(double K, xSizeType sizeX=16, ySizeType sizeY=0);
    ~MetropolisState(void) {}
    void step(unsigned int nSteps);
    void step(void);
    void show(void);
private:
    const double K;
    std::vector<std::vector<bool> > rep; //<! the representation; vector<bool> is overloaded for optimization!
    double getProbablity(xSizeType x, xSizeType y);
    int sumAdjacentSpins(xSizeType x, xSizeType y);
};


#endif //ACP_METROPOLISSTATE_H
