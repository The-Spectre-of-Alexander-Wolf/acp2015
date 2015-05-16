//
// Created by andreas on 16.05.15.
//


#include <iostream>
#include "MCbinning.h"

using namespace std;

bool testDoubleConstructor()
{
    MCbinning<double> mc();
    return true;
}


bool testDoubleBinning()
{
    MCbinning<double> mc;
    for (int i(1); i!=5; ++i) {
        cout << "add " << 1.+1./i << endl;
        mc.add(1.+1./i);
    }
    for (auto it(mc.values2.cbegin()); it!=mc.values2.cend(); ++it) {cout << *it << endl;}
    return mc.evaluate();
}

int main(int, char **)
{
    cout << "testDoubleConstructor: ";
    testDoubleConstructor()?cout << "passed\n":cout << "failed\n";
    testDoubleBinning()?(cout << "probably converging\n"):(cout << "probably not converging\n");
}