//
// Created by andreas on 14.05.15.
//

#include "metropolisState.h"

int main(int, char **)
{
    MetropolisState ms(.4, 16);
    ms.show();
    ms.step();
    ms.show();
    return 0;
}