//
// Created by andreas on 16.05.15.
//

#ifndef ACP_MCBINNING_H
#define ACP_MCBINNING_H

#include <deque>
#include <cmath>

template<typename T>
class MCbinning {
public:
    MCbinning(): values(1, T()), counter(0) {}
    void add(T val)
    {
        ++counter;
        // update all possible bins
        T oldVal(values[0]);
        std::swap(values[0],val); // recycle val!
        // val & oldVal store the same now
        auto i(counter);
        std::cout << "counter=" << i << " " << ((i%2)!=0?"true":"false") << std::endl;
        for (int j(1); i%2==0; ++j) {
            std::cout << "j=" << j << " i=" << i;
            i/=2;
            oldVal=values[j];
            if (j<values.size()) {
                values[j]=(val+values[j-1])/2;
            }
            else {
                values.push_back((val+values[j-1])/2);
            }
            std::cout << " val: " << values[j] << std::endl;
            std::swap(val, oldVal);
        }
    }

    bool evaluate(void) {return false;}
//private:
    std::deque<T> values;
    long long counter;
};


#endif //ACP_MCBINNING_H
