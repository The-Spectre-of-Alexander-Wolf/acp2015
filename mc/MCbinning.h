//
// Created by andreas on 16.05.15.
//

#ifndef ACP_MCBINNING_H
#define ACP_MCBINNING_H

#include <deque>
#include <cmath>
#include <algorithm>

template<typename T>
class MCbinning {
public:
    MCbinning(unsigned int maxLevel=10):
            values(maxLevel, T()),
            values2(maxLevel, T()),
            variances(maxLevel, T()),
            mean(0.),
            counter(0) {}

    void add(const T& val)
    {
        for_each(++values.begin(), values.end(), [&](T& x){x+=val;});
        // update mean value: mean=(mean*counter+val)/(counter+1)
        mean*=counter;
        mean+=val;
        mean/=++counter;
        // update the bins
        auto i(counter);
        values2[0]=(values2[0]*(counter-1)+val*val)/counter;
        for (int j(1); j!=values.size(); ++j) {
            if (i%2==0) {
                i/=2;
                values2[j]=(values2[j]*(i-1)+values[j]*values[j]*pow(2,-2*j))/i;
                values[j]=0;
            }
            else {break;}
        }
    }

    bool evaluate(void)
    {
        transform(values2.cbegin(), values2.cend(), variances.rbegin(), [&](const T& x){return x-(mean*mean); std::cout << x << std::endl;});

        return is_sorted(variances.cbegin(), variances.cend());
    }
//private:
    std::deque<T> values, values2, variances;
    double mean;
    long long counter;
};


#endif //ACP_MCBINNING_H
