//  Modified by Arthur Rodríguez on 8/23/18
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.


#include "cec2013/Header.h"
#include "Function.hpp"

#ifndef CEC2013_H
#define CEC2013_H

class CEC2013: public Function{
    public:
        double evaluate(double *x);
        CEC2013(int);
        // Constructor for EEG problem
        CEC2013(string);
        Benchmarks * getFP();
    private:
        Benchmarks *fp;
};

#endif
