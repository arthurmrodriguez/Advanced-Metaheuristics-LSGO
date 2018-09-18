//
//  Function.hpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#ifndef Function_hpp
#define Function_hpp

#include <stdio.h>
#include <iostream>

class Function {
    
public:
    
    virtual double evaluate (double *d) = 0;
    int getDim();
    double getLb();
    double getUb();
    std::string getName();
    void setDim(int d);
    void setLb(double lowerBound);
    void setUb(double upperBound);
    
protected:
    
    int dim;
    double lb;
    double ub;
    std::string name;
};

#endif /* Function_hpp */
