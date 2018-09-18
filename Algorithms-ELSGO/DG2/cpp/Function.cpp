//
//  Function.cpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#include "Function.hpp"

int Function::getDim(){
    return dim;
}

double Function::getLb(){
    return lb;
}

double Function::getUb(){
    return ub;
}

void Function::setDim(int d){
    dim = d;
}

void Function::setLb(double lowerBound){
    lb = lowerBound;
}

void Function::setUb(double upperBound){
    ub = upperBound;
}

std::string Function::getName(){
    return name;
}


