//
//  main.cpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Function.hpp"
#include "DG2.hpp"
#include "CEC2013.hpp"

int main(int argc, char **argv) {

    char filename[100];
    for (int i = 1 ; i < argc ; i++){
        int funID = atoi(argv[i]);
        cout<<"Function: "<<funID<<endl;
        CEC2013 f(funID);
        DG2 dg2(&f);
        dg2.ism();
        dg2.dsm();
        sprintf(filename, "./results/theta%02d.txt", funID);
        dg2.save(filename);
        cout<<"FEs: "<<dg2.getNumEvaluations()<<endl;
    }
    return 0;
}
