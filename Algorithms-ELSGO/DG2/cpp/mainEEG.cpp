//
//  main.cpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Modified by Arthur Rodríguez on 8/23/18
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Function.hpp"
#include "DG2.hpp"
#include "CEC2013.hpp"

// This main is specifically designed for EEG problem
// Proper changes over CEC2013 were applied
int main(int argc, char **argv) {

    char filename[100];

    if(argc < 2){
        cout<<"Error: wrong use. Try ./dg2 <problem name>" <<endl;
        cout<<"Problem name: D4 D4N D12 D12N D19 D19N" <<endl;
        return -1;
    }
    else if(argc >= 2){
        for(int i = 1; i < argc; i++){
            string problem = argv[i];
            cout<<"Function: "<<problem<<endl;
            CEC2013 f(problem);
            DG2 dg2(&f);
            dg2.ism();
            dg2.dsm();
            sprintf(filename, "./results/theta%02d.txt", f.getFP()->getID());
            dg2.save(filename);
            cout<<"FEs: "<<dg2.getNumEvaluations()<<endl;
        }

    }

}
