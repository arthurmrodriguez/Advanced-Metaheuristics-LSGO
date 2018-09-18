//  Modified by Arthur Rodríguez on 8/23/18
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.

#include "Function.hpp"
#include "CEC2013.hpp"
#include "cec2013/Header.h"

CEC2013::CEC2013(int funcID){
    if (funcID==1){
      fp = new F1();
    }else if (funcID==2){
      fp = new F2();
    }else if (funcID==3){
      fp = new F3();
    }else if (funcID==4){
      fp = new F4();
    }else if (funcID==5){
      fp = new F5();
    }else if (funcID==6){
      fp = new F6();
    }else if (funcID==7){
      fp = new F7();
    }else if (funcID==8){
      fp = new F8();
    }else if (funcID==9){
      fp = new F9();
    }else if (funcID==10){
      fp = new F10();
    }else if (funcID==11){
      fp = new F11();
    }else if (funcID==12){
      fp = new F12();
    }else if (funcID==13){
      fp = new F13();
    }else if (funcID==14){
      fp = new F14();
    }else if (funcID==15){
      fp = new F15();
    }else{
      cerr<<"Fail to locate Specified Function Index"<<endl;
      exit(-1);
    }

    lb = fp->getMinX();
    ub = fp->getMaxX();
    dim = fp->getDimension();
}

// Constructor for EEG problem
CEC2013::CEC2013(string problem){

    fp = new EEG(problem);
    lb = fp->getMinX();
    ub = fp->getMaxX();
    dim = fp->getDimension();

}

double CEC2013::evaluate(double *x){
    return fp->compute(x);
}

Benchmarks* CEC2013::getFP(){
    return fp;
}
