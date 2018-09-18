//
//  DG2.cpp
//  DG2
//
//  Created by Keyhan Kouhkiloui on 2/24/17.
//  Copyright © 2017 Keyhan Kouhkiloui. All rights reserved.
//

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include "Function.hpp"
#include "DG2.hpp"

using namespace std;

DG2::DG2 (Function *f): g(f->getDim()-1){

    this->func = f;
    int dim = func->getDim();

    boost::add_vertex (g);

    delta1    = new double* [dim];
    delta2    = new double* [dim];
    lambda    = new double* [dim];
    theta     = new double* [dim];
    epsilon   = new double* [dim];
    archive_f = new double* [dim];
    for (int i = 0; i < dim; i++){
        delta1[i]    = new double [dim];
        delta2[i]    = new double [dim];
        lambda[i]    = new double [dim];
        theta[i]     = new double [dim];
        epsilon[i]   = new double [dim];
        archive_f[i] = new double [dim];
    }

    archive_fhat = new double [dim];

    eval_count = 0;
    archive_base = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < dim; i++){
        archive_fhat[i] = std::numeric_limits<double>::quiet_NaN();
        for (int j = 0; j < dim; j++){
            delta1[i][j] = std::numeric_limits<double>::quiet_NaN();
            delta2[i][j] = std::numeric_limits<double>::quiet_NaN();
            lambda[i][j] = std::numeric_limits<double>::quiet_NaN();
            epsilon[i][j] = std::numeric_limits<double>::quiet_NaN();
            archive_f[i][j] = std::numeric_limits<double>::quiet_NaN();
            if (i == j){
                theta[i][j] = 1;
            }
            else{
                theta[i][j] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}

void DG2::arrayCopy(double *dest, double *src, int size){
    for (int i = 0; i < size; i++){
        dest[i] = src[i];
    }
}

void DG2::ism (void){

    double fp1, fp2, fp3, fp4, d1, d2;
    double temp = (func->getLb()+func->getUb())/2.0;
    int dim = func->getDim();

    double *p1 = new double [dim];
    double *p2 = new double [dim];
    double *p3 = new double [dim];
    double *p4 = new double [dim];

    for (int i = 0; i < dim; i++){
        p1[i] = func->getLb();
    }

    fp1 = func->evaluate(p1);
    archive_base = fp1;
    eval_count += 1;

    int counter = 0;
    int prev = 0;
    int prog = 0;
    for (int i = 0; i < dim - 1; i++){
        if (!std::isnan(archive_fhat[i])){
            fp2 = archive_fhat[i];
        }else{
            arrayCopy(p2, p1, func->getDim());
            p2[i] = temp;
            fp2 = func->evaluate(p2);
            eval_count += 1;
            archive_fhat[i] = fp2;

        }

        for (int j = i + 1; j < func->getDim(); j++){
            counter++;
            prev = prog;
            prog = (int)(counter/(float)(func->getDim()*(func->getDim()-1))*2*100);
            if (prog % 5 == 0 && prev != prog){
                printf("Progress = %02d\r", prog);
                std::cout.flush();
            }

            if (!std::isnan(archive_fhat[j])){
                fp3 = archive_fhat[j];
            }else {
                arrayCopy(p3, p1, func->getDim());
                p3[j] = temp;
                fp3 = func->evaluate(p3);
                eval_count += 1;
                archive_fhat[j] = fp3;
            }

            arrayCopy(p4, p1, func->getDim());
            p4[i] = temp;
            p4[j] = temp;
            fp4 = func->evaluate(p4);
            eval_count += 1;
            archive_f[i][j] = fp4;
            archive_f[j][i] = fp4;

            d1 = fp2 - fp1;
            d2 = fp4 - fp3;

            delta1[i][j] = d1;
            delta2[i][j] = d2;
            lambda[i][j] = fabs(d1 - d2);
            lambda[j][i] = fabs(d1 - d2);
        }
    }

    delete [] p1;
    delete [] p2;
    delete [] p3;
    delete [] p4;

};

double DG2::gammaFunc(double d){
    double muM = (std::numeric_limits<double>::epsilon())/2.0;
    return (d * muM)/(1 - (d * muM));

}

void DG2::dg(double epsilon){
    for (int i = 0 ; i < func->getDim() - 1; i++){
        for (int j = i + 1; j < func->getDim(); j++){
            if (lambda[i][j] < epsilon){
                theta[i][j] = 0;
            }else if (lambda[i][j] > epsilon){
                theta[i][j] = 1;
            }
        }
    }
}

void DG2::dsm (void){

    double array [4] = {};
    double fMax = 0;
    double eInf = 0, eSup = 0;
    double etha0 = 0, etha1 = 0;
    double eps = 0;

    for (int i = 0 ; i < func->getDim() - 1; i++){
        for (int j = i + 1; j < func->getDim(); j++){
            array[0] = archive_base;
            array[1] = archive_f[i][j];
            array[2] = archive_fhat[i];
            array[3] = archive_fhat[j];
            fMax = *std::max_element(array, array + 4);
            eInf = gammaFunc(2) * std::max(array[0] + array[1], array[2] + array[3]);
            //eSup = gammaFunc(std::sqrt(func->getDim())) * fMax;
            eSup = gammaFunc(pow((double)func->getDim(), (double)0.5)) * fMax;
            if (lambda[i][j] <= eInf){
                theta[i][j] = 0;
                theta[j][i] = 0;
                etha0 += 1;
            }
            else if (lambda[i][j] >= eSup){
                theta[i][j] = 1;
                theta[j][i] = 1;
                boost::add_edge (i, j, g);
                etha1 += 1;
            }
        }
    }

    for (int i = 0 ; i < func->getDim() - 1; i++){
        for (int j = i + 1; j < func->getDim(); j++){
            array[0] = archive_base;
            array[1] = archive_f[i][j];
            array[2] = archive_fhat[i];
            array[3] = archive_fhat[j];
            fMax = *std::max_element(array, array + 4);
            eInf = gammaFunc(2) * std::max(array[0] + array[1], array[2] + array[3]);
            //eSup = gammaFunc(std::sqrt(func->getDim())) * fMax;
            eSup = gammaFunc(pow((double)func->getDim(), (double)0.5)) * fMax;

            if (std::isnan(theta[i][j])){
                eps = (etha0 / (etha0 + etha1)) * eInf + (etha1 / (etha0 + etha1)) * eSup;
                epsilon[i][j] = eps;
                epsilon[j][i] = eps;
                if (lambda[i][j] > eps){
                    theta[i][j] = 1;
                    theta[j][i] = 1;
                    boost::add_edge (i, j, g);
                }
                else{
                    theta[i][j] = 0;
                    theta[j][i] = 0;
                }
            }

        }
    }

    std::vector<int> component (boost::num_vertices (g));
    size_t num_components = boost::connected_components (g, &component[0]);
    std::cout <<"number of components: " << num_components << std::endl;

    for (int j = 0 ; j < num_components ; j++){
        std::vector<int> temp_group;
        for (size_t i = 0; i < boost::num_vertices (g); ++i){
            if (component[i] == j){
                temp_group.push_back(i);
                std::cout << i << " ";
            }
        }
        std::cout << std::endl;
        if(temp_group.size() == 1){
            seps.push_back(temp_group.back());
        }else{
            nonseps.push_back(temp_group);
        }
    }
};

bool DG2::save (std::string filename){
    //TODO
    std::ofstream fh;
    fh.open (filename.c_str());
    for (int i = 0; i < func->getDim(); i++){
        for (int j = 0; j < func->getDim(); j++){
            fh << theta[i][j] << ((j == func->getDim()-1) ? "" : ",");
        }

        fh << std::endl;
    }
    fh.close();
    return true;
};

bool DG2::load (std::string){
    //TODO
    return true;
};


int DG2::getNumEvaluations(){
    return eval_count;
};


DG2::~DG2(void){
    int dim = func->getDim();
    for (int i = 0; i < dim; i++){
        delete [] delta1[i]    ;
        delete [] delta2[i]    ;
        delete [] lambda[i]    ;
        delete [] theta[i]     ;
        delete [] epsilon[i]   ;
        delete [] archive_f[i] ;
    }
    delete [] delta1    ;
    delete [] delta2    ;
    delete [] lambda    ;
    delete [] theta     ;
    delete [] epsilon   ;
    delete [] archive_f ;
    delete [] archive_fhat;
};
