#include <iostream>
#include <stdlib.h>

#include "CEC2013LSGOWrapper.h"

int main(int argc, char **argv) {
    long double sol[1000];
    long double f;

    for (int i = 0; i < 1000; i++)
        sol[i] = 0;

    Initialize_CEC2013_LSGO_Wrapper(argv[1]);

    int nfunc = atoi(argv[2]);
    int D = 1000;

    if (nfunc == 13) D = 905;

    bench_func(sol, &f, nfunc, D);

    std::cout << "Resultado: " << f << std::endl;

    Terminate_CEC2013_LSGO_Wrapper();

    return 0;
}
