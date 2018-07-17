#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include <armadillo>
#include <matio.h>

arma::mat  loadMatrix (mat_t* matfp, char* varname);
arma::uvec loadIndices(mat_t* matfp, char* varname);

#endif
