#include "load_data.h"

#include <stdexcept>

arma::mat loadMatrix(mat_t* matfp, char* varname) {
  matvar_t* var = Mat_VarRead(matfp, varname);
  //Mat_VarPrint(var, 1);

  if (var->rank > 2) {
    std::cerr << "Error: maximum number of dimensions for matrices is 2. Trying to load a " << var->rank << "-D matrix" << std::endl;
    exit(-1);
  }

  arma::mat res=arma::mat(var->dims[0], var->dims[1]);

  int dp = 0;
  int*    data_i=(int*   )var->data;
  long double* data_d=(long double*)var->data;

  for (int j=0; j<var->dims[1]; j++) {
    for (int i=0; i<var->dims[0]; i++) {
      switch(var->data_type) {
        case MAT_T_INT32:
          res(i, j)=data_i[dp++];
          break;
        case MAT_T_DOUBLE:
          res(i, j)=data_d[dp++];
          break;
        default:
          throw std::runtime_error("[Load_Matrix] Error: Unknown type in MAT file.");
          break;
      }
    }
  }

  Mat_VarFree(var);

  return res;
}


arma::uvec loadIndices(mat_t* matfp, char* varname) {
  matvar_t* var = Mat_VarRead(matfp, varname);

  if (var->rank > 2) {
    std::cerr << "Error: maximum number of dimensions for matrices is 2. Trying to load a " << var->rank << "-D matrix" << std::endl;
    exit(-1);
  }

  unsigned dims=(var->dims[0] == 1) ? var->dims[1] : var->dims[0];

  arma::uvec res=arma::uvec(dims);

  int*    data_i=(int*   )var->data;
  long double* data_d=(long double*)var->data;

  bool print_warning=true;

  for (int i=0; i<dims; i++) {
    // We keep this switch in case indices must be of another integral type (UINT, etc.)
    switch(var->data_type) {
      case MAT_T_INT32:
        // Important!!! We must substract 1 to every index as C++ arrays start in zero
        res(i)=data_i[i]-1;
        break;
        case MAT_T_DOUBLE:
        // Important!!! We must substract 1 to every index as C++ arrays start in zero
        res(i)=((int)data_d[i])-1;
        if (print_warning) {
          print_warning=false;
          std::cerr << "[Load_indices] Warning: indices vector contains long doubles instead of ints. Trying to automatically convert..." << std::endl;
        }
        break;
      default:
        throw std::runtime_error("[Load_Indices] Error: Unexpected type in MAT file.");
        break;
    }
  }

  Mat_VarFree(var);

  return res;
}
