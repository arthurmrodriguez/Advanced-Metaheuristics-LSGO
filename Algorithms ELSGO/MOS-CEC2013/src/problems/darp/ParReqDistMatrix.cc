#include "ParReqDistMatrix.h"

static const int SERDATA_STARTPOS = 3;

ParReqDistMatrix::ParReqDistMatrix(vector<Request>& reqs, CostMatrix& costmatrix, DARPEvaluator& eval,
                                   double TWVWeight, double rideVWeight, double loadVWeight, bool useconstrpen,
                                   bool usewaitingpen, long waitingPenThreshold, double waitingPenConstant,
                                   CommManager& comm) :
                                     ReqDistMatrix(costmatrix,eval,TWVWeight,rideVWeight,loadVWeight,useconstrpen,usewaitingpen,waitingPenThreshold,waitingPenConstant) {
  setUpDistMatrix(reqs);
  initializeDistMatrix(reqs,comm);
}

void ParReqDistMatrix::initializeDistMatrix(vector<Request>& reqs, CommManager& comm) {
  assert(comm.getNumIslands() > 1);

  // First, we set the appropriate size of the matrix
  for (int i=0; i<reqs.size(); i++) {
    if (i>0) resize(i,i);
    reqid2matrixpos_[reqs[i].pickup_vert->id_] = i;
  }

  vector<double> computed_data = computeCorrespondingVectorData(reqs, comm);

  parallelfillOfData(comm,computed_data);

//  /*DEBUG*/ if (comm.isIslandMaster()) {
//  /*DEBUG*/   for (int i=0; i<reqs.size(); i++) {
//  /*DEBUG*/     for (int j=0; j<i; j++) {
//  /*DEBUG*/       double exp_value = bestScoreFromAllCombs(reqs[i],reqs[j]);
//  /*DEBUG*/       if (getDist(i,j) != exp_value) {
//  /*DEBUG*/         cout << "diferencias en la fila: " << i << " columna: " << j;
//  /*DEBUG*/         cout << " valor esperado=" << exp_value << " en la matriz=" << getDist(i,j) << endl;
//  /*DEBUG*/         exit(-1);
//  /*DEBUG*/       }
//  /*DEBUG*/     }
//  /*DEBUG*/   }
//  /*DEBUG*/ }
}

inline vector<double> ParReqDistMatrix::computeCorrespondingVectorData(vector<Request>& reqs, CommManager& comm) {
  // Since we are not computing the diagonal the number of elements comes from the smaller triangular matrix of dim -1
  int matrix_nelems = reqs.size()*(reqs.size()-1)/2;

  // Then, we compute the starting row and column as well as the number of elements
  int matrix_firstpos_row,matrix_firstpos_col,vector_nelems;
  computeFirstRowColumnAndNElems(matrix_nelems,comm,matrix_firstpos_row,matrix_firstpos_col,vector_nelems);

  vector<double> computed_data (vector_nelems+3); // we are also storing the row, column, and number of elements
  int datapos = 0;
  computed_data[datapos++] = matrix_firstpos_row;
  computed_data[datapos++] = matrix_firstpos_col;
  computed_data[datapos++] = vector_nelems;       // The number of elements is being saved because the MPI call could
  assert(datapos==SERDATA_STARTPOS);              // use a greater buffer and we need a way to store this number

  int row = matrix_firstpos_row;
  int col = matrix_firstpos_col;
  for (int i=0; i<vector_nelems; i++) {
    computed_data[SERDATA_STARTPOS+i] = bestScoreFromAllCombs(reqs[row],reqs[col]);
    col += 1;
    if (col >= row) {
      row += 1;
      col =  0;
    }
  }

  return computed_data;
}

inline void ParReqDistMatrix::computeFirstRowColumnAndNElems(int matrix_size, CommManager& comm,
                                                             /*out*/ int& matrix_firstpos_row,
                                                             /*out*/ int& matrix_firstpos_col,
                                                             /*out*/ int& nelems){
  int numoppernode = matrix_size / comm.getNumIslands();
  int firstpos = comm.getMyRank() * numoppernode;
  int lastpos  = firstpos + numoppernode - 1;
  assert( comm.getMyRank() != comm.getLastRank() || (matrix_size -1 - lastpos) < comm.getNumIslands() );
  if (comm.getMyRank() == comm.getLastRank() ) lastpos = matrix_size - 1;

  // This formulas (adapted to this problem) have been inspired from
  // http://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients
  // Since we are not using the diagonal, we need to increase the number of rows, i.e.,  firstpos_row = 1 + ( ...
  matrix_firstpos_row = 1 + ( floor(sqrt(8.0*firstpos+1)) - 1.0) / 2;
  matrix_firstpos_col = firstpos - ( (matrix_firstpos_row-1)*matrix_firstpos_row/2);

  nelems = lastpos-firstpos+1;         // could be different than numoppernode if we have the last rank
}


inline void ParReqDistMatrix::parallelfillOfData(CommManager& comm, vector<double>& data) {
  if (comm.isIslandMaster()) {
    fillMatrixWithSerData(data);

    for (int i=1;i<=comm.getLastRank(); i++) {
      // The comm manager (and MPI) needs to know the buffer size of the data that is expecting to receive.
      // Since the last island could have more than nelems (nelems + nislands - 1 at most).
      // For clarity reasons, we always for the worst scenario buffer size
      vector<double> received_data = comm.receiveDoubleVectorFromIsland(i,data.size()+comm.getNumIslands());
      fillMatrixWithSerData(received_data);
    }
  }
  else { //slave node
    comm.sendDoubleVectorToMaster(data);
  }
  GARandomSeed( GAGetRandomNoRankSeed(), comm.getMyRank() ); // Hack to avoid the bug that after calling MPI_Recv the random function gets called
                                                             // and all the random calls loose the synchronization
}

inline void ParReqDistMatrix::fillMatrixWithSerData(vector<double>& data) {
  int row    = data[0];
  int col    = data[1];
  int nelems = data[2]; // We cannot use the data.size because this data has come from a MPI call where the buffer
                        // size (i.e. the data size) could be greater than the actual number of stored values.

  for (int i=0; i<nelems; i++) {
    setDist(row,col,data[SERDATA_STARTPOS+i]);
    col += 1;
    if (col >= row) {
      row += 1;
      col =  0;
    }
  }
}
