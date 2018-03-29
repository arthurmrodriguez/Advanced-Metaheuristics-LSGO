#include "DistMatrix.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

#include <stdlib.h>

// We create a matrix greater than the actual size so that indices are used directly (and not with -1)
DistMatrix::DistMatrix(std::string distFile, int size) : data_(size+1, std::vector<long double>(size+1, 0.0)) {

  // Open distMatrix file as an input stream
  std::ifstream fin (distFile.c_str());

  // Check if file can be correctly open
  if (fin.fail() or fin.bad()) {
    std::stringstream msg;
    msg << "[DistMatrix] Error: file '" << distFile << "' not found or could not be open" << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Process the whole file
  while (!fin.eof()) {

    // Read one line
    std::string buffer;
    getline(fin, buffer);

    // Check if this line was properly read
    if (!fin.fail()) {

      std::istringstream iss(buffer);
      std::string token;
      std::vector<std::string> tokens;

      while (getline(iss, token, '|'))
        tokens.push_back(token);

      int    idOrig = atoi  (tokens[3].c_str());
      int    idDest = atoi  (tokens[2].c_str());
      long double cost   = strtod(tokens[1].c_str(), NULL);

      assert(idOrig >= 0);
      assert(idDest >= 0);
      assert(idOrig < size);
      assert(idDest < size);

      data_[idOrig][idDest] = cost;

    }

  }

  // Close input file
  fin.close();

}
