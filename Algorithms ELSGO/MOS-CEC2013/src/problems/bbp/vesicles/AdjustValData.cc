#include <iostream>

#include "Parser.h"

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Usage: ";
    std::cout << argv[0] << " input_CSV input_image output_CSV" << std::endl;
    return -1;
  }

  Parser parser(argv[1], argv[2], (char*)0);
  parser.parseValData();
  parser.bestFitVesicles(argv[3]);

  return 0;
}
