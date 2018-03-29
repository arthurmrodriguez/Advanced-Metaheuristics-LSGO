#include <iostream>

#include "Parser.h"

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Usage: ";
    std::cout << argv[0] << " input_CSV input_image [marked_image] output_image" << std::endl;
    return -1;
  }
  else if (argc == 4) {
    Parser parser(argv[1], argv[2], (char*)0);
    parser.run();
    parser.drawSelectedVesicles(argv[3]);
    return 0;
  }
  else {
    Parser parser(argv[1], argv[2], argv[3]);
    parser.run();
    parser.drawSelectedVesicles(argv[4]);
    return 0;
  }
}
