#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <sys/stat.h>
#include <vector>

#include <boost/tokenizer.hpp>
#include <Magick++.h>

typedef struct {
  int x;
  int y;
  int inArea;
  std::string rowID;
  std::string line;
} InfoT;

std::vector<InfoT> parse_file(const char* pathCSV, std::string& header, bool isKnime);
void write_file(const char* pathCSV, std::string header, std::vector<InfoT> lines);
bool isInInterestArea(const InfoT& ves, const Magick::Image& imgMarked);

void print_usage(const char* prog) {
  std::cerr << "Usage: " << prog << " fin fout imgMarked [knime]" << std::endl;
}

int main(int argc, char** argv) {
  // Check for the correct number of arguments
  if (argc != 4 && argc != 5) {
    print_usage(argv[0]);
    exit(-1);
  }

  // Configure parser to read knime output files
  bool isKnime = false;

  if (argc == 5 && strcmp(argv[4], "knime") == 0)
    isKnime = true;

  struct stat stFileInfo;

  // Check for input file to exist
  if (stat(argv[1], &stFileInfo) != 0) {
    std::cerr << "Error: input file '" << argv[1] << "' does not exist." << std::endl;
    print_usage(argv[0]);
    exit(-1);
  }

  // Check for marked image to exist
  if (stat(argv[3], &stFileInfo) != 0) {
    std::cerr << "Error: marked image '" << argv[3] << "' does not exist." << std::endl;
    print_usage(argv[0]);
    exit(-1);
  }

  // Read input file
  std::string header;
  std::vector<InfoT> lines=parse_file(argv[1], header, isKnime);

  // Read marked image
  Magick::Image imgMarked;
  imgMarked.read(argv[3]);

  std::vector<InfoT> filteredLines;

  for (std::vector<InfoT>::const_iterator it=lines.begin(); it!=lines.end(); it++) {
    InfoT newLine = *it;

    if (isInInterestArea(*it, imgMarked))
      newLine.inArea = 1;
    else
      newLine.inArea = 0;

    filteredLines.push_back(newLine);
  }

  // Write output file
  write_file(argv[2], header+", \"inArea\"", filteredLines);

  return 0;
}


std::vector<InfoT> parse_file(const char* pathCSV, std::string& header, bool isKnime) {
  std::ifstream csv(pathCSV);
  std::vector<InfoT> lines;

  bool firstLine = true;

  if (csv.is_open()) {
    std::string line;

    while (getline(csv, line)) {
      // Ignore headers
      if (firstLine) {
        header = line;
        firstLine = false;
        continue;
      }

      InfoT info;

      boost::tokenizer<> tok(line);
      boost::tokenizer<>::iterator beg=tok.begin();

      // Knime specific: Ignore first column (Row ID)
      if (isKnime) {
        std::istringstream issID(*beg);
        issID >> info.rowID;
        beg++;
      }

      std::istringstream issX(*beg);
      issX >> info.x;

      beg++;

      std::istringstream issY(*beg);
      issY >> info.y;

      info.line = line;

      lines.push_back(info);
    }

    csv.close();
  }

  return lines;
}


// Checks if a detected vesicle belongs to the interest area or not
bool isInInterestArea(const InfoT& ves, const Magick::Image& imgMarked) {
  int x = ves.x;
  int y = ves.y;

  int maxX = imgMarked.baseColumns();

  // This is equivalent to RGB(200, 50, 200, 0) in 8 bits
  Magick::Color border (51400, 12850, 51400, 0);
  Magick::Color curPxColor;

  int foundBorders = 0;

  while (x <= maxX) {
    curPxColor = imgMarked.pixelColor(x, y);
    int borderPixels = 0;

    if (curPxColor == border) {
      foundBorders++;
      while (curPxColor == border && x <= maxX) {
        x++;
        borderPixels++;
        curPxColor = imgMarked.pixelColor(x, y);
      }
    }
    x++;
  }

  if ((foundBorders % 2) == 0)
    return false;
  else
    return true;
}


void write_file(const char* pathCSV, std::string header, std::vector<InfoT> lines) {
  std::ofstream of(pathCSV);

  if (of.good()) {
    of << "\"row ID\", \"inArea\"" << std::endl;

    for (std::vector<InfoT>::const_iterator it=lines.begin(); it!=lines.end(); it++)
      of << it->rowID << ", " << it->inArea << std::endl;
  }

  of.close();
}
