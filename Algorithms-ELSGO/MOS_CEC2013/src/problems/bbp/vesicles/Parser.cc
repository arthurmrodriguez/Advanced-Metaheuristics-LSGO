#include "Parser.h"

#include "ImageManager.h"
#include "ImagePrinter.h"
#include "Vesicle.h"

#include <fstream>

// Class constructor
Parser::Parser(char* _pathCSV, char* _pathImg, char* _pathImgMarked) : pathCSV(_pathCSV), imgMarked(NULL) {
  img    = new ImagePrinter(_pathImg);
  imgMgr = new ImageManager(_pathImg);

  std::ifstream marked (_pathImgMarked);
  if (!marked.fail()) {
    imgMarked = new Magick::Image();
    imgMarked->read(_pathImgMarked);
  }
}

Parser::~Parser() {
  if (img)
    delete img;
  if (imgMgr)
    delete imgMgr;
  if (imgMarked)
    delete imgMarked;
}

void Parser::run() const {
  std::ifstream csv(pathCSV);

  if (csv.is_open()) {
    int x, y, radius, width, intAux, len, valid;
    double value;
    std::string line, strAux;
    Vesicle* vesicle;

    while (csv.good()) {
      getline (csv, line);
      len = line.length()-1;

      // X
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issX(strAux);
      issX >> x;
      line = line.substr(intAux+1, len);

      // Y
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issY(strAux);
      issY >> y;
      line = line.substr(intAux+1, len);
      len = line.length()-1;

      // Radius
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issR(strAux);
      issR >> radius;
      line = line.substr(intAux+1, len);
      len = line.length()-1;

      // Width
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issW(strAux);
      issW >> width;
      line = line.substr(intAux+1, len);
      len = line.length()-1;

      // Fitness Value
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issV(strAux);
      issV >> value;
      line = line.substr(intAux+1, len);
      len = line.length()-1;

      // Is this vesicle valid?
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issValid(strAux);
      issValid >> valid;

      vesicle = new Vesicle(x, y, radius, width, value, valid);
      img->addVesicle(*vesicle);
      delete vesicle;
    }
    csv.close();
  }
}


void Parser::parseValData() const {
  std::ifstream csv(pathCSV);

  if (csv.is_open()) {
    int x, y, intAux, len;
    std::string line, strAux;
    Vesicle* vesicle;

    while (csv.good()) {
      getline (csv, line);
      len = line.length()-1;

      // X
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issX(strAux);
      issX >> x;
      line = line.substr(intAux+1, len);

      // Y
      intAux = line.find_first_of(",");
      strAux = line.substr (0, intAux);
      std::istringstream issY(strAux);
      issY >> y;
      line = line.substr(intAux+1, len);
      len = line.length()-1;

      double score = imgMgr->objectiveFunction(x, y, 4, 2, (VesicleInfo*)0);

      vesicle = new Vesicle(x, y, 4, 2, score, 0);
      img->addVesicle(*vesicle);
      delete vesicle;
    }
    csv.close();
  }
}


void Parser::drawSelectedVesicles(char* _path) const {
  img->drawSelectedVesicles(_path, imgMarked);
}


void Parser::bestFitVesicles(char* _path) const {
  img->bestFitVesicles(_path);
}
