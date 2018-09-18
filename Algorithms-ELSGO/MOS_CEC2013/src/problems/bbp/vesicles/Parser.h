#ifndef PARSER
#define PARSER

#include <fstream>

#include "Magick++.h"

class ImageManager;
class ImagePrinter;

class Parser {

  public:
    // Class constructor
    Parser(char* _pathCSV, char* _pathImg, char* _pathImgMarked);
    ~Parser();

    void run()          const;
    void parseValData() const;

    void drawSelectedVesicles(char* _path) const;
    void bestFitVesicles     (char* _path) const;

  protected:
    ImagePrinter* img;  // Data structure for image access
    ImageManager* imgMgr;
    char* pathCSV;
    Magick::Image* imgMarked;
};

#endif
