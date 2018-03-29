#ifndef IMAGEMGR
#define IMAGEMGR

#include <vector>

#include "Magick++.h"

#include "ITKConfig.h"
#include "Objective.h"
#include "Vesicle.h"
#include "VesicleInfo.h"

class ImageManager {
  public:
    // Class constructor
    ImageManager(char* _path);
    ~ImageManager();

    void addVesicle(const Vesicle& _vesicle);
    void drawSelectedVesicles (char* _path) const;
    void printSelectedVesicles(char* _path) const;
    long double objectiveFunction(int _x, int _y, int _radius, int _width, VesicleInfo* info=NULL);
    int getDimX() const {return dimX;}
    int getDimY() const {return dimY;}

  protected:
    // Attributes
    int dimX, dimY, dimImg;  // X and Y image dimensions in pixels and total
    char * path;

    std::vector<bool> selected_pixels;      // Image map of selected pixels
    ReaderType::Pointer reader;             // Image reader
    ImageType::ConstPointer image;
    std::vector<Vesicle> selected_vesicles; // Array of currently selected vesicles
};

#endif
