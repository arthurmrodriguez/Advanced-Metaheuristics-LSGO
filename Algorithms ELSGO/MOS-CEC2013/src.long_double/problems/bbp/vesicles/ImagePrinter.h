#ifndef IMAGEPRINTER
#define IMAGEPRINTER

#include <vector>

#include "Magick++.h"

#include "ITKConfig.h"

class Vesicle;

class ImagePrinter {
  public:
    // Class constructor
    ImagePrinter(char* _path);
    ~ImagePrinter();

    void addVesicle(const Vesicle& _vesicle);
    void drawSelectedVesicles(char* _path, const Magick::Image* imgMarked) const;
    void bestFitVesicles(char* pathCSVOut) const;

  private:
    bool isInInterestArea(const Vesicle& ves, const Magick::Image& imgMarked) const;

    std::vector<Vesicle> selected_vesicles; // Array of currently selected vesicles
    ReaderType::Pointer reader;             // Image reader
    ImageType::ConstPointer image;
    char* path;
};

#endif
