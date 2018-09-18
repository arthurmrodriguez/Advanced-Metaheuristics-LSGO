#include "ImagePrinter.h"

#include <iostream>

#include "DrawVesicle.h"
#include "Vesicle.h"
#include "ImageManager.h"

// Class constructor
ImagePrinter::ImagePrinter(char* _path) : path(_path) {
  reader = ReaderType::New();
  reader->SetFileName(_path);
  try {
    reader->Update();
  }
  catch (itk::ExceptionObject &err) {
    std::cerr << "Error: the input image '" << path << "' could not be read." << std::endl;
    std::cout << err << std::endl;
  }
  image = reader->GetOutput();
}

ImagePrinter::~ImagePrinter() {
}

// Adds the given vesicle to the selected vesicles structure
void ImagePrinter::addVesicle(const Vesicle& _vesicle) {
  selected_vesicles.push_back(_vesicle);
}

// Creates and draws onto the given path image the representation
// of the objective function value of every vesicle in the
// selected vesicles structure
/*
bool vesiclesCmpLesser (Vesicle& v1, Vesicle& v2) { return v1.getValue() < v2.getValue(); }
bool vesiclesCmpGreater(Vesicle& v1, Vesicle& v2) { return v1.getValue() > v2.getValue(); }
*/
void ImagePrinter::drawSelectedVesicles(char* _path, const Magick::Image* imgMarked) const {
  int len = selected_vesicles.size();

  if (len>0) {
    /* NOT USED
    double minValue = (std::min_element(selected_vesicles.begin(), selected_vesicles.end(), vesiclesCmpLesser ))->getValue();
    double maxValue = (std::max_element(selected_vesicles.begin(), selected_vesicles.end(), vesiclesCmpGreater))->getValue();
    */

    // Draw the selected vesicles in the output image
    Magick::Image outImage(path);

    int goodVesicles = 0;

    for (int i=0; i<len; i++) {
      const Vesicle& vesicle = selected_vesicles[i];

      if (!imgMarked || isInInterestArea(vesicle, *imgMarked)) {
        goodVesicles++;
        int radius=vesicle.getRadius();

        NeighborhoodIteratorType::RadiusType rt;
        rt.Fill(radius);
        NeighborhoodIteratorType* nit = new NeighborhoodIteratorType(rt, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());

        DrawVesicle draw(vesicle, nit, outImage);
        draw.draw();
        delete nit;
      }
    }
    outImage.write(_path);

    std::cout << "Detected vesicles: " << goodVesicles << " / " << len << std::endl;
  }
}

// Checks if a detected vesicle belongs to the interest area or not
bool ImagePrinter::isInInterestArea(const Vesicle& ves, const Magick::Image& imgMarked) const {
  int x = ves.getX();
  int y = ves.getY();

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


void ImagePrinter::bestFitVesicles(char* pathCSVOut) const {
  ImageManager imgMgr(path);

  std::vector<Vesicle> new_vesicles;
  std::vector<Vesicle>::const_iterator it;

  for (it = selected_vesicles.begin(); it != selected_vesicles.end(); it++) {
    Vesicle bestVesicle = *it;

    for (unsigned x = 0; x <= 2; x++) {
      for (unsigned y = 0; y <= 2; y++) {
        for (unsigned radius = 4; radius <= 6; radius++) {
          VesicleInfo info;
          double newScore = imgMgr.objectiveFunction(it->getX() + x, it->getY() + y, radius, it->getWidth(), &info);

          if (newScore < bestVesicle.getValue()) {
            Vesicle v(it->getX() + x, it->getY() + y, radius, it->getWidth(), newScore, it->isValid());
            bestVesicle = v;
          }
        }
      }
    }

    new_vesicles.push_back(bestVesicle);

    std::cout << "Old vesicle: (" << it->getX() << ", " << it->getY() << "), radius = " << it->getRadius() << " and score = " << it->getValue() << std::endl;
    std::cout << "New vesicle: (" << bestVesicle.getX() << ", " << bestVesicle.getY() << "), radius = " << bestVesicle.getRadius() << " and score = " << bestVesicle.getValue()  << std::endl << std::endl;
  }

  std::ofstream of(pathCSVOut);

  for (it = new_vesicles.begin(); it != new_vesicles.end(); it++)
    of << it->getX() << ", " << it->getY() << ", 0" << std::endl;

  of.close();

  return;
}
