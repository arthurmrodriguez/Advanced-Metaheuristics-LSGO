#include "ImageManager.h"

#include "DrawVesicle.h"

// Class constructor
ImageManager::ImageManager(char* _path) : path(_path) {
  reader = ReaderType::New();
  reader->SetFileName(path);

  try {
    reader->Update();
  }
  catch (itk::ExceptionObject &err) {
    std::cerr << "Error: the input image '" << path << "' could not be read." << std::endl;
    std::cerr << err << std::endl;
  }

  image = reader->GetOutput();

  dimX = image->GetLargestPossibleRegion().GetSize()[0];
  dimY = image->GetLargestPossibleRegion().GetSize()[1];

  dimImg = dimX*dimY;
  selected_pixels.resize(dimImg, false);
}

ImageManager::~ImageManager() {
}

// Adds the given vesicle to the selected vesicles structure
void ImageManager::addVesicle(const Vesicle& _vesicle) {
  selected_vesicles.push_back(_vesicle);

  int _x      = _vesicle.getX();
  int _y      = _vesicle.getY();
  int _radius = _vesicle.getRadius();
  int _width  = _vesicle.getWidth();

  NeighborhoodIteratorType::RadiusType rt;
  rt.Fill(_radius);
  NeighborhoodIteratorType* nit = new NeighborhoodIteratorType(rt, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());

  Objective obj(_x, _y, _radius, _width, nit, selected_pixels, dimX, dimY);
  obj.markVesiclesMap();
  delete nit;
}

// Creates and draws onto the given path image the representation
// of the objective function value of every vesicle in the
// selected vesicles structure
/*
bool vesiclesCmpLesser (Vesicle& v1, Vesicle& v2) { return v1.getValue() < v2.getValue(); }
bool vesiclesCmpGreater(Vesicle& v1, Vesicle& v2) { return v1.getValue() > v2.getValue(); }
*/
void ImageManager::drawSelectedVesicles(char* _path) const {
  int len = selected_vesicles.size();

  if (len>0) {
    /* NOT USED
    double minValue = (std::min_element(selected_vesicles.begin(), selected_vesicles.end(), vesiclesCmpLesser ))->getValue();
    double maxValue = (std::max_element(selected_vesicles.begin(), selected_vesicles.end(), vesiclesCmpGreater))->getValue();
    */

    // Draw the selected vesicles in the output image
    Magick::Image outImage(path);

    for (int i=0; i<len; i++) {
      const Vesicle& vesicle = selected_vesicles[i];
      int radius=vesicle.getRadius();

      NeighborhoodIteratorType::RadiusType rt;
      rt.Fill(radius);
      NeighborhoodIteratorType* nit = new NeighborhoodIteratorType(rt, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());

      DrawVesicle draw(vesicle, nit, outImage);
      draw.draw();
      delete nit;
    }

    // BEGIN: debug trace, output marked pixels in yellow
    int dimY = image->GetLargestPossibleRegion().GetSize()[1];

    for (unsigned i = 0; i < selected_pixels.size(); i++) {
      if (selected_pixels[i]) {
        int absX = i / dimY;
        int absY = i % dimY;
        outImage.pixelColor(absX, absY, Magick::Color(MaxRGB, MaxRGB, 0, 0));
      }
    }
    // END: debug trace, output marked pixels in yellow

    outImage.write(_path);
  }
}

void ImageManager::printSelectedVesicles(char* _path) const {
  std::ofstream fout;
  fout.open (_path);

  for (int i=0; i<selected_vesicles.size(); i++) {
    const Vesicle& vesicle = selected_vesicles[i];
    fout << vesicle.getX     () << "," << vesicle.getY    () << ","
         << vesicle.getRadius() << "," << vesicle.getWidth() << ","
         << vesicle.getValue () << "," << std::endl;
  }

  fout.close();
}

// Calculates the objective function value for the specific given parameters
double ImageManager::objectiveFunction(int _x, int _y, int _radius, int _width, VesicleInfo* info) {
  NeighborhoodIteratorType::RadiusType rt;
  rt.Fill(_radius);
  NeighborhoodIteratorType* nit = new NeighborhoodIteratorType(rt, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());

  Objective obj(_x, _y, _radius, _width, nit, selected_pixels, dimX, dimY);
  double res = obj.objectiveFunction(info);
  delete nit;

  return res;
}
