#include "DrawVesicle.h"

#include "Circle.h"
#include "Vesicle.h"

Magick::Color DrawVesicle::currentColor (188, 234, 23, 0);

DrawVesicle::DrawVesicle (int _x, int _y, int _radius, int _width, NeighborhoodIteratorType* _nit, Magick::Image& _outImage, int _valid)
  :  x(_x), y(_y), radius(_radius), width(_width), valid(_valid), nit(_nit), outImage(_outImage) {
}

DrawVesicle::DrawVesicle (Vesicle _vesicle, NeighborhoodIteratorType* _nit, Magick::Image &_outImage)
  : x(_vesicle.getX()), y(_vesicle.getY()), radius(_vesicle.getRadius()), width(_vesicle.getWidth()), value(_vesicle.getValue()), valid(_vesicle.isValid()), nit(_nit), outImage(_outImage) {
}

DrawVesicle::~DrawVesicle() {
}

void DrawVesicle::drawPixel (int pixelPosition, bool isEdge) const {
  int absX = nit->GetIndex(pixelPosition)[0];
  int absY = nit->GetIndex(pixelPosition)[1];

  if (isEdge)
    outImage.pixelColor(absX, absY, DrawVesicle::currentColor);
}

// Calculates the objective function value for the given pixel
void DrawVesicle::draw() const {
  // Initialize the image region handler
  ImageType::IndexType index;
  index[0] = x;
  index[1] = y;
  nit->SetLocation(index);

  Circle* circle = new Circle(radius, width);
  int cStat = 0, dim = 2*radius, diameter = dim+1, in = dim-width, aux;

  for (int h=0; h<diameter; h++)
    if (circle->isMarked(0, h))
      drawPixel(h, true);

  for (int i=1; i<dim; i++) {
    for (int j=0; j<diameter; j++) {
      aux = circle->isMarked(i, j);

      switch(cStat) {
        case 0:
               if (aux==1)              { drawPixel(i*diameter+j, true ); cStat=1; }
          else if (aux==2)              { drawPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 1:
               if (aux==0)              { drawPixel(i*diameter+j, true ); cStat=2; }
          else if (aux==1)              { drawPixel(i*diameter+j, true );          }
          else                          { drawPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 2:
               if (aux==0)              { drawPixel(i*diameter+j, true );          }
          else if ((i<width)||(i>in))   { drawPixel(i*diameter+j, true ); cStat=7; }
          else                          { drawPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 3:
               if ((i==width)||(i==in)) { drawPixel(i*diameter+j, true ); cStat=5; }
          else if (aux==0)              { drawPixel(i*diameter+j, false); cStat=4; }
          else                          { drawPixel(i*diameter+j, true );          }
          break;
        case 4:
               if (aux==0)              { drawPixel(i*diameter+j, false);          }
          else                          { drawPixel(i*diameter+j, true ); cStat=5; }
          break;
        case 5:
               if (aux==0)              { drawPixel(i*diameter+j, true ); cStat=6; }
          else if (aux==1)              { drawPixel(i*diameter+j, true );          }
          else                          { drawPixel(i*diameter+j, true ); cStat=7; }
          break;
        case 6:
               if (aux==0)              { drawPixel(i*diameter+j, true );          }
          else                          { drawPixel(i*diameter+j, true ); cStat=7; }
          break;
        case 7:
              if (aux==1)               { drawPixel(i*diameter+j, true);           }
          else cStat=-1;
          break;
        default:
          break;
      } // switch
    } // for j
    cStat=0;
  } // for i

  for (int k=0; k<diameter; k++)
    if (circle->isMarked(dim, k))
      drawPixel(dim*diameter+k, true);

  delete(circle);
}


void DrawVesicle::draw2() const {
  int dim = 2*radius, diameter = dim+1;

  // Initialize the image region handler
  ImageType::IndexType index;
  index[0] = x;
  index[1] = y;
  nit->SetLocation(index);

  Magick::Color oldColor = DrawVesicle::currentColor;

  if (valid == 0)
    DrawVesicle::currentColor = Magick::Color(0, MaxRGB, 0, 0);
  else if (valid == 1)
    DrawVesicle::currentColor = Magick::Color(MaxRGB, 0, 0, 0);
  else if (valid == 2)
    DrawVesicle::currentColor = Magick::Color(0, 0, MaxRGB, 0);
  else if (valid == 5)
    DrawVesicle::currentColor = Magick::Color(MaxRGB, MaxRGB, 0, 0);
  else
    DrawVesicle::currentColor = Magick::Color(0, MaxRGB, MaxRGB, 0);

  drawPixel((diameter+1)*radius, true);
  drawPixel(((diameter+1)*radius)-1, true);
  drawPixel(((diameter+1)*radius)+1, true);
  drawPixel(((diameter+1)*radius)-diameter, true);
  drawPixel(((diameter+1)*radius)+diameter, true);

  DrawVesicle::currentColor = oldColor;
}
