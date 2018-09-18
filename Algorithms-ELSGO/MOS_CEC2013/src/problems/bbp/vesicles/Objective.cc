#include "Objective.h"

#include "Circle.h"

Objective::Objective (int _x, int _y, int _radius, int _width, NeighborhoodIteratorType* _nit, std::vector<bool> &_occupation, int _dimX, int _dimY)
  : x(_x), y(_y), radius(_radius), width(_width), nit(_nit), occupation(_occupation), tamX(_dimX), tamY(_dimY) {
  edgePx             = 0;
  centerPx           = 0;
  edgeIntensity      = 0;
  centerIntensity    = 0;
  edgeAvgIntensity   = 0;
  centerAvgIntensity = 0;
  edgeDeviation      = 0;
  centerDeviation    = 0;
  kEdge              = 1;
  kCenter            = 0;
}

Objective::~Objective() {
}

void Objective::handlePixel (int pixelPosition, bool isEdge) {
  int px;
  int absX = nit->GetIndex(pixelPosition)[0];
  int absY = nit->GetIndex(pixelPosition)[1];

  if ((absX>0) && (absX<tamX) && (absY>0) && (absY<tamY)) {
    int pos = absX*tamY+absY;
    if (isEdge) {
      if (occupation[pos])
        px = Objective::white;
      else
        px = nit->GetPixel(pixelPosition);
      edgePx++;
      edgeIntensity += px;
    }
    else {
      if (occupation[pos])
        px = black;
      else
        px = nit->GetPixel(pixelPosition);
      centerPx++;
      centerIntensity += px;
    }
  }
}

void Objective::handlePixelDeviation (int pixelPosition, bool isEdge) {
  int px;
  int absX = nit->GetIndex(pixelPosition)[0];
  int absY = nit->GetIndex(pixelPosition)[1];

  if ((absX > 0) && (absX<tamX) && (absY>0) && (absY<tamY)) {
    int pos = absX*tamY+absY;
    if (isEdge) {
      if (occupation[pos])
        px = white;
      else
        px = nit->GetPixel(pixelPosition);
        edgeDeviation += (px-edgeAvgIntensity)*(px-edgeAvgIntensity);
        //std::cout << "diff: " << px - edgeAvgIntensity << std::endl;
    }
    else {
      if (occupation[pos])
        px = black;
      else
        px = nit->GetPixel(pixelPosition);
        centerDeviation += (px-centerAvgIntensity)*(px-centerAvgIntensity);
    }
  }
}

void Objective::iterateCircle() {
  Circle circle (radius, width);
  int cStat = 0, dim = 2*radius, diameter = dim+1, in = dim-width, aux;

  for (int h=0; h<diameter; h++)
    if (circle.isMarked(0, h))
      handlePixel(h, true);

  for (int i=1; i<dim; i++) {
    for (int j=0; j<diameter; j++) {
      aux = circle.isMarked(i, j);
      switch (cStat) {
        case 0:
               if (aux==1)              { handlePixel(i*diameter+j, true); cStat=1; }
          else if (aux==2)              { handlePixel(i*diameter+j, true); cStat=3; }
          break;
        case 1:
               if (aux==0)              { handlePixel(i*diameter+j, true); cStat=2; }
          else if (aux==1)              { handlePixel(i*diameter+j, true);          }
          else                          { handlePixel(i*diameter+j, true); cStat=3; }
          break;
        case 2:
               if (aux==0)              { handlePixel(i*diameter+j, true);          }
          else if ((i<width)||(i>in))   { handlePixel(i*diameter+j, true); cStat=7; }
          else                          { handlePixel(i*diameter+j, true); cStat=3; }
          break;
        case 3:
               if ((i==width)||(i==in)) { handlePixel(i*diameter+j, true ); cStat=5; }
          else if (aux==0)              { handlePixel(i*diameter+j, false); cStat=4; }
          else                          { handlePixel(i*diameter+j, true );          }
          break;
        case 4:
               if (aux==0)              { handlePixel(i*diameter+j, false);          }
          else                          { handlePixel(i*diameter+j, true ); cStat=5; }
          break;
        case 5:
               if (aux==0)              { handlePixel(i*diameter+j, true ); cStat=6; }
          else if (aux==1)              { handlePixel(i*diameter+j, true );          }
          else                          { handlePixel(i*diameter+j, true ); cStat=7; }
          break;
        case 6:
               if (aux==0)              { handlePixel(i*diameter+j, true );          }
          else                          { handlePixel(i*diameter+j, true ); cStat=7; }
          break;
        case 7:
               if (aux==1)              { handlePixel(i*diameter+j, true );          }
          else                          {                                   cStat=-1;}
          break;
        default:
          break;
      } // switch
    } // for j
    cStat=0;
  } // for i

  for (int k=0; k<diameter; k++)
      if (circle.isMarked(dim, k))
        handlePixel(dim*diameter+k, true);
}

void Objective::iterateCircleDeviation() {
  Circle circle (radius, width);
  int cStat = 0, dim = 2*radius, diameter = dim+1, in = dim-width, aux;

  for (int h=0; h<diameter; h++)
      if (circle.isMarked(0, h)) {
        //std::cout << "(" << 0 <<", " << h << ")" << std::endl;
        handlePixelDeviation(h, true);
      }

  for (int i=1; i<dim; i++) {
    for (int j=0; j<diameter; j++) {
      bool sel = false;
      aux = circle.isMarked(i, j);
      switch(cStat) {
        case 0:
               if (aux==1)              { handlePixelDeviation(i*diameter+j, sel=true ); cStat=1; }
          else if (aux==2)              { handlePixelDeviation(i*diameter+j, sel=true ); cStat=3; }
          break;
        case 1:
               if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=true ); cStat=2; }
          else if (aux==1)              { handlePixelDeviation(i*diameter+j, sel=true );          }
          else                          { handlePixelDeviation(i*diameter+j, sel=true ); cStat=3; }
          break;
        case 2:
               if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=true );          }
          else if ((i<width)||(i>in))   { handlePixelDeviation(i*diameter+j, sel=true ); cStat=7; }
          else                          { handlePixelDeviation(i*diameter+j, sel=true ); cStat=3; }
          break;
        case 3:
               if ((i==width)||(i==in)) { handlePixelDeviation(i*diameter+j, sel=true ); cStat=5; }
          else if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=false); cStat=4; }
          else                          { handlePixelDeviation(i*diameter+j, sel=true );          }
          break;
        case 4:
               if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=false);          }
          else                          { handlePixelDeviation(i*diameter+j, sel=true ); cStat=5; }
          break;
        case 5:
               if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=true ); cStat=6; }
          else if (aux==1)              { handlePixelDeviation(i*diameter+j, sel=true );          }
          else                          { handlePixelDeviation(i*diameter+j, sel=true ); cStat=7; }
          break;
        case 6:
               if (aux==0)              { handlePixelDeviation(i*diameter+j, sel=true );          }
          else                          { handlePixelDeviation(i*diameter+j, sel=true ); cStat=7; }
          break;
        case 7:
               if (aux==1)              { handlePixelDeviation(i*diameter+j, sel=true );          }
          else                          {                                            cStat=-1;}
          break;
        default:
          break;
      } // switch
      //if (sel)
        //std::cout << "(" << i <<", " << j << ")" << std::endl;
    } // for j
    cStat=0;
  } // for i

  for (int k=0; k<diameter; k++)
      if (circle.isMarked(dim, k)) {
        //std::cout << "(" << dim <<", " << k << ")" << std::endl;
        handlePixelDeviation(dim*diameter+k, true);
      }
}

void Objective::iterateCircleMark() {
  Circle circle (radius, width);
  int cStat = 0, dim = 2*radius, diameter = dim+1, in = dim-width, aux;

  for (int h=0; h<diameter; h++)
      if (circle.isMarked(0, h))
        markPixel(h, true);

  for (int i=1; i<dim; i++) {
    for (int j=0; j<diameter; j++) {
      aux = circle.isMarked(i, j);
      switch(cStat) {
        case 0:
               if (aux==1)              { markPixel(i*diameter+j, true ); cStat=1; }
          else if (aux==2)              { markPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 1:
               if (aux==0)              { markPixel(i*diameter+j, true ); cStat=2; }
          else if (aux==1)              { markPixel(i*diameter+j, true );          }
          else                          { markPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 2:
               if (aux==0)              { markPixel(i*diameter+j, true );          }
          else if ((i<width)||(i>in))   { markPixel(i*diameter+j, true ); cStat=7; }
          else                          { markPixel(i*diameter+j, true ); cStat=3; }
          break;
        case 3:
               if ((i==width)||(i==in)) { markPixel(i*diameter+j, true ); cStat=5; }
          else if (aux==0)              { markPixel(i*diameter+j, false); cStat=4; }
          else                          { markPixel(i*diameter+j, true );          }
          break;
        case 4:
               if (aux==0)              { markPixel(i*diameter+j, false);          }
          else                          { markPixel(i*diameter+j, true ); cStat=5; }
          break;
        case 5:
               if (aux==0)              { markPixel(i*diameter+j, true ); cStat=6; }
          else if (aux==1)              { markPixel(i*diameter+j, true );          }
          else                          { markPixel(i*diameter+j, true ); cStat=7; }
          break;
        case 6:
               if (aux==0)              { markPixel(i*diameter+j, true );          }
          else                          { markPixel(i*diameter+j, true ); cStat=7; }
          break;
        case 7:
               if (aux==1)              { markPixel(i*diameter+j, true);           }
          else                          {                                 cStat=-1;}
          break;
        default:
          break;
      } // switch
    } // for j
    cStat=0;
  } // for i

  for (int k=0; k<diameter; k++)
    if (circle.isMarked(dim, k))
      markPixel(dim*diameter+k, true);
}

// Calculates the objective function value for the given pixel
double Objective::objectiveFunction(VesicleInfo* info) {
  // Initialize the image region handler
  ImageType::IndexType index;
  index[0] = x;
  index[1] = y;
  nit->SetLocation(index);

  iterateCircle();

  // Compute the intensity averages
    edgeAvgIntensity = (  edgePx > 0) ? (double)  edgeIntensity /   edgePx : 0;
  centerAvgIntensity = (centerPx > 0) ? (double)centerIntensity / centerPx : 0;

  iterateCircleDeviation();

  // Compute the deviations
    edgeDeviation = (  edgePx > 0) ? sqrt(  edgeDeviation /   edgePx) : 0;
  centerDeviation = (centerPx > 0) ? sqrt(centerDeviation / centerPx) : 0;

  if (info) {
    info->edgePx             = edgePx;
    info->centerPx           = centerPx;
    info->edgeAvgIntensity   = edgeAvgIntensity;
    info->edgeDeviation      = edgeDeviation;
    info->centerAvgIntensity = centerAvgIntensity;
    info->centerDeviation    = centerDeviation;
  }

  // Compute the objective function value
  return (edgeAvgIntensity + kEdge*edgeDeviation - centerAvgIntensity - kCenter*centerDeviation);
}

void Objective::markPixel (int pixelPosition, bool isEdge) {
  int absX = nit->GetIndex(pixelPosition)[0];
  int absY = nit->GetIndex(pixelPosition)[1];

  if ((absX>0) && (absX<tamX) && (absY>0) && (absY<tamY) && !isEdge)
    occupation[absX*tamY+absY] = true;
}

// Marks the pixels in the occupation map for the given vesicle
void Objective::markVesiclesMap() {
  // Initialize the image region handler
  ImageType::IndexType index;
  index[0] = x;
  index[1] = y;
  nit->SetLocation(index);

  iterateCircleMark();
}
