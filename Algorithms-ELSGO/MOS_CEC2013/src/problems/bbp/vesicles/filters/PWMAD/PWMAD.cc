#include <algorithm>
#include <iostream>
#include <sstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <Magick++.h>

using namespace boost::numeric::ublas;

int N = 3;
int K = 1;
const double INVALID = -9999;

// Aux functions
double median(const matrix<double>&);
void createWindowMatrix(const matrix<double>& origMatrix, matrix< matrix<double> >& windowMatrix, int Wsize);

// Main program
int main(int argc, char** argv) {
  // Check for the correct number of parameters
  if (argc < 6) {
    std::cout << "Usage: " << argv[0] << " input_image window_size iterations output_image noise_prefix" << std::endl;
    exit(-1);
  }

  K = atoi(argv[2]);
  N = atoi(argv[3]);

  // Open input image for reading
  Magick::Image inputImg(argv[1]);

  int rows = inputImg.rows();
  int cols = inputImg.columns();

  matrix<double> noisyImage (cols, rows);
  matrix<double> medianImage(cols, rows);
  matrix<double> devImage   (cols, rows);

  matrix< matrix<double> > noiseWindows (cols, rows);
  matrix< matrix<double> > medianWindows(cols, rows);
  matrix< matrix<double> > devWindows   (cols, rows);

  // Initialize noisyImage with the original image pixels
  for (unsigned i = 0; i < cols; i++) {
    for (unsigned j = 0; j < rows; j++) {
      Magick::ColorGray color = inputImg.pixelColor(i, j);
      noisyImage(i, j) = color.shade();
    }
  }

  // Compute noiseWindows
  createWindowMatrix(noisyImage, noiseWindows, K);

  // Compute medianImage
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      medianImage(i, j) = median(noiseWindows(i, j));
    }
  }

  std::vector< matrix<double>* > dij(N);
  std::vector< matrix< matrix<double> >* > dijW(N);

  // Compute dij0
  Magick::Image noise0(Magick::Geometry(cols, rows), Magick::ColorGray(1));

  dij[0] = new matrix<double>(cols, rows);
  matrix<double>& dij0 = *(dij[0]);

  for (int i = 0; i < cols; i++)
    for (int j = 0; j < rows; j++) {
      dij0(i, j) = fabs(noisyImage(i, j) - medianImage(i, j));
      noise0.pixelColor(i, j, Magick::ColorGray(dij0(i, j)));
    }

  noise0.write((std::string(argv[5]) + "_0.tiff").c_str());

  dijW[0] = new matrix< matrix<double> >(cols, rows);
  matrix< matrix<double> >& dijW0 = *(dijW[0]);

  createWindowMatrix(dij0, dijW0, K);

  // Conduct iterative process
  for (unsigned n = 1; n < N; n++) {
    dij[n]  = new matrix<double>          (cols, rows);
    dijW[n] = new matrix< matrix<double> >(cols, rows);

    matrix<double>& dijn_old = *(dij[n-1]);
    matrix<double>& dijn     = *(dij[n  ]);

    matrix< matrix<double> >& dijWn_old = *(dijW[n-1]);
    matrix< matrix<double> >& dijWn     = *(dijW[n  ]);

    Magick::Image noise(Magick::Geometry(cols, rows), Magick::ColorGray(1));

    for (int i = 0; i < cols; i++)
      for (int j = 0; j < rows; j++) {
        dijn(i, j) = fabs(dijn_old(i, j) - median(dijWn_old(i, j)));
        noise.pixelColor(i, j, Magick::ColorGray(dijn(i, j)));
      }

    std::stringstream fname;
    fname << argv[5] << "_" << n << ".tiff";

    noise.write(fname.str());

    createWindowMatrix(dijn, dijWn, K);
  }


  Magick::Image filtImg(Magick::Geometry(cols, rows), Magick::ColorGray(1));
  matrix<double>& dijFinal = *(dij[N-1]);

  filtImg.depth(8);

  for (int i = 0; i < cols; i++)
    for (int j = 0; j < rows; j++) {
      if (dijFinal(i, j) > 0)
        filtImg.pixelColor(i, j, Magick::ColorGray(medianImage(i, j)));
      else
        filtImg.pixelColor(i, j, Magick::ColorGray(noisyImage (i, j)));
    }

  filtImg.write(argv[4]);

  for (int i = 0; i < dij.size(); i++) {
    delete dij[i];
    delete dijW[i];
  }

  return 0;
}


double median(const matrix<double>& m) {
  std::vector<double> v;

  for (int i = 0; i < m.size1(); i++)
    for (int j = 0; j < m.size1(); j++)
      if (m(i, j) != INVALID)
        v.push_back(m(i, j));

  int n = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n, v.end());

  return v[n];
}

void createWindowMatrix(const matrix<double>& origMatrix, matrix< matrix<double> >& windowMatrix, int Wsize) {
  int cols = origMatrix.size1();
  int rows = origMatrix.size2();

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      matrix<double> xij(2*Wsize+1, 2*Wsize+1);

      for (int m = i-Wsize, s=0; m <= i+Wsize; m++, s++) {
        for (int n = j-Wsize, t=0; n <= j+Wsize; n++, t++) {
          xij(s, t) = (m < 0 || m >= cols || n < 0 || n >= rows) ? INVALID : origMatrix(m, n);
        }
      }

      windowMatrix(i, j) = xij;
    }
  }
}
