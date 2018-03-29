#include <iomanip>
#include <vector>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <genomes/MOSGenome.h>
#include <GAEDAConfig.h>

#include "Circle.h"
#include "ImageManager.h"
#include "ITKConfig.h"
#include "VesicleGenome.h"

const int minWidth = 2, maxWidth = 2;   // Vesicle width range of iteration
const int minRadius = 4, maxRadius = 6; // Vesicle radius range of iteration

const double alpha = 0.9;

ImageManager* img = NULL;
std::vector<VesicleGenome<double>*> found_vesicles;

bool getEdgeInfo(const Vesicle& ves, double& distAllOUT, double& distAllIN, double& distNOUT, double& distNIN);

double fitness (VesicleGenome<double>& genome) {
  double res = 9999, aux = 9999;
  VesicleInfo info;

  genome.resetFits();
  genome.resetInfos();

  for (int radius=minRadius; radius<=maxRadius; radius++) {
    for (int width=minWidth; width<=maxWidth; width++) {
      aux = img->objectiveFunction((int)genome.gene(0), (int)genome.gene(1), radius, width, &info);

      // Fixed penalty. Should be done adaptive...
      //aux = (info.edgeAvgIntensity < 80 || info.edgeAvgIntensity > 170) || (info.centerAvgIntensity < 80 || info.centerAvgIntensity > 200) ? 1000 : aux;
      genome.addFit (aux );
      genome.addInfo(info);

      if ((radius==minRadius && width==minWidth) || (aux<res)) {
        res = aux;
        genome.radius  (radius);
        genome.width   (width );
      }
    }
  }

  return res;
}


extern "C" double objective (GAGenome& g) {
  VesicleGenome<double>& genome = DYN_CAST (VesicleGenome<double>&, g);

  double fit = fitness(genome);

  return fit;
}


extern "C" void individualInit (GAGenome& g) {
  return RealUniformInitializer (g);
}


extern "C" GAGenome* defineProblem () {

  GAEDAConfig* cfg = GAEDAConfig::handle();

  if (img == NULL)
    img = new ImageManager(cfg->getProblemData());

  // Definition of the genome and the allele set
  GAAlleleSetArray<double> alleles;
  alleles.add(maxRadius + maxWidth, img->getDimX() - maxWidth - maxRadius);
  alleles.add(maxRadius + maxWidth, img->getDimY() - maxWidth - maxRadius);
  alleles.add(0.0, 0.0           );

  VesicleGenome<double>* genome = new VesicleGenome<double> (alleles, objective);

  // Common operators
  genome->initializer (RealUniformInitializer);
  genome->comparator  (RealEuclideanComparator);

  // Specific stuff for GAs
  genome->crossover   (RealBlendCrossover);
  genome->mutator     (RealGaussianMutator);

  // Specific stuff for DE
  genome->crossover   (RealExponentialCrossover);

  // Specific stuff for MOS
  MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);

  return genome;

}


extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  MOSGenome& g = dynamic_cast<MOSGenome&>(pop->best());
  VesicleGenome<double>& genome = dynamic_cast<VesicleGenome<double>&>(*g.getGenome(GAID::RealEncoding));

  genome.valid(0);

  Vesicle ves ((int)genome.gene(0), (int)genome.gene(1), genome.radius(), genome.width(), genome.score(), genome.valid());
  img->addVesicle(ves);

  VesicleGenome<double>* new_vesicle = dynamic_cast<VesicleGenome<double>*>(genome.clone());
  found_vesicles.push_back(new_vesicle);

  // If we are in the last repetition, we can tidy up things...
  if (repetition == GAEDAConfig::handle()->getRepetitions() - 1) {
    // Print found vesicles
    for (unsigned i = 0; i < found_vesicles.size(); i++) {
      VesicleGenome<double>& vgenome = *found_vesicles[i];

      // Print row header
      if (i == 0) {
        std::cerr << "x, y, radius, width, score, fails, order, distAllOUT, distAllIN, distNOUT, distNIN";

        for (unsigned j = 0; j < vgenome.infos().size(); j++)
          std::cerr << ", edgeAvg" << j << ", edgeDev" << j << ", centerAvg" << j << ", centerDev" << j << ", edgePx" << j << ", centerPx" << j << ", score" << j;
        std::cerr << std::endl;
      }

      int nFails = 0;
      double avgFits = 0.0;
      for (unsigned cnt = 0; cnt <  vgenome.fits().size(); cnt++) {
        avgFits += vgenome.fits()[cnt];
        if (vgenome.fits()[cnt] == 1000)
          nFails++;
      }
      avgFits /= vgenome.fits().size();

      if (nFails == 0 && avgFits > 20)
        nFails = 5;

      std::cerr << setw(3) << right << (int) vgenome.gene (0) << ", "
                << setw(3) << right << (int) vgenome.gene (1) << ", "
                << setw(1) << right <<       vgenome.radius() << ", "
                << setw(1) << right <<       vgenome.width () << ", "
                << setw(9) << fixed << setprecision (4) << right << vgenome.score () << ", "
                << nFails  << ", "
                << setw(3) << right << (int) i                << ", ";

      Vesicle vesAux ((int)vgenome.gene(0), (int)vgenome.gene(1), vgenome.radius(), vgenome.width(), vgenome.score(), vgenome.valid());
      double distAllOUT, distAllIN, distNOUT, distNIN;
      getEdgeInfo(vesAux, distAllOUT, distAllIN, distNOUT, distNIN);

      std::cerr << setw(9) << fixed << setprecision (4) << right << distAllOUT << ", "
                << setw(9) << fixed << setprecision (4) << right << distAllIN  << ", "
                << setw(9) << fixed << setprecision (4) << right << distNOUT   << ", "
                << setw(9) << fixed << setprecision (4) << right << distNIN    << ", ";

      std::vector<VesicleInfo> infos = vgenome.infos();
      std::vector<double> fits = vgenome.fits();

      for (unsigned j = 0; j < infos.size(); j++) {
        std::cerr << setw(8) << fixed << setprecision (4) << right << infos[j].edgeAvgIntensity   << ", "
                  << setw(8) << fixed << setprecision (4) << right << infos[j].edgeDeviation      << ", "
                  << setw(8) << fixed << setprecision (4) << right << infos[j].centerAvgIntensity << ", "
                  << setw(8) << fixed << setprecision (4) << right << infos[j].centerDeviation    << ", "
                  << setw(2) << right << infos[j].edgePx   << ", "
                  << setw(2) << right << infos[j].centerPx << ", "
                  << setw(9) << fixed << setprecision (4) << right << fits[j];
        if (j != fits.size() - 1)
          std::cerr << ", ";
      }
      std::cerr << std::endl;

    }

    img->drawSelectedVesicles((char*)(std::string(GAEDAConfig::handle()->getProblemData()) + ".res.png").c_str());

    // Free image memory
    delete img;

    // Free found vesicles memory
    for (unsigned i = 0; i < found_vesicles.size(); i++)
      delete found_vesicles[i];
    found_vesicles.clear();
  }

  return true;
}


extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}


extern "C" const char* describeProblem (void) {
  return "CajalBBP: vesicles detection problem.";
}


void getCirclePoints (int x, int y, int radius, int width, NeighborhoodIteratorType* nit, std::vector<int>& points_in, std::vector<int>& points_out) {
  // Initialize the image region handler
  ImageType::IndexType index;
  index[0] = x;
  index[1] = y;
  nit->SetLocation(index);

  Circle* circle = new Circle(radius, width);
  int cStat = 0, dim = 2*radius, diameter = dim+1, in = dim-width, aux;

  for (int h=0; h<diameter; h++)
    if (circle->isMarked(0, h))
      points_out.push_back(h);

  for (int i=1; i<dim; i++) {
    for (int j=0; j<diameter; j++) {
      aux = circle->isMarked(i, j);

      switch(cStat) {
        case 0:
               if (aux==1)              { points_out.push_back(i*diameter+j); cStat=1; }
          else if (aux==2)              { points_out.push_back(i*diameter+j); cStat=3; }
          break;
        case 1:
               if (aux==0)              { points_out.push_back(i*diameter+j); cStat=2; }
          else if (aux==1)              { points_out.push_back(i*diameter+j);          }
          else                          { points_out.push_back(i*diameter+j); cStat=3; }
          break;
        case 2:
               if (aux==0)              { points_out.push_back(i*diameter+j);          }
          else if ((i<width)||(i>in))   { points_out.push_back(i*diameter+j); cStat=7; }
          else                          { points_out.push_back(i*diameter+j); cStat=3; }
          break;
        case 3:
               if ((i==width)||(i==in)) { points_out.push_back(i*diameter+j); cStat=5; }
          else if (aux==0)              { points_in .push_back(i*diameter+j); cStat=4; }
          else                          { points_out.push_back(i*diameter+j);          }
          break;
        case 4:
               if (aux==0)              { points_in .push_back(i*diameter+j);          }
          else                          { points_out.push_back(i*diameter+j); cStat=5; }
          break;
        case 5:
               if (aux==0)              { points_out.push_back(i*diameter+j); cStat=6; }
          else if (aux==1)              { points_out.push_back(i*diameter+j);          }
          else                          { points_out.push_back(i*diameter+j); cStat=7; }
          break;
        case 6:
               if (aux==0)              { points_out.push_back(i*diameter+j);          }
          else                          { points_out.push_back(i*diameter+j); cStat=7; }
          break;
        case 7:
              if (aux==1)               { points_out.push_back(i*diameter+j);          }
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
      points_out.push_back(dim*diameter+k);

  delete circle;
}


bool getEdgeInfo(const Vesicle& ves, double& distAllOUT, double& distAllIN, double& distNOUT, double& distNIN) {
  // Open the image again for our own purposes
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(GAEDAConfig::handle()->getProblemData());

  try {
    reader->Update();
  }
  catch (itk::ExceptionObject &err) {
    std::cerr << "Error: the input image '" << GAEDAConfig::handle()->getProblemData() << "' could not be read." << std::endl;
    std::cerr << err << std::endl;
  }

  // Create the Neighborhodd iterator
  NeighborhoodIteratorType::RadiusType rt;
  rt.Fill(ves.getRadius());
  NeighborhoodIteratorType* nit = new NeighborhoodIteratorType(rt, reader->GetOutput(), reader->GetOutput()->GetRequestedRegion());

  // Points in the inside part of circle and parts in the outside part of the circle
  std::vector<int> points_in;
  std::vector<int> points_out;
  getCirclePoints(ves.getX(), ves.getY(), ves.getRadius(), ves.getWidth(), nit, points_in, points_out);

  // Initialize clearest and darkest points and average intensity of the outside part of the circle
  int darkestOUT = points_out[0], clearestOUT = points_out[0];
  double avgAllOUT = 0.0;

  for (std::vector<int>::const_iterator it = points_out.begin(); it != points_out.end(); it++) {
    double px = nit->GetPixel(*it);
    avgAllOUT += px;

    if (px < nit->GetPixel(darkestOUT))
      darkestOUT = *it;

    if (px > nit->GetPixel(clearestOUT))
      clearestOUT = *it;
  }
  avgAllOUT /= (double) points_out.size();

  // Get X and Y coordinates of the clearest point of the outside part of the circle
  int clearestOUTX = nit->GetIndex(clearestOUT)[0];
  int clearestOUTY = nit->GetIndex(clearestOUT)[1];

  // Find a set of closest points to the clearest point in the outside part of the circle
  std::vector<int> closestToClearestOUT;
  int closestToClearestOUTDist = 99999;

  for (std::vector<int>::const_iterator it = points_out.begin(); it != points_out.end(); it++) {
    int absX = nit->GetIndex(*it)[0];
    int absY = nit->GetIndex(*it)[1];

    int dist = 0;
    dist += abs(clearestOUTX - absX);
    dist += abs(clearestOUTY - absY);

    if (dist < closestToClearestOUTDist && dist > 0) {
      closestToClearestOUTDist = dist;
      closestToClearestOUT.clear();
      closestToClearestOUT.push_back(*it);
    }
    else if (dist == closestToClearestOUTDist) {
      closestToClearestOUT.push_back(*it);
    }
  }

  // Find a set of closest points to the clearesr point in the inside part of the circle
  std::vector<int> closestToClearestIN;
  int closestToClearestINDist = 99999;

  // Compute average intensity of the inside at the same time
  double avgAllIN = 0.0;

  for (std::vector<int>::const_iterator it = points_in.begin(); it != points_in.end(); it++) {
    int absX = nit->GetIndex(*it)[0];
    int absY = nit->GetIndex(*it)[1];
    double px = nit->GetPixel(*it);

    avgAllIN += px;

    int dist = 0;
    dist += abs(clearestOUTX - absX);
    dist += abs(clearestOUTY - absY);

    if (dist < closestToClearestINDist) {
      closestToClearestINDist = dist;
      closestToClearestIN.clear();
      closestToClearestIN.push_back(*it);
    }
    else if (dist == closestToClearestINDist) {
      closestToClearestIN.push_back(*it);
    }
  }
  avgAllIN /= (double) points_in.size();

  // Compute average intensity of the closest points to the clearest point in the outside part of the circle
  double avgNeighborhoodOUT = 0.0;
  for (std::vector<int>::const_iterator it = closestToClearestOUT.begin(); it != closestToClearestOUT.end(); it++)
    avgNeighborhoodOUT += nit->GetPixel(*it);
  avgNeighborhoodOUT /= (double) closestToClearestOUT.size();

  // Compute average intensity of the closest points to the clearest point in the inside part of the circle
  double avgNeighborhoodIN = 0.0;
  for (std::vector<int>::const_iterator it = closestToClearestIN.begin(); it != closestToClearestIN.end(); it++)
    avgNeighborhoodIN += nit->GetPixel(*it);
  avgNeighborhoodIN /= (double) closestToClearestIN.size();

  // Compute distances between clearest point intensity and all the other computed average intensities
  distAllOUT = fabs(nit->GetPixel(clearestOUT) - avgAllOUT);
  distAllIN  = fabs(nit->GetPixel(clearestOUT) - avgAllIN);
  distNIN    = fabs(nit->GetPixel(clearestOUT) - avgNeighborhoodIN);
  distNOUT   = fabs(nit->GetPixel(clearestOUT) - avgNeighborhoodOUT);

  std::cout << "Dist to All OUT: " << distAllOUT << std::endl;
  std::cout << "Dist to All IN: "  << distAllIN  << std::endl;
  std::cout << std::endl;
  std::cout << "Dist to Neighborhood OUT: " << distNOUT << std::endl;
  std::cout << "Dist to Neighborhood IN: " << distNIN << std::endl;

  return true;
}
