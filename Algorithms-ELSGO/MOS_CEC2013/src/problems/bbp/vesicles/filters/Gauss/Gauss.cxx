#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 4 ) 
  { 
    std::cerr << "Uso: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  outputImageFile  sigma ";
    std::cerr <<  std::endl;
    return EXIT_FAILURE;
  }
  
 const int dim = 2;
 typedef float tipoPixelEntrada;
 typedef float tipoPixelSalida;
  
 typedef itk::Image <tipoPixelEntrada, dim> imagenEntrada;
 typedef itk::Image <tipoPixelSalida, dim> imagenSalida;

 typedef itk::ImageFileReader <imagenEntrada> lectorImagen;

 typedef itk::RecursiveGaussianImageFilter
    <imagenEntrada, imagenSalida> tipoFiltro;


 lectorImagen::Pointer pLector = lectorImagen::New();
 pLector->SetFileName(argv[1]);

 tipoFiltro::Pointer filtroX = tipoFiltro::New();
 tipoFiltro::Pointer filtroY = tipoFiltro::New();
  
 filtroX->SetDirection(0);   // 0 --> X
 filtroY->SetDirection(1);   // 1 --> Y

 filtroX->SetOrder(tipoFiltro::ZeroOrder);
 filtroY->SetOrder(tipoFiltro::ZeroOrder);

 filtroX->SetNormalizeAcrossScale(false);
 filtroY->SetNormalizeAcrossScale(false);

 filtroX->SetInput(pLector->GetOutput());
 filtroY->SetInput(filtroX->GetOutput() );

 const double sigma = atof(argv[3]);

 filtroX->SetSigma(sigma);
 filtroY->SetSigma(sigma);

 filtroY->Update();

 typedef  unsigned char tipoPixelEscritor;
 typedef itk::Image<tipoPixelEscritor, dim> escritorImagen;
 typedef itk::RescaleIntensityImageFilter
    <imagenSalida, escritorImagen> RescalaTipoFiltro;

 RescalaTipoFiltro::Pointer pResc = RescalaTipoFiltro::New();
 
 pResc->SetOutputMinimum(0);
 pResc->SetOutputMaximum(255);
 typedef itk::ImageFileWriter<escritorImagen>  tipoEscritor;
 tipoEscritor::Pointer pEscritor = tipoEscritor::New();
 pEscritor->SetFileName(argv[2]);

 pResc->SetInput(filtroY->GetOutput());
 pEscritor->SetInput(pResc->GetOutput());
 pEscritor->Update();

 return EXIT_SUCCESS;
}
