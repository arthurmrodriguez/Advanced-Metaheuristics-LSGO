/*=========================================================================
  Programa:  Filtro Otsu
  Modulo:    Filtro Otsu
  Fecha:     11/03/2010

=========================================================================*/

#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
  {
    std::cerr << "Uso: " << argv[0];
    std::cerr << " imagenEntrada ImagenSalida ";  
    std::cerr << " valorInterior valorExterior "; 
    std::cerr << " umbralInferior umbralSuperior "  << std::endl; 
    return EXIT_FAILURE;
  }

  typedef  unsigned char  tipoPixelEntrada;
  typedef  unsigned char  tipoPixelSalida;

  typedef itk::Image <tipoPixelEntrada, 2> imagenEntrada;
  typedef itk::Image <tipoPixelSalida, 2> imagenSalida;

//  typedef itk::OtsuThresholdImageFilter
//   <imagenEntrada, imagenSalida> tipoFiltro;

  typedef itk::BinaryThresholdImageFilter
    <imagenEntrada, imagenSalida> tipoFiltro;

  typedef itk::ImageFileReader <imagenEntrada> lectorImagen;
  typedef itk::ImageFileWriter <imagenSalida> escritorImagen;

  lectorImagen::Pointer pLector = lectorImagen::New();
  tipoFiltro::Pointer pFiltro = tipoFiltro::New();
  escritorImagen::Pointer pEscritor = escritorImagen::New();

  pEscritor->SetInput(pFiltro->GetOutput());
  pLector->SetFileName(argv[1]);
  pFiltro->SetInput(pLector->GetOutput());

  const tipoPixelSalida valorExterior = atoi(argv[3]);
  const tipoPixelSalida valorInterior  = atoi(argv[4]);
  const tipoPixelSalida umbralInferior = atoi(argv[5]);
  const tipoPixelSalida umbralSuperior  = atoi(argv[6]);

  pFiltro->SetOutsideValue(valorExterior);
  pFiltro->SetInsideValue(valorInterior);

  pFiltro->SetLowerThreshold(umbralInferior);
  pFiltro->SetUpperThreshold(umbralSuperior);

//  pFiltro->SetNumberOfHistogramBins(128);
  pFiltro->Update();

//  int threshold = pFiltro->GetThreshold();
//  std::cout << "Umbral = " << threshold << std::endl;

  pEscritor->SetFileName(argv[2]);
  pEscritor->Update();

  return EXIT_SUCCESS;
}

