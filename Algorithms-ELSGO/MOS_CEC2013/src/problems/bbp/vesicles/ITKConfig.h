#ifndef ITKCONFIG
#define ITKCONFIG

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageLinearConstIteratorWithIndex.h"

typedef float PixelType;
typedef itk::Image<PixelType, 2> ImageType;
typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
typedef itk::ImageFileReader<ImageType> ReaderType;

#endif
