set (COMPILE_BBP_VESICLES OFF CACHE BOOL "Compile Cajal BBP Vesicles Segmentation Problem")

if (COMPILE_BBP_VESICLES)
  # Find and use Magick++
  find_package(ImageMagick COMPONENTS Magick++ REQUIRED)
  include_directories(${ImageMagick_Magick++_INCLUDE_DIR})

  # Find and use ITK
  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})
  include_directories(${ITK_INCLUDE_DIRS})
  link_directories(${ITK_LIBRARY_DIRS})

  find_path(TOKENIZER_INCLUDE "boost/tokenizer.hpp")
  if ("${TOKENIZER_INCLUDE}" STREQUAL "TOKENIZER_INCLUDE-NOTFOUND")
    MESSAGE(FATAL_ERROR "Couldn't find Boost TOKENIZER header files.")
  endif()

  include_directories(${GAEDA_INCLUDE} ${TOKENIZER_INCLUDE})

  # Add targets
  add_executable(GenCSV GenCSV.cc Vesicle.cc Circle.cc ImageManager.cc Objective.cc BruteForceAlgorithm.cc DrawVesicle.cc)
  target_link_libraries(GenCSV ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})

  add_executable(DrawCSV DrawCSV.cc Vesicle.cc Circle.cc ImageManager.cc ImagePrinter.cc Objective.cc Parser.cc DrawVesicle.cc)
  target_link_libraries(DrawCSV ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})

  add_executable(AdjustValData AdjustValData.cc Vesicle.cc Circle.cc ImageManager.cc ImagePrinter.cc Objective.cc Parser.cc DrawVesicle.cc)
  target_link_libraries(AdjustValData ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})

  add_library(vesicles SHARED vesicles.cc Circle.cc DrawVesicle.cc ImageManager.cc Objective.cc Vesicle.cc)
  target_link_libraries(vesicles ${GAEDA_LIBRARY} ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})
  set_target_properties(vesicles PROPERTIES SUFFIX ".so")

  add_executable(testObjFunc testObjFunc.cc ImageManager.cc Objective.cc Circle.cc Vesicle.cc DrawVesicle.cc)
  target_link_libraries(testObjFunc ${GAEDA_LIBRARY} ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})

  add_executable(evalVesicle evalVesicle.cc ImageManager.cc Objective.cc Circle.cc Vesicle.cc DrawVesicle.cc)
  target_link_libraries(evalVesicle ${GAEDA_LIBRARY} ${LIBS} ${ITK_LIBRARIES} ${ImageMagick_Magick++_LIBRARY})

  add_executable(FilterInArea FilterInArea.cc)
  target_link_libraries(FilterInArea ${ImageMagick_Magick++_LIBRARY})

  add_subdirectory (gaussfilter)
  add_subdirectory (filters/Canny)
  add_subdirectory (filters/Edge)
  add_subdirectory (filters/Gauss)
  add_subdirectory (filters/Hough)
  add_subdirectory (filters/Otsu)
  add_subdirectory (filters/PWMAD)
endif (COMPILE_BBP_VESICLES)
