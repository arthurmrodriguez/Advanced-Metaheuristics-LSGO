# Find and use Magick++
find_package(ImageMagick COMPONENTS Magick++ REQUIRED)
include_directories(${ImageMagick_Magick++_INCLUDE_DIR})

add_executable(PWMAD PWMAD.cc)

target_link_libraries(PWMAD ${ImageMagick_Magick++_LIBRARY})
