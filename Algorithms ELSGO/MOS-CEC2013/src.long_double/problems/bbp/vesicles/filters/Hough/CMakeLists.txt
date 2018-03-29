# Find and use ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories(${ITK_INCLUDE_DIRS})
link_directories(${ITK_LIBRARY_DIRS})

ADD_EXECUTABLE(Hough Hough.cxx)

TARGET_LINK_LIBRARIES(Hough ${ITK_LIBRARIES})
