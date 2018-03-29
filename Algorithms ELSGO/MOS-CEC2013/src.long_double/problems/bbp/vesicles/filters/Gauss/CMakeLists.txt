# Find and use ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories(${ITK_INCLUDE_DIRS})
link_directories(${ITK_LIBRARY_DIRS})

ADD_EXECUTABLE(Gauss Gauss.cxx)

TARGET_LINK_LIBRARIES(Gauss ${ITK_LIBRARIES})
