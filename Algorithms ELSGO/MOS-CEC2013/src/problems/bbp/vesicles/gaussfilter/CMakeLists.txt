# Find and use ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories(${ITK_INCLUDE_DIRS})
link_directories(${ITK_LIBRARY_DIRS})

add_executable(Image Image.cxx)

target_link_libraries(Image ${ITK_LIBRARIES})
