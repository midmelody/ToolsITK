PROJECT(CropSlice)
# This project is designed to be built outside the Insight source tree.
cmake_minimum_required(VERSION 2.6)
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(CropSlice GetVOI.cxx)
TARGET_LINK_LIBRARIES(CropSlice ${ITK_LIBRARIES})

