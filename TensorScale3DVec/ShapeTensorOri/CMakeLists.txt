PROJECT(ShapeTensorOri)
# This project is designed to be built outside the Insight source tree.
cmake_minimum_required(VERSION 2.8)
# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(ShapeTensorOri main.cxx)
TARGET_LINK_LIBRARIES(ShapeTensorOri ${ITK_LIBRARIES})