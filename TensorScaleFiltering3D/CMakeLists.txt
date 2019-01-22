PROJECT(TensorScaleFiltering3D)
# This project is designed to be built outside the Insight source tree.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

# Find ITK.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(TensorScaleFiltering3D main.cxx)
TARGET_LINK_LIBRARIES(TensorScaleFiltering3D ${ITK_LIBRARIES})