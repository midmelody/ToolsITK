PROJECT(ResampleIsotropicMask)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(ResampleIsotropicMask main.cxx)
TARGET_LINK_LIBRARIES(ResampleIsotropicMask ${ITK_LIBRARIES})
