INCLUDE(${CMAKE_SOURCE_DIR}/BinaryThinningSkeleton/IJMacros.txt)

PROJECT(BinaryThinning3D)

cmake_minimum_required(VERSION 2.8)


ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
  ITK
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#IF(BUILD_TESTING)
  SET(CurrentExe "BinaryThinning3D")
  ADD_EXECUTABLE(${CurrentExe} BinaryThinning3D.cxx) 
  TARGET_LINK_LIBRARIES(${CurrentExe} ${ITK_LIBRARIES})
#ENDIF(BUILD_TESTING)

