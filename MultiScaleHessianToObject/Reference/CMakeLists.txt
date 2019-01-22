cmake_minimum_required(VERSION 2.4)

PROJECT(FOO)

SET(FOO_SRCS
  foo.cxx
  )

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)

IF (USE_ITK_FILE)
  INCLUDE (${USE_ITK_FILE})
ENDIF (USE_ITK_FILE)

ADD_EXECUTABLE(foo ${FOO_SRCS})
TARGET_LINK_LIBRARIES(foo ITKCommon ITKBasicFilters ITKAlgorithms ITKIO ITKStatistics)

