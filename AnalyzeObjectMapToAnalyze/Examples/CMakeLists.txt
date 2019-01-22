cmake_minimum_required(VERSION 2.6)

PROJECT(AnalyzeObjectMapExamples)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(${AnalyzeObjectMap_SOURCE_DIR})
INCLUDE_DIRECTORIES(${AnalyzeObjectMapIO_SOURCE_DIR})

SET(TEST_DATA_ROOT     ${AnalyzeObjectMapExamples_SOURCE_DIR}/Data/Input)
SET(TESTING_OUTPUT_DIR ${AnalyzeObjectMapExamples_BINARY_DIR}/Temporary)

INCLUDE_DIRECTORIES( ${AnalyzeObjectMapExamples_SOURCE_DIR} )

ADD_EXECUTABLE( CreatingObjects CreatingObjects.cxx )
TARGET_LINK_LIBRARIES( CreatingObjects ${ITK_LIBRARIES} AnalyzeObjectMapIO AnalyzeObjectMap)

ADD_EXECUTABLE( DisplayingObjectMaps DisplayingObjectMaps.cxx )
TARGET_LINK_LIBRARIES( DisplayingObjectMaps ${ITK_LIBRARIES} AnalyzeObjectMapIO AnalyzeObjectMap) 

ADD_EXECUTABLE( PickOneObjectEntry PickOneObjectEntry.cxx )
TARGET_LINK_LIBRARIES( PickOneObjectEntry ${ITK_LIBRARIES} AnalyzeObjectMapIO AnalyzeObjectMap)


