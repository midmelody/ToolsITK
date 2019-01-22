cmake_minimum_required(VERSION 2.6)

PROJECT(AnalyzeObjectMapTest)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(${AnalyzeObjectMap_SOURCE_DIR})
INCLUDE_DIRECTORIES(${AnalyzeObjectMapIO_SOURCE_DIR})

INCLUDE_DIRECTORIES( ${AnalyzeObjectMapTest_SOURCE_DIR} )

ADD_EXECUTABLE( AnalyzeObjectMapToAnalyze AnalyzeObjectMapTest.cxx )
TARGET_LINK_LIBRARIES(AnalyzeObjectMapToAnalyze ${ITK_LIBRARIES} AnalyzeObjectMapIO AnalyzeObjectMap)
