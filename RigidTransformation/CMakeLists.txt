PROJECT(RigidTransformation)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(RigidTransformation main.cxx)
TARGET_LINK_LIBRARIES(RigidTransformation ${ITK_LIBRARIES})