cmake_minimum_required(VERSION 2.6)

PROJECT(LungDeformReg)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)




SET( SOURCE_FILES

  ImageType.h
  IMGIO.h
  ConfigReader.h
  ConfigReader.cxx
  ireg-ssd-vm.cxx
  ap.cpp
  )

ADD_EXECUTABLE(LungDeformReg ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(LungDeformReg ${ITK_LIBRARIES} )

SET( SOURCE_FILES

  ImageType.h
  IMGIO.h
  imageWarp.cxx
  )
ADD_EXECUTABLE(imageDeform ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(imageDeform ${ITK_LIBRARIES} )

