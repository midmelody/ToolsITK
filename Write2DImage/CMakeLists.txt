PROJECT(Write2DImage)
# This project is designed to be built outside the Insight source tree.
cmake_minimum_required(VERSION 2.8)
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

ADD_EXECUTABLE(Write2DImage main.cxx)
TARGET_LINK_LIBRARIES(Write2DImage ${ITK_LIBRARIES})