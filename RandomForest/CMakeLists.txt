PROJECT( RandomForest )

CMAKE_MINIMUM_REQUIRED(VERSION 2.8) 

FIND_PACKAGE( OpenCV REQUIRED )

ADD_EXECUTABLE( decisiontree decisiontree.cpp )
TARGET_LINK_LIBRARIES( decisiontree ${OpenCV_LIBS} )
ADD_EXECUTABLE( randomforest randomforest.cpp )
TARGET_LINK_LIBRARIES( randomforest ${OpenCV_LIBS} )
ADD_EXECUTABLE( extremerandomforest extremerandomforest.cpp )
TARGET_LINK_LIBRARIES( extremerandomforest ${OpenCV_LIBS} )
ADD_EXECUTABLE( boosttree boosttree.cpp )
TARGET_LINK_LIBRARIES( boosttree ${OpenCV_LIBS} )