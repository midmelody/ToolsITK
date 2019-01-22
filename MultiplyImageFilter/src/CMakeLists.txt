cmake_minimum_required(VERSION 3.10.2)

project( MultiplyImageByScalar )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

add_executable( MultiplyImageByScalar main.cxx )
target_link_libraries( MultiplyImageByScalar ${ITK_LIBRARIES} )

