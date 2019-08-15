cmake_minimum_required(VERSION 3.12)
project(decimate)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenCV)
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(decimate main.cpp Pair.h Pair.h Kdtree.h)

target_link_libraries(decimate ${OpenCV_LIBS} )
