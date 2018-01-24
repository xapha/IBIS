cmake_minimum_required(VERSION 2.8)

project(IBIS)

#compiler options
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package( OpenCV REQUIRED )

set(CMAKE_CXX_STANDARD 14)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

set(sourceList
    utils.cpp
    ibis.cpp
    main.cpp

)

set(headerList
    ibis.h
    utils.h

)

set(includeDirList ${includeDirList} $ENV{OPENCV_DIR}/include ./ )
include_directories(/usr/local/include/opencv)
include_directories(/usr/local/include/opencv2)

set(libDirList ${libDirList} $ENV{OPENCV_DIR}/x64/vc14/lib )

include_directories(${includeDirList})
link_directories(${libDirList})

add_executable( ${PROJECT_NAME} ${sourceList} ${headerList} )
target_link_libraries( ${PROJECT_NAME} ${OpenCV_LIBS} )
