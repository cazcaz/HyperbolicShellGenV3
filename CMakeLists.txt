cmake_minimum_required(VERSION 3.0.0)
project(ShellGeneration VERSION 0.1.0)

include_directories(LBFGSpp)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_BUILD_TYPE Debug)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fno-omit-frame-pointer" )

add_library(LibsModule
     LineGraphOutputter.cpp
     PolarGraphOutputter.cpp
     circleGen.cpp
     shellGen.cpp
     energyFunction.cpp
     shellParams.h
     shellNaming.cpp
     batchGen.cpp
     surface.cpp
     radialSurface.cpp
     LBFGS.h
     LBFGSB.h
     LBFGSpp/Param.h
     LBFGSpp/BFGSMat.h
     LBFGSpp/Cauchy.h
     LBFGSpp/SubspaceMin.h
     LBFGSpp/LineSearchMoreThuente.h
     LBFGSpp/LineSearchBacktracking.h
     LBFGSpp/LineSearchBracketing.h
     LBFGSpp/LineSearchNocedalWright.h
     LBFGSpp/BKLDLT.h
)

add_executable(ShellGeneration main.cpp)


target_link_libraries(ShellGeneration LibsModule)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
