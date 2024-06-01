cmake_minimum_required(VERSION 3.22)
project(Hericendre VERSION 0.1)
set(TARGET Hericendre)
set(HERICENDRE_SOURCE_LIST src/nuclide.cpp
                           src/chain.cpp
                           src/decay.cpp
                           src/solver.cpp
                           src/source.cpp
                           src/fission.cpp
                           src/utils.cpp
                           src/results.cpp
                           src/input.cpp
                           src/decaysolver.cpp)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_subdirectory(vendor/eigen)
# find_package (Eigen3 3.4 REQUIRED NO_MODULE)
find_package(pugixml)
find_package(fmt)
find_package(yaml-cpp)
include_directories(inc)

add_library(${TARGET} SHARED ${HERICENDRE_SOURCE_LIST})
target_link_libraries(${TARGET} Eigen3::Eigen pugixml fmt yaml-cpp)

add_executable(${TARGET}Exe src/main.cpp)
target_link_libraries(${TARGET}Exe Hericendre)
set_property(TARGET ${TARGET}Exe PROPERTY OUTPUT_NAME ${TARGET})

add_subdirectory(vendor/pybind11)
pybind11_add_module(pyHericendre src/python/pyHericendre.cpp src/python/nuclide.cpp src/python/chain.cpp)
target_link_libraries(pyHericendre PRIVATE Hericendre)
set_property(TARGET pyHericendre PROPERTY OUTPUT_NAME hericendre)

#
# target_include_directories(${TARGET} PUBLIC ${PROJECT_BINARY_DIR}/inc ${})
# target_link_libraries (${TARGET} Eigen3::Eigen pugixml fmt yaml-cpp nuclide)
# target_include_directories(${TARGET}_exe PUBLIC ${PROJECT_BINARY_DIR}/inc ${})
# target_link_libraries (${TARGET}_exe Eigen3::Eigen pugixml fmt yaml-cpp)
#
# target_link_libraries(nuclide PRIVATE ${PROJECT_NAME})
