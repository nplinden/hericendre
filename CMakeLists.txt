cmake_minimum_required(VERSION 3.22)
project(Hericendre VERSION 0.1)
set(CMAKE_CXX_STANDARD 20)
set(TARGET Hericendre)
set(HERICENDRE_SOURCE_LIST src/nuclide.cpp
        src/chain.cpp
        src/decay.cpp
        src/cramsolver.cpp
        src/source.cpp
        src/fission.cpp
        src/utils.cpp
        src/results.cpp
        src/model.cpp
        src/decaysolver.cpp
)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(inc)

find_package(HighFive CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)

add_executable(Hericendre src/main.cpp ${HERICENDRE_SOURCE_LIST})
target_link_libraries(Hericendre PRIVATE fmt::fmt HighFive Eigen3::Eigen pugixml yaml-cpp::yaml-cpp tomlplusplus::tomlplusplus)