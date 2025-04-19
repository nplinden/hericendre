cmake_minimum_required(VERSION 3.22)
project(Hericendre VERSION 0.1)
set(CMAKE_CXX_STANDARD 20)
set(TARGET hericendre)
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

find_package(fmt CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(pugixml CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)
find_package(HighFive CONFIG REQUIRED)

add_executable(hericendre src/main.cpp ${HERICENDRE_SOURCE_LIST})
target_link_libraries(hericendre PRIVATE fmt::fmt HighFive Eigen3::Eigen pugixml tomlplusplus::tomlplusplus)

install(TARGETS hericendre DESTINATION ${SKBUILD_SCRIPTS_DIR})
