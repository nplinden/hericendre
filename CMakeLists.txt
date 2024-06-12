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

add_subdirectory(vendor/eigen)
add_subdirectory(vendor/fmt)
add_subdirectory(vendor/pugixml)
add_subdirectory(vendor/yaml-cpp)

include_directories(inc)

add_library(${TARGET} SHARED ${HERICENDRE_SOURCE_LIST})
target_link_libraries(${TARGET} Eigen3::Eigen pugixml fmt yaml-cpp)

add_executable(${TARGET}Exe src/main.cpp)
target_link_libraries(${TARGET}Exe Hericendre)
set_property(TARGET ${TARGET}Exe PROPERTY OUTPUT_NAME ${TARGET})

add_subdirectory(vendor/pybind11)
pybind11_add_module(pyHericendre
        src/python/pyHericendre.cpp
        src/python/nuclide.cpp
        src/python/chain.cpp
        src/python/model.cpp)
target_link_libraries(pyHericendre PRIVATE Hericendre)
set_property(TARGET pyHericendre PROPERTY OUTPUT_NAME hericendre)
set_property(TARGET HericendreExe PROPERTY OUTPUT_NAME hericendre)
