cmake_minimum_required(VERSION 3.22)

project(Hericendre VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Enable FetchContent module
include(FetchContent)

# Fetch Eigen3.4
message("Fetching Eigen 3.4")
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE
  EXCLUDE_FROM_ALL
)
set(BUILD_TESTING OFF)
set(EIGEN_TEST_NOQT ON)
set(EIGEN_BUILD_CMAKE_PACKAGE OFF)
FetchContent_MakeAvailable(eigen)

# Fetch fmt 10.2
message("Fetching fmt 10.2")
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG        10.2.1
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(fmt)


# Fetch pugixml 1.14
message("Fetching pugixml 1.14")
FetchContent_Declare(
  pugixml
  GIT_REPOSITORY https://github.com/zeux/pugixml.git
  GIT_TAG        v1.14
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(pugixml)

# Fetch toml++ 3.4
message("Fetching toml++ 3.4")
FetchContent_Declare(
  tomlplusplus
  GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
  GIT_TAG        v3.4.0
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(tomlplusplus)


# Get HDF5
message(STATUS "Downloading HDF5 v1.14.3")
set(HDF5_ENABLE_Z_LIB_SUPPORT OFF)
set(HDF5_EXTERNALLY_CONFIGURED  1)
set(HDF5_ENABLE_SZIP_SUPPORT  OFF)
set(HDF5_BUILD_EXAMPLES       OFF)
set(HDF5_BUILD_TOOLS          OFF)
set(HDF5_BUILD_UTILS          OFF)
set(HDF5_BUILD_HL_LIB         OFF)
set(BUILD_SHARED_LIBS         OFF)
FetchContent_Declare(hdf5
  GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
  GIT_TAG        hdf5-1_14_3
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(hdf5)

# Fetch HighFive
message("Fetching Highfive 2.10")
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_TEST_SPAN OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_FIND_HDF5 OFF)
FetchContent_Declare(
  highfive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
  GIT_TAG        v2.10.0
)
FetchContent_MakeAvailable(highfive)

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

add_executable(hericendre src/main.cpp ${HERICENDRE_SOURCE_LIST})
target_link_libraries(hericendre PRIVATE fmt::fmt hdf5-static HighFive Eigen3::Eigen pugixml tomlplusplus::tomlplusplus)

install(TARGETS hericendre DESTINATION ${SKBUILD_SCRIPTS_DIR})
