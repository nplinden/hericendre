cmake_minimum_required(VERSION 3.22)

project(Hericendre VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find system HDF5
find_package(HDF5 REQUIRED COMPONENTS C)

# Try to find system Eigen3, otherwise fetch it
find_package(Eigen3 3.4 QUIET NO_MODULE)

# Enable FetchContent module
include(FetchContent)

if(NOT Eigen3_FOUND)
    # Fetch Eigen3.4
    message("Eigen3 not found, fetching Eigen 3.4")
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
else()
    message("Using system Eigen3: ${Eigen3_DIR}")
endif()

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

# Fetch HighFive
message("Fetching Highfive 3.2.0")
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_TEST_SPAN OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_USE_BOOST OFF)
FetchContent_Declare(
  highfive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
  GIT_TAG        v3.2.0
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
target_link_libraries(hericendre PRIVATE fmt::fmt HDF5::HDF5 HighFive Eigen3::Eigen pugixml tomlplusplus::tomlplusplus)

install(TARGETS hericendre DESTINATION ${SKBUILD_SCRIPTS_DIR})
