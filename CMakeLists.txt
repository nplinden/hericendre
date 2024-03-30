cmake_minimum_required(VERSION 3.22)
project(Hericendre VERSION 0.1)
set(TARGET Hericendre)
find_package (Eigen3 3.4 REQUIRED NO_MODULE)
find_package(pugixml)
find_package(fmt)

include_directories(inc)
add_executable(${TARGET} src/main.cpp src/nuclide.cpp src/chain.cpp src/decay.cpp src/nreaction.cpp src/solver.cpp src/source.cpp src/fission.cpp src/utils.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/inc)
target_link_libraries (${PROJECT_NAME} Eigen3::Eigen pugixml fmt)

