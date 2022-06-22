cmake_minimum_required(VERSION 3.21)

project(discrete_bif LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### INCLUSION OF Numerical Libs
find_package (Eigen3 REQUIRED)
find_package (nlohmann_json REQUIRED)
set(Boost_USE_STATIC_LIBS        ON)  # only find static libs
set(Boost_USE_DEBUG_LIBS        OFF)  # ignore debug libs and
set(Boost_USE_RELEASE_LIBS       ON)  # only find release libs
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED)
include_directories("${Boost_INCLUDE_DIRS}")
link_directories("${Boost_LIBRARY_DIR}")


### INCLUSION OF Project Files
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")

include_directories(${INCLUDE_DIR})
include_directories(${SOURCE_DIR})

file(GLOB_RECURSE SOURCES
    "${INCLUDE_DIR}/*.h"
    "${INCLUDE_DIR}/*.hpp"
    "${SOURCE_DIR}/*.cpp"
)

### SETUP EXECUTABLE
add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${INCLUDE_DIR})
target_include_directories(main PRIVATE ${SOURCE_DIR})

target_link_libraries (main PRIVATE
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
  ${Boost_LIBRARIES}
)

# target_compile_options(main PUBLIC -O3 -mtune=native -march=native -mfpmath=both)
target_compile_options(main PUBLIC -O3 -mtune=native)