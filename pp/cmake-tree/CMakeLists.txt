cmake_minimum_required(VERSION 3.17)

project(discrete_pp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### INCLUSION OF Qt
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
#-- Qt6::Core Qt6::Widgets Qt6::PrintSupport are required by qcustomplot
find_package (Qt6 COMPONENTS Core Widgets PrintSupport REQUIRED)
qt_standard_project_setup()

### INCLUSION OF Numerical Libs
find_package (Eigen3 REQUIRED)
find_package (nlohmann_json REQUIRED)

### INCLUSION OF Project Files
set(FORM_DIR "${CMAKE_SOURCE_DIR}/form")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(LIB_DIR "${CMAKE_SOURCE_DIR}/lib")
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")

include_directories(${FORM_DIR})
include_directories(${INCLUDE_DIR})
include_directories(${LIB_DIR})
include_directories(${SOURCE_DIR})

file(GLOB_RECURSE SOURCES
    "${FORM_DIR}/*.ui"
    "${INCLUDE_DIR}/*.h"
    "${INCLUDE_DIR}/*.hpp"
    "${LIB_DIR}/*.cpp"
    "${SOURCE_DIR}/*.cpp"
)

### SETUP EXECUTABLE
set(CMAKE_AUTOUIC_SEARCH_PATHS ${CMAKE_AUTOUIC_SEARCH_PATHS} ${FORM_DIR})

add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${FORM_DIR})
target_include_directories(main PRIVATE ${INCLUDE_DIR})
target_include_directories(main PRIVATE ${LIB_DIR})
target_include_directories(main PRIVATE ${SOURCE_DIR})

target_link_libraries (main PRIVATE
  Qt6::Core
  Qt6::Widgets
  Qt6::PrintSupport
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
)

# target_compile_options(main PUBLIC -O3 -mtune=native -march=native -mfpmath=both)
target_compile_options(main PUBLIC -O3 -mtune=native)