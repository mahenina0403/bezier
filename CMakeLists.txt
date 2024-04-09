cmake_minimum_required(VERSION 3.5)

project(compare)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(include)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

file(GLOB source_files src/*.cpp)
file(GLOB include_files include/*.h)

add_executable(compare ${source_files} ${include_files})

target_link_libraries(compare ${MPFR_LIBRARIES})