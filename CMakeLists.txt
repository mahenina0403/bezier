cmake_minimum_required(VERSION 3.5)

project(compare)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(include)

file(GLOB source_files src/*.cpp)
file(GLOB include_files include/*.h)

add_executable(compare ${source_files} ${include_files})