cmake_minimum_required(VERSION 3.5)

project(compare)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(include)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
find_package(MPFR REQUIRED )
if (NOT MPFR_FOUND)
	message(STATUS "This project requires the mpfr library, and will not be compiled.")
	return()  
endif()
include_directories(${MPFR_INCLUDE_DIRS})

file(GLOB source_files src/*.cpp)
file(GLOB include_files include/*.h)

add_executable(compare ${source_files} ${include_files})

target_link_libraries(compare ${MPFR_LIBRARIES})