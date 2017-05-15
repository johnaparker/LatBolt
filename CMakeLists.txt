cmake_minimum_required(VERSION 3.5.2)
project(latbolt)
#set(CMAKE_BUILD_TYPE DEBUG)
enable_testing()
 
#Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

if(${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -Og")
endif()
message(STATUS "FLAGS: ${CMAKE_CXX_FLAGS}")

#Dependencies

#HDF5 library
find_library(h5_lib hdf5)
if(NOT h5_lib)
    message(FATAL_ERROR "hdf5 library not found")
else()
    message(STATUS "hdf5 library found: ${h5_lib}")
endif()

#h5cpp library
find_library(h5cpp_lib h5cpp)
if(NOT h5_lib)
    message(FATAL_ERROR "h5cpp library not found")
else()
    message(STATUS "h5cpp library found: ${h5cpp_lib}")
endif()

#Include the headers
include_directories(include)
 
#Define the sources and headers
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB MAIN_HEADERS "include/*.h" )
#file(GLOB SUB_HEADERS "include/sub/*.h")
 
#Generate the shared library from the sources
add_library(latbolt SHARED ${SOURCES})
 
#Build executable
set (PROJECT_LINK_LIBS latbolt ${h5cpp_lib} ${h5_lib})
link_directories( build )

add_executable(main tests/main.cpp)
target_link_libraries(main ${PROJECT_LINK_LIBS} )
add_test(test1 main)


#Set the location for library installation
install(TARGETS latbolt DESTINATION lib)
install(FILES ${MAIN_HEADERS} DESTINATION include/latbolt)
#install(FILES ${SUB_HEADERS} DESTINATION include/latbolt/sub)
install(FILES include/latbolt.h DESTINATION include)
