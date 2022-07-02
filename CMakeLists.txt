cmake_minimum_required(VERSION 3.23)

project(rtree CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(rtree INTERFACE)
target_include_directories(rtree INTERFACE include)

option(RTREE_TESTS OFF)
if(RTREE_TESTS)
  add_subdirectory(test)
endif()

set(SRC_FILE ${CMAKE_BINARY_DIR}/compile_commands.json)
set(DST_FILE ${CMAKE_SOURCE_DIR}/compile_commands.json)
if(EXISTS ${DST_FILE})
  file(REMOVE ${DST_FILE})
endif()
if(EXISTS ${SRC_FILE})
  file(CREATE_LINK ${SRC_FILE} ${DST_FILE} SYMBOLIC)
endif()
