cmake_minimum_required(VERSION 3.20)
project(mdpc-gf4 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
    message("release mode")
    add_definitions(-DNDEBUG)
else()
    message("debug mode")
    add_definitions(-DCANRESIZE)
    add_definitions(-D_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

add_executable(mdpc-gf4 main.c src/gf4_poly.c src/gf4_poly.h src/gf4.h src/gf4.c)