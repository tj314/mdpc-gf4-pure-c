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
elseif (CMAKE_BUILD_TYPE MATCHES Weights)
    message("write weights mode")
    add_definitions(-DWRITE_WEIGHTS)
elseif (CMAKE_BUILD_TYPE MATCHES GJS)
    message("run GJS simulations mode")
    add_definitions(-DGJS)
elseif (CMAKE_BUILD_TYPE MATCHES Iterations)
    message("test iterations mode")
    add_definitions(-DTEST_ITERATIONS)
elseif(CMAKE_BUILD_TYPE MATCHES Testing)
    message("testing mode")
    add_definitions(-D_DEBUG -DRUNTESTS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
else()
    message("debug mode")
    add_definitions(-DCANRESIZE)
    add_definitions(-D_DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

set(SOURCES src/gf4.h src/gf4.c src/gf4_poly.h src/gf4_poly.c src/contexts.h src/contexts.c src/random.c src/random.h src/enc.c src/enc.h src/dec_symbol_flipping.c src/dec.h src/utils.h src/utils.c src/tests.c src/tests.h src/dec_threshold.c src/dec_sf_with_delta.c src/dec_utils.c src/gf4_matrix.c src/gf4_matrix.h src/gf4_array.c src/gf4_array.h)

if(CMAKE_BUILD_TYPE MATCHES GJS)
    add_executable(mdpc-gf4 main-gjs.c ${SOURCES})
else()
    add_executable(mdpc-gf4 main.c ${SOURCES})
endif()

# target_link_libraries(mdpc-gf4 pthread)