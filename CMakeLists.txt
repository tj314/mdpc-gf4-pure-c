cmake_minimum_required(VERSION 3.20)
project(mdpc-gf4 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_FLAGS "-Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# add_definitions(-DWRITE_WEIGHTS)

if(CMAKE_BUILD_TYPE MATCHES Release)
    message("release mode")
    add_definitions(-DNDEBUG)
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

set(SOURCES main.c src/gf4.h src/gf4.c src/gf4_poly.h src/gf4_poly.c src/contexts.h src/contexts.c src/random.c src/random.h src/enc.c src/enc.h src/dec.c src/dec.h src/utils.h src/utils.c src/tests.c src/tests.h)

add_executable(mdpc-gf4 ${SOURCES})
# target_link_libraries(mdpc-gf4 pthread)