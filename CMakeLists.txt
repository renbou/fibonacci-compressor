cmake_minimum_required(VERSION 3.18)
project(fibonacci_compressor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wredundant-move -Wpessimizing-move")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(phicc_sources
        main.cpp)
list(TRANSFORM phicc_sources PREPEND "src/")

add_executable(fibonacci_compressor ${phicc_sources})

add_subdirectory(src/gcem)
target_link_libraries(fibonacci_compressor PRIVATE gcem)