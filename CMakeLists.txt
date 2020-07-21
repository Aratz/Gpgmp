cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Gpgmp VERSION 0.1 LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Qt5Test REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

include_directories(include)

file(
    GLOB SOURCES
    src/*.cpp
    )

enable_testing(true)

macro (register_test name source)
    add_executable(${name} ${source} ${SOURCES})
    target_link_libraries(${name} PRIVATE Qt5::Test)
    target_link_libraries(${name} LINK_PUBLIC ${Boost_LIBRARIES})
    add_test(${name} ${name})
endmacro ()

register_test(SsaTest tests/SsaTest.cpp)
register_test(MultiparticleTest tests/MultiparticleTest.cpp)
register_test(GpgmpTest tests/GpgmpTest.cpp)

register_test(MultiparticleCpuBenchmark tests/MultiparticleCpuBenchmark.cpp)
