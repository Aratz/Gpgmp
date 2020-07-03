cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Gpgmp VERSION 0.1)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_STANDARD 14)

find_package(Qt5Test REQUIRED)

include_directories(include)

file(
    GLOB SOURCES
    src/*.cpp
    )

enable_testing(true)

add_executable(SsaTest tests/SsaTest.cpp ${SOURCES})
target_link_libraries(SsaTest PRIVATE Qt5::Test)
add_test(SsaTest SsaTest)
