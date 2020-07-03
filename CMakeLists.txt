cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Gpgmp VERSION 0.1)

include_directories(include)

file(
    GLOB SOURCES
    src/*.cpp
    )

add_executable(gpgmp ${SOURCES})
