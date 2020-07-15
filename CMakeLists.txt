cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Gpgmp VERSION 0.1 LANGUAGES CUDA CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_STANDARD 14)

enable_language(CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-10.2/bin/nvcc")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -lcnpy -O3 -std=c++14 -arch=sm_75 --expt-relaxed-constexpr -rdc=true -lcudadevrt")

find_package(Qt5Test REQUIRED)
find_package(Boost REQUIRED)
find_package(CUDA REQUIRED)

include_directories(${Boost_INCLUDE_DIR})

include_directories(include)

file(
    GLOB SOURCES
    src/*.cpp
    )

file(
    GLOB CUDA_SOURCES
    src/multiparticle.cu
    )

enable_testing(true)

function(CUDA_CONVERT_FLAGS EXISTING_TARGET)
    get_property(old_flags TARGET ${EXISTING_TARGET} PROPERTY INTERFACE_COMPILE_OPTIONS)
    if(NOT "${old_flags}" STREQUAL "")
        string(REPLACE ";" "," CUDA_flags "${old_flags}")
        set_property(TARGET ${EXISTING_TARGET} PROPERTY INTERFACE_COMPILE_OPTIONS
            "$<$<BUILD_INTERFACE:$<COMPILE_LANGUAGE:CXX>>:${old_flags}>$<$<BUILD_INTERFACE:$<COMPILE_LANGUAGE:CUDA>>:-Xcompiler=${CUDA_flags}>"
            )
    endif()
endfunction()

macro (register_test name source)
    add_executable(${name} ${source} ${SOURCES})
    target_link_libraries(${name} PRIVATE Qt5::Test)
    target_link_libraries(${name} LINK_PUBLIC ${Boost_LIBRARIES})
    add_test(${name} ${name})
endmacro ()

macro (cuda_register_test name source)
    cuda_add_executable(${name} ${source} ${CUDA_SOURCES})
    #target_link_libraries(${name} PRIVATE Qt5::Test)
    qt5_use_modules(${name} Test)
    #target_link_libraries(${name} LINK_PUBLIC ${Boost_LIBRARIES})
    add_test(${name} ${name})
endmacro ()

register_test(SsaTest tests/SsaTest.cpp)
register_test(MultiparticleTest tests/MultiparticleTest.cpp)
register_test(GpgmpTest tests/GpgmpTest.cpp)

cuda_register_test(CudaMultiparticleTest tests/MultiparticleGpuTest.cpp)
