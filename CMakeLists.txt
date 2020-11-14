set(WIN32 1)

cmake_minimum_required(VERSION 3.2)
project(CSC417)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(APPLE)
    set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
    set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
    set(OpenMP_C "${CMAKE_C_COMPILER}")
    set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})

    include_directories("${OPENMP_INCLUDES}")
    link_directories("${OPENMP_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")
endif()
if(UNIX AND NOT APPLE)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
endif()

### libIGL options: choose between header only and compiled static library
### Replace libigl/extern/glad/include/glad/glad.h with
### https://raw.githubusercontent.com/alecjacobson/computer-graphics-shader-pipeline/master/glad/include/glad/glad.h
set(LIBIGL ${PROJECT_SOURCE_DIR}/libigl)
set(LIBIGL_INCLUDE_DIR ${LIBIGL}/include)
option(LIBIGL_USE_STATIC_LIBRARY        "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL               "Use OpenGL"            ON)
option(LIBIGL_WITH_PNG                  "Use PNG"               ON)
option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"              ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use ImGui"             ON)

find_package(LIBIGL REQUIRED QUIET)

set(EIGENROOT "${PROJECT_SOURCE_DIR}/libigl/external/eigen")


add_subdirectory(A1)
add_subdirectory(A2)
add_subdirectory(A3)
add_subdirectory(A4)
add_subdirectory(A5)