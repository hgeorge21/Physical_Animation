cmake_minimum_required(VERSION 3.2)

project(a3-finite-elements-3d)

#Important CMAKE stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#include directories
include_directories(${PROJECT_SOURCE_DIR}/shared/include)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${LIBIGL_INCLUDE_DIR})

#include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB H_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB CPP_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/src/*.cpp)

#build executable 
add_executable(${PROJECT_NAME} main.cpp assignment_setup.h ${H_FILES} ${H_SHARED_FILES} ${CPP_FILES} ${CPP_SHARED_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)


