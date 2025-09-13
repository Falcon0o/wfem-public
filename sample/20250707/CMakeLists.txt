cmake_minimum_required(VERSION 3.10) 

project(samples VERSION 0.0)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE "Debug")
add_executable(demo main.cpp)
target_link_directories(demo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib)
target_link_libraries(demo PUBLIC wfem)
# target_link_libraries(demo PUBLIC wfem dl)

target_include_directories(demo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
