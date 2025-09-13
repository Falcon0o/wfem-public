cmake_minimum_required(VERSION 3.10) 
# set(CMAKE_CUDA_ARCHITECTURES "native")
enable_language(CUDA)
### Configuration

include(cmake/options.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/find)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(wfem VERSION 0.0)


### 设置输出颜色
string(ASCII 27 Esc)

# 显示方式  
# 0                终端默认设置
# 1                高亮显示
# 4                使用下划线
# 5                闪烁
# 7                反白显示
# 8                不可见

# 前景色            背景色           颜色
# ---------------------------------------
# 30                40              黑色
# 31                41              红色
# 32                42              绿色
# 33                43              黃色
# 34                44              蓝色
# 35                45              紫红色
# 36                46              青蓝色
# 37                47              白色
set(E "${Esc}[m")
set(G "${Esc}[1;32m")
set(B "${Esc}[1;36m")


# set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/thirdparty)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes)

add_subdirectory(src)

# set(CUDA_NVCC_FLAGS --extended-lambda; -Wext-lambda-captures-this)
# add_library(${PROJECT_NAME} SHARED ${WFEM_SOURCE_FILES})
add_executable(${PROJECT_NAME} main.cpp ${WFEM_SOURCE_FILES})

# set_property(TARGET ${PROJECT_NAME} PROPERTY CUDA_SEPARABLE_COMPILATION ON)


# -expt-relaxed-constexpr 常量表达式函数
target_compile_options(
    ${PROJECT_NAME}
    PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-extended-lambda;-Wext-lambda-captures-this;-arch=sm_75;-expt-relaxed-constexpr>
)

# target_compile_options(
#     ${PROJECT_NAME}
#     PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe;"--diag_suppress=unrecognized_gcc_pragma">
# )


target_compile_options(
    ${PROJECT_NAME}
    PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native;-fconcepts>)


find_package(CUDAToolkit REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE CUDA::cudart cudss)

find_package(ZLIB)
target_link_libraries(${PROJECT_NAME} PRIVATE ZLIB::ZLIB)

find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
# 调用 backtrace_symbol
# target_link_options(${PROJECT_NAME} PRIVATE -rdynamic)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# 使用 NVTX
if (${WFEM_ENABLE_NVTX})
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_ENABLE_NVTX)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/NVTX)
    target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/contrib/NVTX/c/include)
    target_link_libraries(${PROJECT_NAME} PRIVATE nvtx3-cpp)
endif()


# 支持 MPI
if (${WFEM_ENABLE_MPI})
    find_package(MPI REQUIRED)
    target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_ENABLE_MPI)
endif()

message("WFEM_ENABLE_MPI: ${WFEM_ENABLE_MPI}" )
message("MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}" )

# target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third/include)
# target_link_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third/lib)

target_link_libraries(${PROJECT_NAME} PRIVATE metis)
# target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third/lib/libGKlib.a)



# include(spdlog)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/spdlog/)
# target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/contrib/spdlog/include)

# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen3)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/contrib/eigen3)

if (${WFEM_ENABLE_USER_DEFINED_BACKTRACE})
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_ENABLE_USER_DEFINED_BACKTRACE)
endif()

if (${WFEM_CHECK_SLOW})
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_CHECK_SLOW)
endif()

target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_WORKING_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/contrib/rapidjson/include)

# find_package(PETSc REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE PETSc::PETSc)


if (${WFEM_ENABLE_MKL})
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_ENABLE_MKL)
    set(MKL_DIR /opt/intel/oneapi/mkl/latest/lib/cmake/mkl)
    find_package(MKL CONFIG REQUIRED)
    target_link_directories(${PROJECT_NAME} PRIVATE /opt/intel/oneapi/mkl/latest/lib)
    target_link_libraries(${PROJECT_NAME} PRIVATE MKL::MKL)
    
    target_include_directories(${PROJECT_NAME} PRIVATE /opt/intel/oneapi/mkl/latest/include)
endif()

if (DEFINED WFEM_MUMPS_DIR AND EXISTS ${WFEM_MUMPS_DIR})
    target_compile_definitions(${PROJECT_NAME} PRIVATE WFEM_ENABLE_MUMPS)
    # target_link_libraries(${PROJECT_NAME} PRIVATE ${WFEM_MUMPS_DIR}/lib/libmumps_common.so)
    # target_link_directories(${PROJECT_NAME} PRIVATE ${WFEM_MUMPS_DIR}/lib)
    target_link_directories(${PROJECT_NAME} PRIVATE ${WFEM_MUMPS_DIR}/lib)
    target_link_libraries(${PROJECT_NAME} PRIVATE dmumps)

    target_include_directories(${PROJECT_NAME} PRIVATE ${WFEM_MUMPS_DIR}/include)
endif()

target_include_directories(${PROJECT_NAME} PRIVATE /usr/lib/x86_64-linux-gnu/hdf5/serial/include)
target_link_directories(${PROJECT_NAME} PRIVATE /usr/lib/x86_64-linux-gnu/hdf5/serial)
target_link_libraries(${PROJECT_NAME} PRIVATE hdf5)
# MUMPS https://coin-or-tools.github.io/ThirdParty-Mumps/

# AMGX
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/contrib/AMGX/include)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/AMGX/Release/libamgxsh.so)

install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# install(FILES ${WFEM_HEADER_FILES} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)


