#  CMakeLists.txt

cmake_minimum_required(VERSION 3.21.0)

################## compiler options ##################
# compiler
set(CMAKE_CXX_COMPILER          ${CXX})
set(CMAKE_Fortran_COMPILER          ${FC})

# compiler flags
set(CMAKE_CXX_FLAGS             ${CXXFLAGS}) 

project(SWDTTI LANGUAGES CXX Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/build/obj)

set(CMAKE_CXX_FLAGS  "-g -Wall")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") 
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O0 -D_GLIBCXX_DEBUG")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -xHost")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -xHost -extend-source")
else ()
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -march=native -ffixed-line-length-none")
  set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
  #set(CMAKE_CXX_FLAGS_RELEASE  "-O0 -D_GLIBCXX_DEBUG")
endif()

set(CMAKE_VERBOSE_MAKEFILE on)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
include_directories(${FFTW_INCLUDE_DIRS})

find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/RF)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/SWD)