# - Main CMake buildscript for costheta
# Comments in a CMake Script are lines begining with a '#'

# - Set CMake requirements then declare project
cmake_minimum_required(VERSION 4.3)

# Declare project
project(FitModel VERSION 0.1)

# - When Makefiles are generated, output all command lines by default
#   Do this to begin with so it's easy to see what compiler command/flags
#   are used. This can also be done by removing the 'set' command and
#   running make as 'make VERBOSE=1'.
set(CMAKE_VERBOSE_MAKEFILE ON)

# - Don't allow C++ Compiler Vendor Extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# - Set C++ Standard directly to c++17
set(CMAKE_CXX_STANDARD 17)

# - Use our standard set of flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfatal-errors -Wshadow -pedantic")

# - ROOT installation
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

# - Find ROOT
find_package(ROOT COMPONENTS XMLIO MathCore MathMore Minuit2 GenVector RooFit)

# - Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

# - Debug messages to check ROOT dependencies
message(STATUS "ROOT inc dir =${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT lib path=${ROOT_LIBRARIES}")

# - Include dir
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# - Sources
file(GLOB FitModelSrcs ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# ==============================================================================
# Generate ROOT Dictionary for PyROOT
# ==============================================================================
# Requires a LinkDef.h file (see instructions below)
ROOT_GENERATE_DICTIONARY(
    G__FitModel 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/FeldmanCousins.h
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h    
)

# ==============================================================================
# Build the Shared Library (.so)
# ==============================================================================
add_library(FitModel SHARED ${FitModelSrcs} G__FitModel.cxx)

# Link against ROOT and RooFit libraries
target_link_libraries(FitModel PUBLIC ${ROOT_LIBRARIES} ROOT::RooFit ROOT::RooFitCore)
