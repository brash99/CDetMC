cmake_minimum_required(VERSION 3.10)
project(AnaBarMC)

# Set G4TARGET and G4EXLIB
set(G4TARGET AnaBarMC)
set(G4EXLIB true)

# Check for G4INSTALL and set default value
if(NOT DEFINED G4INSTALL)
    set(G4INSTALL ../../..)
endif()

# Set C++ standard and compiler flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Define custom target for cleanup
add_custom_target(cleanup
        COMMAND rm -f g4*.prim g4*.eps g4*.wrl
        COMMAND rm -f .DAWN_* *~ src/*~ include/*~ macros/*~
        COMMAND rm -f *.so *.d
        )

########################### ROOT #################################

# Check for ROOTSYS and set necessary configurations
if(DEFINED ENV{ROOTSYS})
    if(NOT DEFINED G4UI_USE_ROOT)
        execute_process(
                COMMAND $ENV{ROOTSYS}/bin/root-config --cflags
                OUTPUT_VARIABLE ROOTCPPFLAGS
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4ANALYSIS_USE_ROOT ${ROOTCPPFLAGS}")
        execute_process(
                COMMAND $ENV{ROOTSYS}/bin/root-config --nonew --glibs
                OUTPUT_VARIABLE ROOTLIBS
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        string(REPLACE "-lNew" "" ROOTLIBS ${ROOTLIBS})
        set(ROOTLIBS "${ROOTLIBS} -lMinuit")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOTLIBS}")
    endif()
endif()

# Add subdirectories for 'lib' and 'bin' (assumed in the original Makefile)
add_subdirectory(lib)
add_subdirectory(bin)
