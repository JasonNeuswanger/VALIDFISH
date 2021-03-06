cmake_minimum_required(VERSION 3.9)

# GENERAL PROJECT SETTINGS FOR BOTH EXECUTABLE AND PYTHON MODUlE VERSION OF THE MODEL

set(CMAKE_OSX_SYSROOT "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(THREADS_PREFER_PTHREAD_FLAG ON)

project(validfish)

set(SOURCES
    Forager.cpp Forager.h
    EnvironmentInhabitant.cpp
    Geometry.cpp
    Analysis.cpp
    PreyType.cpp PreyType.h
    Swimmer.cpp Swimmer.h
    ManeuverInterpolation.cpp ManeuverInterpolation.h
    Wolf.h Wolf.cpp
    Optimizer.cpp Optimizer.h
    utility.h utility.cpp
    lib/gsl_function_pp.h
    lib/fast-cpp-csv-parser/csv.h
    lib/flat_hash_map.hpp
    ForagerGetters.cpp
    Printing.cpp
    DetectionSubmodel.cpp
    DiscriminationSubmodel.cpp
    Energetics.cpp Bounds.cpp)

find_package(GSL REQUIRED)
find_package(Threads REQUIRED)

set(LIBRARIES_TO_LINK GSL::gsl GSL::gslcblas stdc++fs Threads::Threads)

include_directories("/Users/Jason/Dropbox/Drift Model Project/Calculations/VALIDFISH/VALIDFISH/lib/eigen")

# GENERATE THE EXECUTABLE BINARY VERSION OF THE MODEL

add_executable(validfish main.cpp ${SOURCES} lib/ExecutionTimer.h)
target_link_libraries(validfish ${LIBRARIES_TO_LINK})

# GENERATE THE PYTHON LIBRARY VERSION OF THE MODEL

set(PYBIND11_CPP_STANDARD -std=c++1z) # use experimental c++17
include_directories(/Users/Jason/anaconda/envs/driftmodelenv/include/python3.5m) # for Python.h
add_subdirectory(lib/pybind11)
pybind11_add_module(pyvalidfish python_bindings.cpp ${SOURCES}) # should include all the source files eventually
target_link_libraries(pyvalidfish PRIVATE ${LIBRARIES_TO_LINK})
#target_compile_options(pyvalidfish PRIVATE -g) # USE THIS VERSION TO DEBUG
target_compile_options(pyvalidfish PRIVATE -O3 -DNDEBUG) # USE THIS VERSION FOR RELEASE

# The inclusion of /usr/local/gcc-7.3/lib/libstdc++fs.a in target_link_libraries is to make experimental::filesystem work.

# Note: I used the following tutorial to install the latest gcc compiler, because the version of clang that comes
# with xcode does not support the <filesystem> feature of c++17.
# https://solarianprogrammer.com/2017/05/21/compiling-gcc-macos/

# The "make -j 4" command for compiling GCC itself kept getting stuck with "no such file or directory" for things
# like <bits/stl_tempbuf.h>, but kept trying, got a  different error every time, and eventually it worked.
# Now and up-to-date GCC compiler is installed at:
# /usr/local/gcc-7.3/bin/gcc-7.3
# /usr/local/gcc-7.3/bin/g++-7.3

# Also installed valgrind following directions here:
# https://stackoverflow.com/questions/46050202/building-valgrind-on-os-x
# using ./configure --prefix=/usr/local/gcc-7.3

# Original valgrind options were --leak-check=full --leak-resolution=med  --track-origins=yes
# For me that's "/Users/Jason/Dropbox/Drift Model Project/Calculations/VALIDFISH/VALIDFISH/cmake-build-debug/validfish"
# To profile with callgrind, go into a terminal window to the folder where this program is output
# and use: valgrind --tool callgrind ./validfish
# Pipe output from callgrind_annotate to a text file to see sensible results:
# callgrind_annotate callgrind.out.6943 > callgrind.out.6943.txt

# Configuring Eigen following thse instructions
# https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html
#
