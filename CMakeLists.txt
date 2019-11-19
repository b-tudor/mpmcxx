# CMakeLists.txt : CMake project for MPMC++

project(mpmc++)
cmake_minimum_required (VERSION 3.5.1)

execute_process(COMMAND bash "-c" "git rev-list HEAD| wc -l |sed 's: ::g'" VERBATIM OUTPUT_VARIABLE REV)
add_definitions(-DVERSION=${REV})
message("Compiling For ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ggdb -Wall")

option(MPI "Use MPI to parallelize the calculations (requires MPI)" OFF)

set(LIB m)
set(INCLUDE src)

set(SRC 
	src/Atom.cpp
	src/Fugacity.cpp
	src/Molecule.cpp
	src/Output.cpp
	src/PeriodicBoundary.cpp
	src/Quaternion.cpp
	src/Rando.cpp
	src/SafeOps.cpp
	src/SimulationControl.Gibbs.cpp
	src/SimulationControl.PathIntegral.cpp
	src/SimulationControl.cpp
	src/System.Averages.cpp
	src/System.Cavity.cpp
	src/System.Energy.cpp
	src/System.Histogram.cpp
	src/System.MonteCarlo.cpp
	src/System.MPI.cpp
	src/System.Output.cpp
	src/System.Pairs.cpp
	src/System.cpp
	src/Vector3D.cpp
	src/main.cpp
)


if(MPI)
    message("-- MPI Enabled")
    find_package(MPI REQUIRED)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "-- MPI not found! Exiting ...")
    endif()
    include_directories(${MPI_INCLUDE_PATH})
    add_definitions(" -D_MPI ")
    set(LIB ${LIB} ${MPI_CXX_LIBRARIES})
else()
    message("-- MPI Disabled")
endif()


# Add source to this project's executable
add_executable(mpmc++ ${SRC})

# Add required libraries to the linker list
target_link_libraries(${PROJECT_NAME} ${LIB})
