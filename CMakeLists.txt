cmake_minimum_required(VERSION 3.10)

# set the project name
project(huerto)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Boost REQUIRED)
find_package(Kokkos)
find_package(Schnek REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# add the executable
add_executable(huerto_test
    io/table_data_source.cpp
    maths/random.cpp
    tables/table_lookup.cpp
    tests/main.cpp
    tests/maths/runge_kutta1d.cpp
    tests/maths/test_interpolate1d.cpp
    tests/maths/test_interpolate2d.cpp
    tests/maths/vector1d.cpp
    tests/maths/vector2d.cpp
    tests/maths/vector3d.cpp
    tests/io/test_table_data_source.cpp
    tests/tables/test_table_lookup.cpp
)

target_compile_definitions(huerto_test PRIVATE HUERTO_ONE_DIM)

function(setoptions target)
    target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
    target_include_directories(${target} PUBLIC ${HDF5_INCLUDE_DIRS})

    target_link_libraries(${target} ${MPI_C_LIBRARIES})
    target_link_libraries(${target} ${HDF5_LIBRARIES})
    target_link_libraries(${target} ${Schnek_LIBRARIES})
    target_link_libraries(${target} schnek)
endfunction()

setoptions(huerto_test)

message("MPI ${MPI_INCLUDE_PATH} | ${MPI_C_LIBRARIES}")
message("HDF5 ${HDF5_INCLUDE_DIRS} | ${HDF5_LIBRARIES}")
