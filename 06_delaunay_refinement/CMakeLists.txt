cmake_minimum_required(VERSION 3.16)
# SET(CMAKE_BUILD_TYPE release)
project(main)

set(DFETCHCONTENT_SOURCE_DIR_LIBIGL ${CMAKE_CURRENT_SOURCE_DIR}/../externals/libigl)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# include polyscope
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../externals/polyscope" "polyscope")

# Include project src files
include_directories("../src/")

# Add your project files
file(GLOB SRC_FILES *.cpp ../src/*.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../externals/polyscope/deps/args/") # grab args.hxx from polyscope
target_link_libraries(${PROJECT_NAME} PUBLIC polyscope igl::core)
target_compile_options(${PROJECT_NAME} PUBLIC "-fno-omit-frame-pointer")
