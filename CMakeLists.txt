cmake_minimum_required(VERSION 3.14)
project(galoiscpp VERSION 0.9)

set(CMAKE_CXX_STANDARD 17)
include_directories("src" "include")

add_subdirectory(src/Algorithms)
add_subdirectory(src/Arith)
add_subdirectory(src/Aux)
add_subdirectory(src/GaloisField)
add_subdirectory(src/GFelement)
add_subdirectory(src/GFpoly)


add_executable(galoiscpp main.cpp)
target_link_libraries(galoiscpp Algorithms Aux GFpoly GFelement GaloisField Arith)