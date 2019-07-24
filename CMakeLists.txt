cmake_minimum_required(VERSION 3.14)
project(GaloisCPP)

set(CMAKE_CXX_STANDARD 11)

include_directories("src" "src/GaloisField" "src/GFelement" "src/Arith" "src/GFpoly" "src/Algorithms" "src/Codec")
add_subdirectory(src/GaloisField)
add_subdirectory(src/GFelement)
add_subdirectory(src/Arith)
add_subdirectory(src/GFpoly)
add_subdirectory(src/Algorithms)
add_subdirectory(src/Codec)

add_executable(GaloisCPP main.cpp)
target_link_libraries(GaloisCPP Codec Algorithms GFpoly GFelement GaloisField Arith)