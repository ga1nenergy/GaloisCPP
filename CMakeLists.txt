cmake_minimum_required(VERSION 3.14)
project(GaloisCPP VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
#include(GNUInstallDirs)
#include_directories("src" "include")
#
#add_library(galoiscpp SHARED src/Algorithms/gfalgorithms.cpp
#        src/GaloisField/galoisfield.cpp src/GaloisField/irreducible_poly.hpp
#        src/GFelement/gfelement.cpp
#        src/Arith/modular_arith.cpp
#        src/Arith/polynomial_arith.cpp
#        src/GFpoly/gfpoly.cpp
#        src/Algorithms/gfalgorithms.cpp)
#
#set_target_properties(galoiscpp PROPERTIES
#        VERSION ${PROJECT_VERSION}
#        PUBLIC_HEADER include/galoiscpp.h)
#target_include_directories(galoiscpp PRIVATE .)
#
#install(TARGETS galoiscpp
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include_directories("src" "include" "api")
add_subdirectory(src/GaloisField)
add_subdirectory(src/GFelement)
add_subdirectory(src/Arith)
add_subdirectory(src/GFpoly)
add_subdirectory(src/Algorithms)
add_subdirectory(src/Aux)

add_executable(GaloisCPP main.cpp)
target_link_libraries(GaloisCPP Aux Algorithms GFpoly GFelement GaloisField Arith)