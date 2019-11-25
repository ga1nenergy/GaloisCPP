add_subdirectory(Algorithms)
add_subdirectory(Arith)
add_subdirectory(Auxillary)
add_subdirectory(GaloisField)
add_subdirectory(GFelement)
add_subdirectory(GFpoly)

set(GALOISCPP_LIBS Algorithms Auxillary GFpoly GFelement GaloisField Arith)