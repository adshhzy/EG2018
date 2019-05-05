project(fcm)


SET(EIGEN_INCLUDE_DIRS "/opt/local/include/")

SET(GUROBI_INCLUDE_DIRS "/Library/gurobi811/mac64/include/")
SET(GUROBI_LIB_DIRS "/Library/gurobi811/mac64/lib/")


SET(SUITESPARSE_LIB_DIR "/usr/local/Cellar/suite-sparse/5.2.0/lib/")


SET(GUROBI_LIB gurobi_c++ gurobi81)
SET(SUITESPARSE_LIB cholmod colamd)
SET(SRC_DIR "src")




SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 ")
include_directories(${EIGEN_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIRS} ${SRC_DIR}) 
LINK_DIRECTORIES(${GUROBI_LIB_DIRS} ${SUITESPARSE_LIB_DIR})

add_subdirectory(lib/triangle) 
add_subdirectory(src/constrainOpt) 
add_subdirectory(src/mymesh) 
add_executable(fcm src/main.cpp)
#target_link_libraries(fcm Opt Geo ${GUROBI_LIB} ${TRIANGLE_LIB} ${SUITESPARSE_LIB}) 
target_link_libraries(fcm Opt Geo Triangle ${GUROBI_LIB} ${SUITESPARSE_LIB}) 