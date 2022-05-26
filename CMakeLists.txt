cmake_minimum_required(VERSION 3.12)
project(PCSaftSuperAncillary)
set(CMAKE_CXX_STANDARD 17)

set(CHEBTOOLS_STATIC_LIBRARY true CACHE INTERNAL "Build ChebTools as a static library")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(fitmain main.cpp)
target_link_libraries(fitmain ChebTools teqpinterface autodiff)

# target_include_directories(fastcheb PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/boost_multiprecision")

# target_include_directories(fastcheb
#     PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools/externals/Eigen"
#     PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals"
# )