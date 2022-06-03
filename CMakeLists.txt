cmake_minimum_required(VERSION 3.12)
project(PCSaftSuperAncillary)
set(CMAKE_CXX_STANDARD 17)

find_program(CCACHE_FOUND "ccache")
set(CCACHE_SUPPORT ON CACHE BOOL "Enable ccache support")
if (CCACHE_FOUND AND CCACHE_SUPPORT)
  message(STATUS "Enabling ccache")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" # GNU is GNU GCC
      OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # without this compiler messages in `make` backend would be uncolored
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fdiagnostics-color=auto")
  endif()
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "ccache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "ccache")
endif()

set(CHEBTOOLS_STATIC_LIBRARY true CACHE INTERNAL "Build ChebTools as a static library")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(fitmain main.cpp)
target_link_libraries(fitmain ChebTools teqpinterface autodiff)

add_executable(fitbench bench.cpp)
target_link_libraries(fitbench ChebTools teqpinterface autodiff Catch2WithMain)

add_executable(fittest test.cpp)
target_link_libraries(fittest ChebTools teqpinterface autodiff) 