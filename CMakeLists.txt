cmake_minimum_required(VERSION 3.12)
project(PCSaftSuperAncillary)
set(CMAKE_CXX_STANDARD 17)

if (MSVC)
  add_compile_options(/bigobj)
else ()
endif ()

macro(add_nix_libraries target)
    # See https://stackoverflow.com/a/29871891
    # Add DL and pthreads
    FIND_PACKAGE(Threads REQUIRED)
    find_package(Threads REQUIRED)
    if(THREADS_HAVE_PTHREAD_ARG)
      target_compile_options(${target} PRIVATE "-pthread")
    endif()
    if(CMAKE_THREAD_LIBS_INIT)
      target_link_libraries(${target} PRIVATE "${CMAKE_THREAD_LIBS_INIT}")
    endif()
endmacro()

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
set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp")
# set_property(TARGET ChebTools PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(fitmain main.cpp)
target_link_libraries(fitmain PRIVATE ChebTools PRIVATE teqpinterface PRIVATE autodiff)
if (NOT MSVC)
  ADD_NIX_LIBRARIES(fitmain)
endif()

add_executable(fitbench bench.cpp)
target_link_libraries(fitbench ChebTools teqpinterface autodiff Catch2WithMain)

add_executable(fittest test.cpp)
target_link_libraries(fittest ChebTools teqpinterface autodiff) 

pybind11_add_module(PCSAFTsuperanc "${CMAKE_CURRENT_SOURCE_DIR}/implementation/interface.cpp")
target_compile_definitions(PCSAFTsuperanc PUBLIC -DPYBIND11)
target_link_libraries(PCSAFTsuperanc PUBLIC ChebTools PUBLIC teqpinterface)
# if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
#   include(CheckIPOSupported)
#   check_ipo_supported(RESULT iposupported OUTPUT error)
#   if( iposupported )
#       message(STATUS "IPO / LTO enabled")
#       set_property(TARGET PCSAFTsuperanc PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
#   else()
#       message(STATUS "IPO / LTO not supported: <${error}>")
#   endif()
#   #set_property(TARGET PCSAFTsuperanc PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
# endif()

add_executable(PCSAFTsuperanctest "${CMAKE_CURRENT_SOURCE_DIR}/implementation/interface.cpp")
target_link_libraries(PCSAFTsuperanctest PUBLIC ChebTools)
target_include_directories(PCSAFTsuperanctest PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/Eigen")