cmake_minimum_required(VERSION 2.6)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

if(NOT(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/contrib/nsbas-playground/CMakeLists.txt"))
  message(FATAL_ERROR "nsbas-playground submodule not found. You probably want to run:\n$ git submodule init\n$ git submodule update")
endif()
add_subdirectory(contrib/nsbas-playground nsbas-playground EXCLUDE_FROM_ALL)
add_custom_target(nsbas-playground ALL DEPENDS nsb_preview_int nsb_preview_unw nsb_preview_slc)
