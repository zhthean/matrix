cmake_minimum_required (VERSION 3.29)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project("supplement_math")
set(LIBRARY_NAME "SupMath")

#set(CMAKE_CXX_CLANG_TIDY "clang-tidy")

add_library(${LIBRARY_NAME} "src/matrix.cpp" "src/vector.cpp" "src/utils/matrix_funcs.cpp")

target_include_directories(${LIBRARY_NAME}
	PUBLIC ${PROJECT_SOURCE_DIR}/include
	PRIVATE ${PROJECT_SOURCE_DIR}/src
)

option(ENABLE_TESTS "Enable or disable test cases" OFF)

if (ENABLE_TESTS)
  enable_testing()
  add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif()
