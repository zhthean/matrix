cmake_minimum_required (VERSION 3.29)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project("supplement_math")

add_library(MyMatrix "src/matrix.cpp" "src/utils/matrix_func.cpp")

target_include_directories(MyMatrix
	PUBLIC ${PROJECT_SOURCE_DIR}/include
	PRIVATE ${PROJECT_SOURCE_DIR}/src/utils
)

# Add source to this project's executable.
add_executable (matrix "src/main.cpp")

target_link_libraries(matrix PRIVATE MyMatrix)

option(ENABLE_TESTS "Enable or disable test cases" OFF)

if (ENABLE_TESTS)
  enable_testing()
  add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
endif()
