find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include(${PROJECT_SOURCE_DIR}/cmake/add_gtest.cmake)

ADD_GTEST(
    NAME test_matrix
    FILES test_matrix_values.cpp test_matrix_single_type.cpp test_matrix_dual_types.cpp
    LIBRARY MyMatrix
    ADDITIONAL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/test_data
)
