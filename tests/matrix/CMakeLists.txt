CREATE_GTEST(
    NAME test_matrix
    FILES test_matrix_values.cpp test_matrix_single_type.cpp test_matrix_dual_types.cpp
    LIBRARY ${LIBRARY_NAME}
    ADDITIONAL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}
)
