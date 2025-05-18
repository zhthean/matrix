CREATE_GTEST(
    NAME test_vector
    FILES test_vector_values.cpp test_vector_single_type.cpp test_vector_dual_types.cpp
    LIBRARY ${LIBRARY_NAME}
    ADDITIONAL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}
)
