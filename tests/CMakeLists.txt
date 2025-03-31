find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test_data)

add_executable(test_matrix "test_matrix_values.cpp" "test_matrix_single_type.cpp" "test_matrix_dual_types.cpp" "../src/matrix.cpp" )

target_link_libraries(
  test_matrix
  GTest::gtest
  GTest::gtest_main
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET test_matrix PROPERTY CXX_STANDARD 20)
endif()

add_test(NAME test_matrix
         COMMAND test_matrix)
