#include <tuple>

#include <gtest/gtest.h>

#include "exceptions/matrix_exceptions.h"
#include "matrix.h"
#include "matrix_test_data.h"

class MatrixAdditionPositiveTest
    : public ::testing::TestWithParam<
          std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};

TEST_P(MatrixAdditionPositiveTest, MatrixAddition) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = first_test_martix + second_test_matrix;

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MAdditionPositive, MatrixAdditionPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x2, expected_addition_2x2),
        std::make_tuple(sample_double_3x3, sample_int_3x3, expected_addition_3x3),
        std::make_tuple(sample_double_4x4, sample_int_4x4, expected_addition_4x4),
        std::make_tuple(sample_double_7x7, sample_int_7x7, expected_addition_7x7),
        std::make_tuple(sample_double_2x3, sample_int_2x3, expected_addition_2x3),
        std::make_tuple(sample_double_3x2, sample_int_3x2, expected_addition_3x2)
    )
);

class MatrixAdditionNegativeTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>>> {};
TEST_P(MatrixAdditionNegativeTest, MatrixAddition) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());

  EXPECT_THROW(first_test_martix + second_test_matrix, supp_math::MatrixMismatchException);
}

INSTANTIATE_TEST_SUITE_P(
    MAdditionNegative, MatrixAdditionNegativeTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x3), std::make_tuple(sample_double_3x3, sample_int_4x4),
        std::make_tuple(sample_double_2x3, sample_int_3x2)
    )
);

class MatrixSubtractionPositiveTest
    : public ::testing::TestWithParam<
          std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};

TEST_P(MatrixSubtractionPositiveTest, MatrixSubtraction) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = first_test_martix - second_test_matrix;

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MSubtractionPositive, MatrixSubtractionPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x2, expected_subtraction_2x2),
        std::make_tuple(sample_double_3x3, sample_int_3x3, expected_subtraction_3x3),
        std::make_tuple(sample_double_4x4, sample_int_4x4, expected_subtraction_4x4),
        std::make_tuple(sample_double_7x7, sample_int_7x7, expected_subtraction_7x7),
        std::make_tuple(sample_double_2x3, sample_int_2x3, expected_subtraction_2x3),
        std::make_tuple(sample_double_3x2, sample_int_3x2, expected_subtraction_3x2)
    )
);

class MatrixSubtractionNegativeTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>>> {};
TEST_P(MatrixSubtractionNegativeTest, MatrixSubtraction) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());

  EXPECT_THROW(first_test_martix - second_test_matrix, supp_math::MatrixMismatchException);
}

INSTANTIATE_TEST_SUITE_P(
    MSubtractionNegative, MatrixSubtractionNegativeTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x3), std::make_tuple(sample_double_3x3, sample_int_4x4),
        std::make_tuple(sample_double_2x3, sample_int_3x2)
    )
);

class MatrixMultiplicationPositiveTest
    : public ::testing::TestWithParam<
          std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};

TEST_P(MatrixMultiplicationPositiveTest, MatrixMultiplication) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = first_test_martix * second_test_matrix;

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MMultiplicationPositive, MatrixMultiplicationPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x2, expected_multiplication_2x2),
        std::make_tuple(sample_double_3x3, sample_int_3x3, expected_multiplication_3x3),
        std::make_tuple(sample_double_4x4, sample_int_4x4, expected_multiplication_4x4),
        std::make_tuple(sample_double_7x7, sample_int_7x7, expected_multiplication_7x7),
        std::make_tuple(sample_double_2x3, sample_int_3x2, expected_multiplication_2x3_3x2),
        std::make_tuple(sample_double_3x2, sample_int_2x3, expected_multiplication_3x2_2x3),
        std::make_tuple(sample_double_3x2, sample_int_2x2, expected_multiplication_3x2_2x2),
        std::make_tuple(sample_double_3x3, sample_int_3x2, expected_multiplication_3x3_3x2)
    )
);

class MatrixMultiplicationNegativeTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>>> {};
TEST_P(MatrixMultiplicationNegativeTest, MatrixMultiplication) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());

  EXPECT_THROW(first_test_martix * second_test_matrix, supp_math::MatrixMismatchException);
}

INSTANTIATE_TEST_SUITE_P(
    MMultiplicationNegative, MatrixMultiplicationNegativeTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_3x3), std::make_tuple(sample_double_3x2, sample_int_3x2),
        std::make_tuple(sample_double_2x3, sample_int_4x4)
    )
);

class MatrixLeftScalarMultiPositiveTest
    : public ::testing::TestWithParam<std::tuple<double, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};
TEST_P(MatrixLeftScalarMultiPositiveTest, MatrixLeftScalarMulti) {
  auto test_scalar = std::get<0>(GetParam());
  auto test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = test_scalar * test_matrix;

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MLeftScalarPositive, MatrixLeftScalarMultiPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_scalar, sample_int_3x3, expected_scalar_multi_3x3),
        std::make_tuple(sample_scalar, sample_int_2x3, expected_scalar_multi_2x3)
    )
);

class MatrixRightScalarMultiPositiveTest
    : public ::testing::TestWithParam<std::tuple<double, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};
TEST_P(MatrixRightScalarMultiPositiveTest, MatrixRightScalarMulti) {
  auto test_scalar = std::get<0>(GetParam());
  auto test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = test_matrix * test_scalar;

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MRightScalarPositive, MatrixRightScalarMultiPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_scalar, sample_int_3x3, expected_scalar_multi_3x3),
        std::make_tuple(sample_scalar, sample_int_2x3, expected_scalar_multi_2x3)
    )
);

class MatrixElementWisePositiveTest
    : public ::testing::TestWithParam<
          std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>, supp_math::Matrix<double>>> {};

TEST_P(MatrixElementWisePositiveTest, MatrixElementWise) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result = first_test_martix.element_wise_product(second_test_matrix);

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MElementWisePositive, MatrixElementWisePositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x2, expected_element_wise_2x2),
        std::make_tuple(sample_double_3x3, sample_int_3x3, expected_element_wise_3x3),
        std::make_tuple(sample_double_4x4, sample_int_4x4, expected_element_wise_4x4),
        std::make_tuple(sample_double_7x7, sample_int_7x7, expected_element_wise_7x7),
        std::make_tuple(sample_double_2x3, sample_int_2x3, expected_element_wise_2x3),
        std::make_tuple(sample_double_3x2, sample_int_3x2, expected_element_wise_3x2)
    )
);

class MatrixElementWiseNegativeTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<int>>> {};
TEST_P(MatrixElementWiseNegativeTest, MatrixElementWise) {
  auto first_test_martix = std::get<0>(GetParam());
  auto second_test_matrix = std::get<1>(GetParam());

  EXPECT_THROW(first_test_martix.element_wise_product(second_test_matrix), supp_math::MatrixMismatchException);
}

INSTANTIATE_TEST_SUITE_P(
    MElementWiseNegative, MatrixElementWiseNegativeTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, sample_int_2x3), std::make_tuple(sample_double_3x3, sample_int_4x4),
        std::make_tuple(sample_double_2x3, sample_int_3x2)
    )
);

class MatrixTransposePositiveTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<double>>> {};

TEST_P(MatrixTransposePositiveTest, MatrixTranspose) {
  auto test_martix = std::get<0>(GetParam());
  auto expected_result = std::get<1>(GetParam());

  auto result = test_martix.transpose();

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MTransposePositive, MatrixTransposePositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_3x3, expected_transposed_double_3x3),
        std::make_tuple(sample_double_3x2, expected_transposed_double_3x2),
        std::make_tuple(sample_double_2x3, expected_transposed_double_2x3)
    )
);

class MatrixDeterminantPositiveTest : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, double>> {};

TEST_P(MatrixDeterminantPositiveTest, MatrixDeterminant) {
  auto test_martix = std::get<0>(GetParam());
  auto expected_result = std::get<1>(GetParam());

  auto result = test_martix.determinant();

  EXPECT_DOUBLE_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MDeterminantPositive, MatrixDeterminantPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, expected_determinant_double_2x2),
        std::make_tuple(sample_double_3x3, expected_determinant_double_3x3),
        std::make_tuple(sample_double_4x4, expected_determinant_double_4x4),
        std::make_tuple(sample_double_7x7, expected_determinant_double_7x7),
        std::make_tuple(sample_singular_double_3x3, expected_determinant_singular_double_3x3)
    )
);

class MatrixDeterminantNegativeTest : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>>> {};

TEST_P(MatrixDeterminantNegativeTest, MatrixDeterminant) {
  auto test_martix = std::get<0>(GetParam());

  EXPECT_THROW(test_martix.determinant(), supp_math::MatrixNonSquareException);
}

INSTANTIATE_TEST_SUITE_P(
    MDeterminantNegative, MatrixDeterminantNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_3x2), std::make_tuple(sample_double_2x3))
);

class MatrixInversePositiveTest
    : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>, supp_math::Matrix<double>>> {};

TEST_P(MatrixInversePositiveTest, MatrixInverse) {
  auto test_martix = std::get<0>(GetParam());
  auto expected_result = std::get<1>(GetParam());

  auto result = test_martix.inverse();

  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MInversePositive, MatrixInversePositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2x2, expected_inversed_2x2),
        std::make_tuple(sample_double_3x3, expected_inversed_3x3),
        std::make_tuple(sample_double_4x4, expected_inversed_4x4),
        std::make_tuple(sample_double_7x7, expected_inversed_7x7)
    )
);

class MatrixInverseNegativeTest : public ::testing::TestWithParam<std::tuple<supp_math::Matrix<double>>> {};

TEST_P(MatrixInverseNegativeTest, MatrixInverse) {
  auto test_martix = std::get<0>(GetParam());

  EXPECT_THROW(test_martix.inverse(), supp_math::MatrixNonSquareException);
}

INSTANTIATE_TEST_SUITE_P(
    MInverseNegative, MatrixInverseNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_3x2), std::make_tuple(sample_double_2x3))
);

TEST(MatrixInverseSingularTest, MatrixInverse) {
  EXPECT_THROW(sample_singular_double_3x3.inverse(), supp_math::MatrixNonInvertibleException);
}

class MatrixSwapPositiveTest
    : public ::testing::TestWithParam<
          std::tuple<supp_math::Matrix<double>, unsigned int, unsigned int, char, supp_math::Matrix<double>>> {};

TEST_P(MatrixSwapPositiveTest, MatrixSwap) {
  auto test_martix = std::get<0>(GetParam());
  auto chosen_index = std::get<1>(GetParam());
  auto target_index = std::get<2>(GetParam());
  auto axis = std::get<3>(GetParam());
  auto expected_result = std::get<4>(GetParam());

  test_martix.swap(chosen_index, target_index, axis);

  EXPECT_EQ(test_martix, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    MatrixSwapPositiveTest, MatrixSwapPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_3x3, sample_chosen_index, sample_target_index, 'r', expected_swap_row_3x3),
        std::make_tuple(sample_double_3x2, sample_chosen_index, sample_target_index, 'r', expected_swap_row_3x2),
        std::make_tuple(sample_double_2x3, sample_chosen_index, sample_target_index, 'r', expected_swap_row_2x3),
        std::make_tuple(sample_double_3x3, sample_chosen_index, sample_target_index, 'c', expected_swap_col_3x3),
        std::make_tuple(sample_double_3x2, sample_chosen_index, sample_target_index, 'c', expected_swap_col_3x2),
        std::make_tuple(sample_double_2x3, sample_chosen_index, sample_target_index, 'c', expected_swap_col_2x3)
    )
);
