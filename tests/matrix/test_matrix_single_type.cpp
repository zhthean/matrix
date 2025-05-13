#include <gtest/gtest.h>

#include "matrix.h"

using namespace supmath;

using MatrixType = ::testing::Types<int, float, double>;

template<typename T>
class MatrixSingleTypeTest : public ::testing::Test
{
protected:
	const double       test_scalar           = 3;
	const unsigned int chosen_index          = 0;
	const unsigned int target_index          = 2;

	Matrix<T>          test_matrix           = Matrix<T>({
    {3,  5, -1},
    {0, -5, -4},
    {8,  1,  3}
  });

	double             expected_determinant  = -233;

	Matrix<T>          expected_scalar_multi = Matrix<T>({
    { 9,  15,  -3},
    { 0, -15, -12},
    {24,   3,   9}
  });

	Matrix<T>          expected_transposed   = Matrix<T>({
    { 3,  0, 8},
    { 5, -5, 1},
    {-1, -4, 3}
  });

	Matrix<double>     expected_inversed     = Matrix<double>({
    { 11 / 233.0,  16 / 233.0,  25 / 233.0},
    { 32 / 233.0, -17 / 233.0, -12 / 233.0},
    {-40 / 233.0, -37 / 233.0,  15 / 233.0}
  });

	Matrix<T>          expected_row_swapped  = Matrix<T>({
    {8,  1,  3},
    {0, -5, -4},
    {3,  5, -1}
  });

	Matrix<T>          expected_col_swapped  = Matrix<T>({
    {-1,  5, 3},
    {-4, -5, 0},
    { 3,  1, 8}
  });
};
TYPED_TEST_SUITE_P(MatrixSingleTypeTest);

TYPED_TEST_P(MatrixSingleTypeTest, MatrixRightScalarMulti)
{
	auto result = this->test_matrix * this->test_scalar;

	ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixLeftScalarMulti)
{
	auto result = this->test_scalar * this->test_matrix;

	ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixTranspose)
{
	auto result = this->test_matrix.transpose();

	ASSERT_EQ(result, this->expected_transposed);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixDeterminant)
{
	auto result = this->test_matrix.determinant();

	ASSERT_EQ(result, this->expected_determinant);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixInverse)
{
	auto result = this->test_matrix.inverse();

	ASSERT_EQ(result, this->expected_inversed);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixRowSwap)
{
	this->test_matrix.swap(this->chosen_index, this->target_index, 'r');

	ASSERT_EQ(this->test_matrix, this->expected_row_swapped);
}

TYPED_TEST_P(MatrixSingleTypeTest, MatrixColSwap)
{
	this->test_matrix.swap(this->chosen_index, this->target_index, 'c');

	ASSERT_EQ(this->test_matrix, this->expected_col_swapped);
}

REGISTER_TYPED_TEST_SUITE_P(
    MatrixSingleTypeTest, MatrixRightScalarMulti, MatrixLeftScalarMulti, MatrixTranspose, MatrixDeterminant,
    MatrixInverse, MatrixRowSwap, MatrixColSwap
);
INSTANTIATE_TYPED_TEST_SUITE_P(MSingleTypeTest, MatrixSingleTypeTest, MatrixType);
