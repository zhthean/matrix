#include <tuple>
#include <type_traits>

#include <gtest/gtest.h>

#include "matrix.h"

using namespace supmath;

using MatrixTypes = ::testing::Types<
    std::tuple<short, long long>, std::tuple<short, long double>, std::tuple<int, int>, std::tuple<double, int>,
    std::tuple<double, float>, std::tuple<double, double>, std::tuple<long double, double>,
    std::tuple<long double, long double>>;

template<typename TupleType>
class MatrixDualTypesTest : public ::testing::Test
{
protected:
	using T                                                      = typename std::tuple_element<0, TupleType>::type;
	using U                                                      = typename std::tuple_element<1, TupleType>::type;

	Matrix<T>                        first_test_matrix           = Matrix<T>({
    {3,  5, -1},
    {0, -5, -4},
    {8,  1,  3}
  });

	Matrix<U>                        second_test_matrix          = Matrix<U>({
    {  7, -9, 1},
    { 11, -7, 2},
    {-12,  6, 3}
  });

	Matrix<std::common_type_t<T, U>> expected_addition           = Matrix<std::common_type_t<T, U>>({
    {10,  -4,  0},
    {11, -12, -2},
    {-4,   7,  6}
  });

	Matrix<std::common_type_t<T, U>> expected_subtraction        = Matrix<std::common_type_t<T, U>>({
    { -4, 14, -2},
    {-11,  2, -6},
    { 20, -5,  0}
  });

	Matrix<std::common_type_t<T, U>> expected_multiplication     = Matrix<std::common_type_t<T, U>>({
    {88, -68,  10},
    {-7,  11, -22},
    {31, -61,  19}
  });

	Matrix<std::common_type_t<T, U>> expected_element_wise_multi = Matrix<std::common_type_t<T, U>>({
	  { 21, -45, -1},
    {  0,  35, -8},
    {-96,   6,  9}
	});
};
TYPED_TEST_SUITE_P(MatrixDualTypesTest);

TYPED_TEST_P(MatrixDualTypesTest, MatrixAddition)
{
	auto result = this->first_test_matrix + this->second_test_matrix;

	ASSERT_EQ(result, this->expected_addition);
}

TYPED_TEST_P(MatrixDualTypesTest, MatrixSubtraction)
{
	auto result = this->first_test_matrix - this->second_test_matrix;

	ASSERT_EQ(result, this->expected_subtraction);
}

TYPED_TEST_P(MatrixDualTypesTest, MatrixMultiplication)
{
	auto result = this->first_test_matrix * this->second_test_matrix;

	ASSERT_EQ(result, this->expected_multiplication);
}

TYPED_TEST_P(MatrixDualTypesTest, MatrixElementWiseMulti)
{
	auto result = this->first_test_matrix.element_wise_product(this->second_test_matrix);

	ASSERT_EQ(result, this->expected_element_wise_multi);
}

REGISTER_TYPED_TEST_SUITE_P(
    MatrixDualTypesTest, MatrixAddition, MatrixSubtraction, MatrixMultiplication, MatrixElementWiseMulti
);
INSTANTIATE_TYPED_TEST_SUITE_P(MDualTypesTest, MatrixDualTypesTest, MatrixTypes);
