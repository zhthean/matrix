#include <gtest/gtest.h>
#include "exceptions/vector_exceptions.h"
#include "vector.h"

using namespace supmath;

// Define test types for single-type tests
using VectorTypes = ::testing::Types<int, float, double>;

template<typename T>
class VectorSingleTypeTest : public ::testing::Test
{
protected:
	const double   test_scalar           = 2.0;

	Vector<T>      test_vector           = Vector<T>({3, -4, 5});
	Vector<T>      expected_scalar_multi = Vector<T>({6, -8, 10});
	Vector<double> expected_unit_vector  = Vector<double>({0.424264, -0.565685, 0.707107});
	double         expected_magnitude    = 7.0710678118654755;
};

TYPED_TEST_SUITE_P(VectorSingleTypeTest);

TYPED_TEST_P(VectorSingleTypeTest, ScalarMultiplicationRight)
{
	auto result = this->test_vector * this->test_scalar;
	ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(VectorSingleTypeTest, ScalarMultiplicationLeft)
{
	auto result = this->test_scalar * this->test_vector;
	ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(VectorSingleTypeTest, Magnitude)
{
	auto result = this->test_vector.magnitude();
	EXPECT_NEAR(result, this->expected_magnitude, 1e-6);
}

TYPED_TEST_P(VectorSingleTypeTest, UnitVector)
{
	auto result = this->test_vector.unitVector();
	for ( size_t i = 0; i < result.getSize(); ++i ) {
		EXPECT_NEAR(result[i], this->expected_unit_vector[i], 1e-6);
	}
}

TYPED_TEST_P(VectorSingleTypeTest, OutOfBoundsAccess)
{
	EXPECT_THROW(this->test_vector[10], std::out_of_range);
	EXPECT_THROW(this->test_vector[-10], std::out_of_range);
}

// Define test types for dual-type tests
using VectorDualTypes = ::testing::Types<std::tuple<int, double>, std::tuple<float, double>>;

template<typename TupleType>
class VectorDualTypeTest : public ::testing::Test
{
protected:
	using T                                                 = typename std::tuple_element<0, TupleType>::type;
	using U                                                 = typename std::tuple_element<1, TupleType>::type;

	Vector<T>                        first_vector           = Vector<T>({1, 2, 3});
	Vector<U>                        second_vector          = Vector<U>({4, 5, 6});

	Vector<std::common_type_t<T, U>> expected_addition      = Vector<std::common_type_t<T, U>>({5, 7, 9});
	Vector<std::common_type_t<T, U>> expected_subtraction   = Vector<std::common_type_t<T, U>>({-3, -3, -3});
	double                           expected_dot_product   = 32.0;
	Vector<std::common_type_t<T, U>> expected_cross_product = Vector<std::common_type_t<T, U>>({-3, 6, -3});
};

TYPED_TEST_SUITE(VectorDualTypeTest, VectorDualTypes);

TYPED_TEST(VectorDualTypeTest, Addition)
{
	auto result = this->first_vector + this->second_vector;
	ASSERT_EQ(result, this->expected_addition);
}

TYPED_TEST(VectorDualTypeTest, Subtraction)
{
	auto result = this->first_vector - this->second_vector;
	ASSERT_EQ(result, this->expected_subtraction);
}

TYPED_TEST(VectorDualTypeTest, DotProduct)
{
	auto result = this->first_vector.dotProduct(this->second_vector);
	EXPECT_NEAR(result, this->expected_dot_product, 1e-6);
}

TYPED_TEST(VectorDualTypeTest, CrossProduct)
{
	auto result = this->first_vector.crossProduct(this->second_vector);
	ASSERT_EQ(result, this->expected_cross_product);
}

TYPED_TEST(VectorDualTypeTest, DimensionMismatch)
{
	Vector<T> mismatched_vector({1, 2});
	EXPECT_THROW(this->first_vector + mismatched_vector, VectorInconsistentDimensionException);
	EXPECT_THROW(this->first_vector - mismatched_vector, VectorInconsistentDimensionException);
	EXPECT_THROW(this->first_vector.dotProduct(mismatched_vector), VectorInconsistentDimensionException);
	EXPECT_THROW(this->first_vector.crossProduct(mismatched_vector), UnsupportVectorSizeException);
}

TYPED_TEST(VectorDualTypeTest, ZeroMagnitudeException)
{
	Vector<T> zero_vector({0, 0, 0});
	EXPECT_THROW(zero_vector.unitVector(), ZeroMagnitudeException);
	EXPECT_THROW(zero_vector.angleBetween(this->second_vector), ZeroMagnitudeException);
}
