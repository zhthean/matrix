#include <gtest/gtest.h>

#include "vector.h"

using namespace supmath;

using VectorTypes = ::testing::Types<int, float, double>;

template<typename T>
class VectorSingleTypeTest : public ::testing::Test
{
protected:
  const double   test_scalar           = 2.0;

  Vector<T>      test_vector           = Vector<T>({3, -4, 5});
  Vector<T>      expected_scalar_multi = Vector<T>({6, -8, 10});
  Vector<double> expected_unit_vector  = Vector<double>({0.424264069, -0.565685425, 0.707106781});
  double         expected_magnitude    = 7.0710678118654755;
};

TYPED_TEST_SUITE_P(VectorSingleTypeTest);

TYPED_TEST_P(VectorSingleTypeTest, VectorRigthScalarMulti)
{
  auto result = this->test_vector * this->test_scalar;
  ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(VectorSingleTypeTest, VectorLeftScalarMulti)
{
  auto result = this->test_scalar * this->test_vector;
  ASSERT_EQ(result, this->expected_scalar_multi);
}

TYPED_TEST_P(VectorSingleTypeTest, VectorMagnitude)
{
  auto result = this->test_vector.magnitude();
  EXPECT_NEAR(result, this->expected_magnitude, 1e-6);
}

TYPED_TEST_P(VectorSingleTypeTest, VectorUnitVector)
{
  auto result = this->test_vector.unitVector();
  ASSERT_EQ(result, this->expected_unit_vector);
}

REGISTER_TYPED_TEST_SUITE_P(
    VectorSingleTypeTest, VectorRigthScalarMulti, VectorLeftScalarMulti, VectorMagnitude, VectorUnitVector
);
INSTANTIATE_TYPED_TEST_SUITE_P(VSingleTypeTest, VectorSingleTypeTest, VectorTypes);
