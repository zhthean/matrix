#include <tuple>
#include <type_traits>

#include <gtest/gtest.h>

#include "exceptions/vector_exceptions.h"
#include "vector.h"

using namespace supmath;

using VectorTypes = ::testing::Types<std::tuple<int, double>, std::tuple<double, short>, std::tuple<double, double>>;

template<typename TupleType>
class VectorDualTypesTest : public ::testing::Test
{
protected:
  using T                                                    = typename std::tuple_element<0, TupleType>::type;
  using U                                                    = typename std::tuple_element<1, TupleType>::type;

  Vector<T>                        first_test_vector         = Vector<T>({1, 2, 3});
  Vector<U>                        second_test_vector        = Vector<U>({4, 5, 6});

  Vector<std::common_type_t<T, U>> expected_addition         = Vector<std::common_type_t<T, U>>({5, 7, 9});
  Vector<std::common_type_t<T, U>> expected_subtraction      = Vector<std::common_type_t<T, U>>({-3, -3, -3});
  double                           expected_dot_product      = 32.0;
  Vector<std::common_type_t<T, U>> expected_cross_product    = Vector<std::common_type_t<T, U>>({-3, 6, -3});
  double                           expected_angle_between    = 0.225726129;
  double                           expected_distance_between = 5.196152422;
  Vector<double>                   expected_projection_first_onto_second =
      Vector<std::common_type_t<T, U>>({1.662337662, 2.077922078, 2.493506493});
};

TYPED_TEST_SUITE_P(VectorDualTypesTest);

TYPED_TEST_P(VectorDualTypesTest, VectorAddition)
{
  auto result = this->first_test_vector + this->second_test_vector;
  ASSERT_EQ(result, this->expected_addition);
}

TYPED_TEST_P(VectorDualTypesTest, VectorSubtraction)
{
  auto result = this->first_test_vector - this->second_test_vector;
  ASSERT_EQ(result, this->expected_subtraction);
}

TYPED_TEST_P(VectorDualTypesTest, VectorDotProduct)
{
  auto result = this->first_test_vector.dotProduct(this->second_test_vector);
  EXPECT_NEAR(result, this->expected_dot_product, 1e-6);
}

TYPED_TEST_P(VectorDualTypesTest, VectorCrossProduct)
{
  auto result = this->first_test_vector.crossProduct(this->second_test_vector);
  ASSERT_EQ(result, this->expected_cross_product);
}

TYPED_TEST_P(VectorDualTypesTest, VectorAngleBetween)
{
  auto result = this->first_test_vector.angleBetween(this->second_test_vector);
  EXPECT_NEAR(result, this->expected_angle_between, 1e-6);
}

TYPED_TEST_P(VectorDualTypesTest, VectorDistanceBetween)
{
  auto result = this->first_test_vector.distanceBetween(this->second_test_vector);
  EXPECT_NEAR(result, this->expected_distance_between, 1e-6);
}

TYPED_TEST_P(VectorDualTypesTest, VectorProjectionOnto)
{
  auto result = this->first_test_vector.projectionOnto(this->second_test_vector);
  ASSERT_EQ(result, this->expected_projection_first_onto_second);
}

REGISTER_TYPED_TEST_SUITE_P(
    VectorDualTypesTest, VectorAddition, VectorSubtraction, VectorDotProduct, VectorCrossProduct, VectorAngleBetween,
    VectorDistanceBetween, VectorProjectionOnto
);
INSTANTIATE_TYPED_TEST_SUITE_P(VDualTypesTest, VectorDualTypesTest, VectorTypes);
