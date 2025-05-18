#include <gtest/gtest.h>

#include "exceptions/vector_exceptions.h"
#include "matrix.h"
#include "vector.h"
#include "vector_test_data.h"

using namespace supmath;

class VectorAdditionPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, Vector<double>>>
{};

TEST_P(VectorAdditionPositiveTest, VectorAddition)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector + second_test_vector;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VAdditionPositive, VectorAdditionPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_addition_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_addition_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_addition_4)
    )
);

class VectorAdditionNegativeTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>>>
{};

TEST_P(VectorAdditionNegativeTest, VectorAddition)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());

  EXPECT_THROW(first_test_vector + second_test_vector, VectorInconsistentDimensionException);
}

INSTANTIATE_TEST_SUITE_P(
    VAdditionNegative, VectorAdditionNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_2, sample_int_3), std::make_tuple(sample_double_3, sample_double_4))
);

class VectorSubtractionPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, Vector<double>>>
{};

TEST_P(VectorSubtractionPositiveTest, VectorSubtraction)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector - second_test_vector;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VSubtractionPositive, VectorSubtractionPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_subtraction_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_subtraction_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_subtraction_4)
    )
);

class VectorSubtractionNegativeTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>>>
{};

TEST_P(VectorSubtractionNegativeTest, VectorSubtraction)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());

  EXPECT_THROW(first_test_vector - second_test_vector, VectorInconsistentDimensionException);
}

INSTANTIATE_TEST_SUITE_P(
    VSubtractionNegative, VectorSubtractionNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_2, sample_int_3), std::make_tuple(sample_double_3, sample_double_4))
);

class VectorLeftScalarMultiTest : public ::testing::TestWithParam<std::tuple<double, Vector<int>, Vector<double>>>
{};

TEST_P(VectorLeftScalarMultiTest, VectorLeftScalarMulti)
{
  auto test_scalar     = std::get<0>(GetParam());
  auto test_vector     = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result          = test_scalar * test_vector;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VLeftScalarMultiTest, VectorLeftScalarMultiTest,
    ::testing::Values(
        std::make_tuple(sample_scalar, sample_int_2, expected_scalar_multi_2),
        std::make_tuple(sample_scalar, sample_int_3, expected_scalar_multi_3),
        std::make_tuple(sample_scalar, sample_int_4, expected_scalar_multi_4)
    )
);

class VectorRightScalarMultiTest : public ::testing::TestWithParam<std::tuple<Vector<int>, double, Vector<double>>>
{};

TEST_P(VectorRightScalarMultiTest, VectorRightScalarMulti)
{
  auto test_scalar     = std::get<0>(GetParam());
  auto test_vector     = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result          = test_vector * test_scalar;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VRightScalarMultiTest, VectorRightScalarMultiTest,
    ::testing::Values(
        std::make_tuple(sample_int_2, sample_scalar, expected_scalar_multi_2),
        std::make_tuple(sample_int_3, sample_scalar, expected_scalar_multi_3),
        std::make_tuple(sample_int_4, sample_scalar, expected_scalar_multi_4)
    )
);

TEST(VectorCrossProductPositiveTest, VectorCrossProductPositive)
{
  auto first_test_vector  = sample_int_3;
  auto second_test_vector = sample_double_3;
  auto expected_result    = expected_cross_product_3;

  auto result             = first_test_vector.crossProduct(second_test_vector);

  ASSERT_EQ(result, expected_result);
}

TEST(VectorCrossProductNegativeTest, VectorCrossProductNegative)
{
  auto first_test_vector  = sample_int_2;
  auto second_test_vector = sample_double_2;

  EXPECT_THROW(first_test_vector.crossProduct(second_test_vector), UnsupportVectorSizeException);
}

class VectorDotProductPositiveTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, double>>
{};

TEST_P(VectorDotProductPositiveTest, VectorDotProduct)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector.dotProduct(second_test_vector);

  EXPECT_NEAR(result, expected_result, 1e-6);
}

INSTANTIATE_TEST_SUITE_P(
    VDotProductPositive, VectorDotProductPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_dot_product_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_dot_product_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_dot_product_4)
    )
);

class VectorDotProductNegativeTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>>>
{};

TEST_P(VectorDotProductNegativeTest, VectorDotProduct)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());

  EXPECT_THROW(first_test_vector.dotProduct(second_test_vector), VectorInconsistentDimensionException);
}

INSTANTIATE_TEST_SUITE_P(
    VDotProductNegative, VectorDotProductNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_2, sample_int_3), std::make_tuple(sample_double_3, sample_double_4))
);

class VectorAngleBetweenPositiveTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, double>>
{};

TEST_P(VectorAngleBetweenPositiveTest, VectorAngleBetween)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector.angleBetween(second_test_vector);

  EXPECT_NEAR(result, expected_result, 1e-6);
}

INSTANTIATE_TEST_SUITE_P(
    VAngleBetweenPositive, VectorAngleBetweenPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_angle_between_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_angle_between_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_angle_between_4)
    )
);

class VectorAngleBetweenNegativeTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>>>
{};

TEST_P(VectorAngleBetweenNegativeTest, VectorAngleBetween)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());

  EXPECT_THROW(first_test_vector.angleBetween(second_test_vector), VectorInconsistentDimensionException);
}

INSTANTIATE_TEST_SUITE_P(
    VAngleBetweenNegative, VectorAngleBetweenNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_2, sample_int_3), std::make_tuple(sample_double_3, sample_double_4))
);

class VectorDistanceBetweenPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, double>>
{};

TEST_P(VectorDistanceBetweenPositiveTest, VectorDistanceBetween)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector.distanceBetween(second_test_vector);

  EXPECT_NEAR(result, expected_result, 1e-6);
}

INSTANTIATE_TEST_SUITE_P(
    VDistanceBetweenPositive, VectorDistanceBetweenPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_distance_between_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_distance_between_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_distance_between_4)
    )
);

class VectorDistanceBetweenNegativeTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>>>
{};

TEST_P(VectorDistanceBetweenNegativeTest, VectorDistanceBetween)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());

  EXPECT_THROW(first_test_vector.distanceBetween(second_test_vector), VectorInconsistentDimensionException);
}

INSTANTIATE_TEST_SUITE_P(
    VDistanceBetweenNegative, VectorDistanceBetweenNegativeTest,
    ::testing::Values(std::make_tuple(sample_double_2, sample_int_3), std::make_tuple(sample_double_3, sample_double_4))
);

class VectorProjectionOntoPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<int>, Vector<double>>>
{};

TEST_P(VectorProjectionOntoPositiveTest, VectorProjectionOnto)
{
  auto first_test_vector  = std::get<0>(GetParam());
  auto second_test_vector = std::get<1>(GetParam());
  auto expected_result    = std::get<2>(GetParam());

  auto result             = first_test_vector.projectionOnto(second_test_vector);

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VProjectionOntoPositive, VectorProjectionOntoPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, sample_int_2, expected_projection_onto_2),
        std::make_tuple(sample_double_3, sample_int_3, expected_projection_onto_3),
        std::make_tuple(sample_double_4, sample_int_4, expected_projection_onto_4)
    )
);

TEST(VectorProjectionOntoNegativeTest, VectorProjectionOntoNegative)
{
  auto first_test_vector  = sample_int_2;
  auto second_test_vector = sample_double_3;

  EXPECT_THROW(first_test_vector.projectionOnto(second_test_vector), VectorInconsistentDimensionException);
}

TEST(VectorProjectionOntoZeroMagnitudeTest, VectorProjectionOntoZeroMagnitude)
{
  auto first_test_vector  = sample_int_2;
  auto second_test_vector = sample_zero_2;

  EXPECT_THROW(first_test_vector.projectionOnto(second_test_vector), ZeroMagnitudeException);
}

class VectorMagnitudePositiveTest : public ::testing::TestWithParam<std::tuple<Vector<double>, double>>
{};

TEST_P(VectorMagnitudePositiveTest, VectorMagnitude)
{
  auto test_vector     = std::get<0>(GetParam());
  auto expected_result = std::get<1>(GetParam());

  auto result          = test_vector.magnitude();

  EXPECT_NEAR(result, expected_result, 1e-6);
}

INSTANTIATE_TEST_SUITE_P(
    VMagnitudePositive, VectorMagnitudePositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, expected_magnitude_2), std::make_tuple(sample_double_3, expected_magnitude_3),
        std::make_tuple(sample_double_4, expected_magnitude_4)
    )
);

class VectorUnitVectorPositiveTest : public ::testing::TestWithParam<std::tuple<Vector<double>, Vector<double>>>
{};

TEST_P(VectorUnitVectorPositiveTest, VectorUnitVector)
{
  auto test_vector     = std::get<0>(GetParam());
  auto expected_result = std::get<1>(GetParam());

  auto result          = test_vector.unitVector();

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VUnitVectorPositive, VectorUnitVectorPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_2, expected_unit_vector_2),
        std::make_tuple(sample_double_3, expected_unit_vector_3),
        std::make_tuple(sample_double_4, expected_unit_vector_4)
    )
);

class VectorVectorMatrixPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Matrix<double>, Vector<double>>>
{};

TEST_P(VectorVectorMatrixPositiveTest, VectorVectorMatrix)
{
  auto test_vector     = std::get<0>(GetParam());
  auto test_matrix     = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result          = test_vector * test_matrix;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VVectorMatrixPositive, VectorVectorMatrixPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_3, sample_double_3x3, expected_vector_matrix_3),
        std::make_tuple(sample_double_4, sample_double_4x4, expected_vector_matrix_4)
    )
);

TEST(VectorVectorMatrixNegativeTest, VectorVectorMatrixNegative)
{
  auto test_vector = sample_double_3;
  auto test_matrix = sample_double_4x4;

  EXPECT_THROW(test_vector * test_matrix, MatrixVectorSizeMismatchException);
}

class VectorMatrixVectorPositiveTest :
    public ::testing::TestWithParam<std::tuple<Vector<double>, Matrix<double>, Vector<double>>>
{};

TEST_P(VectorMatrixVectorPositiveTest, VectorMatrixVector)
{
  auto test_vector     = std::get<0>(GetParam());
  auto test_matrix     = std::get<1>(GetParam());
  auto expected_result = std::get<2>(GetParam());

  auto result          = test_matrix * test_vector;

  ASSERT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(
    VMatrixVectorPositive, VectorMatrixVectorPositiveTest,
    ::testing::Values(
        std::make_tuple(sample_double_3, sample_double_3x3, expected_matrix_vector_3),
        std::make_tuple(sample_double_4, sample_double_4x4, expected_matrix_vector_4)
    )
);

TEST(VectorMatrixVectorNegativeTest, VectorMatrixVectorNegative)
{
  auto test_vector = sample_double_3;
  auto test_matrix = sample_double_4x4;

  EXPECT_THROW(test_matrix * test_vector, MatrixVectorSizeMismatchException);
}
