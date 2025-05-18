#ifndef VECTOR_TEST_DATA_H
#define VECTOR_TEST_DATA_H

#include "matrix.h"
#include "vector.h"

inline double                  sample_scalar            = 3.0;

inline supmath::Vector<int>    sample_int_2             = supmath::Vector<int>({1, 2});
inline supmath::Vector<double> sample_double_2          = supmath::Vector<double>({3.5, -1.5});
inline supmath::Vector<int>    sample_int_3             = supmath::Vector<int>({-1, 0, 2});
inline supmath::Vector<double> sample_double_3          = supmath::Vector<double>({2.0, 4.0, -3.0});
inline supmath::Vector<int>    sample_int_4             = supmath::Vector<int>({-1, 0, -6, 6});
inline supmath::Vector<double> sample_double_4          = supmath::Vector<double>({2.0, -4.0, -2.0, 9.0});
inline supmath::Vector<double> sample_zero_2            = supmath::Vector<double>({0.0, 0.0});

supmath::Matrix<double>        sample_double_3x3        = supmath::Matrix<double>({
  {  6,   -4, 2.11},
  { -7, 6.88,  7.1},
  {9.6,  -11,    8}
});
supmath::Matrix<double>        sample_double_4x4        = supmath::Matrix<double>({
  {  6,   -4, 2.11, 6},
  { -7, 6.88,  7.1, 1},
  {  6,   -4,    8, 3},
  {9.6,  -11,    8, 9},
});

inline supmath::Vector<double> expected_addition_2      = supmath::Vector<double>({4.5, 0.5});
inline supmath::Vector<double> expected_addition_3      = supmath::Vector<double>({1.0, 4.0, -1.0});
inline supmath::Vector<double> expected_addition_4      = supmath::Vector<double>({1.0, -4.0, -8.0, 15.0});

inline supmath::Vector<double> expected_subtraction_2   = supmath::Vector<double>({2.5, -3.5});
inline supmath::Vector<double> expected_subtraction_3   = supmath::Vector<double>({3.0, 4.0, -5.0});
inline supmath::Vector<double> expected_subtraction_4   = supmath::Vector<double>({3.0, -4.0, 4.0, 3.0});

inline supmath::Vector<double> expected_scalar_multi_2  = supmath::Vector<double>({3, 6});
inline supmath::Vector<double> expected_scalar_multi_3  = supmath::Vector<double>({-3, 0, 6});
inline supmath::Vector<double> expected_scalar_multi_4  = supmath::Vector<double>({-3, 0, -18.0, 18.0});

inline supmath::Vector<double> expected_cross_product_3 = supmath::Vector<double>({-8, 1, -4});

inline supmath::Vector<double> expected_proj            = supmath::Vector<double>({3.5, -1.5});
inline supmath::Vector<double> expected_unit_2          = supmath::Vector<double>({3.5 / 3.8078866, -1.5 / 3.8078866});

inline double                  expected_dot_product_2   = 0.5;
inline double                  expected_dot_product_3   = -8;
inline double                  expected_dot_product_4   = 64;

inline double                  expected_angle_between_2 = 1.512040504;
inline double                  expected_angle_between_3 = 2.297438667;
inline double                  expected_angle_between_4 = 0.750993983;

inline double                  expected_distance_between_2 = 4.301162634;
inline double                  expected_distance_between_3 = 7.071067812;
inline double                  expected_distance_between_4 = 7.071067812;

inline supmath::Vector<double> expected_projection_onto_2  = supmath::Vector<double>({0.1, 0.2});
inline supmath::Vector<double> expected_projection_onto_3  = supmath::Vector<double>({1.6, 0, -3.2});
inline supmath::Vector<double> expected_projection_onto_4 =
    supmath::Vector<double>({-0.876712329, 0, -5.260273973, 5.260273973});

inline double                  expected_magnitude_2   = 3.807886553;
inline double                  expected_magnitude_3   = 5.385164807;
inline double                  expected_magnitude_4   = 10.246950766;

inline supmath::Vector<double> expected_unit_vector_2 = supmath::Vector<double>({0.919145030, -0.393919299});
inline supmath::Vector<double> expected_unit_vector_3 =
    supmath::Vector<double>({0.371390676, 0.742781353, -0.557086015});
inline supmath::Vector<double> expected_unit_vector_4 =
    supmath::Vector<double>({0.195180015, -0.390360029, -0.195180015, 0.878310066});

inline supmath::Vector<double> expected_vector_matrix_3 = supmath::Vector<double>({-44.8, 52.52, 8.62});
inline supmath::Vector<double> expected_vector_matrix_4 = supmath::Vector<double>({114.4, -126.52, 31.82, 83});

inline supmath::Vector<double> expected_matrix_vector_3 = supmath::Vector<double>({-10.33, -7.78, -48.8});
inline supmath::Vector<double> expected_matrix_vector_4 = supmath::Vector<double>({77.78, -46.72, 39, 128.2});

#endif    // VECTOR_TEST_DATA_H
