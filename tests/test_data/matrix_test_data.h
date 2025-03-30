#include "matrix.h"

double sample_scalar = 3.2;

unsigned int sample_chosen_index = 0;
unsigned int sample_target_index = 1;

Matrix<int> sample_int_2x2 = Matrix<int>({
    {3, -7},
    {5,  7}
});
Matrix<double> sample_double_2x2 = Matrix<double>({
    {1.2,    6},
    {  9, -1.6}
});
Matrix<int> sample_int_3x3 = Matrix<int>({
    {  5,  1, -8},
    { -2, -1, 11},
    {-15,  6,  8}
});
Matrix<double> sample_double_3x3 = Matrix<double>({
    {  6,   -4, 2.11},
    { -7, 6.88,  7.1},
    {9.6,  -11,    8}
});
Matrix<int> sample_int_4x4 = Matrix<int>({
    {  5,  1, -8, 1},
    { -2, -1, 11, 7},
    { -2, -1, 11, 8},
    {-15,  6,  8, 4}
});
Matrix<double> sample_double_4x4 = Matrix<double>({
    {  6,   -4, 2.11, 6},
    { -7, 6.88,  7.1, 1},
    {  6,   -4,    8, 3},
    {9.6,  -11,    8, 9},
});
Matrix<int> sample_int_7x7 = Matrix<int>({
    {  5,  1, -8, 1,  9,   4,  -1},
    { -2, -1, 11, 7,  6,   3, -10},
    { -2, -1, 11, 8,  7,  -9, -11},
    {-15,  6,  8, 4,  8,  -3, -13},
    { -2, -1, 11, 8, -4,   6,   3},
    {  5,  1, -8, 1,  4, -21,   5},
    {  5,  1, -8, 1, -5,  -6,  -1},
});
Matrix<double> sample_double_7x7 = Matrix<double>({
    {    6,    -4,  2.11,     6,   9.3,     4,     5},
    { 8.88,   9.3,  10.1,   2.5,  -6.3, -15.4, -0.23},
    {    6,    -4,     8,     3,    -9,   5.5,  0.13},
    {  9.6,   -11,     8,     9, -0.99,  0.99,     2},
    { -6.2,     1,  3.33,  8.12,     6,  7.89,     3},
    {-8.21, -9.63,  7.77, -1.23,     6,   9.6,    -4},
    { 7.13,  6.54, -8.96, -4.32,  0.15,  9.88,  -4.2},
});
Matrix<int> sample_int_2x3 = Matrix<int>({
    {3, 4, 9},
    {5, 8, 8}
});
Matrix<double> sample_double_2x3 = Matrix<double>({
    {  6.3, 8, -4.6},
    {-7.77, 9,    1}
});
Matrix<int> sample_int_3x2 = Matrix<int>({
    { 6, -7},
    {-1,  3},
    { 9,  3}
});
Matrix<double> sample_double_3x2 = Matrix<double>({
    {-9.16, 6.33},
    {    1,  6.2},
    { -7.6,    8}
});

Matrix<double> sample_singular_double_3x3 = Matrix<double>({
    {2.1, 1, 2.1},
    {  1, 0,   1},
    {4.2, 1, 4.2}
});

Matrix<double> expected_addition_2x2 = Matrix<double>({
    {4.2,  -1},
    { 14, 5.4}
});
Matrix<double> expected_addition_3x3 = Matrix<double>({
    {  11,   -3, -5.89},
    {  -9, 5.88,  18.1},
    {-5.4,   -5,    16}
});
Matrix<double> expected_addition_4x4 = Matrix<double>({
    {  11,   -3, -5.89,  7},
    {  -9, 5.88,  18.1,  8},
    {   4,   -5,    19, 11},
    {-5.4,   -5,    16, 13}
});
Matrix<double> expected_addition_7x7 = Matrix<double>({
    {   11,    -3,  -5.89,     7,  18.3,     8,      4},
    { 6.88,   8.3,   21.1,   9.5,  -0.3, -12.4, -10.23},
    {    4,    -5,     19,    11,    -2,  -3.5, -10.87},
    { -5.4,    -5,     16,    13,  7.01, -2.01,    -11},
    { -8.2,     0,  14.33, 16.12,     2, 13.89,      6},
    {-3.21, -8.63,  -0.23, -0.23,    10, -11.4,      1},
    {12.13,  7.54, -16.96, -3.32, -4.85,  3.88,   -5.2}
});
Matrix<double> expected_addition_2x3 = Matrix<double>({
    {  9.3, 12, 4.4},
    {-2.77, 17,   9}
});
Matrix<double> expected_addition_3x2 = Matrix<double>({
    {-3.16, -0.67},
    {    0,   9.2},
    {  1.4,    11}
});

Matrix<double> expected_subtraction_2x2 = Matrix<double>({
    {-1.8,   13},
    {   4, -8.6}
});
Matrix<double> expected_subtraction_3x3 = Matrix<double>({
    {   1,   -5, 10.11},
    {  -5, 7.88,  -3.9},
    {24.6,  -17,     0}
});
Matrix<double> expected_subtraction_4x4 = Matrix<double>({
    {   1,   -5, 10.11,  5},
    {  -5, 7.88,  -3.9, -6},
    {   8,   -3,    -3, -5},
    {24.6,  -17,     0,  5}
});
Matrix<double> expected_subtraction_7x7 = Matrix<double>({
    {     1,     -5, 10.11,     5,   0.3,     0,     6},
    { 10.88,   10.3,  -0.9,  -4.5, -12.3, -18.4,  9.77},
    {     8,     -3,    -3,    -5,   -16,  14.5, 11.13},
    {  24.6,    -17,     0,     5, -8.99,  3.99,    15},
    {  -4.2,      2, -7.67,  0.12,    10,  1.89,     0},
    {-13.21, -10.63, 15.77, -2.23,     2,  30.6,    -9},
    {  2.13,   5.54, -0.96, -5.32,  5.15, 15.88,  -3.2}
});
Matrix<double> expected_subtraction_2x3 = Matrix<double>({
    {   3.3, 4, -13.6},
    {-12.77, 1,    -7}
});
Matrix<double> expected_subtraction_3x2 = Matrix<double>({
    {-15.16, 13.33},
    {     2,   3.2},
    { -16.6,     5}
});

Matrix<double> expected_multiplication_2x2 = Matrix<double>({
    {33.6,  33.6},
    {  19, -74.2}
});
Matrix<double> expected_multiplication_3x3 = Matrix<double>({
    {   6.35, 22.66, -75.12},
    {-155.26, 28.72, 188.48},
    {    -50,  68.6, -133.8}
});
Matrix<double> expected_multiplication_4x4 = Matrix<double>({
    {-56.22,  43.89, -20.79,  18.88},
    {-77.96, -14.98, 217.78, 101.96},
    {   -23,     20,     20,     54},
    {   -81,   66.6,  -37.8,   32.6}
});
Matrix<double> expected_multiplication_7x7 = Matrix<double>({
    { -29.82,  43.59,    9.51, 102.28,  46.57,  -83.19,  -24.31},
    { -97.45,  -4.85,   218.1,  98.75, 191.17,     252, -341.15},
    {  23.15,  34.63, -124.04, -12.37, 167.35, -239.28,  -92.63},
    { -64.07,  70.58,  -72.61,  27.67, 146.32, -132.33, -104.62},
    {-119.01,  43.08,  141.07, 118.81,  31.03, -223.82,  -91.54},
    {  -2.88, -14.13,   56.58,  35.22, -52.72, -269.57,  105.03},
    { 133.39, -10.84, -162.01, -29.17,  66.05,  -39.64,  136.24}
});
Matrix<double> expected_multiplication_2x3_3x2 = Matrix<double>({
    { -11.6, -33.9},
    {-46.62, 84.39}
});
Matrix<double> expected_multiplication_3x2_2x3 = Matrix<double>({
    {4.17,   14, -31.8},
    {  34, 53.6,  58.6},
    {17.2, 33.6,  -4.4}
});
Matrix<double> expected_multiplication_3x2_2x2 = Matrix<double>({
    {4.17, 108.43},
    {  34,   36.4},
    {17.2,  109.2}
});
Matrix<double> expected_multiplication_3x3_3x2 = Matrix<double>({
    {58.99, -47.67},
    {15.02,  90.94},
    {140.6,  -76.2}
});

Matrix<double> expected_scalar_multi_3x3 = Matrix<double>({
    {  16,  3.2, -25.6},
    {-6.4, -3.2,  35.2},
    { -48, 19.2,  25.6}
});
Matrix<double> expected_scalar_multi_2x3 = Matrix<double>({
    {9.6, 12.8, 28.8},
    { 16, 25.6, 25.6}
});

Matrix<double> expected_element_wise_2x2 = Matrix<double>({
    {3.6,   -42},
    { 45, -11.2}
});
Matrix<double> expected_element_wise_3x3 = Matrix<double>({
    {  30,    -4, -16.88},
    {  14, -6.88,   78.1},
    {-144,   -66,     64}
});
Matrix<double> expected_element_wise_4x4 = Matrix<double>({
    {  30,    -4, -16.88,  6},
    {  14, -6.88,   78.1,  7},
    { -12,     4,     88, 24},
    {-144,   -66,     64, 36}
});
Matrix<double> expected_element_wise_7x7 = Matrix<double>({
    {    30,    -4, -16.88,     6,  83.7,     16,    -5},
    {-17.76,  -9.3,  111.1,  17.5, -37.8,  -46.2,   2.3},
    {   -12,     4,     88,    24,   -63,  -49.5, -1.43},
    {  -144,   -66,     64,    36, -7.92,  -2.97,   -26},
    {  12.4,    -1,  36.63, 64.96,   -24,  47.34,     9},
    {-41.05, -9.63, -62.16, -1.23,    24, -201.6,   -20},
    { 35.65,  6.54,  71.68, -4.32, -0.75, -59.28,   4.2}
});
Matrix<double> expected_element_wise_2x3 = Matrix<double>({
    {  18.9, 32, -41.4},
    {-38.85, 72,     8}
});
Matrix<double> expected_element_wise_3x2 = Matrix<double>({
    {-54.96, -44.31},
    {    -1,   18.6},
    { -68.4,     24}
});

Matrix<double> expected_transposed_double_3x3 = Matrix({
    {   6,   -7, 9.6},
    {  -4, 6.88, -11},
    {2.11,  7.1,   8}
});
Matrix<double> expected_transposed_double_3x2 = Matrix({
    {-9.16,   1, -7.6},
    { 6.33, 6.2,    8}
});
Matrix<double> expected_transposed_double_2x3 = Matrix({
    { 6.3, -7.77},
    {   8,     9},
    {-4.6,     1}
});

double expected_determinant_double_2x2 = -55.92;
double expected_determinant_double_3x3 = 325.30872;
double expected_determinant_double_7x7 = 48387659.4170776556976;
double expected_determinant_singular_double_3x3 = 0;

Matrix<double> expected_inversed_2x2 = Matrix<double>({
    {0.0286123033,  0.1072961373},
    {0.1609442060, -0.0214592275}
});
Matrix<double> expected_inversed_3x3 = Matrix<double>({
    {0.4092727671, 0.0270204869, -0.1319263744},
    {0.3816682197, 0.0852851408, -0.1763555554},
    {0.0336664815, 0.0848424844,  0.0408227606},
});
Matrix<double> expected_inversed_4x4 = Matrix<double>({
    { 0.1827938041, -0.0674274736,  0.1890009173, -0.1773709003},
    { 0.2438745820,  0.0278214018,  0.1149913597, -0.2040047748},
    {-0.0807245858,  0.0371382612,  0.0954614784,  0.0178694242},
    { 0.1748440633,  0.0729147862, -0.1459106307,  0.0350836363}
});
Matrix<double> expected_inversed_7x7 = Matrix<double>({
    { 0.0425653670,  0.0158893179,  0.0079832641,  0.0111730474, -0.0352851395, -0.0009661208,  0.0310869704},
    { 0.0117618017,  0.0361712588,  0.0142194659, -0.0471518148,  0.0386910430, -0.0059762736,  0.0233363674},
    { 0.0512991958,  0.0410066200,  0.0545130072, -0.0568808708, -0.0104621119,  0.0468562743, -0.0186719112},
    {-0.1077104765,  0.0045598922, -0.0813216688,  0.1285436828,  0.1014566215, -0.0313027636,  0.0324988833},
    { 0.0502963588,  0.0227282011, -0.0545687781, -0.0028830245, -0.0069736425,  0.0368288530,  0.0155138043},
    { 0.0249537001, -0.0149473927,  0.0560779212, -0.0350805871,  0.0161375022,  0.0073151857,  0.0201159809},
    { 0.1524211247, -0.0432233787,  0.1330127595, -0.1479509054, -0.0439769594, -0.0601855205, -0.0947028447}
});

Matrix<double> expected_swap_row_3x3 = Matrix<double>({
    { -7, 6.88,  7.1},
    {  6,   -4, 2.11},
    {9.6,  -11,    8}
});
Matrix<double> expected_swap_row_2x3 = Matrix<double>({
    {-7.77, 9,    1},
    {  6.3, 8, -4.6}
});
Matrix<double> expected_swap_row_3x2 = Matrix<double>({
    {    1,  6.2},
    {-9.16, 6.33},
    { -7.6,    8}
});

Matrix<double> expected_swap_col_3x3 = Matrix<double>({
    {  -4,   6, 2.11},
    {6.88,  -7,  7.1},
    { -11, 9.6,    8}
});
Matrix<double> expected_swap_col_2x3 = Matrix<double>({
    {8,   6.3, -4.6},
    {9, -7.77,    1}
});
Matrix<double> expected_swap_col_3x2 = Matrix<double>({
    {6.33, -9.16},
    { 6.2,     1},
    {   8,  -7.6}
});
