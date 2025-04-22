#ifndef MATRIX_FUNC_H
#define MATRIX_FUNC_H

#include <tuple>

#include "matrix.h"

std::tuple<supp_math::Matrix<double>, int> gaussian_elimination(supp_math::Matrix<double> matrix);

supp_math::Matrix<double> gaussian_elimination_inverse(supp_math::Matrix<double> matrix);

#endif // MATRIX_FUNC_H
