#ifndef MATRIX_FUNC_H
#define MATRIX_FUNC_H

#include <tuple>

#include "matrix.h"

namespace supmath
{
std::tuple<Matrix<double>, int> gaussian_elimination(Matrix<double> matrix);

Matrix<double>                  gaussian_elimination_inverse(Matrix<double> matrix);
}    // namespace supmath

#endif    // MATRIX_FUNC_H
