#ifndef MATRIX_FUNCS_H
#define MATRIX_FUNCS_H

#include <tuple>

#include "matrix.h"

namespace supmath
{
namespace matrix_funcs
{
std::tuple<Matrix<double>, int> gaussianElimination(Matrix<double> matrix);

Matrix<double>                  gaussianEliminationInverse(Matrix<double> matrix);
}    // namespace matrix_funcs
}    // namespace supmath

#endif    // MATRIX_FUNCS_H
