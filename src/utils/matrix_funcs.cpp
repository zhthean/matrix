#include "matrix_funcs.h"
#include "exceptions/matrix_exceptions.h"
#include "matrix_enums.h"

using namespace supmath;

std::tuple<Matrix<double>, int> supmath::matrix_funcs::gaussianElimination(Matrix<double> matrix)
{
  int num_swaps = 0;

  for ( unsigned int pivot = 0; pivot < matrix.getRowSize(); pivot++ ) {
    unsigned int max_value_index = pivot;

    for ( unsigned int row = pivot + 1; row < matrix.getRowSize(); row++ ) {
      if ( std::abs(matrix[row][pivot]) > std::abs(matrix[max_value_index][pivot]) ) {
        max_value_index = row;
      }
    }

    if ( pivot != max_value_index ) {
      matrix.swap(pivot, max_value_index, MatrixOrder::Row);
      num_swaps++;
    }

    if ( matrix[pivot][pivot] == 0 ) {
      throw MatrixSingularException("The matrix is singular.");
    }

    for ( unsigned int row = pivot + 1; row < matrix.getRowSize(); row++ ) {
      double factor = static_cast<double>(matrix[row][pivot]) / matrix[pivot][pivot];

      for ( unsigned int col = pivot + 1; col < matrix.getColSize(); col++ ) {
        matrix[row][col] -= matrix[pivot][col] * factor;
      }

      matrix[row][pivot] = 0;
    }
  }

  return std::tuple<Matrix<double>, int>(matrix, num_swaps);
}

Matrix<double> supmath::matrix_funcs::gaussianEliminationInverse(Matrix<double> matrix)
{
  if ( matrix.getRowSize() != matrix.getColSize() ) {
    throw InvalidMatrixDimensionException("The matrix must be square.");
  }

  Matrix<double> augmented_matrix(matrix.getRowSize(), matrix.getColSize() * 2);

  for ( unsigned int i = 0; i < matrix.getRowSize(); i++ ) {
    for ( unsigned int j = 0; j < matrix.getColSize(); j++ ) {
      augmented_matrix[i][j] = matrix[i][j];
      if ( i == j ) {
        augmented_matrix[i][j + matrix.getColSize()] = 1;
      }
    }
  }

  for ( unsigned int pivot = 0; pivot < augmented_matrix.getRowSize(); pivot++ ) {
    unsigned int max_value_index = pivot;

    for ( unsigned int row = pivot + 1; row < augmented_matrix.getRowSize(); row++ ) {
      if ( std::abs(augmented_matrix[row][pivot]) > std::abs(augmented_matrix[max_value_index][pivot]) ) {
        max_value_index = row;
      }
    }

    if ( pivot != max_value_index ) {
      augmented_matrix.swap(pivot, max_value_index, MatrixOrder::Row);
    }

    if ( augmented_matrix[pivot][pivot] == 0 ) {
      throw MatrixSingularException("The matrix is singular.");
    }

    for ( unsigned row = 0; row < augmented_matrix.getRowSize(); row++ ) {
      if ( row == pivot ) {
        continue;
      }

      double factor = static_cast<double>(augmented_matrix[row][pivot]) / augmented_matrix[pivot][pivot];
      for ( unsigned int col = 0; col < augmented_matrix.getColSize(); col++ ) {
        augmented_matrix[row][col] -= augmented_matrix[pivot][col] * factor;
      }
    }
  }

  Matrix<double> inverse_matrix(matrix.getRowSize(), matrix.getColSize());

  for ( unsigned int row = 0; row < inverse_matrix.getRowSize(); row++ ) {
    double pivot_value = augmented_matrix[row][row];

    for ( unsigned int col = 0; col < inverse_matrix.getColSize(); col++ ) {
      inverse_matrix[row][col] = augmented_matrix[row][col + inverse_matrix.getColSize()] / pivot_value;
    }
  }

  return inverse_matrix;
}
