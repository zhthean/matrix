#include "matrix_func.h"
#include "exceptions/matrix_exceptions.h"

std::tuple<supp_math::Matrix<double>, int> gaussian_elimination(supp_math::Matrix<double> matrix) {
  int num_swaps = 0;

  for (unsigned int pivot = 0; pivot < matrix.get_row_size(); pivot++) {
    unsigned int max_value_index = pivot;

    for (unsigned int row = pivot + 1; row < matrix.get_row_size(); row++) {
      if (std::abs(matrix[row][pivot]) > std::abs(matrix[max_value_index][pivot])) {
        max_value_index = row;
      }
    }

    if (pivot != max_value_index) {
      matrix.swap(pivot, max_value_index, 'r');
      num_swaps++;
    }

    if (matrix[pivot][pivot] == 0) {
      throw supp_math::MatrixSingularException("The matrix is singular.");
    }

    for (unsigned int row = pivot + 1; row < matrix.get_row_size(); row++) {
      double factor = static_cast<double>(matrix[row][pivot]) / matrix[pivot][pivot];

      for (unsigned int col = pivot + 1; col < matrix.get_col_size(); col++) {
        matrix[row][col] -= matrix[pivot][col] * factor;
      }

      matrix[row][pivot] = 0;
    }
  }

  return std::tuple<supp_math::Matrix<double>, int>(matrix, num_swaps);
}

supp_math::Matrix<double> gaussian_elimination_inverse(supp_math::Matrix<double> matrix) {
  if (matrix.get_row_size() != matrix.get_col_size()) {
    throw supp_math::InvalidMatrixDimensionException("The matrix must be square.");
  }

  supp_math::Matrix<double> augmented_matrix(matrix.get_row_size(), matrix.get_col_size() * 2);

  for (unsigned int i = 0; i < matrix.get_row_size(); i++) {
    for (unsigned int j = 0; j < matrix.get_col_size(); j++) {
      augmented_matrix[i][j] = matrix[i][j];
      if (i == j) {
        augmented_matrix[i][j + matrix.get_col_size()] = 1;
      }
    }
  }

  for (unsigned int pivot = 0; pivot < augmented_matrix.get_row_size(); pivot++) {
    unsigned int max_value_index = pivot;

    for (unsigned int row = pivot + 1; row < augmented_matrix.get_row_size(); row++) {
      if (std::abs(augmented_matrix[row][pivot]) > std::abs(augmented_matrix[max_value_index][pivot])) {
        max_value_index = row;
      }
    }

    if (pivot != max_value_index) {
      augmented_matrix.swap(pivot, max_value_index, 'r');
    }

    if (augmented_matrix[pivot][pivot] == 0) {
      throw supp_math::MatrixSingularException("The matrix is singular.");
    }

    for (unsigned row = 0; row < augmented_matrix.get_row_size(); row++) {
      if (row == pivot) {
        continue;
      }

      double factor = static_cast<double>(augmented_matrix[row][pivot]) / augmented_matrix[pivot][pivot];
      for (unsigned int col = 0; col < augmented_matrix.get_col_size(); col++) {
        augmented_matrix[row][col] -= augmented_matrix[pivot][col] * factor;
      }
    }
  }

  supp_math::Matrix<double> inverse_matrix(matrix.get_row_size(), matrix.get_col_size());

  for (unsigned int row = 0; row < inverse_matrix.get_row_size(); row++) {
    double pivot_value = augmented_matrix[row][row];

    for (unsigned int col = 0; col < inverse_matrix.get_col_size(); col++) {
      inverse_matrix[row][col] = augmented_matrix[row][col + inverse_matrix.get_col_size()] / pivot_value;
    }
  }

  return inverse_matrix;
}
