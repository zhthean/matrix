#ifndef MATRIX_EXCEPTIONS_H
#define MATRIX_EXCEPTIONS_H

#include <stdexcept>
#include <string>

class MatrixExceptions : public std::runtime_error {
public:
  explicit MatrixExceptions(const std::string &msg) : std::runtime_error(msg) {}
};

class InvalidMatrixDimensionException : public MatrixExceptions {
public:
  explicit InvalidMatrixDimensionException(const std::string &msg) : MatrixExceptions(msg) {}
};

class MatrixMismatchException : public MatrixExceptions {
public:
  explicit MatrixMismatchException(const std::string &msg) : MatrixExceptions(msg) {}
};

#endif // MATRIX_EXCEPTIONS_H
