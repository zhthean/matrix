#ifndef MATRIX_EXCEPTIONS_H
#define MATRIX_EXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace supmath
{
class MatrixExceptions : public std::runtime_error
{
public:
	explicit MatrixExceptions(const std::string &msg) : std::runtime_error(msg)
	{
	}
};

class InvalidMatrixDimensionException : public MatrixExceptions
{
public:
	explicit InvalidMatrixDimensionException(const std::string &msg) : MatrixExceptions(msg)
	{
	}
};

class MatrixMismatchException : public MatrixExceptions
{
public:
	explicit MatrixMismatchException(const std::string &msg) : MatrixExceptions(msg)
	{
	}
};

class MatrixNonSquareException : public MatrixExceptions
{
public:
	explicit MatrixNonSquareException(const std::string &msg) : MatrixExceptions(msg)
	{
	}
};

class MatrixSingularException : public MatrixExceptions
{
public:
	explicit MatrixSingularException(const std::string &msg) : MatrixExceptions(msg)
	{
	}
};

class MatrixNonInvertibleException : public MatrixExceptions
{
public:
	explicit MatrixNonInvertibleException(const std::string &msg) : MatrixExceptions(msg)
	{
	}
};
}    // namespace supmath

#endif    // MATRIX_EXCEPTIONS_H
