#ifndef VECTOR_EXCEPTIONS_H
#define VECTOR_EXCEPTIONS_H

#include <stdexcept>
#include <string>

class VectorException : public std::runtime_error {
public:
  explicit VectorException(const std::string &msg) : std::runtime_error(msg) {}
};

class VectorInconsistentDimensionException : public VectorException {
public:
  explicit VectorInconsistentDimensionException(const std::string &msg) : VectorException(msg) {}
};

class UnsupportVectorSizeException : public VectorException {
public:
  explicit UnsupportVectorSizeException(const std::string &msg) : VectorException(msg) {}
};

class ZeroMagnitudeException : public VectorException {
public:
  explicit ZeroMagnitudeException(const std::string &msg) : VectorException(msg) {}
};

#endif // VECTOR_EXCEPTIONS_H
