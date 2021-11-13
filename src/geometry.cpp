#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

template<class T, size_t H_, size_t W_>
struct Matrix {
  static constexpr size_t H = H_, W = W_;

  T values[H*W];

  Matrix() {}
  Matrix(const T &value) {
    std::fill(values, values+H*W, value);
  }
  Matrix(const T (&src)[H][W]) {
    for (size_t i = 0; i < H; ++i)
      for (size_t j = 0; j < W; ++j)
        values[i*W + j] = src[i][j];
  }

  T &operator()(size_t i, size_t j) {
    return values[i*W + j];
  }
  const T &operator()(size_t i, size_t j) const {
    return values[i*W + j];
  }

  Matrix &operator*=(const T &s) {
    for (size_t i = 0; i < H*H; ++i)
      values[i] *= s;
    return *this;
  }
  Matrix operator*(const T &s) const {
    Matrix res = *this;
    return res *= s;
  }

  template<size_t OW> 
  Matrix &operator*=(const Matrix<T, W, OW> &rhs) {
    return *this = *this * rhs;
  }
};

template<class T, size_t W, size_t M, size_t H>
Matrix<T, H, W> operator*(const Matrix<T, H, M> &lhs, const Matrix<T, M, W> &rhs) {
  Matrix<T, H, W> result(0);
  for (size_t i = 0; i < H; ++i)
    for (size_t j = 0; j < M; ++j)
      for (size_t k = 0; k < W; ++k)
        result(i, k) += lhs(i, j) * rhs(j, k);
  return result;
}

template<class T, size_t N>
auto identity(const T &additive_identity = 0, const T &multiplicative_identity = 1) {
  Matrix<T, N, N> res(additive_identity);
  for (size_t i = 0; i < N; ++i)
    res(i, i) = multiplicative_identity;
  return res;
}

using Mat3 = Matrix<double, 3, 3>;

Mat3 rot_x(double theta) {
  double st = std::sin(theta);
  double ct = std::cos(theta);
  return Mat3({
    {  1,  0,  0},
    {  0, ct,-st},
    {  0, st, ct}
  });
}

Mat3 rot_y(double theta) {
  double st = std::sin(theta);
  double ct = std::cos(theta);
  return Mat3({
    { ct,  0, st},
    {  0,  1,  0},
    {-st,  0, ct}
  });
}

Mat3 rot_z(double theta) {
  double st = std::sin(theta);
  double ct = std::cos(theta);
  return Mat3({
    { ct,-st,  0},
    { st, ct,  0},
    {  0,  0,  1}
  });
}

struct V3 {
  double x, y, z;

  V3() = default;
  V3(double x_, double y_, double z_)
    : x(x_), y(y_), z(z_) {}
  V3(const Matrix<double, 3, 1> &m)
    : x(m(0, 0)), y(m(1, 0)), z(m(2, 0)) {}
  V3(const Matrix<double, 1, 3> &m)
    : x(m(0, 0)), y(m(0, 1)), z(m(0, 2)) {}

  // [3][1]
  auto as_col() const {
    return Matrix<double, 3, 1>({{x}, {y}, {z}});
  }
  // [1][3]
  auto as_row() const {
    return Matrix<double, 1, 3>({{x, y, z}});
  }

  template<size_t OW>
  friend Matrix<double, 3, OW> operator*(const V3 &lhs, const Matrix<double, 1, OW> &rhs) {
    return lhs.as_col() * rhs;
  }

  template<size_t OW>
  friend Matrix<double, 1, OW> operator*(const V3 &lhs, const Matrix<double, 3, OW> &rhs) {
    return lhs.as_row() * rhs;
  }

  template<size_t OH>
  friend Matrix<double, OH, 1> operator*(const Matrix<double, OH, 3> &lhs, const V3 &rhs) {
    return lhs * rhs.as_col();
  }
  template<size_t OH>
  friend Matrix<double, OH, 3> operator*(const Matrix<double, OH, 1> &lhs, const V3 &rhs) {
    return lhs * rhs.as_row();
  }

  friend V3 operator+(const V3 &lhs, const V3 &rhs) {
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
  }
};

template<size_t N>
auto operator+(const std::array<V3, N> &lhs, const V3 &rhs) {
  std::array<V3, N> res;
  for (size_t i = 0; i < N; ++i)
    res[i] = lhs[i] + rhs;
  return res;
}

template<size_t N>
auto operator*(const std::array<V3, N> &lhs, const Mat3 &rhs) {
  std::array<V3, N> res;
  for (size_t i = 0; i < N; ++i)
    res[i] = lhs[i] * rhs;
  return res;
}
