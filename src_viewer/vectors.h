#pragma once
#include <array>
#include <cmath>
#include <iostream>

template <size_t dim>
using vec = std::array<float, dim>;

template <size_t dim>
vec<dim> operator+(const vec<dim>& a, const vec<dim>& b) {
  vec<dim> res;
  for (int i = 0; i < dim; i++)
    res[i] = a[i] + b[i];
  return res;
}

template <size_t dim>
vec<dim> operator-(const vec<dim>& a, const vec<dim>& b) {
  vec<dim> res;
  for (int i = 0; i < dim; i++)
    res[i] = a[i] - b[i];
  return res;
}

template <size_t dim>
float operator*(const vec<dim>& a, const vec<dim>& b) {
  float res = .0f;
  for (int i = 0; i < dim; i++)
    res += a[i] * b[i];
  return res;
}

template <size_t dim>
vec<dim> operator*(const float a, const vec<dim>& b) {
  vec<dim> res;
  for (int i = 0; i < dim; i++)
    res[i] = a * b[i];
  return res;
}

template <size_t dim>
vec<dim> operator*(const vec<dim>& b, const float a) {
  vec<dim> res;
  for (int i = 0; i < dim; i++)
    res[i] = a * b[i];
  return res;
}

template <size_t dim>
vec<dim> operator/(const vec<dim>& b, const float a) {
  vec<dim> res;
  for (int i = 0; i < dim; i++)
    res[i] = b[i] / a;
  return res;
}

vec<3> cross_prod(const vec<3>& a, const vec<3>& b);

float cross_prod(const vec<2>& a, const vec<2>& b);

template <size_t dim>
void normalize(vec<dim>& a) {
  a = a / std::sqrt(a * a);
}

template <size_t dim>
std::ostream& operator<<(std::ostream& o, const vec<dim>& v) {
  o << '(';
  int i=0;
  for(;i<dim-1;i++)
    o << v[i] << ", ";
  o << v[i] << ')';
  return o;
}

typedef vec<2> vec2f;
typedef vec<3> vec3f;
typedef vec<4> vec4f;