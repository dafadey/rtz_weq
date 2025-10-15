#include "vectors.h"

vec<3> cross_prod(const vec<3>& a, const vec<3>& b) {
  return vec<3>{ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

float cross_prod(const vec<2>& a, const vec<2>& b) {
  return a[0] * b[1] - a[1] * b[0];
}