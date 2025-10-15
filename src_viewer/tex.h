#pragma once
#include <array>

struct tex_manager {
  int nx=0;
  int ny=0;
    
  void init(int nx, int ny);
  
  std::array<float, 4> (*value_to_color)(float);
  std::array<float, 4> (*value_to_color_cplx)(float, float);
  
  float* tex_buff = nullptr;
  unsigned int tex = 0;
  
  void update_on_host(int Dx, int Dy, float* data);
  void update_on_host(int Dx, int Dy, float* data_re, float* data_im);
  void update_on_device();
};
