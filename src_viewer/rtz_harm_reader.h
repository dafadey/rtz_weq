#pragma once

#include <string>
#include <vector>
#include <array>

struct energy_record {
  std::vector<float> ew; // spec energy
  std::vector<float> ew2;
  std::vector<float> ew3;
  std::vector<float> rw; // spec aperture
  std::vector<float> rw2;
  std::vector<float> rw3;
  std::vector<float> sw; // spec width
  std::vector<float> sw2;
  std::vector<float> sw3;
  std::vector<float> total; // total energy
};

struct src_record : public std::array<std::vector<float>, 3> {
  int nx, ny;
  float maxval() const;
};

float maxval(const std::vector<float>& vec);


float RefIndexAir_nm(float lam);

struct rtz_harm_reader {

  float k0;
  float x0;
  float t0;
  float z0;
  float dz;
  int save_interval{-1};
  
  enum src_rec_type {WWW, W2W2_W, PW, W3};
  
  std::string path;
  std::string refpath;// contains w3 field envelope
  int records_count;
  float max_src_www = .0f;
  float max_src_2w2w_w = .0f;
  float max_src_pw = .0f;
  float max_w3 = .0f;
  
  void init();
  
  energy_record e_rec;
	std::vector<int> iterations; // if save interval is dynamic we need intertion ids
  std::vector<float> www_src_integral;
  std::vector<float> www_src_instant;
  std::vector<float> w2w2w_src_integral;
  std::vector<float> w2w2w_src_instant;
  std::vector<float> pw_src_integral;
  std::vector<float> pw_src_instant;
  
  src_record get_src(int record_id, src_rec_type t, double fac = 1.) const;
};
