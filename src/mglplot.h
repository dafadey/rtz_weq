#pragma once

#include "dfs.h"

#include <string>

struct postprocess
{
  postprocess(int Nt, int Nx, double DT, double DX);
  ~postprocess();
  FL_DBL* field{};
  FL_DBL* spec{};
  FL_DBL* filtered_w{};
  FL_DBL* filtered_2w{};
  FL_DBL* filtered_3w{};
  FL_DBL* tmp{};
  int Nt, Nx;
  double DT;
  double DX;
  double w;
  std::string pictureFilename;
  std::string energyFilename;
  std::string signalFilename;
  
  void prepare_spectrum();
  void filter();  
  void w2w3w_plot();
};


