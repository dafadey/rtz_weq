#pragma once

#include "dfs.h"

#include <string>
#include <vector>

struct postprocess
{
  bool drawPics = true; 
  bool save1D = true; 
  bool saveEnergy = true; 
  bool saveSource = true; 
  bool saveFieldEnvelope = false; 
  bool saveW3envelope = false; 
  const int CLEAN=1;
  const int SPEC_READY=2;
  const int FILTERS_READY=4;
  postprocess(int Nt, int Nx, double DT, double DX);
  ~postprocess();
  FL_DBL* plasma{};
  FL_DBL* field{};
  FL_DBL* spec{};
  FL_DBL* filtered_w{};
  FL_DBL* filtered_2w{};
  FL_DBL* filtered_3w{};
  std::vector<FL_DBL> tmp;
  FL_DBL* srcWWW{};
  FL_DBL* src2W2W_W{};
  FL_DBL* srcPW{};
	
  FL_DBL* spec_on_aperture{};
  int Nt, Nx;
  double dimt,dimx,dimz;
  double DT;
  double DX;
  double w0;
  double gap;
  double Kerr;
  std::string pictureFilename;
  std::string energyFilename;
  std::string signalFilename;
  std::string plasmaProfileFilename;
  std::string sourcePrefix;
  
  int status{CLEAN};
  
  void reset();
  void prepare_spectrum();
  void filter(FL_DBL** field, double freq, double gap);
  void w2w3w_filter();
  void plot();
  void saveField();
  void savePlasma();
  void saveEnergies();
  void saveSources();
};


