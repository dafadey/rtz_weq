#include "dfs.h"

extern "C" 
void cuda_load(FL_DBL*, int, int, FL_DBL, FL_DBL, FL_DBL, FL_DBL=.0);
extern "C" 
void get_spec(FL_DBL*, int, int);
extern "C" 
void get_field(FL_DBL*, int, int);
extern "C" 
void get_density(FL_DBL*, int, int);
extern "C" 
void get_THz_source(FL_DBL*, int, int);
extern "C" 
FL_DBL step(int Nt, int Nx, FL_DBL z, FL_DBL T, bool save, FL_DBL airDensity, FL_DBL Kerr, FL_DBL tunExp, FL_DBL tunFac);

extern "C"
FL_DBL LaplassianA(int i, FL_DBL dr);

extern "C"
FL_DBL LaplassianB(int i, FL_DBL dr);

extern "C"
FL_DBL LaplassianC(int i, FL_DBL dr);

extern "C"
FL_DBL LaplassianMetric(int i, FL_DBL dr);
