struct descr {
  real dr;
  real dz=0.5;
  int save_interval=333;
  real dt;
  real dw;
  real T;
  real Lx;
  real ampW;
  real r_0;
  real z_shift;
  real k0;
  real airDensity=0;
  int Nt;
  int Nx;
};

descr read_descr(string filename) {
  descr d;
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  for(string s : raw) {
    string[] items = split(s, '=');
    string item = items[0];
    real val = (real) items[1];
    if(item == "k0") d.k0 = val;
    if(item == "z_shift") d.z_shift = -val;
    if(item == "dz") d.dz = -val;
    if(item == "save_interval") d.save_interval = (int) val;
    if(item == "dw") d.dw = val;
    if(item == "r_0") d.r_0 = val;
    if(item == "ampW") d.ampW = val;
    if(item == "Lx") d.Lx = val;
    if(item == "Nx") d.Nx = (int) val;
    if(item == "T") d.T = val;
    if(item == "Nt") d.Nt = (int) val;
    if(item == "airDensity") d.airDensity = val;
  }
  d.dr = d.Lx/d.Nx;
  d.dt = d.T/d.Nt;
  return d;
}

string number(int i)
{
  string num;
  string _num = string(i);
  for(int j=0; j<7-length(_num);++j)
    num += '0';
  num += _num;
  return num;
}
