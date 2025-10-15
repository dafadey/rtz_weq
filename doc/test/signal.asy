import dims;

string number(int i)
{
  string num;
  string _num = string(i);
  for(int j=0; j<7-length(_num);++j)
    num += '0';
  num += _num;
  return num;
}

struct descr {
  real dr;
  real dz;
  int save_interval;
  real dt;
  real dw;
  real z_shift;
  real k0;
  real airDensity=0;
  void print() {
    if(abs(this.k0 - k0)/k0 > 1.e-5)
      write("ACHTUNG!!!");
    write("dr = "+string(dr)+" : "+string(dr*x0)+" cm");
    write("dt = "+string(dt)+" : "+string(dt*t0)+" s");
    write("save_interval = "+string(save_interval)+" : "+string(save_interval*dz*z0)+" cm");
    write("ne = "+string(airDensity)+" : "+string(airDensity*n0)+" cm^-3");
  }
};

descr read_descr(string filename) {
  descr d;
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  int Nx;
  int Nt;
  real Lx;
  real T;
  for(string s : raw) {
    string[] items = split(s, '=');
    string item = items[0];
    real val = (real) items[1];
    if(item == "k0") d.k0 = val;
    if(item == "z_shift") d.z_shift = -val;
    if(item == "dz") d.dz = -val;
    if(item == "save_interval") d.save_interval = (int) val;
    if(item == "dw") d.dw = val;
    if(item == "Lx") Lx = val;
    if(item == "Nx") Nx = (int) val;
    if(item == "T") T = val;
    if(item == "Nt") Nt = (int) val;
    if(item == "airDensity") d.airDensity = val;
  }
  d.dr = Lx/Nx;
  d.dt = T/Nt;
  return d;
}

real[] read_signal(string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  string[] sinput = split(raw[4], '\t');
  int n = (int) raw[2];
  real[] res=new real[n];
  for(int i=0;i<n;++i)
    res[i] = (real) raw[i+4];
  return res;
}

pair getPeriodAndPhase(real[] signal, int hint) {
  int n = signal.length;
  //write("n=",n);
  real pt1=0;
  real pt2=0;
  hint = (hint+n)%n;
  for(int i = hint;i<n-1;++i) {
    if(signal[i]< .0 && signal[i+1] > .0 && pt1==0) {
      pt1 = i + abs(signal[i])/(abs(signal[i])+abs(signal[i+1]));
      continue;
    }
    if(signal[i] < .0 && signal[i+1] > .0 && pt2==0) {
      pt2 = i + abs(signal[i])/(abs(signal[i])+abs(signal[i+1]));
      break;
    }
  }
  //write("phase:", (pt1, pt2));
  return (pt1, pt2);
}
