import palette;
import graph;
import fft;

defaultpen(linewidth(0.1pt));

int fileCount=49;

real beamDia_mm=5.;
real focalLength_mm = 85.;

real k0=80553.65779;
real x0 = 1. / (2.^.5 * k0);
real t0 = 1. / k0;
real z0 = 2/k0;
real dz = .5*z0;
real d0 = .1 * beamDia_mm / x0;// in cm
real z_foc = .1 * focalLength_mm / z0;
real dw = (d0^2/2. - (d0^4/4.-1024.*z_foc^2)^.5)^.5;
real z_start = -.8/z0*160./dw*z0;
int save_interval=2000;

int set=7;

void read_point(real[][] data, string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  string[] sinput = split(raw[4], '\t');
  real[] input;
  for(string v : sinput)
    input.push((real) v);
  data.push(input);
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

real[][] filter(real[] data) {  
  int n=data.length;
  //write("n=", n);
  qqFFT_freal(n, data);
  real[] data_w = new real[n];
  real[] data_2w = new real[n];
  real[] data_3w = new real[n];
  for(int i=0;i<n;++i) {
    data_w[i] = data[i];
    data_2w[i] = data[i];
    data_3w[i] = data[i];
  }
  int n_2 = floor(n/2);
  for(int i=0;i<n_2;++i) {
    real fw = exp(-((i-40)/30)^4);
    real f2w = exp(-((i-80)/30)^4);
    real f3w = exp(-((i-120)/30)^4);
    data_w[i*2] *= fw;
    data_w[i*2+1] *= fw;
    data_2w[i*2] *= f2w;
    data_2w[i*2+1] *= f2w;
    data_3w[i*2] *= f3w;
    data_3w[i*2+1] *= f3w;
  }
  qqFFT_freal_1(n, data_w);
  qqFFT_freal_1(n, data_2w);
  qqFFT_freal_1(n, data_3w);
  real[][] res;
  res.push(data_w);
  res.push(data_2w);
  res.push(data_3w);
  return res;
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

void plot_energy(picture pic, string data_path, pen p=defaultpen)
{
  real[][] data;

  for(int i=0;i<=fileCount;++i)
    read_point(data, data_path+"/"+"field_energy"+number(i)+".dat");
  
  write("read enrgies from" + data_path);
  
  guide g;
  for(int i=0;i<data.length;++i)
  {
    real z=z_start+i*save_interval*dz;
    //g=g--(z,log10(data[i][2]/data[i][0]^3+1e-13));
    g=g--(z,log10(data[i][2]+1e-13));
  }
  draw(pic, g, p);
}

void plot_radius(picture pic, string data_path, pen p=defaultpen)
{
  real[][] data;

  for(int i=0;i<=fileCount;++i)
    read_point(data, data_path+"/"+"field_energy"+number(i)+".dat");

  write("read radii from" + data_path);

  if(data[0].length>3)
  {
    guide g,g2w;
    for(int i=0;i<data.length;++i)
    {
      real z=z_start+i*save_interval*dz;
      g=g--(z,data[i][3]*x0*10.);
      g2w=g2w--(z,data[i][4]*x0*10.);
    }
    draw(pic, g, p);
    draw(pic, g2w, p+dotted);
  }
}

void plot_spec_width(picture pic, string data_path, pen p=defaultpen)
{
  real[][] data;

  for(int i=0;i<=fileCount;++i)
    read_point(data, data_path+"/"+"field_energy"+number(i)+".dat");

  write("read spec widths from" + data_path);

  if(data[0].length>6)
  {
    guide gw;
    guide g3w;
    for(int i=0;i<data.length;++i)
    {
      real z=z_start+i*save_interval*dz;
      gw = gw -- (z, data[i][6]);
      g3w = g3w -- (z, data[i][8]);
    }
    draw(pic, gw, p);
    draw(pic, g3w, p + linetype(new real[] {3,11}));
  }
}


void plot_total_energy(picture pic, string data_path, pen p=defaultpen)
{
  real[][] data;

  for(int i=0;i<=fileCount;++i)
    read_point(data, data_path+"/"+"field_energy"+number(i)+".dat");

  write("read total energy from" + data_path);

  if(data[0].length>9)
  {
    guide g;
    for(int i=0;i<data.length;++i)
    {
      real z = z_start + i * save_interval * dz;
      g = g -- (z, data[i][9]);
    }
    draw(pic, g, p);
  }
}


void plot_phase(picture pic, string data_path, pen p=defaultpen) {

  real[][] phase;

  real[] pp0 = new real[3];

  pair[] pp = new pair[3];
  real[] cm;
  for(int i=0;i<=fileCount;++i) {

    real addition=0;

    real cntr=0;
    real norm=0;

    real[] signal = read_signal(data_path+"/"+"field_axis"+number(i)+".dat");
    real[][] sigs = filter(signal);
    real[] item = new real[3];
    for(int j=0;j<sigs[0].length;++j)
    {
      cntr += j*sigs[0][j]^2;
      norm += sigs[0][j]^2;
    }
    for(int j=0;j<3;++j)
    {
      pp[j] = getPeriodAndPhase(sigs[j], i==0 ? floor(signal.length/2) : floor(pp[j].x)-3);
      if(i==0)
        pp0[j] = pp[j].x;
      if(abs(pp[j].y-pp[j].x)!=0)
        item[j] = (pp[j].x-pp0[j]-addition)/abs(pp[j].y-pp[j].x)*2*pi;
      else
        item[j] = 0;
    }
    phase.push(item);
    cm.push(cntr/norm);
  }

  write("done phases from" + data_path);

  guide gw, g2w, g3w, a;
  for(int i=0;i<phase.length;++i) {
    real z=z_start+i*save_interval*dz;
    //gw=gw--(z, cm[i]);
    gw=gw--(z, phase[i][0]);
    g2w=g2w--(z, phase[i][1]);
    g3w=g3w--(z, phase[i][2]);
    //a = a--(z,0.00031*(z_start-z)*2./z0*(0.5+0.5*tanh(5.)));
  }
  draw(pic, gw, p);
  //draw(pic,a,dotted);
  //draw(pic, g2w, p+dashed);
  draw(pic, g3w, p+linetype(new real[] {3,11}));
}

picture pic, picp, picr, picw, pic_total;

if(set==7)
{
  //0.00011 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0011 0.0011744
  pen[] colors = Rainbow(16);
  for(int i=0; i<colors.length; ++i)
    colors[i] *= .7;
  colors.push(black);
  z_start = -13402./1.7*z0;//-.4/z0*160./dw*z0;
  fileCount=92;
  save_interval = 333;
  void plotall(string name, pen p)
  {
    plot_total_energy(pic_total, name, p);
    plot_energy(pic, name, p);
    plot_radius(picr, name, p);
    plot_phase(picp, name, p);
    plot_spec_width(picw, name, p);
  }

/* 
	string[] phases = {"0.6","0.7","0.8","0.9","1.0","1.1","1.2","1.3"};
  for(int i=0;i<8;++i)
   plotall("harm_ampW0=0.0005_ampW2=0.0003_phase2w="+phases[i]+"_/pics/res", colors[i*2]);
*/

	string[] amps = {"0.0001","0.0002","0.0003","0.0004","0.0005","0.0006","0.00065","0.0007"};
  for(int i=0;i<7;++i)
    plotall("harm_ampW0="+amps[i]+"_ampW2=0.0003_phase2w=1.2_/pics/res", colors[i*2]);

/*
  for(int i=0;i<5;++i)
    plotall("harm_ampW0=0.000"+string(i+1)+"_ampW2=0.0005_phase2w=1.0_/pics/res", colors[i*3]);
*/
  
  /*
  plotall("w_harm_ampW0=0.00011744_ratio2w=0_/pics/res", colors[0]);
  plotall("w_harm_ampW0=0.0002_ratio2w=0_/pics/res", colors[1]);
  plotall("w_harm_ampW0=0.0003_ratio2w=0_/pics/res", colors[2]);
  plotall("w_harm_ampW0=0.0004_ratio2w=0_/pics/res", colors[3]);
  plotall("w_harm_ampW0=0.0005_ratio2w=0_/pics/res", colors[4]);
  plotall("w_harm_ampW0=0.0006_ratio2w=0_/pics/res", colors[5]);
  plotall("w_harm_ampW0=0.0007_ratio2w=0_/pics/res", colors[6]);
  plotall("w_harm_ampW0=0.0008_ratio2w=0_/pics/res", colors[7]);
  plotall("w_harm_ampW0=0.0009_ratio2w=0_/pics/res", colors[8]);
  plotall("w_harm_ampW0=0.0010_ratio2w=0_/pics/res", colors[9]);
  plotall("w_harm_ampW0=0.0011_ratio2w=0_/pics/res", colors[10]);
  //plotall("w_harm_ampW0=0.0011_ratio2w=0_dz=-.25_save_interval=1000_draw_interval=500_/pics/res", colors[11]);
  plotall("w_harm_ampW0=0.0011744_ratio2w=0_/pics/res", colors[11]);
  plotall("test_harm_ampW0=0.0011744_/pics/res", colors[11]);
  */
  
  /*
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_/pics/res", colors[11]);
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_phase2w=0.3_/pics/res", colors[11]);
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_phase2w=0.5_/pics/res", colors[11]);
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_phase2w=0.7_/pics/res", colors[11]);
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_phase2w=0.8_/pics/res", colors[11]);
  plotall("w2w_harm_ampW0=0.0010_ampW2=0.00025_phase2w=1.1_/pics/res", colors[11]);
  */
}

if(set==6)
{
  z_start = -.4/z0*160./dw*z0;
  fileCount=99;
  save_interval = 500;
  
  plot_energy(pic, "test_with_Kerr/pics/res",.7*red);
  plot_radius(picr, "test_with_Kerr/pics/res",.7*red);
  plot_phase(picp, "test_with_Kerr/pics/res",.7*red);

  plot_energy(pic, "test_with_Kerr_dimmed/pics/res", .7*blue);
  plot_radius(picr, "test_with_Kerr_dimmed/pics/res", .7*blue);
  plot_phase(picp, "test_with_Kerr_dimmed/pics/res", .7*blue);

  plot_energy(pic, "test_with_Kerr_dimmed_more/pics/res", .5*green);
  plot_radius(picr, "test_with_Kerr_dimmed_more/pics/res", .5*green);
  plot_phase(picp, "test_with_Kerr_dimmed_more/pics/res",.5*green);
}

if(set==5)
{
  plot_energy(pic, "test_nl_koeff/pics/res");
  plot_energy(pic, "test_lin_koeff/pics/res", red);

  plot_phase(picp, "test_nl_koeff/pics/res", dotted);
  plot_phase(picp, "test_lin_koeff/pics/res", red+dashed);

  //plot_phase(picp, "test_nl_koeff_hi/pics/res", gray);
  //plot_phase(picp, "test_lin_koeff_hi/pics/res", red+0.7*green);
}
  
if(set==4)
{
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=0_/pics/res",0.8*red);
  //plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=0_Nt=2048_/pics/res",0.6*green);
  //plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=0_Nx=2048_/pics/res",0.8*blue);
  //plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=0_dz=-0.25_save_interval=4000_/pics/res");
  plot_energy(pic, "harm_ampW0=0.000744_ratio2w=0_/pics/res",0.5*green);
  plot_energy(pic, "harm_ampW0=0.000144_ratio2w=0_/pics/res",0.8*blue);

  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=0_/pics/res",0.8*red);
  //plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=0_Nt=2048_/pics/res",0.6*green);
  //plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=0_Nx=2048_/pics/res",0.8*blue);
  //plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=0_dz=-0.25_save_interval=4000_/pics/res");
  plot_phase(picp, "harm_ampW0=0.000744_ratio2w=0_/pics/res",0.5*green);
  plot_phase(picp, "harm_ampW0=0.000144_ratio2w=0_/pics/res",0.8*blue);
}

if(set==3)
{
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.033_/pics/res",0.8*red);
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.333_/pics/res",0.5*green);
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.733_/pics/res",0.8*blue);
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1033_/pics/res",0.6*red+0.5*green);
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1333_/pics/res",0.6*red+0.6*blue);
  plot_energy(pic, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1533_/pics/res",0.5*green+0.6*blue);

  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.033_/pics/res",0.8*red);
  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.333_/pics/res",0.5*green);
  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.733_/pics/res",0.8*blue);
  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1033_/pics/res",0.6*red+0.5*green);
  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1333_/pics/res",0.6*red+0.6*blue);
  plot_phase(picp, "harm_ampW0=0.0011744_ratio2w=1_phase2w=0.1533_/pics/res",0.5*green+0.6*blue);
}

if(set==2)
{
  plot_energy(pic, "harm_w2amp=1_w2pahse=0_w2curv=1/pics/res", .8*red);
  plot_energy(pic, "harm_w2amp=1_w2pahse=0.3_w2curv=1/pics/res", .5*green);
  plot_energy(pic, "harm_w2amp=1_w2pahse=0.7_w2curv=1/pics/res", .8*blue);
  plot_energy(pic, "harm_w2amp=1_w2pahse=0_w2curv=1/pics/res");

  plot_phase(picp, "harm_w2amp=1_w2pahse=0_w2curv=1/pics/res", .8*red);
  plot_phase(picp, "harm_w2amp=1_w2pahse=0.3_w2curv=1/pics/res", .5*green);
  plot_phase(picp, "harm_w2amp=1_w2pahse=0.7_w2curv=1/pics/res", .8*blue);
}

if(set==1)
{
  plot_energy(pic, "harm_w2amp=0_w2pahse=0_w2curv=1/pics/res");
  plot_phase(picp, "harm_w2amp=0_w2pahse=0_w2curv=1/pics/res", .8*red);
}

if(set==0)
{
  plot_energy(pic, "test_vacuum/pics/res");
  plot_phase(picp, "test_vacuum/pics/res", .8*red);
}

scale(pic, Linear, Log);
xaxis(pic, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(pic, "$W_{3\omega}$", LeftRight, RightTicks);
size(pic, 8cm, 5cm, point(pic, SW), point(pic, NE));
pic = shift(-point(pic,SW))*pic;
add(pic.fit());

xaxis(picp, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(picp, "$\varphi$ [rad]", LeftRight, RightTicks);
size(picp, 8cm, 5cm, point(picp, SW), point(picp, NE));
picp = shift(-point(picp,SW))*picp;
add(picp.fit(), (0cm, 6.5cm));

ylimits(picr,0);
xaxis(picr, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(picr, "$r$ [mm]", LeftRight, RightTicks);
size(picr, 8cm, 5cm, point(picr, SW), point(picr, NE));
picr = shift(-point(picr,SW))*picr;
add(picr.fit(), (0cm, -6.3cm));

ylimits(picw,0);
xaxis(picw, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(picw, "spec width [$\mathrm{fs}^{-1}$]", LeftRight, RightTicks);
size(picw, 8cm, 5cm, point(picw, SW), point(picw, NE));
picw = shift(-point(picw,SW))*picw;
add(picw.fit(), (0cm, -13cm));

ylimits(pic_total,0);
xaxis(pic_total, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(pic_total, "total energy", LeftRight, RightTicks);
size(pic_total, 8cm, 5cm, point(pic_total, SW), point(pic_total, NE));
pic_total = shift(-point(pic_total,SW))*pic_total;
add(pic_total.fit(), (0cm, -20.5cm));
