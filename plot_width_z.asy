import graph;
import patterns;

add("crosshatch", crosshatch(.7mm, .8*cyan+.1*blue));
add("crosshatchdense", crosshatch(.3mm, .4*cyan+.1*blue+.5*red));

defaultpen(linewidth(0.5pt));
defaultpen(fontsize(8pt));

int lastFile=92;

real k0=80553.65779; // cm^-1
real x0 = 1. / (2.^.5 * k0);
real t0 = 1. / k0;
real z0 = 2/k0;
real dz = .5*z0;
real n0 = 1.83485883*10.^21;
real Ea_Vcm = 5.14e9;


struct descr {
  real dr;
  real dz=0.5;
  int save_interval=333;
  real dt;
  real dw;
  real T;
  real ampW;
  real r_0;
  real z_shift;
  real k0;
  real airDensity=0;
};

descr read_descr(string filename) {
  descr d;
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  int Nx;
  int Nt;
  real Lx;
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
    if(item == "Lx") Lx = val;
    if(item == "Nx") Nx = (int) val;
    if(item == "T") d.T = val;
    if(item == "Nt") Nt = (int) val;
    if(item == "airDensity") d.airDensity = val;
  }
  d.dr = Lx/Nx;
  d.dt = d.T/Nt;
  return d;
}

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

real[][] read_1ddata(string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  int n = (int) raw[2];
  int items_count = (int) split(raw[1], '\t')[1];
  real[][] data;
  for(int i = 0; i < n; ++i) {
    string[] sitems = split(raw[i+4], '\t');
    real[] items;
    for(int j = 0; j < sitems.length; ++j)
      items.push((real) sitems[j]);
    data.push(items);
  }
  return data;
}

pair[][] read_2dcomplexdata(string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  pair[][] res;
  int lineid=0;
	while(true) {
		int n = (int) split(raw[lineid], '\t')[1];
		lineid += 1;
		pair[] data;
		for(int i = 0; i < n; ++i) {
			string[] sitems = split(raw[lineid], '\t');
			lineid += 1;
			data.push(((real) sitems[0], (real) sitems[1]));
		}
		res.push(data);
		if(lineid == raw.length)
			break;
	}
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
//ew, e2w, e3w, rw, r2w, r3w, w, w2, w3, total
real get_width(string data_path, int zi)
{
  string filename = data_path+"/"+"field_energy"+number(zi)+".dat";
  write("reading file "+filename);
  real[][] data;
  read_point(data, filename);

  return data[0][3];
}


real get_a(string data_path, int zi)
{
  string filename = data_path+"/"+"field_axis"+number(zi)+".dat";
  write("reading file "+filename);
  real[][] data = read_1ddata(filename);
	real maxA = 0;
	for(int i=0; i<data.length; ++i)
		maxA = max(abs(data[i][0]), maxA);
  return maxA;
}


real get_src(string filepath, int zi, string postfix) {
  real[] res;
	pair[][] data = read_2dcomplexdata(filepath+"src"+number(zi)+postfix+".dat");
	descr desc = read_descr(filepath+"/../../descr.txt");
	//write("read cplx data ",i);

	real maxA3 = .0;
	for(int tt=0; tt<data[0].length; ++tt)
		maxA3 = max(abs(data[0][tt]), maxA3);// * 2. * pi * (3/8+xx) * (desc.dr * x0)^2;
  return maxA3;
}

real mylog10(real x)
{
  return x > 1e-5 ? log10(x) : -5;
}

picture pic, pica, picsa;

void plot_width(string pth, pen p) {
	guide g;
	descr d = read_descr(pth + "/descr.txt");
	for(int zi=0; zi<lastFile; ++zi) {
		real z = d.dz * d.save_interval * zi + d.z_shift;
		real w = get_width(pth+'/pics/res/', zi);
		g=g--(z*z0*10, w*x0*10^4);
	}
	draw(pic, g, p);
}

void plot_a(string pth, pen p) {
	guide g;
	descr d = read_descr(pth + "/descr.txt");
	real a0 = get_a(pth+'/pics/res/', 0);
	for(int zi=0; zi<lastFile; ++zi) {
		real z = d.dz * d.save_interval * zi + d.z_shift;
		real a = get_a(pth+'/pics/res/', zi);
		g=g--(z*z0*10, a/a0);
	}
	draw(pica, g, p);
}

void plot_sa(string pth, pen p) {
	guide g;
	descr d = read_descr(pth + "/descr.txt");
	real s0 = get_src(pth+'/pics/res/', 0, "_src_www");
	for(int zi=0; zi<lastFile; zi+=7) {
		write(zi);
		real z = d.dz * d.save_interval * zi + d.z_shift;
		real s = get_src(pth+'/pics/res/', zi, "_src_www");
		g=g--(z*z0*10, (s/s0)^(1/3));
	}
	draw(picsa, g, p);
}




plot_width("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=0.8_",.8*blue*.8+.3*green);
plot_width("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1._",.7*green);
plot_width("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1.25_",.8*red);

plot_a("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=0.8_",.8*blue*.8+.3*green);
plot_a("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1._",.7*green);
plot_a("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1.25_",.8*red);

plot_sa("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=0.8_",.8*blue*.8+.3*green);
plot_sa("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1._",.7*green);
plot_sa("harm_Kerr=0.00754_collimated=1_townes_mode=1_townes_mode_factor=1.25_",.8*red);


ylimits(pic, 0);
xaxis(pic, "$z$ [mm]", BottomTop, LeftTicks);
yaxis(pic, "$R$ [um] (half maximum by intensity)", LeftRight, RightTicks);
size(pic, 7.5cm, 5cm, point(pic, SW), point(pic, NE));

ylimits(pica, 0);
xaxis(pica, "$z$ [mm]", BottomTop, LeftTicks);
yaxis(pica, "$E_{max}/E_{max}(z=0)$", LeftRight, RightTicks);
size(pica, 7.5cm, 5cm, point(pica, SW), point(pica, NE));

ylimits(picsa, 0);
xaxis(picsa, "$z$ [mm]", BottomTop, LeftTicks);
yaxis(picsa, "$A_{\omega}/A_{\omega}(z=0)$", LeftRight, RightTicks);
size(picsa, 7.5cm, 5cm, point(picsa, SW), point(picsa, NE));



void myaddpic(picture p, pair pos=(.0, .0)) {
  p = shift(-point(p,SW))*p;
  add(p.fit(), pos);
}

myaddpic(pic);

myaddpic(pica, (0, 6.2cm));

myaddpic(picsa, (0, 12.4cm));

write("writing output to "+defaultfilename);
