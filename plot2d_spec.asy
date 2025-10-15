string path="/home/dan/rtz_weq/harm_w2amp=1_w2pahse=0_w2curv=1/pics/res";

import graph;
import fft;

defaultpen(fontsize(8pt));

real[][] read(int id)
{
  file inf = input(path + "/field" + string(id) + ".dat").line();
  string[] raw = inf;
  write(raw[2]);
  int nt = (int) split(raw[2],'\t')[0];
  int nr = (int) split(raw[2],'\t')[1];
  write("nt=",nt);
  write("nr=",nr);
  real[][] res = new real[nr][nt];
  for(int j = 0; j < nr; ++j)
  {
    for(int i = 0; i < nt; ++i)
      res[j][i] = ((real) raw[j*nt + i + 4]);
  }
  return res;
}

picture pic;

void plot(picture pic, int id, pen p = defaultpen)
{
  real[][] data = read(id);

  int nr = data.length;
  int nt = data[0].length;
  for(int i=0;i<nr;++i)
    qqFFT_freal(data[i].length, data[i]);

  guide g;

  int n = floor(nt/2);
  real[] spec_aperture = new real[n];
  for(int i=0;i<n;++i)
  {
    spec_aperture[i]=0;
    for(int j=0;j<nr;++j)
      spec_aperture[i] += (j+0.5)*(data[j][i*2]*data[j][i*2]+data[j][i*2+1]*data[j][i*2+1]);
  }


  for(int i=0;i<floor(n/4);++i)
    g = g -- (i, spec_aperture[i]);


  draw(pic, g, p);
}

for(int i=0;i<80;i=i+16)
  plot(pic, i, i/80 * red);

plot(pic, 38, blue);
plot(pic, 39, .5*green);

xaxis(pic,"$\omega_i$",BottomTop,LeftTicks);
yaxis(pic,"$E$",LeftRight,RightTicks);
size(pic, 8cm, 5cm, point(pic, SW), point(pic, NE));
add(pic.fit());
