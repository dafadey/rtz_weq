string path="/home/dan/rtz_weq/harm_w2amp=0_w2pahse=0_w2curv=1/pics/res";

import graph;
import fft;

defaultpen(fontsize(8pt));

real[] read(int id)
{
  file inf = input(path + "/field_axis" + string(id) + ".dat").line();
  string[] raw = inf;
  real[] res;
  int n = (int) raw[2];
  for(int i = 0; i < n; ++i)
    res.push((real) raw[i + 4]);
  return res;
}

picture pic;

void plot(picture pic, int id, pen p = defaultpen)
{
  real[] data = read(id);

  qqFFT_freal(data.length, data);

  guide g;
/*
  for(int i=0;i<data.length;++i)
    g = g -- (i, data[i]);
*/
  for(int i=0;i<floor(data.length/8);++i)
    g = g -- (i, (data[i*2]*data[i*2]+data[i*2+1]*data[i*2+1])^.5);


  draw(pic, g, p);
}

for(int i=0;i<80;i=i+4)
  plot(pic, i, i/80 * red);

plot(pic, 38, blue);
plot(pic, 39, .5*green);

xaxis(pic,"$\omega_i$",BottomTop,LeftTicks);
yaxis(pic,"$E$",LeftRight,RightTicks);
size(pic, 8cm, 5cm, point(pic, SW), point(pic, NE));
add(pic.fit());
