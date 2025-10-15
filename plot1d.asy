string path="/home/dan/rtz_weq/harm_w2amp=0_w2pahse=0_w2curv=1/pics/res";

import graph;
import fft;

defaultpen(fontsize(8pt));
defaultpen(linewidth(.1pt));

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

  //qqFFT_freal(data.length, data);

  guide g;
  real absmax = max(max(data), abs(min(data)));

  for(int i=0;i<data.length;++i)
    g = g -- (i, data[i]/absmax);

  draw(pic, g, p);
}

for(int i=0;i<118;i=i+58)
  plot(pic, i, i/118 * red);

plot(pic, 77, blue);
plot(pic, 79, .5*green);

xaxis(pic,"$\tau_i$",BottomTop,LeftTicks);
yaxis(pic,"$E$",LeftRight,RightTicks);
size(pic, 8cm, 1cm, point(pic, SW), point(pic, NE));
add(pic.fit());
