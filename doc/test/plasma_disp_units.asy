import graph;
import signal;

string wpath_weak="../../_WEAK_PLASMA_DISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";

string wpath="../../_PLASMA_DISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";

//we need this to get pure plasma dispersion
string wpathFIX="../../_FIXDISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";

int filecount = 91;

real slope;
real shft;

guide doplot(string p) {
  descr d = read_descr(p+"descr.txt");
  guide g;
  int hint=-1;
  //y=ax+b
  // \sum_i (yi-a*xi-b)^2
  // sum_i xi * (yi-a*xi-b) = 0
  // sum_i (yi-a*xi-b) = 0
  // xi yi - a(xi xi) - b xi = 0
  // yi - a xi - n b = 0
  // n sum(xi yi) - sum(yi) sum(xi) = a(n sum(xi xi)-sum(xi)^2)
  real sumxiyi=0;
  real sumxi=0;
  real sumyi=0;
  real sumxixi=0;
  real n = 0;
  
  for(int i=0;i<filecount;++i) {
    real[] signal = read_signal(p+"pics/res/"+"field_axis"+number(i)+".dat");
    if(hint == -1)
      hint = floor(signal.length/2)-4;
    pair pp = getPeriodAndPhase(signal, hint);
    hint = floor(pp.x - 7);
    real z = i * d.dz * d.save_interval; 
    real xi = z;
    real yi = pp.x * d.dt;
    if(i<filecount/2) {
      sumxiyi += xi*yi;
      sumxi += xi;
      sumyi += yi;
      sumxixi += xi^2;
      n+=1;
    }
    g = g -- (xi, yi);
  }
  slope = (n * sumxiyi - sumyi * sumxi) / (n * sumxixi - sumxi^2);
  shft = (sumyi - slope * sumxi) / n;
  return g;
}

doplot(wpathFIX);
real air_w_slope = slope;

guide gw_weak=doplot(wpath_weak);
real weak_w_slope = slope;

guide gw=doplot(wpath);
real w_slope = slope;

write("--------");


picture p;
draw(p, gw, 0.8*red+0.2*blue);

draw(p, gw_weak, dashed+0.8*red+0.2*blue);

//exp(i(wt - kz))
// wt = kz
// w \tau + k0 z = k z = k0 * (1+dn) z
// w \tau = k0 dn z
// c \tau = dn z
//label(p, format("%.3g",wslope), (.5*(min(gwFIX)+max(gwFIX))), SE, 0.8*red+0.2*blue);
//label(p, format("%.3g",w2slope), (.5*(min(gwFIX)+max(gwFIX))), NW, 0.8*blue+0.3*green);

//draw(p, (min(gwFIX).x, min(gwFIX).x * wslope + wshft) -- (max(gwFIX).x, max(gwFIX).x * wslope + wshft), linewidth(0.07pt));
//draw(p, (min(g2wFIX).x, min(g2wFIX).x * w2slope + w2shft) -- (max(g2wFIX).x, max(g2wFIX).x * w2slope + w2shft), linewidth(0.07pt));

xaxis(p, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(p, "$c \cdot \tau$ [cm]", LeftRight, RightTicks);
size(p, 7cm, 5cm, point(p,SW), point(p,NE));
add(p.fit());

real[] plasma_weak = read_signal(wpath_weak+"pics/res/plasma0000000.dat");
real[] plasma = read_signal(wpath+"pics/res/plasma0000000.dat");

real w = 1.;

real n_weak = plasma_weak[floor(plasma_weak.length/2)];
real wp2_weak = n_weak;
write("weak_w_slope: ", weak_w_slope - air_w_slope);
write("a_weak_w_slope: ", - wp2_weak/w^2);

write("----------------------------------------------");
write("----------------------------------------------");

real n = plasma[floor(plasma.length/2)];
real wp2 = n;
write("w_slope: ", w_slope - air_w_slope);
write("a_w_slope: ", - wp2/w^2);


