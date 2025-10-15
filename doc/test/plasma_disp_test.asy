import graph;
import signal;

real e0=4.8032 * 10.^(-10);
real m0=9.1094 * 10.^(-28);

string wpath_weak="../../_WEAK_PLASMA_DISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";
string w2path_weak="../../_WEAK_PLASMA_DISP__harm_ampW0=0_ampW2=0.00001_phase2w=0.0_/";

string wpath="../../_PLASMA_DISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";
string w2path="../../_PLASMA_DISP__harm_ampW0=0_ampW2=0.00001_phase2w=0.0_/";

//we need this to get pure plasma dispersion
string wpathFIX="../../_FIXDISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";
string w2pathFIX="../../_FIXDISP__harm_ampW0=0_ampW2=0.00001_phase2w=0.0_/";

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
    real z = i * d.dz * z0 * d.save_interval; 
    real xi = z;
    real yi = pp.x * d.dt * t0 * c0;
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

picture p;

void draw_approx(guide g, real a, real b) {
  draw(p, (min(g).x, min(g).x * a + b) -- (max(g).x, max(g).x * a + b), opacity(.3)+linewidth(0.07pt));
}

doplot(wpathFIX);
real air_w_slope = slope;
doplot(w2pathFIX);
real air_w2_slope = slope;

guide gw_weak=doplot(wpath_weak);
draw_approx(gw_weak, slope, shft);
real weak_w_slope = slope;
guide g2w_weak=doplot(w2path_weak);
draw_approx(g2w_weak, slope, shft);
real weak_w2_slope = slope;

guide gw=doplot(wpath);
draw_approx(gw, slope, shft);
real w_slope = slope;
guide g2w=doplot(w2path);
draw_approx(g2w, slope, shft);
real w2_slope = slope;

write("--------");

void add_label(guide g, string txt, pair pos, pen pn) {
  label(p, txt, point(g, (int)(filecount/2)), pos, pn+fontsize(8pt));
}

draw(p, gw, 0.8*red+0.2*blue);
add_label(gw, format("%.3g", w_slope), SW, 0.8*red+0.2*blue);
draw(p, g2w, 0.8*blue+0.3*green);
add_label(g2w, format("%.3g", w2_slope), NE, 0.8*blue+0.3*green);

draw(p, gw_weak, dashed+0.8*red+0.2*blue);
add_label(gw_weak, format("%.3g", weak_w_slope), NE, 0.8*red+0.2*blue);
draw(p, g2w_weak, dashed+0.8*blue+0.3*green);
add_label(g2w_weak, format("%.3g", weak_w2_slope), N, 0.8*blue+0.3*green);


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

real w = 2*pi/(2*pi*t0);
real ww = 2*pi/(.5*2*pi*t0);

real n_weak = plasma_weak[floor(plasma_weak.length/2)] * n0;
real wp2_weak = 4*pi*n_weak*e0^2/m0;
write("weak: n=", n_weak);
write("weak_w_slope: ", weak_w_slope - air_w_slope);
write("a_weak_w_slope: ", - wp2_weak/2/w^2);
write("A_weak_w_slope: ", sqrt(1 - wp2_weak/w^2) - 1);
write("----------------------------------------------");
write("weak_w2_slope: ", weak_w2_slope - air_w2_slope);
write("a_weak_w2_slope: ", - wp2_weak/2/ww^2);
write("A_weak_w2_slope: ", sqrt(1 - wp2_weak/ww^2) - 1);

write("----------------------------------------------");
write("----------------------------------------------");

real n = plasma[floor(plasma.length/2)] * n0;
real wp2 = 4*pi*n*e0^2/m0;
write("n=", n);
write("w_slope: ", w_slope - air_w_slope);
write("a_w_slope: ", - wp2/2/w^2);
write("A_w_slope: ", sqrt(1 - wp2/w^2) - 1);
write("----------------------------------------------");
write("w2_slope: ", w2_slope - air_w2_slope);
write("a_w2_slope: ", - wp2/2/ww^2);
write("A_w2_slope: ", sqrt(1 - wp2/ww^2) - 1);


write("n0=",(2*pi/(780*10.^(-9)*100/c0))^2*m0/4/pi/e0^2);

