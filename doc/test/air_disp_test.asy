import graph;
import signal;

string wpathNO="../../_NODISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";
string w2pathNO="../../_NODISP__harm_ampW0=0_ampW2=0.00001_phase2w=0.0_/";

string wpathBUG="../../_DISP__harm_ampW0=0.00001_ampW2=0_phase2w=0.0_/";
string w2pathBUG="../../_DISP__harm_ampW0=0_ampW2=0.00001_phase2w=0.0_/";

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
    hint = floor(pp.x - 3);
    real z = i * d.dz * z0 * d.save_interval; 
    real xi = z;
    real yi = pp.x * d.dt * t0 * c0;
    sumxiyi += xi*yi;
    sumxi += xi;
    sumyi += yi;
    sumxixi += xi^2;
    n+=1;
    g = g -- (xi, yi);
  }
  slope = (n * sumxiyi - sumyi * sumxi) / (n * sumxixi - sumxi^2);
  shft = (sumyi - slope * sumxi) / n;
  return g;
}

guide gwNO=doplot(wpathNO);
write(slope);
guide g2wNO=doplot(w2pathNO);
write(slope);

write("--------");

guide gwBUG=doplot(wpathBUG);
guide g2wBUG=doplot(w2pathBUG);

guide gwFIX=doplot(wpathFIX);
real wslope = slope;
real wshft = shft;
guide g2wFIX=doplot(w2pathFIX);
real w2slope = slope;
real w2shft = shft;

picture p;
draw(p, gwNO, 0.8*red+0.2*blue);
draw(p, g2wNO, 0.8*blue+0.3*green);

draw(p, gwFIX, dotted+0.8*red+0.2*blue);
draw(p, g2wFIX, dotted+0.8*blue+0.3*green);

//exp(i(wt - kz))
// wt = kz
// w \tau + k0 z = k z = k0 * (1+dn) z
// w \tau = k0 dn z
// c \tau = dn z

void add_label(guide g, string txt, pair pos, pen pn) {
  label(p, txt, (.5*(min(g)+max(g))), pos, pn);
}

add_label(gwFIX, format("%.3g",wslope), SE, 0.8*red+0.2*blue);
add_label(gwFIX, format("%.3g",w2slope), NW, 0.8*blue+0.3*green);
  
//label(p, format("%.3g",wslope), (.5*(min(gwFIX)+max(gwFIX))), SE, 0.8*red+0.2*blue);
//label(p, format("%.3g",w2slope), (.5*(min(gwFIX)+max(gwFIX))), NW, 0.8*blue+0.3*green);

draw(p, gwBUG, dashed+0.8*red+0.2*blue);
draw(p, g2wBUG, dashed+0.8*blue+0.3*green);

void draw_approx(guide g, real a, real b) {
  draw(p, (min(g).x, min(g).x * a + b) -- (max(g).x, max(g).x * a + b), opacity(.3)+linewidth(0.07pt));
}

draw_approx(gwFIX, wslope, wshft);
draw_approx(g2wFIX, w2slope, w2shft);

xaxis(p, "$z$ [cm]", BottomTop, LeftTicks);
yaxis(p, "$c \cdot \tau$ [cm]", LeftRight, RightTicks);
size(p, 7cm, 5cm, point(p,SW), point(p,NE));
add(p.fit());

real refIndex_um(real lam) {
  return 0.05792105/(238.0185-1/lam^2) + 0.00167917/(57.362-1/lam^2);
}

write(wslope);
write(refIndex_um(.78));
write("diff:" , wslope - refIndex_um(.78));
write("--------");
write(w2slope);
write(refIndex_um(.78/2));
write("diff:" , w2slope - refIndex_um(.78/2));
