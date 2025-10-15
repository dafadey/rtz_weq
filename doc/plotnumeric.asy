import graph;

file inf=input("SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_ampW0=0.000085_airDensity=0.0122_/pics/res/src0000000src_pw.dat").line();
string[] raw = inf;

write(raw[0]);
int n = (int) split(raw[0], '\t')[1];

pair[] src = new pair[n];
for(int i = 0; i < n; ++i) {
  string[] items = split(raw[i+1], '\t');
  src[i] = ((real) items[0], (real) items[1]);
}

real n0 = 1.83e21; //n0 in cm-3

guide gr, gi, ga;
for(int i=0; i<n; ++i) {
  gr = gr -- (i, src[i].x*n0);
  gi = gi -- (i, src[i].y*n0);
  ga = ga -- (i, abs(src[i])*n0);
}

picture pic;

draw(pic, gr, red);
draw(pic, gi, blue);
draw(pic, ga);

xaxis(pic, "$i$", BottomTop, LeftTicks);
yaxis(pic, "src", LeftRight, RightTicks);
size(pic, 7cm, 5cm, point(pic, SW), point(pic, NE));
add(pic.fit());

