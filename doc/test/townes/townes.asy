//laplassian A + 3/8 K A^3 - alpha A = 0
import graph;

real eps0 = 8.85*1e-12; //SI epsilon0 dielectric constant of vacuum
real c = 3*1e8; //SI speed of light
real E_unit = 5.14 * 1e11; // SI atomic field
real lambda = 780*1e-9; // SI wavelength
real k0 = 2*pi/lambda; // SI wave number
real z_unit = 2/k0;
real tau_unit = (lambda/c)/(2*pi);
real r_unit = 1/sqrt(2)/k0;

string str(real v) {return format("%.3g",v);}

write("z_unit = "+str(z_unit*1e6)+" um");
write("tau_unit = "+str(tau_unit*1e15)+" fs");
write("r_unit = "+str(r_unit*1e6)+" um");

struct townsend
{
  real K = 0.00725;//0.0377;//0.035;//0.00943;
  int n=2048;
  real dx = 2500/n;
 
  real LaplassianFactor=2;
  
  real A(int i) {
  //  return ((i==0 ? 1:2)+1/(i+.5))/dx/dx;
    return -(2+1/(i+.5))/dx/dx;
  }
  
  real B(int i) {
    return 1/dx/dx;
  }
  
  real C(int i) {
  //  return (-1-1/(i+.5))/dx/dx;
    return ((i==0 ? 2 : 1)+1/(i+.5))/dx/dx;
  }
  
  real metric(int i) {
    return (i==0 ? 3/8 : i+.5) * dx * 2 * pi;
  }
  
  real[] solve(real f_0, real alpha) {
    real[] res = new real[n];
    res[0] = f_0;
    for(int i=0;i<n-1;++i)
    {
      real f_i_1 = i>0 ? res[i-1] : .0;
      real f_i = res[i];
      real f_i1 = -(LaplassianFactor * (f_i * A(i) + B(i) * f_i_1) + 3/4*K*f_i^3 - alpha * f_i) / (LaplassianFactor*C(i));
      res[i+1] = f_i1;
    }
    return res;
  }
  
  bool trend(real[] f) {
    for(int i=0;i<f.length;++i) {
      if(f[i]<0)
        return false;
    }
    return true;
  }
  
  real f0 = 0.0611344;
  real r0 = 395.98;
  real alpha0;
  
  void upd_r0(real f0) {
    r0 = (4*LaplassianFactor/(2*3/4*f0^2*K))^.5;
  }
  
  void upd_alpha0(real f0) {
    alpha0 = 3/4*f0^2*K-4*LaplassianFactor/r0^2;
  }

  real[] solg;
  
  void gauss()  {
    solg = new real[n];
    for(int i=0;i<n;++i)
      solg[i] = f0 * exp(-(i*dx)^2/r0^2);
  }

  real[] solution;
  real alpha;
  
  void find()
  {
    alpha = alpha0;
    real dalpha = alpha0/2;
    //write(alpha);
    //write(f0);
    //write("----");
    while(abs(dalpha) > abs(alpha0*1e-13)) {
      real[] s = solve(f0, alpha);
      guide g;
      for(int i=0;i<n;++i)
        g=g--(i*dx, s[i]);
      if((!trend(s) && dalpha<0) || (trend(s) && dalpha>0))
        dalpha = -dalpha/2;
      alpha+=dalpha;
      //write((alpha,dalpha));
    }
    solution = solve(f0, alpha);
  }
  
  void find_by_f0() {
    alpha = alpha0;
    real f00=f0;
    real df0 = f0/2; 
    while(abs(df0) > abs(f00*1e-13)) {
      real[] s = solve(f0, alpha);
      guide g;
      for(int i=0;i<n;++i)
        g=g--(i*dx, s[i]);
      if((!trend(s) && df0>0) || (trend(s) && df0<0))
        df0 = -df0/2;
      f0+=df0;
      //write((f0,df0));
    }
    solution = solve(f0, alpha);
  }
  
  void upd() {
    upd_r0(f0);
    upd_alpha0(f0);
    gauss();
    find();
  }
  
  void plotgauss(picture pic)
  {
    guide g;
    for(int i=0;i<n;++i)
      g=g--(i*dx, solg[i]);
    draw(pic,g,blue);
    guide ga;
    real a = (-2 * solg[0] + 2 * solg[1])/dx/dx/2;
    for(int i=0;i<n;++i) {
      real r = i*dx;
      real v = f0+a*r^2;
      if(v>0)
        ga=ga--(i*dx, v);
    }
    draw(pic,ga,dotted);
  }
  
  void plotgauss_check(picture pic)
  {
    real[] term1=new real[n];
    for(int i=0;i<n-1;++i)
      term1[i] = LaplassianFactor*((i==0?0:solg[i-1]*B(i)) + solg[i]*A(i) + solg[i+1]*C(i));
    term1[n-1] = 0;
    real[] term2=new real[n];
    for(int i=0;i<n;++i)
      term2[i] = K*solg[i]^3*3/4-alpha0*solg[i];
    guide gt1,gt2;
    real factor = .5 * solg[0]/max(abs(term1[0]),abs(term2[0]));
    for(int i=0;i<n;++i) {
      gt1=gt1--(i*dx, -factor*term1[i]);
      gt2=gt2--(i*dx, factor*term2[i]);
    }
    draw(pic,gt1,.4*green+dashed);
    draw(pic,gt2,.5*blue+dotted);
  }

  void plot(picture pic)
  {
    guide g;
    for(int i=0;i<n;++i)
      g=g--(i*dx, solution[i]);
    draw(pic,g,red);
  }
  
  void plot_check(picture pic)
  {
    guide ga;
    real a = (-2 * solution[0] + 2 * solution[1])/dx/dx/2;
    for(int i=0;i<n;++i) {
      real r = i*dx;
      real v = f0+a*r^2;
      if(v>0)
        ga=ga--(i*dx, v);
    }
    draw(pic,ga,red+dotted);
    real[] term1=new real[n];
    for(int i=0;i<n-1;++i)
      term1[i] = LaplassianFactor*((i==0?0:solution[i-1]*B(i)) + solution[i]*A(i) + solution[i+1]*C(i));
    term1[n-1] = 0;
    real[] term2=new real[n];
    for(int i=0;i<n;++i)
      term2[i] = K*solution[i]^3*3/4-alpha*solution[i];
    guide gt1,gt2;
    real factor = .5 * solution[0]/max(abs(term1[0]),abs(term2[0]));
    for(int i=0;i<n;++i) {
      gt1=gt1--(i*dx, -factor*term1[i]);
      gt2=gt2--(i*dx, factor*term2[i]);
    }
    draw(pic,gt1,0.5*red+.4*green+dashed);
    draw(pic,gt2,0.5*red+.5*blue+dotted);
  }
  
  real en(real[] s) {
    real res=0;
    for(int i=0;i<n;++i)
      res += s[i]^2 * metric(i) * dx;
    return res;
  }
  
  real Energy() {
    return en(solution);
  }

  real EnergyG() {
    return en(solg);
  }
  
};

townsend T;
T.f0=0.08;
T.upd();

//write(2*3/4*T.f0^2*T.K);
//write(8/T.r0^2);
real energyF = .5 * c * eps0 * (E_unit *  r_unit)^2;

write("Classic Pcr Sasha : " + str(1.86*lambda^2/4/pi/((7.4*.8+9.2*.2)*1e-24)/1e9) + " GW");

write("Classic Pcr : " + str(1.86*lambda^2/4/pi/(4*1e-23)/1e9) + " GW");

write("Townsend power : " + str(T.Energy()*energyF/1e9)+" GW");
write("Gauss power : " + str(T.EnergyG()*energyF/1e9)+" GW");
write("ratio = " + string(T.Energy()/T.EnergyG()));
picture pic;
//gauss:
T.plotgauss(pic);
T.plotgauss_check(pic);

T.plot(pic);
T.plot_check(pic);

xaxis(pic,"$r$",BottomTop,LeftTicks);
yaxis(pic,"$f$",LeftRight,RightTicks);
size(pic,13cm,4cm,point(pic,SW),point(pic,NE));
add(pic.fit());

/*
picture pic_test;
T.LaplassianFactor=1.;
T.alpha0=1.;
T.f0=2.2;
T.dx=11/T.n;
T.K=4/3;
T.find_by_f0();
write(T.Energy()/2/pi);
T.plot(pic_test);
T.plot_check(pic_test);
xaxis(pic_test,"$r$",BottomTop,LeftTicks);
yaxis(pic_test,"$f$",LeftRight,RightTicks);
size(pic_test,13cm,4cm,point(pic_test,SW),point(pic_test,NE));
add(pic_test.fit(), (0,5cm));
*/