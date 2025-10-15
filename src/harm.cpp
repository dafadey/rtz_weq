#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>
#include <cufft.h>
#include <math_constants.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>
#include <unistd.h>

#include "getGPUtemp.h"
#include "getTemps.h"

#include "dfs.h"
#include "queue.h"
#include "postprocess.h"
#include "fft_real_3.cpp"
#include "host_methods.h"

int save_interval = 333;
double phase2w=1.57;
double ratio2w=0.0; //ratio of second harmonic to first harmonic (by amplitude)
double ampW2=0.0;
double ampW3=0.0;
double curvativity2w=1.0;
double ampW0=0.0007;//0.000987;//0.0011744335174163965;//0.015;//0.003;//0.01;
double airDensity = .0;
double Kerr = 0.0377;//0.00943;//0.014156; // for lambda 780nm, and E measured in E_a
double focalLength = 85.; //mm
double beamDia = 5.; //mm FWHM by I
double zmax=0;
double k0;
double tunExp=1.;
double tunFac=1.;
double E_a = 5.14220675112 * 1.e11; // V/m
double e0 = 8.8541878188 * 1.e-12; // F/m
double c = 2.99792458 * 1.e8; //m/s
//numerical conf:
double Lx=2500.0F;
double T=250.0F;
double TAU=35.0F;
int Nx=1024*2;
int Nt=1024;
bool townes_mode=false;
double townes_mode_factor = 1.;
bool collimated=false;
bool skipPics=false;
bool skip1D=false;
bool skipEnergies=false;
bool skipSources=false;
bool saveW3envelope=false;
bool saveFieldEnvelope=false;
bool singleStepInFocus=false;
bool dynamic_save_interval=false;

std::string folder_name;
FL_DBL dz;

bool cool_dev(int dev)
{
	int t=getGPUtemp(dev);
	if(t < 75)
		return true;
	while(getGPUtemp(dev)>75)
		sleep(7);
	printf("temp is %d\n",t);
	return true;
}

std::vector<double> solveFlatEikonal(double amp0, double alpha) {
	std::vector<double> f(Nx);
	double dr = Lx / (double) Nx;
	f[0] = amp0;
	for(int i=0;i<Nx-1;i++) {
		f[i+1] = -(2. * (f[i] * LaplassianA(i, dr) + (i>0 ? f[i-1] * LaplassianB(i, dr) : .0)) + 3./4.*Kerr*f[i]*f[i]*f[i] - alpha * f[i]) / (2. * LaplassianC(i,dr));
	}
	return f;
}

std::vector<double> getTownsend(double amp, double r) {
	double alpha0 = 3./4.*amp*amp*Kerr-8./(r*r);
	double alpha = alpha0;
	double dalpha = alpha*.5;
	std::cout << "alpha0=" << alpha0 << ", amp=" << amp << '\n';
	auto trend = [](const std::vector<double> in) -> bool {
		for(const auto& v : in) {
			if(v < 0)
				return false;
		}
		return true;
	};
	int count=0;
	while(std::abs(dalpha) > std::abs(alpha0*1.e-13)) {
		auto sol = solveFlatEikonal(amp, alpha); //reallocates every time. bad, but ok
		if((!trend(sol) && dalpha < .0) || (trend(sol) && dalpha > .0))
			dalpha *= -.5;
		alpha += dalpha;
		count++;
		std::cout << "alpha=" << alpha << " dalpha=" << dalpha << '\n';
		if(count>113)
			break;
	}
	return solveFlatEikonal(amp, alpha);
}

int setinitdata_harm(FL_DBL* arr, int Nt, int Nx, FL_DBL dt, FL_DBL dx)
{
	
	Nt=Nt-2;
	
	k0 = 80553.65779; //k0[cm^-1] for lambda = 780nm
	double z0 = 2. / k0;
  double x0 = 1. / (sqrt(2.) * k0); // .1 to convert to cm
  double t0 = 1. / k0;
  double d0 = .1 * beamDia / x0; // .1 to convert to cm
  d0 = d0 * 1.7; // convert from FWHM to A/e notation for beam dia/R.
  double d02 = d0 * d0;
  double z_foc = .1 * focalLength / z0; // .1 to convert to cm
  
  // r^2 = r_0^2+64 z^2 / r_0^2
  // r_0^2 * r^2 = r_0^4 + 64 * z^2
  // r_0^4 - r_0^2 * r^2 + r^4/4 - r^4/4 + 64 * z^2 = 0
  // (r_0^2 - r^2/2)^2 = r^4/4 - 64*z^2
  // (d_0^2/4 - d^2/8)^2 = d^4/4/16 - 64 * z^2
  // (d_0^2 - d^2/2)^2/16 = d^4/4/16 - 64 * z^2
  // (d_0^2 - d^2/2)^2 = d^4/4 - 1024 * z^2
  // d_0^2 = d^2/2 - sqrt(d^4/4-1024 * z^2) 
  
  double dw = std::sqrt(d02/2. - std::sqrt(d02*d02/4.-1024.*z_foc*z_foc));

	double z_shift = 13402./1.7 ;//.4/z0*160./dw; //some value for z-shift

  zmax = singleStepInFocus ? std::abs(.5*dz) : 2. * z_shift;

  //z_shift = .0;
  //zmax *= 0.5;

	double aaa=dw/2.; //focal radus of the beam

	double r_0 = (singleStepInFocus || collimated) ? aaa : sqrt(pow(z_shift/aaa,2.0)*64.0+pow(aaa,2.0));
	double V = (singleStepInFocus || collimated) ? 0 : 8.0*z_shift/(64.0*pow(z_shift,2.0)+pow(aaa,4.0));
	
  double ampW_init = ampW0 * d0 / (2. * r_0);
  
  double ampW2_init = ampW2 * d0 / (2. * r_0);

  double ampW3_init = ampW3 * d0 / (2. * r_0);

  if(collimated) {
	  ampW_init /= 7.;
	  ampW2_init /= 7.;
	  ampW3_init /= 7.;
	  r_0 *= 7.;
	}
  
  std::cout << "r_0=" << r_0 << '\n';
  std::cout << "ampW_init=" << ampW_init << '\n';
  std::cout << "8/r0^2=" << 8. / std::pow(r_0, 2.) << " vs A_init^2/2*Kerr=" << std::pow(ampW_init, 2.) / 2. * Kerr << '\n';
  std::cout << "dw = " << dw * x0 << " cm\n";
  std::cout << "d0 = " << d0 * x0 << " vs " << 2.*sqrt(pow(z_foc/aaa, 2.)*64. + pow(aaa,2.)) * x0 << " cm\n";
  std::cout << "d_start = " << r_0 * 2. * x0 << " cm\n";
  std::cout << "z_shift = " << z_shift * z0 << " cm\n";
  
  double tau_0=(TAU/2.66*2.0*Pi);
	
	std::vector<double> townes;
	if(townes_mode)
	  townes = getTownsend(ampW_init, r_0);
  
  auto envR=[&](int i) {
		if(townes_mode)
			return townes[i]/townes[0]*townes_mode_factor;
		else
			return exp(-pow(i*dx/r_0,2.0));
	};
  auto envT=[&](double arg) {return std::abs(arg) < M_PI ? .5+.5*cos(arg) : .0;};
	
	for(int j=0;j<Nx;j++)
	{
		double r=double(j)*dx;
		for(int i=0;i<Nt;i++)
		{
			double t=double(i-Nt/2)*dt;
			double arg = (t+r*r*V)/tau_0*M_PI;
			double env = envT(arg);
			arr[j*Nt+i]=ampW_init*env*envR(j)*sin(t+r*r*V) 
                  + (ampW2 == .0 ? ratio2w * ampW_init : ampW2_init)*std::pow(env * envR(j), 2.) * sin(2.0*t+r*r*V*curvativity2w*2.0+phase2w)
                  + ampW3_init * std::pow(env * envR(j), 3.) * sin(3.0 * t+r*r*V*3.);
		}
	}
  
  double P = .0;
  for(int i=0;i<Nx;i++)
    P += pow(ampW_init * envR(i), 2.) * LaplassianMetric(i,dx) * dx;

  P *= .5 * c * e0 * std::pow(E_a * x0/1.e2, 2.);
  
  double I = .5 * e0 * c * std::pow(ampW_init * E_a, 2.); // W/m^2
  double I_lens = .5 * e0 * c * std::pow(ampW0 * E_a, 2.); // W/m^2
  
  std::cout << "P_init = " << P * 1.e-9  << " GW\n";
  std::cout << "I_init = " << I << " [W/m^2]\n";
  std::cout << "I_lens = " << I_lens << " [W/m^2]\n";
  
  std::ofstream descr(folder_name+"/descr.txt");
	if(!dynamic_save_interval)
		descr << "save_interval=" << save_interval << '\n';
  descr << "Lx=" << Lx << '\n';
  descr << "T=" << T << '\n';
  
  descr << "Nx=" << Nx << '\n';
  descr << "Nt=" << Nt << '\n';
  descr << "dz=" << dz << '\n';
  
  descr << "k0=" << k0 << '\n';
  descr << "dw=" << dw << '\n';
  descr << "TAU=" << TAU << '\n';
  descr << "z_shift=" << z_shift << '\n';
  descr << "r_0=" << r_0 << '\n';
  descr << "ampW=" << ampW_init << '\n';
  if(ampW2 > 0)
		descr << "ampW2=" << ampW2_init << '\n';
  else
		descr << "ratio2w=" << ratio2w << '\n';
  if(ampW3 > 0)
		descr << "ampW3=" << ampW3_init << '\n';
  descr << "phase2w=" << phase2w << '\n';
  descr << "curvativity2w=" << curvativity2w << '\n';
  descr << "airDensity=" << airDensity << '\n';
  descr << "tunExp=" << tunExp << '\n';
  descr << "tunFac=" << tunFac << '\n';
  descr << "Kerr=" << Kerr << '\n';
  descr << "P_init[GW]=" << P * 1.e-9 << '\n';
  descr << "I_init[W/m^2]=" << I << '\n';
  descr << "I_lens[W/m^2]=" << I_lens << '\n';
  if(townes_mode) {
		descr << "townes_mode=1\n";
		descr << "townes_mode_factor=" << townes_mode_factor << '\n';
	} else
		descr << "townes_mode=0\n";

  descr.close();
  
	std::ofstream ofiterations(folder_name+"/iterations.txt");
	ofiterations.close();
  
	Nt=Nt+2;
	return 0;
}

#ifdef DRAW
	//#include <cutil_gl_inline.h>
	//#include <cutil_gl_error.h>
	//#include <cuda_gl_interop.h>
	#include "simpledraw2D.h"
#endif


std::string number(int i) {
  std::stringstream ss;
  ss << std::setw(7) << std::setfill('0') << i;
  return ss.str();
}

//using namespace std;
bool askforloadmore=false;
bool esc=false;
int adr;

float norma(float x, float y, float z)
{return sqrt(x*x+y*y+z*z);}

int main(int argc, char **argv)
{
	printf("Hello world\n");
	int GPU2use=0;
/*
	if(argc>1) sscanf(argv[1],"%d",&GPU2use);
	if(argc>2) sscanf(argv[2],"%lf",&ratio2w);
	if(argc>3) sscanf(argv[3],"%lf",&phase2w);
	if(argc>4) sscanf(argv[4],"%lf",&curvativity2w);
*/
  std::string prefix="";
  std::string fullname="";
  
  dz=-0.5;

  int draw_interval = 333;

  folder_name="harm_";
	#define PARSE(ARG) if(name == #ARG) { sval >> ARG; folder_name = folder_name + #ARG + "=" + sval.str() + "_"; continue;}
	#define PARSE_SIMPLE(ARG) if(name == #ARG) { sval >> ARG; continue;}
	for(int i=1;i<argc;i++)
	{
		std::string inp = std::string(argv[i]);
		std::cout << "---" << inp << '\n';
		size_t pos = inp.find("=");
		if(pos == std::string::npos)
			printf("you specified parameter wrong way use <name>=<value> format. NOTE: no \'-\' and spaces\n");
		else
		{
			std::string name = inp.substr(0,pos);
			std::stringstream sval;
			sval << inp.substr(pos+1,std::string::npos);
			printf("parameter[%d] has name %s and value %s\n",i-1,name.c_str(), sval.str().c_str());
			PARSE_SIMPLE(prefix)
			PARSE_SIMPLE(fullname)
			PARSE_SIMPLE(GPU2use)
			PARSE(TAU)
			PARSE(ampW0)
			PARSE(ampW2)
			PARSE(ampW3)
			PARSE(ratio2w)
			PARSE(curvativity2w)
			PARSE(phase2w)
			PARSE(Kerr)
			PARSE(airDensity)
			PARSE(Nx) 
			PARSE(Nt) 
			PARSE(Lx) 
			PARSE(T) 
			PARSE(dz)
      PARSE(tunExp)
      PARSE(tunFac)
      PARSE_SIMPLE(save_interval)     
      PARSE_SIMPLE(dynamic_save_interval)     
      PARSE_SIMPLE(draw_interval)
      PARSE(townes_mode)     
      PARSE(townes_mode_factor)
      PARSE(collimated)
      PARSE_SIMPLE(skipPics)
      PARSE_SIMPLE(skip1D)
      PARSE_SIMPLE(skipEnergies)
      PARSE_SIMPLE(skipSources)
      PARSE_SIMPLE(saveW3envelope)
      PARSE_SIMPLE(saveFieldEnvelope)
      PARSE_SIMPLE(singleStepInFocus)
    }
  }
  
  if(singleStepInFocus) {
		save_interval=1;
	}
  
  if(!prefix.empty())
    folder_name = prefix+'_'+folder_name;
  
  if(!fullname.empty())
    folder_name = fullname;
  
  std::cout << "output folder is " << folder_name << '\n';
  std::filesystem::create_directory(folder_name);
  std::filesystem::create_directory(folder_name+"/pics");
  std::filesystem::create_directory(folder_name+"/pics/jpg");
  std::filesystem::create_directory(folder_name+"/pics/res");
  
  int Ntau=Nt+2;
  FL_DBL dx=Lx/FL_DBL(Nx);
  FL_DBL dtau=T/FL_DBL(Ntau-2);
  FL_DBL* E_THZ;
  FL_DBL* plasma;
  FL_DBL* Ew_THZ;
  FL_DBL* E_THZs=new FL_DBL[Nx*Ntau/SCR/SCR];
  int i,j;
  FL_DBL x,tau;
  FL_DBL z=0.0;
  int compute=0;
  
  printf("dz/dx/dx=%f\n",dz/dx/dx);
  
	//printf("you've chosen %d GPU and second harmonic ratio=%g; second harmonic phase=%g; second harmonic curvativity=%g\n",GPU2use,ratio2w,phase2w,curvativity2w);
	cudaSetDevice(GPU2use);
	E_THZ=new FL_DBL[Ntau*Nx];
	Ew_THZ=new FL_DBL[Ntau*Nx];
	plasma=new FL_DBL[Ntau*Nx];
	printf("INIT DATA...\n");
	setinitdata_harm(E_THZ,Ntau,Nx,dtau,dx);
	printf("INIT DATA COMPLETE\n");
	cuda_load(E_THZ, Ntau, Nx, dtau, dx, dz, airDensity);

#ifdef DRAW
	if(draw_interval)
		fadey_init((Ntau-2)/SCR,Nx/SCR,2);
#endif
	get_field(Ew_THZ,Ntau,Nx);
	//get_spec(Ew_THZ,Ntau,Nx);
	for(i=0;i<Nx/SCR;i++)
	{
		for(j=0;j<(Ntau-2)/SCR;j++) E_THZs[j+i*(Ntau-2)/SCR]=Ew_THZ[j*SCR+i*(Ntau-2)*SCR];
	}
	#ifdef DRAW
	if(draw_interval)
		fadey_draw(E_THZs,(Ntau-2)/SCR,Nx/SCR,0);
	#endif
int flnum=-1;

#ifdef SAVE
	FILE* fp;
  std::string filename=folder_name+"/pics/res/field"+number(flnum)+".dat";
	fp=fopen(filename.c_str(),"w");
	fprintf(fp,"head\n2d\t1\n%d\t%d\ndata\n",(Ntau-2)/SCR,Nx/SCR);
	for(i=0;i<Nx/SCR;i++)
	{
		for(j=0;j<(Ntau-2)/SCR;j++) fprintf(fp,"%g\n",E_THZs[j+i*(Ntau-2)/SCR]);
	}
	fclose(fp);

	FL_DBL* dens=new FL_DBL[Nx/SCR];
	sprintf(filename,"pics/res/plasma.dat");
	fp=fopen(filename,"w");
	fprintf(fp,"head\n2d\t1\n%d\tx\ndata\n",Nx/SCR);
	fclose(fp);
#endif

  double e0 = .0;
  int iter=0;
	
	flnum = iter;
  
  busyQueue q;
      
	while(-z<zmax)
	{
		int current_save_interval = save_interval;
		if(dynamic_save_interval && std::abs(std::abs(z)-std::abs(zmax*.49)) < .07 * std::abs(zmax))
			current_save_interval = save_interval / 10;
		//std::cout << "in loop z=" << z << " zmax=" << zmax << '\n';
    FL_DBL stptime;

    stptime=step(Ntau,Nx,z,T,(iter-1) % current_save_interval == 0, airDensity, Kerr, tunExp, tunFac);
    z+=dz;
    iter++;
    //std::cout << "iter=" << iter << " current_save_interval=" << current_save_interval << " criteria=" << (iter % current_save_interval) << '\n';
    if((iter % current_save_interval) && (!draw_interval || (iter % draw_interval)))
      continue;
		
		cool_dev(GPU2use);
		SystemTemperatures::coolDown();
		
    get_field(Ew_THZ, Ntau, Nx);
    get_density(plasma, Ntau, Nx);
		//get_spec(Ew_THZ,Ntau,Nx);
		for(i=0;i<Nx/SCR;i++)
		{
			for(j=0;j<(Ntau-2)/SCR;j++)
        E_THZs[j+i*(Ntau-2)/SCR]=Ew_THZ[j*SCR+i*(Ntau-2)*SCR];
		}
		#ifdef DRAW
		if(draw_interval) {
			fadey_draw(E_THZs,(Ntau-2)/SCR,Nx/SCR,0);
			fadey_draw(plasma,(Ntau-2),Nx,1);
		}
		#endif
    FL_DBL energy = .0;
    FL_DBL norm = .0;
		for(j=0;j<Nx;j++)
		{
      FL_DBL r = j==0 ? 3./8. : ((FL_DBL) j + .5); // should be proper metric
      //FL_DBL r = ((FL_DBL) j + .5);
			for(i=0;i<Ntau-2;i++)
      {
        norm += Ew_THZ[j*(Ntau-2)+i] * 2. * M_PI * r * dx * dtau;
        energy += pow(Ew_THZ[j*(Ntau-2)+i], 2.0) * 2. * M_PI * r * dx * dtau;
      }
    }
    
    e0 = e0 == .0 ? energy : e0;
      
    printf("z=%f/%f, norm=%g, E0=%g, E=%g, GPUtime_per_stp=%f(ms)\n", z, zmax, norm, e0, energy, stptime);
    
    if(iter % current_save_interval)
      continue;

		std::cout << "SAVING FILE at iter " << iter << " z=" << z << " zmax=" << zmax << '\n';

#ifdef SAVE
    FILE* fp;
    char filename[100];		sprintf(filename,"pics/res/field%d.dat",flnum);
		fp=fopen(filename,"w");
		fprintf(fp,"head\n2d\t1\n%d\t%d\ndata\n",(Ntau-2)/SCR,Nx/SCR/2);
		for(i=0;i<Nx/SCR/2;i++)
		{
			for(j=0;j<(Ntau-2)/SCR;j++) fprintf(fp,"%g\n",E_THZs[j+i*(Ntau-2)/SCR]);
		}
		fclose(fp);
    
#endif
		std::ofstream ofiterations(folder_name+"/iterations.txt", std::ios_base::app);
		ofiterations << flnum << '\n';
		ofiterations.close();
    static std::vector<postprocess*> pp;
    int slot = q.busyFindSlot();
    std::cout << "got slot " << slot << '\n';
    while(pp.size() <= slot)
      pp.push_back(new postprocess(Ntau-2, Nx, T, Lx));
    pp[slot]->drawPics = !skipPics;
    pp[slot]->save1D = !skip1D;
    pp[slot]->saveEnergy = !skipEnergies;
    pp[slot]->saveSource = !skipSources;
    pp[slot]->saveW3envelope = saveW3envelope;
    pp[slot]->saveFieldEnvelope = saveFieldEnvelope;
    pp[slot]->pictureFilename = folder_name+"/pics/jpg/field"+number(flnum)+".jpg";
    pp[slot]->signalFilename = folder_name+"/pics/res/field_axis"+number(flnum)+".dat";
    pp[slot]->plasmaProfileFilename = folder_name+"/pics/res/plasma"+number(flnum)+".dat";
    pp[slot]->energyFilename = folder_name+"/pics/res/field_energy"+number(flnum)+".dat";
    pp[slot]->sourcePrefix = folder_name+"/pics/res/src"+number(flnum);
    pp[slot]->Kerr = Kerr;
    pp[slot]->w0 = 1.0;
    pp[slot]->dimz=2./k0;
    pp[slot]->dimx=1./std::sqrt(2.)/k0;
    pp[slot]->dimt=1./k0/3./std::pow(10.,10.);
    
    for(j=0;j<Nx;j++) {
			for(i=0;i<Ntau-2;i++) {
        pp[slot]->field[j*(Ntau-2)+i] = Ew_THZ[j*(Ntau-2)+i];
				pp[slot]->plasma[j*(Ntau-2)+i] = plasma[j*(Ntau-2)+i];
			}
    }
    
    pp[slot]->reset();

    q.addTask(slot, [](void* args, int id) {
              postprocess* pp = (postprocess*) args;
              if(pp->drawPics)
                pp->plot();
              if(pp->save1D) {
                pp->saveField();
                pp->savePlasma();
              }
              if(pp->saveEnergy)
                pp->saveEnergies();
              pp->saveSources();
              }, (void*) pp[slot]);

		flnum=iter;

	}
	#ifdef DRAW
	if(draw_interval)
		fadey_close();
	#endif
	delete[] E_THZ;
	delete[] Ew_THZ;
  
  while(!q.empty())
  {
    std::cout << "waiting for postprocess workers\n";
    std::this_thread::sleep_for(std::chrono::milliseconds(1333));
	}
  
  return 0;
}
