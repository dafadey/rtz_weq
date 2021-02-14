#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cutil_inline.h>
#include <cufft.h>
#include <math_constants.h>
#include "dfs.h"
#include "fft_real_3.cpp"
double phase2w=1.57;
double ratio2w=0.0; //ratio of second harmonic to first harmonic (by amplitude)
double curvativity2w=1.0;
#include "thz_funcs.cpp"

#ifdef DRAW
	#include <cutil_gl_inline.h>
	#include <cutil_gl_error.h>
	#include <cuda_gl_interop.h>
	#include "simpledraw2D.h"
#endif

#define	DT	0.1			//	dt 
#define	STP	5			//	steps calculated on GPU before copy data 2 host and visualisation


//using namespace std;
bool askforloadmore=false;
bool esc=false;
int adr;

extern "C" 
void cuda_load(FL_DBL*, int, int, FL_DBL, FL_DBL, FL_DBL);
extern "C" 
void get_spec(FL_DBL*, int, int);
extern "C" 
void get_field(FL_DBL*, int, int);
extern "C" 
void get_density(FL_DBL*, int, int);
extern "C" 
void get_THz_source(FL_DBL*, int, int);
extern "C" 
FL_DBL step(int, int, FL_DBL, bool);

/*
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT-INIT
*/
int Ntau=NNT;
int Nx=NNX;
FL_DBL dx=Lx/FL_DBL(Nx);
FL_DBL dtau=T/FL_DBL(Ntau-2);
FL_DBL dz=-0.5;
//FL_DBL dz=0.0000001;
FL_DBL* E_THZ;
FL_DBL* plasma;
FL_DBL* Ew_THZ;
FL_DBL* E_THZs=new FL_DBL[Nx*Ntau/SCR/SCR];
int i,j;
FL_DBL x,tau;
FL_DBL z=0.0;
int compute=0;

float norma(float x, float y, float z)
{return sqrt(x*x+y*y+z*z);}

int main(int argc, char **argv)
{
//Generate particles
	printf("Hello world\n");
	printf("dz/dx/dx=%f\n",dz/dx/dx);
	int GPU2use=0;
	if(argc>1) sscanf(argv[1],"%d",&GPU2use);
	if(argc>2) sscanf(argv[2],"%lf",&ratio2w);
	if(argc>3) sscanf(argv[3],"%lf",&phase2w);
	if(argc>4) sscanf(argv[4],"%lf",&curvativity2w);
	printf("you've chosen %d GPU and second harmonic ratio=%g; second harmonic phase=%g; second harmonic curvativity=%g\n",GPU2use,ratio2w,phase2w,curvativity2w);
	cudaSetDevice(GPU2use);
	E_THZ=new FL_DBL[Ntau*Nx];
	Ew_THZ=new FL_DBL[Ntau*Nx];
	plasma=new FL_DBL[Ntau*Nx];
	printf("INIT DATA...\n");
	setinitdata(E_THZ,Ntau,Nx,dtau,dx);
	printf("INIT DATA COMPLETE\n");
//fadey DEBUG
//	find_spec(E_THZ, Ew_THZ,Ntau,Nx);

//	cuda_load(Ew_THZ, Ntau, Nx, dtau, dx, dz);
	cuda_load(E_THZ, Ntau, Nx, dtau, dx, dz);

#ifdef DRAW
	fadey_init((Ntau-2)/SCR,Nx/SCR,1);
#endif
	get_field(Ew_THZ,Ntau,Nx);
	//get_spec(Ew_THZ,Ntau,Nx);
	for(i=0;i<Nx/SCR;i++)
	{
		for(j=0;j<(Ntau-2)/SCR;j++) E_THZs[j+i*(Ntau-2)/SCR]=Ew_THZ[j*SCR+i*(Ntau-2)*SCR];
	}
	#ifdef DRAW
	fadey_draw(E_THZs,(Ntau-2)/SCR,Nx/SCR,0);
	#endif
int flnum=-1;

#ifdef SAVE
	FILE* fp;
	char* filename=new char[100];
	sprintf(filename,"pics/res/field%d.dat",flnum);
	fp=fopen(filename,"w");
	fprintf(fp,"head\n2d\t1\n%d\t%d\ndata\n",(Ntau-2)/SCR,Nx/SCR);
	//fprintf(fp,"%d\t%d\n",NNX/SCR,(NNT-2)/SCR);
	for(i=0;i<Nx/SCR;i++)
	{
		for(j=0;j<(Ntau-2)/SCR;j++) fprintf(fp,"%g\n",E_THZs[j+i*(Ntau-2)/SCR]);
	}
	fclose(fp);
	flnum++;

	FL_DBL* dens=new FL_DBL[Nx/SCR];
	sprintf(filename,"pics/res/plasma.dat");
	fp=fopen(filename,"w");
	fprintf(fp,"head\n2d\t1\n%d\tx\ndata\n",Nx/SCR);
	fclose(fp);
#endif

double e0=0;
	while(-z<60000)
	{
	    FL_DBL stptime;
	    //calc energy;
	    double energy=0;
		
		for(i=0;i<REP-1;i++) stptime=step(Ntau,Nx,z,false);
		stptime=step(Ntau,Nx,z,true);

		get_field(Ew_THZ,Ntau,Nx);
		//get_spec(Ew_THZ,Ntau,Nx);
		for(i=0;i<Nx/SCR;i++)
		{
//			for(j=0;j<Ntau/SCR;j++) E_THZs[j+i*Ntau/SCR]=Ew_THZ[j*SCR+i*Ntau*SCR];
			for(j=0;j<(Ntau-2)/SCR;j++) E_THZs[j+i*(Ntau-2)/SCR]=Ew_THZ[j*SCR+i*(Ntau-2)*SCR];
		}
		#ifdef DRAW			
		fadey_draw(E_THZs,(Ntau-2)/SCR,Nx/SCR,0);
		#endif
		FL_DBL norm=0.0;
		for(j=0;j<Nx;j++)
		{			
			for(i=0;i<Ntau-2;i++)	norm+=Ew_THZ[j*(Ntau-2)+i]*Ew_THZ[j*(Ntau-2)+i];
		}

		printf("norm=%g\n",double(norm));

	    printf("z=%f, E0=%g, E=%g, GPUtime_per_stp=%f(ms)\n", z, e0, energy, stptime);
	    z = z + dz*FL_DBL(REP);// cutGetTimerValue(timer);

#ifdef SAVE
		sprintf(filename,"pics/res/field%d.dat",flnum);
		fp=fopen(filename,"w");
		fprintf(fp,"head\n2d\t1\n%d\t%d\ndata\n",(Ntau-2)/SCR,Nx/SCR);
		//fprintf(fp,"%d\t%d\n",NNX/SCR,(NNT-2)/SCR);
		for(i=0;i<Nx/SCR;i++)
		{
			for(j=0;j<(Ntau-2)/SCR;j++) fprintf(fp,"%g\n",E_THZs[j+i*(Ntau-2)/SCR]);
		}
		fclose(fp);

		get_THz_source(plasma,Ntau,Nx);
		sprintf(filename,"pics/res/source%d.dat",flnum);
		fp=fopen(filename,"w");
		fprintf(fp,"head\n2d\t1\n%d\t%d\ndata\n",(Ntau-2)/SCR,Nx/SCR);
		//fprintf(fp,"%d\t%d\n",NNX/SCR,(NNT-2)/SCR);
		for(i=0;i<Nx/SCR;i++)
		{
			for(j=0;j<(Ntau-2)/SCR;j++) fprintf(fp,"%g\n",plasma[j+i*(Ntau-2)/SCR]);
		}
		fclose(fp);

//
		get_density(plasma,Ntau,Nx);
		for(i=0;i<Nx/SCR;i++)	dens[i]=0.0;
		for(i=0;i<Nx/SCR;i++)
		{
			for(j=0;j<(Ntau-2)/SCR;j++)
				dens[i]=(dens[i]>fabs(plasma[j+i*(Ntau-2)/SCR]))?dens[i]:fabs(plasma[j+i*(Ntau-2)/SCR]);
		}
		sprintf(filename,"pics/res/plasma.dat");
		fp=fopen(filename,"a");
		for(i=0;i<Nx/SCR;i++)
			fprintf(fp,"%g\n",dens[i]);
		fclose(fp);
		
		


		flnum++;
#endif
	}
	#ifdef DRAW
	fadey_close();
	#endif
	delete[] E_THZ;
	delete[] Ew_THZ;
	return 0;
}
