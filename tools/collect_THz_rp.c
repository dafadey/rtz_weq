// collects data from single run
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <mgl2/mgl.h>
#include <omp.h>
#include <sys/types.h>
#include <dirent.h>
#define FL_DBL double
#include "../src/fft_real_3.cpp"
#define THS 32 //threads num

#define PATH "../"
#define DIR_PREFIX "test"
#define DATA_PREFIX "plasma"
#define DATA_PREFIX_SKIP "plasma.dat"
int	Nx=1024;
#define	Ny	512
#define	stepx	1
#define	stepy	1

#define TEF 8 // TIME_EXPANSION_FACTOR used to create longer signals then calculated source on fs pulse duration

int smaller(int x, int y) {return (x<y)?x:y;}

bool find_prefix(const char * str, const char * prefix )
{
	//finds prefix string in the begining of str
	bool cmp=true;
	int i=0;
	while(cmp==true)
	{
		cmp=(str[i]!='\0')&(prefix[i]!='\0')&(str[i]==prefix[i]);
		if(prefix[i]=='\0') return true;
		i++;
	}
	return false;
}

int main(int argc, char* const argv[])
{
	char samp[10];
	char sphase[10];
	char scurv[10];
	char resfname[100];
	char picfname[100];
	char appnd=true;
	if (argc>=4)
	{
		  sscanf(argv[1],"%s",&samp);
		  sscanf(argv[2],"%s",&sphase);
		  sscanf(argv[3],"%s",&scurv);
			sprintf(resfname,"src.dat");
			sprintf(picfname,"src.png");
			if(argc==5)
			{
				char resprefix[100];
				sscanf(argv[4],"%s",&resprefix);
				sprintf(resfname,"%s_%s_%s_%s.dat",resprefix,samp,sphase,scurv);
				sprintf(picfname,"%s_%s_%s_%s.png",resprefix,samp,sphase,scurv);
				appnd=false;
			}
			else
			{
				printf("more than 4 parameters is not supported yet. please specify configuration: <second harmonic amp> <second harmonic phase> <second harmonc curvativity> (example <exec> 0.1 1.57 1.0)\n");
				return 0;
			}
	}
	else
	{
		printf("less than 3 parameters is not supported. please specify configuration: <second harmonic amp> <second harmonic phase> <second harmonc curvativity> (example <exec> 0.1 1.57 1.0)\n");
		return 0;
	}
	//generate directory name:
	char path[100];
	//test_w2amp=0.1_w2pahse=0.0_w2curv=1.0
	sprintf(path,"%s%s_w2amp=%s_w2pahse=%s_w2curv=%s/pics/res/",PATH,DIR_PREFIX,samp,sphase,scurv);
	printf("trying path: %s\n",path);
	//detect number of files
	DIR *dp;
  int Nfl=0;
  struct dirent *ep;     
  dp = opendir (path);
  if (dp != NULL)
  {
    while (ep = readdir (dp))
		{
			if((find_prefix(ep->d_name,DATA_PREFIX))&&(!find_prefix(ep->d_name,DATA_PREFIX_SKIP))) Nfl++;
		}
    (void) closedir (dp);
  }
  else	perror ("Couldn't open the directory");
	printf("found %d files with prefix \"%s\"\n",Nfl,DATA_PREFIX);
	
	double norms[Nfl];
	double** src;
	src=new double*[Nfl];
	#pragma omp parallel for num_threads(THS)
	for(int i=0;i<Nfl;i++)
	{
		char filename[100];
		sprintf(filename,"%s%s%d.dat",path,DATA_PREFIX,i);
		printf("%d-th trhread is opening %s\n",omp_get_thread_num(),filename);
		FILE* fp;
		fp=fopen(filename,"r");
		int ifEOF;
		float val;
		src[i]=new double[Nx*TEF];
		for(int j=0;j<Nx;j++)
		{
			fscanf(fp,"%f\n",&val);
			src[i][j]=double(val);
		}
		for(int j=Nx;j<Nx*TEF;j++)	src[i][j]=0.0;
		fclose(fp);
	}


Nx=Nx*TEF;
printf("Nx now eq%d\n",Nx);
printf("test: %d\n",(-3)%10+10);
printf("test: %d\n",(-13)%10+10);
printf("test: %d\n",(-23)%10+10);
double MAX=-1000000000.0;
double min=100000000.0;

double L=2.0;
double dz=L/double(Nfl);

for(int i=0;i<Nfl;i++)
{
	src[i][0]=0.0;
	for(int j=1;j<Nx;j++) src[i][j]=src[i][j-1]*exp(-3.0/double(Nx))+src[i][j];
}


for(int i=0;i<Nfl;i++)
{
	double* srctmp=new double[Nx];
	for(int j=1;j<Nx-1;j++) srctmp[j]=-src[i][j]*2.0+src[i][j-1]+src[i][j+1];
	src[i][0]=0.0;
	src[i][Nx-1]=0.0;
	for(int j=1;j<Nx-1;j++) src[i][j]=srctmp[j];
}

/*
for(int i=0;i<Nfl;i++)
{
	double* srctmp=new double[Nx];
	for(int j=1;j<Nx-1;j++) srctmp[j]=-src[i][j]+src[i][j+1];
	src[i][0]=0.0;
	src[i][Nx-1]=0.0;
	for(int j=1;j<Nx-1;j++) src[i][j]=srctmp[j];
}
*/

char filename[100];
sprintf(filename,"%s%s",PATH,resfname);
FILE* fp=fopen(filename,(appnd==true)?"a":"w");
for(int i=0;i<Nfl;i++)
{
	for(int j=0;j<Nx;j++)
	{
		double v=src[i][j];
		fprintf(fp,"%g\n",v);
		MAX=(MAX>v)?MAX:v;
		min=(min<v)?min:v;
	}
}
fclose(fp);
MAX=(MAX>fabs(min))?MAX:fabs(min);

printf("src is ready\nintegrating radiation pattern...\n");

int Ntheta=Nfl;
double** rp=new double*[Ntheta];
double thetamax=0.1;
double dtheta=thetamax/double(Ntheta);
double dt=0.00819/double(Nx); // 0.00819 is a ful interval in cm 
for(int i=0;i<Ntheta;i++)
{
	rp[i]=new double[Nx];
	double theta=dtheta*double(i);
	for(int k=0;k<Nx;k++)
	{
		double v=0.0;
		for(int j=0;j<Nfl;j++)
		{
			int indx=floor((double(k)*dt-(1.0-cos(theta))*dz*double(j-Nfl/3))/dt);
/*			indx=(indx>0)?indx:0;
			indx=(indx<Nx-1)?indx:Nx-1;
*/
			if(indx>Nx-1) indx=indx%Nx;
			if(indx<0) indx=Nx+indx%Nx;
			v+=src[j][(indx>0)?indx:0];
		}
		rp[i][k]=v;
	}
}
// let Ntheta=Nfl !
printf("data is ready\napplying filter...\n");
//apply fft filter

fp=fopen("rp_new30w.dat","w");
double* tmp=new double[Nx];
for(int i=0;i<Ntheta;i++)
{
	qqFFT_freal(Nx,rp[i],tmp);
	//filter
	for(int j=0;j<30;j++) fprintf(fp,"%g\n",rp[i][j*2]*rp[i][j*2]+rp[i][j*2+1]*rp[i][j*2+1]);
	for(int j=0;j<Nx/2;j++)
	{
			double f=exp(-j*j);
			rp[i][j*2]*=f;
			rp[i][j*2+1]*=f;
	}
	qqFFT_freal_1(Nx,rp[i],tmp);
	double res=0.0;
	for(int j=0;j<Nx;j++) res+=rp[i][j]*rp[i][j];
	printf("W[%d]=%g\n",i,res);
}
fclose(fp);


/*
for(int i=0;i<Ntheta;i++)
{
	int sig_start=-1;
	int sig_end=-1;
	//filter
	double MAX=0.0;
	for(int j=0;j<Nx;j++) MAX=(MAX>fabs(rp[i][j]))?MAX:fabs(rp[i][j]);
	double ths_up=0.5*MAX;
	double ths_down=0.1*MAX;
	for(int j=0;j<Nx;j++)
	{
		sig_start=((sig_start==-1)&(fabs(rp[i][j])>ths_up))?j:sig_start;
		sig_end=(fabs(rp[i][j])>ths_down)?j:sig_end;
	}
	//printf("sig_start=%d, sig_end=%d\n",sig_start,sig_end);
	if((sig_start!=-1)&(sig_end!=-1))
	{
		//for(int j=0;j<Nx;j++) rp[i][j]=((j-sig_start)*(sig_end-j)>0)?1.0:0.0;
		//for(int j=0;j<smaller(sig_end+1024,Nx-1);j++) rp[i][j]=0.0;
		//for(int j=Nx-512;j<Nx;j++) rp[i][j]=0.0;
	}
}
*/


MAX=0.0;
min=0.0;
for(int i=0;i<Ntheta;i++)
{
	for(int k=0;k<Nx;k++)
	{
		double v=rp[i][k];
		MAX=(MAX>v)?MAX:v;
		min=(min<v)?min:v;
	}
}
MAX=(MAX>fabs(min))?MAX:fabs(min);

printf("filtering is ready\nplotting...(MAX=%g)\n",MAX);



HMDT a = mgl_create_data_size(Nx,Nfl,1);

for(int i=0;i<Nfl;i++)
{
	for(int j=0;j<Nx;j++)	mgl_data_put_val(a,(MAX>0)?rp[i][j]/MAX:0.0,j,i,0);
}

HMGL gr = mgl_create_graph(Nx, Nfl);
mgl_rotate(gr,0,0,0);
mgl_surf(gr,a,"BbcyrR","#");
mgl_write_png(gr,picfname,0);
mgl_delete_graph(gr);
printf("picture written size is %d X %d\n",Nfl,Nx);

	return 0;
}
