#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <mgl2/mgl.h>
#define FL_DBL double
#include "fft_real_3.cpp"
//#define NAN 0
int main(int argc, char* const argv[])
{
	int first_file;
	int last_file;
	char filename[100];
	printf("start args is %d\n",argc);
	if (argc==4)
	{
	    sscanf(argv[1],"%d",&first_file);
	    sscanf(argv[2],"%d",&last_file);
	    sscanf(argv[3],"%s",&filename);
	}
	else
	{
	    first_file=0;
	    last_file=0;
		sprintf(filename,"field");
	}
	printf("first file is %s%d, last is %s%d\n",filename,first_file,filename,last_file);
	int Nx=256*4;
	int Ny=128*4;
	int stepx=1;
	int stepy=1;
	int i,j,k;
	double* value=new double[Nx*Ny];
	char in_fl_name[100];
	char out_fl_name[100];
	HMDT a = mgl_create_data_size(Nx/stepx,Ny/stepy,1);
	FILE *fp;
	for(k=first_file;k<last_file+1;k++)
{
	sprintf(in_fl_name,"res/%s%d.dat",filename,k);
	sprintf(out_fl_name,"jpg/%s%03d.jpg",filename,k);
	printf("opening %s\n",in_fl_name);
	fp=fopen(in_fl_name,"r");
	int ifEOF;
	float val;
	i=0;
	char char_a;
	ifEOF=1;
	double MAX=0.0;
	double min=0.0;
	while(ifEOF!=0)
	{
	    i++;
	    if(ifEOF!=0) ifEOF=fscanf(fp,"%f\n",&val);
	    if(i<Nx*Ny) value[i]=double(val);
	    else ifEOF=0;
	}
	fclose(fp);
	printf("file read ok, min=%g, MAX=%g\n",min,MAX);
	double* sig=new double[Nx];
	double* tmpsig=new double[Nx];
	for(int i=0;i<Ny;i++)
	{
		for(int j=0;j<Nx;j++) sig[j]=value[i*Nx+j];
		qqFFT_freal(Nx,sig,tmpsig);		
		for(int j=0;j<Nx/2;j++)
		{
			double freq=double(j)/10.0;
			sig[j*2]=sig[j*2]*exp(-freq*freq);
			sig[j*2+1]=sig[j*2+1]*exp(-freq*freq);
		}
		qqFFT_freal_1(Nx,sig,tmpsig);
		for(int j=0;j<Nx;j++)
		{
			value[i*Nx+j]=sig[j];
			if(sig[j]>MAX) MAX=sig[j]; 
			if(sig[j]<min) min=sig[j];
		}
	}
	printf("applied filter\n");
	MAX=(-min>MAX)?-min:MAX;
	for(i=0;i<Ny/stepy;i++)
	{
//	    for(j=0;j<Nx;j++) mgl_data_put_val(a,exp(-double((i-50)*(i-50)+(j-50)*(j-50))/50.0/50.0),i,j,0);
	    for(j=0;j<Nx/stepx;j++) mgl_data_put_val(a,value[i*Nx*stepy+j*stepx]/MAX,j,i,0);
	}

	HMGL gr = mgl_create_graph(Nx/stepx*2, Ny/stepy*2);

//	mgl_set_alpha(gr, true);
//	mgl_set_light(gr, true);
	mgl_rotate(gr,0,0,0);
	mgl_surf(gr,a,"BbcyrR","#");
	printf("saving %s\n",in_fl_name);
	mgl_write_png(gr,out_fl_name,0);
	mgl_delete_graph(gr);
}
	delete(value);
	return 0;
}
