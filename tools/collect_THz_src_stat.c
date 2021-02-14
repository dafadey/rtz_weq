// collects data from single run
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <mgl2/mgl.h>
#include <omp.h>
#include <sys/types.h>
#include <dirent.h>

#define THS 32 //threads num

#define PATH "../"
#define DIR_PREFIX "test"
#define DATA_PREFIX "plasma"
#define DATA_PREFIX_SKIP "plasma.dat"
#define	Nx	1024
#define	Ny	512
#define	stepx	1
#define	stepy	1


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
	//printf("got %d parameters\n",argc);
	char samp[10];
	char sphase[10];
	char scurv[10];
	char resfname[100];
	char appnd=true;
	if (argc>=4)
	{
		  sscanf(argv[1],"%s",&samp);
		  sscanf(argv[2],"%s",&sphase);
		  sscanf(argv[3],"%s",&scurv);
			sprintf(resfname,"result.dat");
			if(argc==5)
			{
				char resprefix[100];
				sscanf(argv[4],"%s",&resprefix);
				sprintf(resfname,"%s_%s_%s_%s.dat",resprefix,samp,sphase,scurv);
				appnd=false;
			}
			else if(argc>5)
			{
				printf("more than 4 parameters is not supported yet. please specify configuration: <second harmonic amp> <second harmonic phase> <second harmonc curvativity> (example <exec> 0.1 1.57 1.0)");
				return 0;
			}
	}
	else
	{
		printf("less than 3 parameters is not supported. please specify configuration: <second harmonic amp> <second harmonic phase> <second harmonc curvativity> (example <exec> 0.1 1.57 1.0)");
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
		int j=0;
		char char_a;
		ifEOF=1;
		double norm=0.0;
		while(ifEOF!=0)
		{
				j++;
				if(ifEOF!=0) ifEOF=fscanf(fp,"%f\n",&val);
				if(j<Nx*Ny) norm+=double(val)*double(j/Ny);
				else ifEOF=0;
		}
		norms[i]=norm;
		fclose(fp);
	}

char filename[100];
sprintf(filename,"%s%s",PATH,resfname);
FILE* fp=fopen(filename,(appnd==true)?"a":"w");
fprintf(fp,"%d\n",Nfl);
fprintf(fp,"%s\n",samp);
fprintf(fp,"%s\n",sphase);
fprintf(fp,"%s\n",scurv);
for(int i=0;i<Nfl;i++) fprintf(fp,"%g\n",norms[i]);
fclose(fp);

/*
	printf("first file is %s%d, last is %s%d\n",filename,first_file,filename,last_file);
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
	    if(val>MAX) MAX=val; 
	    if(val<min) min=val; 
	}
	fclose(fp);
	printf("file read ok, min=%g, MAX=%g\n",min,MAX);
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
*/
	return 0;
}
