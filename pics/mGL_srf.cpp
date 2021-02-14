#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <mgl2/mgl.h>
#include <vector>

#define FL_DBL double

#include "../src/fft_real_3.cpp"
//#define NAN 0


double* filter(double* in, int nx, int ny, int freq, int gap)
{
	double* out=new double[nx*ny];
	double* tmp=new double[nx];
	double* sig=new double[nx];
	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
			sig[i]=in[i+j*nx];
		qqFFT_freal(nx,sig,tmp);

		sig[0]=0;		
		sig[1]=0;		
		for(int i=0;i<nx/2;i++)
		{
			double f=double(i-freq);
			double w=double(gap);
			sig[i*2]*=exp(-pow(f/w,2));
			sig[i*2+1]*=exp(-pow(f/w,2));
		}

		qqFFT_freal_1(nx,sig,tmp);
		for(int i=0;i<nx;i++)
			out[i+j*nx]=sig[i];
	}
	delete[] tmp;
	delete[] sig;
	return out;
}


int main(int argc, char* const argv[])
{
	int first_file;
	int last_file;
	char filename[100];
	double globalMAX=0.0;
	printf("start args is %d\n",argc);
	if (argc==2)
	{
	    first_file=0;
	    last_file=0;
	    sscanf(argv[1],"%s",(char*) &filename);
	}
	else if (argc==4)
	{
	    sscanf(argv[2],"%d",&first_file);
	    sscanf(argv[3],"%d",&last_file);
	    sscanf(argv[1],"%s",(char*) &filename);
	}
	else if (argc==5) // with global max
	{
	    sscanf(argv[2],"%d",&first_file);
	    sscanf(argv[3],"%d",&last_file);
		float gm;
	    sscanf(argv[4],"%f",&gm);
		globalMAX=double(gm);
	    sscanf(argv[1],"%s",(char*) &filename);
	}
	else
	{
		printf("something wrong with input parameters using \"field 0 0\"\n");
	    first_file=0;
	    last_file=0;
		sprintf(filename,"field");
	}

	printf("first file is %s%d, last is %s%d\n",filename,first_file,filename,last_file);
	int i,j,k;
	#define NFILES 300
	double* Int=new double[512*NFILES];
	for(i=0;i<512*NFILES;i++) Int[i]=0.0;
	for(k=first_file;k<last_file+1;k++)
	{
		FILE *fp;
		char in_fl_name[100];
		char out_fl_name[100];
		if(strcmp(filename,"plasma")==0)
		{
			sprintf(in_fl_name,"res/%s.dat",filename);
			sprintf(out_fl_name,"jpg/%s.png",filename);
		}
		else if(strcmp(filename,"source")==0)
		{
			sprintf(in_fl_name,"res/plasma%d.dat",k);
			sprintf(out_fl_name,"jpg/plasma%03d.png",k);
		}
		else if(strcmp(filename,"field")==0)
		{
			sprintf(in_fl_name,"res/%s%d.dat",filename,k);
			sprintf(out_fl_name,"jpg/%s%03d.png",filename,k);
		}
		else
		{
			printf("%s is not recognized. please plot plasma/source/field\n",filename);
			return -1;
		}
		printf("opening %s\n",in_fl_name);
		fp=fopen(in_fl_name,"r");
		//try reading header
		int Nx,Ny;
		bool matchNy=false;		
		int ifEOF=1;
		char str[80]; 
		if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
		if(strcmp(str,"head")==0)
		{
			if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
			if(strcmp(str,"2d")!=0)
			{
				printf("sorry data is not 2d\n");
				return -1;
			}
			if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
			if(strcmp(str,"1")!=0)
			{
				printf("sorry data is not scalar\n");
				return -1;
			}
			if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
			if(strcmp(str,"x")==0)
			{
				printf("sorry can't match Nx\n");
				return -1;
			}
			else
				Nx=atoi(str);
			if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
			if(strcmp(str,"x")==0)
			{
				matchNy=true;
				printf("match Ny mode is on\n");
			}
			else
				Ny=atoi(str);
			if(ifEOF!=0) ifEOF=fscanf(fp,"%s",str);
			if(strcmp(str,"data")!=0)
			{
				printf("sorry header is broken. expected data section here\n");
				return -1;
			}
			
			printf("header read ok. Nx=%d, Ny=%d\n",Nx,Ny);
		}
		else
		{
			fseek(fp,0,SEEK_SET);
			Nx=512;//256*4;
			Ny=1024;//128*4;
			printf("no header found. using default geometry for old data style. Nx=%d, Ny=%d\n",Nx,Ny);
		}
		int stepx=1;
		int stepy=1;
		std::vector<double> value;
//		double* value=new double[Nx*Ny];
		float val;
		i=0;
		char char_a;
		double MAX=0.0;
		double min=0.0;
		while(!feof(fp))
		{
		    i++;
		   	ifEOF=fscanf(fp,"%f\n",&val);
			if(!feof(fp))
			{
				value.push_back(double(val));
				if(val>MAX) MAX=val; 
				if(val<min) min=val; 
			}
		}
		fclose(fp);
		if(matchNy==true)
			Ny=value.size()/Nx;
		printf("file read ok, min=%g, MAX=%g geom=%dx%d datasize is %ld\n",min,MAX,Nx,Ny,(long int) value.size());
		if(globalMAX!=0.0) printf("global maximum is specified (%g) and will be used for field and source plots\n",globalMAX);
		HMDT a;
		MAX=(-min>MAX)?-min:MAX;
		HMGL gr;
		double* value_w1;
		double* value_w2;
		HMDT b,c;
		if(strcmp(filename,"field")==0)
		{
			a = mgl_create_data_size(Nx/stepx,Ny/stepy,1);
			if(globalMAX!=0.0) MAX=globalMAX;
			printf("ok doing field...\n");
			printf("filtering...\n");
			value_w1=filter(&value[0],Nx,Ny,40,20);
			value_w2=filter(&value[0],Nx,Ny,80,20);
			printf("filtering is done\n");

			b = mgl_create_data_size(Nx/stepx,Ny/stepy,1);
			c = mgl_create_data_size(Nx/stepx,Ny/stepy,1);
			for(i=0;i<Ny/stepy;i++)
			{
				for(j=0;j<Nx/stepx;j++)
				{
					mgl_data_put_val(a,fabs(value_w1[i*Nx*stepy+j*stepx]),j,i,0);
					mgl_data_put_val(b,fabs(value_w2[i*Nx*stepy+j*stepx]),j,i,0);
					mgl_data_put_val(c,value[i*Nx*stepy+j*stepx],j,i,0);
				}
			}

			mgl_set_color('W',1,1,1);

			//gr = mgl_create_graph(Nx/stepx, Ny/stepy*3);
			gr = mgl_create_graph(1024, 768);
		//	mgl_set_alpha(gr, true);
		//	mgl_set_light(gr, true);
		//	mgl_rotate(gr,45,45+double(k),0);



			mgl_subplot(gr,1,2,0,"E_\\omega");
			mgl_title(gr,"E_\\omega","",6);
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',0,MAX);
			mgl_set_range_val(gr,'c',0,MAX);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,100);
			mgl_surf(gr,a,"WRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"E_\\omega",0,"");
			mgl_colorbar(gr,"WRry<");
			mgl_box(gr);		

			mgl_subplot(gr,1,2,1,"E_{2\\omega}");
			mgl_title(gr,"E_{2\\omega}","",6);
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',0,MAX);
			mgl_set_range_val(gr,'c',0,MAX);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,100);
			mgl_surf(gr,b,"WRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"E_\\omega",0,"");
			mgl_colorbar(gr,"WRry<");
			mgl_box(gr);		
/*
			mgl_subplot(gr,1,3,2,"E_{2\\omega}");
			mgl_title(gr,"E","",6);
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',-MAX,MAX);
			mgl_set_range_val(gr,'c',-MAX,MAX);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,100);
			mgl_surf(gr,c,"cbBWRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"E_\\omega",0,"");
			mgl_colorbar(gr,"cbBWRry<");
			mgl_box(gr);		
*/
		}
		else if(strcmp(filename,"source")==0)
		{
			printf("ok doing source...\n");
			a = mgl_create_data_size(Nx/stepx,Ny/stepy,1);
			b = mgl_create_data_size(NFILES,Ny,1);
			if(globalMAX!=0.0) MAX=globalMAX;
			double SMAX=0.0;
			for(i=0;i<Ny/stepy;i++)
			{
				double INT=0.0;
				for(j=0;j<Nx/stepx;j++)
				{
					mgl_data_put_val(a,value[i*Nx*stepy+j*stepx],j,i,0);
					INT+=value[i*Nx*stepy+j*stepx];
				}
				Int[i*NFILES+k]=INT;
				SMAX=(SMAX>fabs(INT))?SMAX:fabs(INT);
			}
			printf("source MAX=%g\n",MAX);
			printf("source int MAX=%g\n",SMAX);
			for(i=0;i<Ny;i++)
			{
				for(j=0;j<NFILES;j++)
					mgl_data_put_val(b,Int[i*NFILES+j],j,i,0);
			}
			mgl_set_color('W',1,1,1);

			gr = mgl_create_graph(1024, 768);
		//	mgl_set_alpha(gr, true);
		//	mgl_set_light(gr, true);
		//	mgl_rotate(gr,45,45+double(k),0);



			mgl_subplot(gr,1,2,0,"src_\\{THz}");
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',-MAX,MAX);
			mgl_set_range_val(gr,'c',-MAX,MAX);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,100);
			mgl_surf(gr,a,"cbBWRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"S_{THz}",0,"");
			mgl_colorbar(gr,"cbBWRry");
			mgl_box(gr);

			mgl_subplot(gr,1,2,1,"\\int_0^\\infty src_\\{THz} d\\tau");
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',-0.4,0.4);
			mgl_set_range_val(gr,'c',-0.4,0.4);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,3);
			mgl_surf(gr,b,"cbBWRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"z \\ [cm]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"S_{THz}",0,"");
			mgl_colorbar(gr,"cbBWRry");
			mgl_box(gr);
		}
		else if(strcmp(filename,"plasma")==0)
		{
			printf("ok doing plasma...\n");
			a = mgl_create_data_size(Ny/stepy,Nx/stepx,1);
			for(i=0;i<Ny/stepy;i++)
			{
				for(j=0;j<Nx/stepx;j++)
				{
					mgl_data_put_val(a,value[i*Nx*stepy+j*stepx],i,j,0);
				}
			}

			mgl_set_color('W',1,1,1);

			gr = mgl_create_graph(2*2*Nx/stepx, 2*2.25*Ny/stepy);
		//	mgl_set_alpha(gr, true);
		//	mgl_set_light(gr, true);
		//	mgl_rotate(gr,45,45+double(k),0);



			mgl_subplot(gr,2,3,0,"plasma density");
			//mgl_rotate(gr,63,30,0);
			mgl_rotate(gr,0,0,0);
			mgl_set_range_val(gr,'z',0,MAX);
			mgl_set_range_val(gr,'c',0,MAX);
			mgl_set_range_val(gr,'y',0,400);
			mgl_set_range_val(gr,'x',0,3);
			mgl_surf(gr,a,"WRry","");
			mgl_axis(gr,"","","");
			mgl_label(gr,'x',"z \\ [cm]",0,"");
			mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'z',"n",0,"");
			mgl_colorbar(gr,"WRry<");
			mgl_box(gr);

//add 1d plots
			b = mgl_create_data_size(Ny/stepx,2,1);
			double lumMAX=0.0;
			int lumMAXi=0.0;			
			for(i=0;i<Ny/stepy;i++)
			{
				double luminosity=0.0;
				for(j=0;j<Nx/stepy;j++)
					luminosity+=fabs(value[i*Nx*stepy+j*stepx]*double(j+0.5));
				if(luminosity>lumMAX)
				{
					lumMAX=luminosity;
					lumMAXi=i;
				}				
			}
			int densMAXi=0;
			double densMAX=0.0;
			for(i=0;i<Ny/stepy;i++)
			{
				mgl_data_put_val(b,fabs(value[i*Nx*stepy]),i,0,0);
				if(fabs(value[i*Nx*stepy])>densMAX)
				{
					densMAX=fabs(value[i*Nx*stepy]);
					densMAXi=i;
				}
				double luminosity=0.0;
				for(j=0;j<Nx/stepy;j++)
					luminosity+=fabs(value[i*Nx*stepy+j*stepx]*double(j+0.5));
				mgl_data_put_val(b,luminosity/lumMAX*MAX,i,1,0);
			}
			mgl_set_range_val(gr,'y',0,MAX);
			mgl_set_range_val(gr,'x',0,3);
			mgl_subplot(gr,2,3,1,"plasma density");
			mgl_plot(gr,b,"rb","");
			mgl_label(gr,'x',"z \\ [cm]",0,"");
			mgl_label(gr,'y',"n/n_{cr} ",0,"");
			mgl_axis(gr,"U","","");
			mgl_box(gr);
			
			c = mgl_create_data_size(Nx/stepx,2,1);
			for(j=0;j<Nx/stepx;j++)
			{
				mgl_data_put_val(c,fabs(value[densMAXi*Nx*stepy+j*stepx]),j,0,0);
				mgl_data_put_val(c,fabs(value[lumMAXi*Nx*stepy+j*stepx]),j,1,0);
			}
			mgl_set_range_val(gr,'y',0,MAX);
			mgl_set_range_val(gr,'x',0,400);
			mgl_subplot(gr,2,3,3,"plasma density");
			mgl_plot(gr,c,"rb","");
			mgl_label(gr,'x',"r \\ [\\mu m]",0,"");
			mgl_label(gr,'y',"n/n_{cr} ",0,"");
			mgl_axis(gr,"U","","");
			mgl_box(gr);
		}

		printf("saving %s\n",in_fl_name);
		mgl_write_png(gr,out_fl_name,0);
		mgl_delete_graph(gr);

		mgl_delete_data(a);
		mgl_delete_data(b);
		value.clear();
		if(strcmp(filename,"field")==0)
		{
			delete[] value_w1;
			delete[] value_w2;
			mgl_delete_data(c);
		}
		if(strcmp(filename,"plasma")==0)
		{
			mgl_delete_data(c);
		}
	}
	delete[] Int;
	return 0;
}
