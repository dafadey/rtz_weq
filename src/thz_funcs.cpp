FL_DBL mod(FL_DBL x, FL_DBL y) {return x*x+y*y;}

int setinitdata(FL_DBL* arr, int Nt, int Nx, FL_DBL dt, FL_DBL dx)
{
/*
initializes arrays of field Nt X Nx
*/
	int i,j;
	Nt=Nt-2;
	
/*
	for(j=0;j<Nx;j++)
	{
		FL_DBL x=FL_DBL(j-0*Nx/2)*dx;
		for(i=0;i<Nt;i++)
		{
			FL_DBL t=FL_DBL(i-Nt/3)*dt;
			FL_DBL tau=t+0.0*0.3*x*x;
			arr[j*Nt+i]=exp(-x*x*0.03-(tau+x*x/43.0)*(tau+x*x/43.0)*1.5)*(sin((tau+x*x/43.0)*7.0));
		}
	}
*/
	double k0=80553.65779; //k0[cm^-1] for lambda = 780nm
	double z0=2.0/k0;

	double aaa=160.0; //focal radius fo 1cm diam beam focused into with 400mm focal length lens
	double z_shift=2.0/z0; //2cm shift from focus
	//z_step=z_shift/dz;
	double r_0=sqrt(pow(z_shift/aaa,2.0)*16.0+pow(aaa,2.0));
	double V=4.0*z_shift/(16.0*pow(z_shift,2.0)+pow(aaa,4.0));
	//double V=2.0*4.0*z_shift/(16.0*pow(z_shift,2.0)+pow(aaa,4.0));
	double tau_0=(35.0/2.66*2.0*Pi);
	for(j=0;j<Nx;j++)
	{
		double r=double(j)*dx;
		for(i=0;i<Nt;i++)
		{
			double t=double(i-Nt/2)*dt;
			arr[j*Nt+i]=0.055*exp(-pow(t/(tau_0/sqrt(2.0)),2.0))*exp(-pow(r/r_0,2.0))*(sin(t+r*r*V)+ratio2w*sin(2.0*t+r*r*V*curvativity2w*2.0+phase2w));
		}
	}
	Nt=Nt+2;
return 0;
}

FL_DBL ffMAX(FL_DBL* arr, int Nt, int Nx)
{
	int i,j;
	FL_DBL MAXarr=0.0;
	FL_DBL minarr=0.0;
	for(j=0;j<Nx;j++)
	{
		for(i=0;i<Nt;i++)
		{
			if(arr[j*Nt+i]>MAXarr) MAXarr=arr[j*Nt+i];
			if(arr[j*Nt+i]<minarr) minarr=arr[j*Nt+i];
		}
	}
	if (MAXarr<-minarr) MAXarr=-minarr;
	return MAXarr;
}

int find_spec(FL_DBL* arr, FL_DBL* spec, int Nt, int Nx)
{
	FL_DBL* tmp=new FL_DBL[Nt];
	int i,j;
	for(i=0;i<Nx;i++)
	{
		for(j=0;j<Nt;j++) tmp[j]=arr[i*Nt+j];
		qqFFT_freal(Nt, tmp, &spec[Nt*i]);
	}
	delete[] tmp;
}

int find_field(FL_DBL* spec, FL_DBL* arr, int Nt, int Nx)
{
	FL_DBL* tmp=new FL_DBL[Nt];
	int i,j;
	for(i=0;i<Nx;i++)
	{
		for(j=0;j<Nt;j++) tmp[j]=spec[i*Nt+j];
		qqFFT_freal_1(Nt, tmp, &arr[Nt*i]);
		for(j=0;j<Nt;j++) arr[Nt*i]=arr[Nt*i]/double(Nt);
	}
	delete[] tmp;
}



