pair cm(pair a, pair b) {return (a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);}
pair cd(pair a, pair b) {return ((a.x*b.x+a.y*b.y)/(b.x^2+b.y^2),(-a.x*b.y+a.y*b.x)/(b.x^2+b.y^2));}
real cmod(pair a) {return sqrt(a.x^2+a.y^2);}
real phase(pair a) {return atan2(a.y,a.x);}
pair csqrt(pair a)
{
  real re=sqrt((a.x+sqrt(a.x^2+a.y^2))/2.0);
	return (re,a.y/2/re);
}
pair cc(pair a)	{return (a.x,-a.y);}

void lin_solver(pair[][] A, pair[] Y, pair[] X, int N)
{
	int i,j,k,q;
	int[] mainstr= new int[N];
	int[] usedstr= new int[N];
	for(i=0;i<N;++i) usedstr[i]=0;
	int k=0;
	for(j=0;j<N;++j)
	{
		i=0;
//		write("---------------------");
//		write(A);
		while(((cmod(A[i][j])<0.00000000000000000000001)||(usedstr[i]==1))&(i<N)) ++i;
		if(i==N) write("achtung!");
		mainstr[N-1-j]=i;
		usedstr[i]=1;
		for(k=0;k<N;++k)
		{
			if(usedstr[k]==0)
			{
			// i last used string, j - number of zero create in others, k string to operate
				for(q=j+1;q<N;++q)
					A[k][q]=A[k][q]-cm(cd(A[k][j],A[i][j]),A[i][q]);
				Y[k]=Y[k]-cm(cd(A[k][j],A[i][j]),Y[i]);
				A[k][j]=0;
			}
		}	
	}
	// now we have triangled matrix A
	for(i=0;i<N;++i)
	{
	  k=N-1-i;
		int ii=mainstr[i];
		X[k]=Y[ii];
		
		for(j=0;j<i;++j)    X[k]=X[k]-cm(A[ii][N-1-j],X[N-1-j]);
		
		X[k]=cd(X[k],A[ii][k]);
	}
	mainstr.delete();
	usedstr.delete();
}

void layer_solver(pair[][] A, pair[] Et, pair[] E, int N, int lay=0)
{
//	write("layer_solver N=",N);
	int i,j;
	pair[] Y=new pair[N];
	pair[] x=new pair[N];
	pair[][] a=new pair[N][N];
	for(i=0;i<N;++i) Y[i]=-A[i][N]*Et[0]-A[i][N+1]*Et[1];
	for(i=0;i<N;++i)
	{
		for(j=0;j<N;++j)	a[i][j]=A[i][j];
	}
	lin_solver(a,Y,x,N);
	E[0]=x[0];
	E[1]=x[1];
	//test solver
/*
								pair[] xx=new pair[N+2];
								for(i=0;i<N;++i) xx[i]=x[i];
								xx[N]=Et[0];
								xx[N+1]=Et[1];
								for(i=0;i<N;++i)
								{
									pair Diff=0;
									for(j=0;j<N+2;++j)
									{
										Diff=Diff+cm(A[i][j],xx[j]);
									}
	//								write("layer_solver Diff=",cmod(Diff));
								}
								xx.delete();
*/
	x.delete();
	Y.delete();
	a.delete();
}

void 	create_A(pair[][] A, pair[] eps, real[] D, real kx, real k0, int N, int type=0)
{
	int i,j;
	pair[] kz=new pair[N];

	for(i=0;i<(N+1)*2;++i)
	{
		for(j=0;j<(N+1)*2+2;++j) A[i][j]=0;
	}

	pair kz0=csqrt((k0^2,0)-(kx^2,0)); // vacuum
	for(i=0;i<N;++i) kz[i]=cc(csqrt(cm(eps[i],(k0^2,0))-(kx^2,0)));

	A[0][0]=1;
	A[0][1]=1;
	A[1][0]=kz0;
	A[1][1]=-kz0;
	real z=0;
	pair factorT=(1,0);
	pair factorR=(1,0);
	for(i=0;i<N;++i)
	{
		
		A[i*2][(i+1)*2]=-1;
		A[i*2][(i+1)*2+1]=-1;
		A[i*2+1][(i+1)*2]=-kz[i];
		A[i*2+1][(i+1)*2+1]=kz[i];

//		z=z+D[i];
		pair arg=-cc(cm(kz[i],(D[i],0)));
		factorT=(cos(arg.x)*exp(arg.y),sin(arg.x)*exp(arg.y));
		factorR=(cos(-arg.x)*exp(-arg.y),sin(-arg.x)*exp(-arg.y));

//write(factor);

		A[(i+1)*2][(i+1)*2]=factorT;
		A[(i+1)*2][(i+1)*2+1]=factorR;
		A[(i+1)*2+1][(i+1)*2]=cm(kz[i],factorT);
		A[(i+1)*2+1][(i+1)*2+1]=-cm(kz[i],factorR);

	}
	if(type==1)
	{
		//stay in last medium
		A[i*2][(i+1)*2]=factorT;
		A[i*2][(i+1)*2+1]=factorR;
		A[i*2+1][(i+1)*2]=cm(kz[i-1],factorT);
		A[i*2+1][(i+1)*2+1]=-cm(kz[i-1],factorR);
	}
	else if(type==0)
	{
		//back 2 vacuum
		A[i*2][(i+1)*2]=-1;
		A[i*2][(i+1)*2+1]=-1;
		A[i*2+1][(i+1)*2]=-kz0;
		A[i*2+1][(i+1)*2+1]=kz0;
	}

	kz.delete();
}




void 	create_debug_A(pair[][] A, pair[] eps, real[] D, real kx, real k0)
{
	int i,j;

	for(i=0;i<(1+1)*2;++i)
	{
		for(j=0;j<(1+1)*2+2;++j) A[i][j]=0;
	}

	pair kz0=csqrt((k0^2,0)-(kx^2,0)); // vacuum
	pair kz1=cc(csqrt(cm(eps[0],(k0^2,0))-(kx^2,0)));

	A[0][0]=1;
	A[0][1]=1;
	A[1][0]=kz0;
	A[1][1]=-kz0;
	A[0][2]=-1;
	A[0][3]=-1;
	A[1][2]=-kz1;
	A[1][3]=kz1;

}





// it is well known that exp in fft needs to be calculayed with method of
// multiplication on itself but in case of big lengths quality is missed
// and function not works.
// so here exp is calculated with  method of multiplication so you can't use this algorythm
// to work with arrays of any length.
// ! but it is Quite Fast Fourier Transformation (qFFT)
//p_o_* - is out array - CONTAINST OUTPUT DATA
//p_i_* - is in array - CONTAINS DATA TO BE TRANSFORMED
//pp_in_* - is internal array - CONTAINS NO DATA ( MUST BE SIZE OF INPUT=OUTPUT DATA)
int qqFFT(int N, real[] p_o_r, real[] p_o_i)
{
real[] pp_in_r=new real[N];
real[] pp_in_i=new real[N];
real[]	tmp_r;
real[]	tmp_i;
real[] p_o_r_0 = p_o_r;
real[] p_o_i_0 = p_o_i;
real[] pp_in_r_0 = pp_in_r;
real[] pp_in_i_0 = pp_in_i;
real b1;
real Xb_r,Xb_i,X_r,X_i,arg;
int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	
	//adr; //address in global array
	//k1;
	//dim; //dim of currrent FT
	//p_dim; //previous dim = dim/2;
	//N_l; //level 1,2,3... for FT in dim=4 N_l=2, for dim=8 N_l=3, etc.
	//cFT; //number of current FT=[0..N/dim]
	
	
/*	for (adr=0;adr<N;adr++) 
		{
		p_o_r[adr]=p_i_r[adr];
		p_o_i[adr]=p_i_i[adr];
		}
*/
	//FT
	a1=floor(N/2); //p_dim*N/dim
	dim=1;
	real N_r=N;
	int max=floor(log(N_r)/log(2.0)-1.0+0.1);
	for (N_l=0;N_l<max;++N_l)
	{
//		for (adr=0;adr<N;adr++) 
//		{
			tmp_r=p_o_r;
			tmp_i=p_o_i;
			p_o_r=pp_in_r;
			p_o_i=pp_in_i;
			pp_in_r=tmp_r;
			pp_in_i=tmp_i;
			//pp_in_r[adr]=p_o_r[adr];
			//pp_in_i[adr]=p_o_i[adr];
//		}
	p_dim=dim;
	dim=dim*2;
	a0=floor(N/dim);
	
	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] + X_r*pp_in_r[a4] - X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] + X_r*pp_in_i[a4] + X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				++a2;
				++a3;
				++a4;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

//		for (adr=0;adr<N;adr++) 
//		{
			tmp_r=p_o_r;
			tmp_i=p_o_i;
			p_o_r=pp_in_r;
			p_o_i=pp_in_i;
			pp_in_r=tmp_r;
			pp_in_i=tmp_i;
			//pp_in_r[adr]=p_o_r[adr];
			//pp_in_i[adr]=p_o_i[adr];
//		}

	//last step
	p_dim=dim;
	dim=dim*2;
	a0=floor(N/dim);

	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(k1*arg);
			//X_i=-sin(k1*arg);

			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] - X_r*pp_in_r[a4] + X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] - X_r*pp_in_i[a4] - X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				++a2;
				++a3;
				++a4;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(pi*(-N/2))/b1;
	for (k1=0;k1<N;++k1) 
		{
		p_o_r[k1]=p_o_r[k1]*Xb_r;
		p_o_i[k1]=p_o_i[k1]*Xb_r;
		Xb_r=-Xb_r;
		pp_in_r[k1]=p_o_r[k1];// because of exchanging we must have result in both arrays
		pp_in_i[k1]=p_o_i[k1];// because of exchanging we must have result in both arrays
		}
		
		// p_o_i must poit to beginning offset
		p_o_r=p_o_r_0;
		p_o_i=p_o_i_0;
		pp_in_r=pp_in_r_0;
		pp_in_i=pp_in_i_0;
	return 0;

}








// fft^-1
// it is well known that exp in fft needs to be calculayed with method of
// multiplication on itself but in case of big lengths quality is missed
// and function not works.
// so here exp is calculated with  method of multiplication so you can't use this algorythm
// to work with arrays of any length.
// ! but it is Quite Fast Fourier Transformation (qFFT)
//p_o_* - is out array - CONTAINST OUTPUT DATA
//p_i_* - is in array - CONTAINS DATA TO BE TRANSFORMED
//pp_in_* - is internal array - CONTAINS NO DATA ( MUST BE SIZE OF INPUT=OUTPUT DATA)
int qqFFT_1(int N, real[] p_o_r, real[] p_o_i)
{
real[] pp_in_r=new real[N];
real[] pp_in_i=new real[N];
real[]	tmp_r;
real[]	tmp_i;
real[] p_o_r_0 = p_o_r;
real[] p_o_i_0 = p_o_i;
real[] pp_in_r_0 = pp_in_r;
real[] pp_in_i_0 = pp_in_i;
real b1;
real Xb_r,Xb_i,X_r,X_i,arg;
int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	
	//adr; //address in global array
	//k1;
	//dim; //dim of currrent FT
	//p_dim; //previous dim = dim/2;
	//N_l; //level 1,2,3... for FT in dim=4 N_l=2, for dim=8 N_l=3, etc.
	//cFT; //number of current FT=[0..N/dim]
	
	
	//FT
	a1=floor(N/2); //p_dim*N/dim
	dim=1;
	real N_r=N;
	int max=floor(log(N_r)/log(2.0)-1.0+0.1);
	for (N_l=0;N_l<max;++N_l)
	{
//		for (adr=0;adr<N;adr++) 
//		{
			tmp_r=p_o_r;
			tmp_i=p_o_i;
			p_o_r=pp_in_r;
			p_o_i=pp_in_i;
			pp_in_r=tmp_r;
			pp_in_i=tmp_i;
			//pp_in_r[adr]=p_o_r[adr];
			//pp_in_i[adr]=p_o_i[adr];
//		}
	p_dim=dim;
	dim=dim*2;
	a0=floor(N/dim);

	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im
	
	X_r=1;//re
	X_i=0;//im
	
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(arg*k1);
			//X_i=sin(arg*k1);

			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] + X_r*pp_in_r[a4] - X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] + X_r*pp_in_i[a4] + X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				++a2;
				++a3;
				++a4;
			}


			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

//		for (adr=0;adr<N;adr++) 
//		{
			tmp_r=p_o_r;
			tmp_i=p_o_i;
			p_o_r=pp_in_r;
			p_o_i=pp_in_i;
			pp_in_r=tmp_r;
			pp_in_i=tmp_i;
			//pp_in_r[adr]=p_o_r[adr];
			//pp_in_i[adr]=p_o_i[adr];
//		}

	//last step
	p_dim=dim;
	dim=dim*2;
	a0=floor(N/dim);

	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] - X_r*pp_in_r[a4] + X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] - X_r*pp_in_i[a4] - X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				++a2;
				++a3;
				++a4;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(pi*(-N/2))/b1;
	for (k1=0;k1<N;++k1) 
		{
		p_o_r[k1]=p_o_r[k1]*Xb_r;
		p_o_i[k1]=p_o_i[k1]*Xb_r;
		Xb_r=-Xb_r;
		pp_in_r[k1]=p_o_r[k1];// because of exchanging we must have result in both arrays
		pp_in_i[k1]=p_o_i[k1];// because of exchanging we must have result in both arrays
		}

		// p_o_i must poit to beginning offset
		p_o_r=p_o_r_0;
		p_o_i=p_o_i_0;
		pp_in_r=pp_in_r_0;
		pp_in_i=pp_in_i_0;

	return 0;
}








//FFT FOR REAL FUNCTION
//f_0,f_1,f_2,f_3...f_{N-1} => A_c_0,(A_s_0=0),A_c_1,A_s_1,...A_c_{N/2-1},A_s_{N/2-1}.
int qqFFT_freal(int N, real[] p_o_r)
{
real[] pp_in_r=new real[N];
real[] pp_in_i;
real[] p_o_i;
real[]	tmp_r;
real[]	tmp_i;
real[] p_o_r_0 = p_o_r;
real[] p_o_i_0 = p_o_i;
real[] pp_in_r_0 = pp_in_r;
real[] pp_in_i_0 = pp_in_i;

int n=floor(N/2);
pp_in_i=pp_in_r;
p_o_i=p_o_r;

real b1;
real Xb_r,Xb_i,X_r,X_i,arg;
int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	
	//FT
	a1=floor(n/2); //p_dim*n/dim
	dim=1;
	real N_r=n;
	int max=floor(log(N_r)/log(2.0)-1.0+0.1);
//write("FFT max=",max);
	for (N_l=0;N_l<max+1;++N_l)
	{
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=pp_in_r;
	p_o_i=pp_in_i;
	pp_in_r=tmp_r;
	pp_in_i=tmp_i;

	p_dim=dim;
	dim=dim*2;
	a0=floor(n/dim);
	
	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<n/dim;cFT++
		{
			a2=k1*a0; //=k1*n/dim
			a3=a2+a2; //=k1*2*n/dim
			a4=a3+a0; //=(k1*2+1)*n/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2*2]=pp_in_r[a3*2] + X_r*pp_in_r[a4*2] - X_i*pp_in_i[a4*2+1];
				p_o_i[a2*2+1]=pp_in_i[a3*2+1] + X_r*pp_in_i[a4*2+1] + X_i*pp_in_r[a4*2];
				
				p_o_r[(a2+a1)*2]=2.0*pp_in_r[a3*2] - p_o_r[a2*2];
				p_o_i[(a2+a1)*2+1]=2.0*pp_in_i[a3*2+1] - p_o_i[a2*2+1];
				
				++a2;
				++a3;
				++a4;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	real b_norm=sqrt(n);
	real X_w,Y_w,X_n,Y_n;
	X_i=0.0;
	X_r=1.0;
	Xb_i=sin(2.0*pi/N);
	Xb_r=cos(2.0*pi/N);
	
	X_w=p_o_r[0];
	Y_w=p_o_r[1];
	X_n=0;
	Y_n=0;

	p_o_r[0]=(X_w+Y_w)/1.0/b_norm;    //A
	p_o_r[1]=(-Y_w+X_w)/1.0/b_norm;    //B
	pp_in_r[0]=p_o_r[0];
	pp_in_r[1]=p_o_r[1];

	b1=X_r;
	X_r=X_r*Xb_r-X_i*Xb_i;
	X_i=X_i*Xb_r + b1*Xb_i;

	
	for (k1=1;k1<n/2+1;++k1) 
	{
//		X_i=sin(2.0*pi/real(N)*real(k1));
//		X_r=cos(2.0*pi/real(N)*real(k1));

		X_w=p_o_r[k1*2];
		Y_w=p_o_r[k1*2+1];
		X_n=p_o_r[(n-k1)*2];
		Y_n=p_o_r[(n-k1)*2+1];

		p_o_r[k1*2]=(X_w+X_n-(X_w-X_n)*X_i+(Y_w+Y_n)*X_r)/2.0/b_norm;    //A
		p_o_r[k1*2+1]=(-(Y_w-Y_n)+(X_w-X_n)*X_r+(Y_w+Y_n)*X_i)/2.0/b_norm;    //B

		p_o_r[(n-k1)*2]=(X_w+X_n+(X_w-X_n)*X_i-(Y_w+Y_n)*X_r)/2.0/b_norm;    //C
		p_o_r[(n-k1)*2+1]=(Y_w-Y_n+(X_w-X_n)*X_r+(Y_w+Y_n)*X_i)/2.0/b_norm;    //D

		pp_in_r[k1*2]=p_o_r[k1*2];
		pp_in_r[k1*2+1]=p_o_r[k1*2+1];
		pp_in_r[(n-k1)*2]=p_o_r[(n-k1)*2];
		pp_in_r[(n-k1)*2+1]=p_o_r[(n-k1)*2+1];

		b1=X_r;
		X_r=X_r*Xb_r-X_i*Xb_i;
		X_i=X_i*Xb_r + b1*Xb_i;
	}
	
	// p_o_i must poit to beginning offset
	p_o_r=p_o_r_0;
	p_o_i=p_o_i_0;
	pp_in_r=pp_in_r_0;
	pp_in_i=pp_in_i_0;

	return 0;

}


//INVERT
int qqFFT_freal_1(int N, real[] p_o_r)
{
real[] pp_in_r=new real[N];
real[] pp_in_i;
real[] p_o_i;
real[]	tmp_r;
real[]	tmp_i;
real[] p_o_r_0 = p_o_r;
real[] p_o_i_0 = p_o_i;
real[] pp_in_r_0 = pp_in_r;
real[] pp_in_i_0 = pp_in_i;

int n=floor(N/2);
pp_in_i=pp_in_r;
p_o_i=p_o_r;

real b1;
real Xb_r,Xb_i,X_r,X_i,arg;
int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	real b_norm=sqrt(n);
	real A,B,C,D;
	X_i=0.0;
	X_r=1.0;
	Xb_i=sin(2.0*pi/N);
	Xb_r=cos(2.0*pi/N);
	
	A=p_o_r[0];
	B=p_o_r[1];
	C=0;
	D=0;

	p_o_r[0]=(A+B)/2.0/b_norm;    //X_w
	p_o_r[1]=(A-B)/2.0/b_norm;    //Y_w

	b1=X_r;
	X_r=X_r*Xb_r-X_i*Xb_i;
	X_i=X_i*Xb_r + b1*Xb_i;

	
	for (k1=1;k1<n/2+1;++k1) 
	{
//		X_i=sin(2.0*pi/real(N)*real(k1));
//		X_r=cos(2.0*pi/real(N)*real(k1));

		A=p_o_r[k1*2];
		B=p_o_r[k1*2+1];
		C=p_o_r[(n-k1)*2];
		D=p_o_r[(n-k1)*2+1];

		p_o_r[k1*2]=(A+C+X_i*(C-A)+(B+D)*X_r)/2.0/b_norm;    //X_w
		p_o_r[k1*2+1]=((A-C)*X_r+D-B+X_i*(D+B))/2.0/b_norm;    //Y_w

		p_o_r[(n-k1)*2]=(A+C+X_i*(A-C)-X_r*(D+B))/2.0/b_norm;    //X_n
		p_o_r[(n-k1)*2+1]=((A-C)*X_r+B-D+X_i*(B+D))/2.0/b_norm;    //Y_n

		b1=X_r;
		X_r=X_r*Xb_r-X_i*Xb_i;
		X_i=X_i*Xb_r + b1*Xb_i;
	}

	//FT
	a1=floor(n/2); //p_dim*n/dim
	dim=1;
	real N_r=n;
	int max=floor(log(N_r)/log(2.0)-1.0+0.1);
	for (N_l=0;N_l<max+1;++N_l)
	{
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=pp_in_r;
	p_o_i=pp_in_i;
	pp_in_r=tmp_r;
	pp_in_i=tmp_i;

	p_dim=dim;
	dim=dim*2;
	a0=floor(n/dim);
	
	arg=pi;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;++k1)     //cFT=0;cFT<n/dim;cFT++
		{
			a2=k1*a0; //=k1*n/dim
			a3=a2+a2; //=k1*2*n/dim
			a4=a3+a0; //=(k1*2+1)*n/dim
			for (cFT=0;cFT<a0;++cFT)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2*2]=pp_in_r[a3*2] + X_r*pp_in_r[a4*2] - X_i*pp_in_i[a4*2+1];
				p_o_i[a2*2+1]=pp_in_i[a3*2+1] + X_r*pp_in_i[a4*2+1] + X_i*pp_in_r[a4*2];
				
				p_o_r[(a2+a1)*2]=2.0*pp_in_r[a3*2] - p_o_r[a2*2];
				p_o_i[(a2+a1)*2+1]=2.0*pp_in_i[a3*2+1] - p_o_i[a2*2+1];
				
				++a2;
				++a3;
				++a4;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	for (k1=0;k1<n;++k1) 
	{
		pp_in_r[k1*2]=p_o_r[k1*2]; // because of exchanging we must have result in both arrays
		pp_in_r[k1*2+1]=p_o_r[k1*2+1]; // because of exchanging we must have result in both arrays
	}

	// p_o_i must poit to beginning offset
	p_o_r=p_o_r_0;
	p_o_i=p_o_i_0;
	pp_in_r=pp_in_r_0;
	pp_in_i=pp_in_i_0;

	return 0;

}

