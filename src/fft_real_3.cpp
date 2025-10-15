// it is well known that exp in fft needs to be calculayed with method of
// multiplication on itself but in case of big lengths quality is missed
// and function not works.
// so here exp is calculated with  method of multiplication so you can't use this algorythm
// to work with arrays of any length.
// ! but it is Quite Fast Fourier Transformation (qFFT)
//p_o_* - is out array - CONTAINST OUTPUT DATA
//p_i_* - is in array - CONTAINS DATA TO BE TRANSFORMED
//pp_in_* - is internal array - CONTAINS NO DATA ( MUST BE SIZE OF INPUT=OUTPUT DATA)
static int qqFFT(long int N, FL_DBL* p_o_r, FL_DBL* p_o_i, FL_DBL* pp_in_r, FL_DBL* pp_in_i)
{
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0 = p_o_r;
FL_DBL* p_o_i_0 = p_o_i;
FL_DBL* pp_in_r_0 = pp_in_r;
FL_DBL* pp_in_i_0 = pp_in_i;
FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,k1,dim,p_dim,N_l,cFT;

	
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
	a1=N/2; //p_dim*N/dim
	dim=1;
	FL_DBL N_r=N;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max;N_l++)
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
	a0=N/dim;
	
	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] + X_r*pp_in_r[a4] - X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] + X_r*pp_in_i[a4] + X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				a2++;
				a3++;
				a4++;
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
	a0=N/dim;

	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(k1*arg);
			//X_i=-sin(k1*arg);

			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] - X_r*pp_in_r[a4] + X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] - X_r*pp_in_i[a4] - X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(M_PI*(-N/2))/b1;
	for (k1=0;k1<N;k1++) 
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
static int qqFFT_1(long int N, FL_DBL* p_o_r, FL_DBL* p_o_i, FL_DBL* pp_in_r, FL_DBL* pp_in_i)
{
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0 = p_o_r;
FL_DBL* p_o_i_0 = p_o_i;
FL_DBL* pp_in_r_0 = pp_in_r;
FL_DBL* pp_in_i_0 = pp_in_i;
FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,k1,dim,p_dim,N_l,cFT;

	
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
	a1=N/2; //p_dim*N/dim
	dim=1;
	FL_DBL N_r=N;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max;N_l++)
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
	a0=N/dim;

	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im
	
	X_r=1;//re
	X_i=0;//im
	
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(arg*k1);
			//X_i=sin(arg*k1);

			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] + X_r*pp_in_r[a4] - X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] + X_r*pp_in_i[a4] + X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				a2++;
				a3++;
				a4++;
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
	a0=N/dim;

	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2]=pp_in_r[a3] - X_r*pp_in_r[a4] + X_i*pp_in_i[a4];
				p_o_i[a2]=pp_in_i[a3] - X_r*pp_in_i[a4] - X_i*pp_in_r[a4];
				
				p_o_r[a2+a1]=2.0*pp_in_r[a3] - p_o_r[a2];
				p_o_i[a2+a1]=2.0*pp_in_i[a3] - p_o_i[a2];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(M_PI*(-N/2))/b1;
	for (k1=0;k1<N;k1++) 
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


/*
void ft(long int N)
{
	FL_DBL x,y,cs,sn,arg;
	int i;
	for(i=0;i<N;i++)
	{
		p_o_r[i]=0;
		p_o_i[i]=0;
	}

	for(i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			arg=2*M_PI;
			arg=arg/N;
			arg=arg*(i-N/2)*(j-N/2);
			cs=cos(arg);
			sn=sin(arg);
			p_o_r[i]=p_o_r[i]+p_i_r[j]*cs-p_i_i[j]*sn;
			p_o_i[i]=p_o_i[i]+p_i_r[j]*sn+p_i_i[j]*cs;
		}
	}
	for(i=0;i<N;i++)
	{
		p_o_r[i]=p_o_r[i]/sqrt(N);
		p_o_i[i]=p_o_i[i]/sqrt(N);
	}
}

void ft_1(long int N)
{
	FL_DBL x,y,cs,sn,arg;
	int i;
	for(i=0;i<N;i++)
	{
		p_o_r[i]=0;
		p_o_i[i]=0;
	}

	for(i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			arg=-2*M_PI;
			arg=arg/N;
			arg=arg*(i-N/2)*(j-N/2);
			cs=cos(arg);
			sn=sin(arg);
			p_o_r[i]=p_o_r[i]+p_i_r[j]*cs-p_i_i[j]*sn;
			p_o_i[i]=p_o_i[i]+p_i_r[j]*sn+p_i_i[j]*cs;
		}
	}
	for(i=0;i<N;i++)
	{
		p_o_r[i]=p_o_r[i]/sqrt(N);
		p_o_i[i]=p_o_i[i]/sqrt(N);
	}
}

*/

//FFT FOR REAL FUNCTION
//f_0,f_1,f_2,f_3...f_{N-1} => A_c_0,(A_s_0=0),A_c_1,A_s_1,...A_c_{N/2-1},A_s_{N/2-1}.
static int qqFFT_freal(long int N, FL_DBL* p_o_r, FL_DBL* pp_in_r)
{
FL_DBL* pp_in_i;
FL_DBL* p_o_i;
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0 = p_o_r;
FL_DBL* p_o_i_0 = p_o_i;
FL_DBL* pp_in_r_0 = pp_in_r;
FL_DBL* pp_in_i_0 = pp_in_i;

int n=N/2;
pp_in_i=&pp_in_r[1];
p_o_i=&p_o_r[1];

FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,k1,dim,p_dim,N_l,cFT;

	
	//FT
	a1=n/2; //p_dim*n/dim
	dim=1;
	FL_DBL N_r=n;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max+1;N_l++)
	{
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=pp_in_r;
	p_o_i=pp_in_i;
	pp_in_r=tmp_r;
	pp_in_i=tmp_i;

	p_dim=dim;
	dim=dim*2;
	a0=n/dim;
	
	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<n/dim;cFT++
		{
			a2=k1*a0; //=k1*n/dim
			a3=a2+a2; //=k1*2*n/dim
			a4=a3+a0; //=(k1*2+1)*n/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2*2]=pp_in_r[a3*2] + X_r*pp_in_r[a4*2] - X_i*pp_in_i[a4*2];
				p_o_i[a2*2]=pp_in_i[a3*2] + X_r*pp_in_i[a4*2] + X_i*pp_in_r[a4*2];
				
				p_o_r[(a2+a1)*2]=2.0*pp_in_r[a3*2] - p_o_r[a2*2];
				p_o_i[(a2+a1)*2]=2.0*pp_in_i[a3*2] - p_o_i[a2*2];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	FL_DBL b_norm=sqrt(FL_DBL(n));
	FL_DBL X_w,Y_w,X_n,Y_n;
	X_i=0.0;
	X_r=1.0;
	Xb_i=sin(2.0*M_PI/FL_DBL(N));
	Xb_r=cos(2.0*M_PI/FL_DBL(N));
	
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

	
	for (k1=1;k1<n/2+1;k1++) 
	{
//		X_i=sin(2.0*M_PI/FL_DBL(N)*FL_DBL(k1));
//		X_r=cos(2.0*M_PI/FL_DBL(N)*FL_DBL(k1));

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
static int qqFFT_freal_1(long int N, FL_DBL* p_o_r, FL_DBL* pp_in_r)
{
FL_DBL* pp_in_i;
FL_DBL* p_o_i;
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0 = p_o_r;
FL_DBL* p_o_i_0 = p_o_i;
FL_DBL* pp_in_r_0 = pp_in_r;
FL_DBL* pp_in_i_0 = pp_in_i;

int n=N/2;
pp_in_i=&pp_in_r[1];
p_o_i=&p_o_r[1];

FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,k1,dim,p_dim,N_l,cFT;

	FL_DBL b_norm=sqrt(FL_DBL(n));
	FL_DBL A,B,C,D;
	X_i=0.0;
	X_r=1.0;
	Xb_i=sin(2.0*M_PI/FL_DBL(N));
	Xb_r=cos(2.0*M_PI/FL_DBL(N));
	
	A=p_o_r[0];
	B=p_o_r[1];
	C=0;
	D=0;

	p_o_r[0]=(A+B)/2.0/b_norm;    //X_w
	p_o_r[1]=(A-B)/2.0/b_norm;    //Y_w

	b1=X_r;
	X_r=X_r*Xb_r-X_i*Xb_i;
	X_i=X_i*Xb_r + b1*Xb_i;

	
	for (k1=1;k1<n/2+1;k1++) 
	{
//		X_i=sin(2.0*M_PI/FL_DBL(N)*FL_DBL(k1));
//		X_r=cos(2.0*M_PI/FL_DBL(N)*FL_DBL(k1));

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
	a1=n/2; //p_dim*n/dim
	dim=1;
	FL_DBL N_r=n;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max+1;N_l++)
	{
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=pp_in_r;
	p_o_i=pp_in_i;
	pp_in_r=tmp_r;
	pp_in_i=tmp_i;

	p_dim=dim;
	dim=dim*2;
	a0=n/dim;
	
	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im
	
	X_r=1.0;//re
	X_i=0.0;//im
	
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<n/dim;cFT++
		{
			a2=k1*a0; //=k1*n/dim
			a3=a2+a2; //=k1*2*n/dim
			a4=a3+a0; //=(k1*2+1)*n/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				p_o_r[a2*2]=pp_in_r[a3*2] + X_r*pp_in_r[a4*2] - X_i*pp_in_i[a4*2];
				p_o_i[a2*2]=pp_in_i[a3*2] + X_r*pp_in_i[a4*2] + X_i*pp_in_r[a4*2];
				
				p_o_r[(a2+a1)*2]=2.0*pp_in_r[a3*2] - p_o_r[a2*2];
				p_o_i[(a2+a1)*2]=2.0*pp_in_i[a3*2] - p_o_i[a2*2];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	for (k1=0;k1<n;k1++) 
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







/*
int qqFFT_2d(long int N, int step, int N_opp, int begin)
{
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0=p_o_r;
FL_DBL* p_o_i_0=p_o_i;
FL_DBL* p_i_r_0=p_i_r;
FL_DBL* p_i_i_0=p_i_i;
FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	//FT
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=p_i_r;
	p_o_i=p_i_i;
	p_i_r=tmp_r;
	p_i_i=tmp_i;

	a1=N/2; //p_dim*N/dim
	dim=1;
	FL_DBL N_r=N;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max;N_l++)
	{
		tmp_r=p_o_r;
		tmp_i=p_o_i;
		p_o_r=p_i_r;
		p_o_i=p_i_i;
		p_i_r=tmp_r;
		p_i_i=tmp_i;

		p_dim=dim;
		dim=dim*2;
		a0=N/dim;
		
		arg=M_PI;
		arg=arg/p_dim;
		Xb_r=cos(arg);//re
		Xb_i=-sin(arg);//im
		
		X_r=1.0;//re
		X_i=0.0;//im
		
		for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				for(int i=begin;i<N_opp;i++)	p_o_r[a2*step+i]=p_i_r[a3*step+i] + X_r*p_i_r[a4*step+i] - X_i*p_i_i[a4*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[a2*step+i]=p_i_i[a3*step+i] + X_r*p_i_i[a4*step+i] + X_i*p_i_r[a4*step+i];
				
				for(i=begin;i<N_opp;i++)	p_o_r[(a2+a1)*step+i]=2.0*p_i_r[a3*step+i]  - p_o_r[a2*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[(a2+a1)*step+i]=2.0*p_i_i[a3*step+i]  - p_o_i[a2*step+i];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=p_i_r;
	p_o_i=p_i_i;
	p_i_r=tmp_r;
	p_i_i=tmp_i;

	//last step
	p_dim=dim;
	dim=dim*2;
	a0=N/dim;

	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=-sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(k1*arg);
			//X_i=-sin(k1*arg);

			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				for(int i=begin;i<N_opp;i++)	p_o_r[a2*step+i]=p_i_r[a3*step+i] - X_r*p_i_r[a4*step+i] + X_i*p_i_i[a4*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[a2*step+i]=p_i_i[a3*step+i] - X_r*p_i_i[a4*step+i] - X_i*p_i_r[a4*step+i];
				
				for(i=begin;i<N_opp;i++)	p_o_r[(a2+a1)*step+i]=2.0*p_i_r[a3*step+i] - p_o_r[a2*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[(a2+a1)*step+i]=2.0*p_i_i[a3*step+i] - p_o_i[a2*step+i];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(M_PI*(-N/2))/b1;
	for (k1=0;k1<N;k1++) 
	{
		for(int i=begin;i<N_opp;i++)
		{
			p_o_r[k1*step+i]=p_o_r[k1*step+i]*Xb_r;
			p_o_i[k1*step+i]=p_o_i[k1*step+i]*Xb_r;
			p_i_r[k1*step+i]=p_o_r[k1*step+i];// because of exchanging we must have result in both arrays
			p_i_i[k1*step+i]=p_o_i[k1*step+i];// because of exchanging we must have result in both arrays
		}
		Xb_r=-Xb_r;

	}
		
		// p_o_i must poit to beginning offset
		p_o_r=p_o_r_0;
		p_o_i=p_o_i_0;
		p_i_r=p_i_r_0;
		p_i_i=p_i_i_0;
	return 0;

}


int qqFFT_2d_1(long int N, int step, int N_opp, int begin)
{
FL_DBL*	tmp_r;
FL_DBL*	tmp_i;
FL_DBL* p_o_r_0 = p_o_r;
FL_DBL* p_o_i_0 = p_o_i;
FL_DBL* p_i_r_0 = p_i_r;
FL_DBL* p_i_i_0 = p_i_i;
FL_DBL b1;
FL_DBL Xb_r,Xb_i,X_r,X_i,arg;
long int a1,a2,a3,a4,a0,adr,k1,dim,p_dim,N_l,cFT;

	//FT
	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=p_i_r;
	p_o_i=p_i_i;
	p_i_r=tmp_r;
	p_i_i=tmp_i;

	a1=N/2; //p_dim*N/dim
	dim=1;
	FL_DBL N_r=N;
	int max=log(N_r)/log(2.0)-1.0+0.1;
	for (N_l=0;N_l<max;N_l++)
	{
		tmp_r=p_o_r;
		tmp_i=p_o_i;
		p_o_r=p_i_r;
		p_o_i=p_i_i;
		p_i_r=tmp_r;
		p_i_i=tmp_i;

		p_dim=dim;
		dim=dim*2;
		a0=N/dim;
		
		arg=M_PI;
		arg=arg/p_dim;
		Xb_r=cos(arg);//re
		Xb_i=sin(arg);//im
		
		X_r=1.0;//re
		X_i=0.0;//im
		
		for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			a2=k1*a0; //=k1*N/dim
			a3=a2+a2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				for(int	i=begin;i<N_opp;i++)	p_o_r[a2*step+i]=p_i_r[a3*step+i] + X_r*p_i_r[a4*step+i] - X_i*p_i_i[a4*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[a2*step+i]=p_i_i[a3*step+i] + X_r*p_i_i[a4*step+i] + X_i*p_i_r[a4*step+i];
				
				for(i=begin;i<N_opp;i++)	p_o_r[(a2+a1)*step+i]=2.0*p_i_r[a3*step+i] - p_o_r[a2*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[(a2+a1)*step+i]=2.0*p_i_i[a3*step+i] - p_o_i[a2*step+i];
				
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}
	}

	tmp_r=p_o_r;
	tmp_i=p_o_i;
	p_o_r=p_i_r;
	p_o_i=p_i_i;
	p_i_r=tmp_r;
	p_i_i=tmp_i;

	//last step
	p_dim=dim;
	dim=dim*2;
	a0=N/dim;

	arg=M_PI;
	arg=arg/p_dim;
	Xb_r=cos(arg);//re
	Xb_i=sin(arg);//im

	X_r=1;//re
	X_i=0;//im
	for (k1=0;k1<p_dim;k1++)     //cFT=0;cFT<N/dim;cFT++
		{
			//X_r=cos(k1*arg);
			//X_i=-sin(k1*arg);

			a2=k1*a0; //=k1*N/dim
			a3=a2*2; //=k1*2*N/dim
			a4=a3+a0; //=(k1*2+1)*N/dim
			for (cFT=0;cFT<a0;cFT++)      //k1=0;k1<p_dim;k1++
			{
				for(int i=begin;i<N_opp;i++)	p_o_r[a2*step+i]=p_i_r[a3*step+i] - X_r*p_i_r[a4*step+i] + X_i*p_i_i[a4*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[a2*step+i]=p_i_i[a3*step+i] - X_r*p_i_i[a4*step+i] - X_i*p_i_r[a4*step+i];
				
				for(i=begin;i<N_opp;i++)	p_o_r[(a2+a1)*step+i]=2.0*p_i_r[a3*step+i] - p_o_r[a2*step+i];
				for(i=begin;i<N_opp;i++)	p_o_i[(a2+a1)*step+i]=2.0*p_i_i[a3*step+i] - p_o_i[a2*step+i];
								
				a2++;
				a3++;
				a4++;
			}

			b1=X_r;
			X_r=X_r*Xb_r-X_i*Xb_i;
			X_i=X_i*Xb_r + b1*Xb_i;
		}

	b1=sqrt(N);
	Xb_r=cos(M_PI*(-N/2))/b1;
	for (k1=0;k1<N;k1++) 
	{
		for(int i=begin;i<N_opp;i++)
		{
			p_o_r[k1*step+i]=p_o_r[k1*step+i]*Xb_r;
			p_o_i[k1*step+i]=p_o_i[k1*step+i]*Xb_r;
			p_i_r[k1*step+i]=p_o_r[k1*step+i];// because of exchanging we must have result in both arrays
			p_i_i[k1*step+i]=p_o_i[k1*step+i];// because of exchanging we must have result in both arrays
		}
		Xb_r=-Xb_r;
	}
		
		// p_o_i must poit to beginning offset
		p_o_r=p_o_r_0;
		p_o_i=p_o_i_0;
		p_i_r=p_i_r_0;
		p_i_i=p_i_i_0;
	return 0;

}

*/
