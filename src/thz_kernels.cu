#include <math_constants.h>
#include <cutil_inline.h>
#include "dfs.h"
#include "sweep_kernel.cu"
#include "nl_kernel.cu"
#define PRINT
FL_DBL* Aw;
FL_DBL* Ax;
FL_DBL* BCw;
FL_DBL* Bx;
FL_DBL* Cx;
FL_DBL* dev_Aw;
FL_DBL* dev_Ax;
FL_DBL* dev_BCw;
FL_DBL* dev_Bx;
FL_DBL* dev_Cx;
FL_DBL* dev_spec;
FL_DBL* dev_plasma1; // to collect plasma distribution on first step
FL_DBL* dev_plasma2; // to collect plasma distribution on second step
FL_DBL* dev_THz_src; // to collect THz source distribution
FL_DBL* dev_spec_zero_freq; // small array containing zero frequency
FL_DBL* dev_Nl_spec;
FL_DBL* dev_Nl_spec_2;
FL_DBL* dev_spec_tmp;
FL_DBL* dev_bb;
FL_DBL* dev_cc;
FL_DBL*	dev_a2;
FL_DBL*	dev_b2;
FL_DBL*	dev_c2;
FL_DBL*	dev_x2;
FL_DBL*	dev_y2;
int Bs, BsX, smem;
cufftHandle planF,planI;
FL_DBL _dt,_dz;

extern "C"
void get_field(FL_DBL*, int, int);

extern "C"
FL_DBL step(int, int, FL_DBL, bool);

void dev_alloc(void** pparr, int sz) {cutilSafeCall( cudaMalloc( pparr, sz*sizeof(FL_DBL)));}

void dev_h2d(FL_DBL* host_arr, FL_DBL* dev_arr, int sz) {cutilSafeCall( cudaMemcpy( dev_arr, host_arr, sz*sizeof(FL_DBL), cudaMemcpyHostToDevice));}

void dev_d2h(FL_DBL* dev_arr, FL_DBL* host_arr, int sz) {cutilSafeCall( cudaMemcpy( host_arr, dev_arr, sz*sizeof(FL_DBL), cudaMemcpyDeviceToHost));}

__global__ void
store_zero_freq_Kernel(FL_DBL* spc, FL_DBL* spc_zerofreq, int Nt, int Nx)
{
//because the method of invertion of A matrix we need to specify some value for omega when i=0 (zero frequency)
//that is why we need to restore the value of zero frequency on each step
	const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
	if(tidx<Nx)
	{
		spc_zerofreq[tidx*2]=spc[tidx*Nt];
		spc_zerofreq[tidx*2+1]=spc[tidx*Nt+1];
		spc_zerofreq[tidx*2+Nx*2]=spc[tidx*Nt+Nt-2];
		spc_zerofreq[tidx*2+1+Nx*2]=spc[tidx*Nt+Nt-1];
	}
}


__global__ void
restore_zero_freq_Kernel(FL_DBL* spc_zerofreq, FL_DBL* spc, int Nt, int Nx)
{
//because the method of invertion of A matrix we need to specify some value for omega when i=0 (zero frequency)
//that is why we need to restore the value of zero frequency on each step
	const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
	if(tidx<Nx)
	{
		spc[tidx*Nt]=spc_zerofreq[tidx*2];
		spc[tidx*Nt+1]=spc_zerofreq[tidx*2+1];
		spc[tidx*Nt+Nt-2]=spc_zerofreq[tidx*2+Nx*2];
		spc[tidx*Nt+Nt-1]=spc_zerofreq[tidx*2+1+Nx*2];
	}
}

extern "C" 
void cuda_load(FL_DBL* spec, int Nt, int Nx, FL_DBL dt, FL_DBL dx, FL_DBL dz)
{
	printf("NOFA=%d\n",NOFA);
	_dt=dt;
	_dz=dz;
	int i,j;
//init matrix
	Ax=new FL_DBL[NOFA*Nx+Nx+Nx+Nt*NOFA+Nt];
	Bx=&Ax[Nx*NOFA];
	Cx=&Bx[Nx];
	Aw=&Cx[Nx];
	BCw=&Aw[Nt*NOFA];
	for(i=0;i<Nx;i++) Ax[i]=1.0; // Identity
	for(i=0;i<Nx;i++) Ax[i+Nx]=(((i==0)||(i==Nx-1))?1.0:2.0)*dz/dx/dx+1.0/(double(i)+0.5)/dx/dx*dz; // Laplassian

	for(i=0;i<Nx;i++) Bx[i]=-1.0*dz/dx/dx; // Laplassian
	//Bx[Nx-1]=-dz/dx/dx;  //non hermit
	Bx[0]=0.0F;   // BECAUSE FAST SWEEP HAS A BUG THIS SHOULD BE 0.0! ALWAYS!!!
	for(i=0;i<Nx;i++) Cx[i]=-1.0*dz/dx/dx-1.0/(double(i)+0.5)/dx/dx*dz; // Laplassian
	//Cx[0]=-2.0F*dz/dx/dx-1.0F/(0.5F)/dx/dx*dz;  // non hermit
	//Ax[0]=0.5*Ax[Nx-1]; // hermit   USE HERMIT MATRIX TO AVOID PROBLEMS WITH ENERGY CONSERVATION CONCERNED WITH BOUNDARY EFFECTS
	Cx[Nx-1]=0.0;   // BECAUSE FAST SWEEP HAS A BUG THIS SHOULD BE 0.0! ALWAYS!!!
	for(i=0;i<Nt/2;i++)
	{
		//NOFA - N of A. here we have 3 different A
		Aw[i*2]=1.0; // Identity
		Aw[i*2+1]=0.0; // Identity

		FL_DBL omega=(i==0)?1.0:1.0/(2.0*Pi/T)/(double(i))*(0.5-0.5*tanh(-(2.0*Pi/250.0*double(i)-0.5)/0.1)); // be carefull with i=0. FAST SWEEP will not work with omega=1 <=> A=I YES! FAST SWEEP can't invert identity matrix.

		Aw[i*2+Nt]=0.0; // Laplassian
		Aw[i*2+1+Nt]=omega; // Laplassian

		BCw[i*2]=0.0; // Laplassian
		BCw[i*2+1]=omega; // Laplassian
	}

//load to device
	dev_alloc((void**) &dev_Ax,	Nx*NOFA);
	dev_h2d(Ax,dev_Ax,Nx*NOFA);

	dev_alloc((void**) &dev_Bx, Nx);
	dev_h2d(Bx,dev_Bx,Nx);

	dev_alloc((void**) &dev_Cx, Nx);
	dev_h2d(Cx,dev_Cx,Nx);

	dev_alloc((void**) &dev_Aw, Nt*NOFA);
	dev_h2d(Aw,dev_Aw,Nt*NOFA);

	dev_alloc((void**) &dev_BCw, Nt);
	dev_h2d(BCw,dev_BCw,Nt);


	int gpu_time;

	Bs=32; //!!!!!
	BsX=32;   //!!!!!!
	// doesnot work for Bs=16 x BsX=128 !!!!!!!! DEBUG THIS !!!!!
	smem=sizeof(FL_DBL)*(BsX*2+BsX*NOFA+Bs*4+Bs*2*NOFA);
	dim3 grid(-floor(-FL_DBL(Nt/2)/FL_DBL(Bs)), Nx/BsX, 1);
	dim3 threads(Bs, 1, 1);
	dim3 simple_grid(-floor(-FL_DBL(Nt/2)/FL_DBL(Bs)), 1, 1);
	dim3 simple_threads(Bs, 1, 1);
	//dim3 forward_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(32)), Nx/16, 1);
	//dim3 forward_threads(32, 16, 1);
	dim3 forward_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(16)), Nx/16, 1);
	dim3 forward_threads(16, 16, 1);

	cufftPlan1d(&planF, (Nt-2), cuFFT_FP, Nx);
	cufftPlan1d(&planI, (Nt-2), cuFFT_BP, Nx);

//set spectrum
	dev_alloc((void**) &dev_spec,Nx*Nt);
	dev_alloc((void**) &dev_plasma1,Nx*Nt);
	dev_alloc((void**) &dev_plasma2,Nx*Nt);
	dev_alloc((void**) &dev_THz_src,Nx*Nt);
	dev_alloc((void**) &dev_spec_zero_freq,Nx*4);
	dev_alloc((void**) &dev_Nl_spec,Nx*Nt);
	dev_alloc((void**) &dev_Nl_spec_2,Nx*Nt);
	dev_alloc((void**) &dev_spec_tmp,Nx*Nt);
//goto Fourier (using GPU)
	cutilSafeCall( cudaMemcpy( dev_spec_tmp, spec, Nx*Nt*sizeof(FL_DBL), cudaMemcpyHostToDevice));
	cudaThreadSynchronize();
	cuFFT_F(planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_spec);
	cudaThreadSynchronize();
	cuFFT_F(planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_Nl_spec);
	cudaThreadSynchronize();

//allocate service arrays
	dev_alloc((void**) &dev_bb,Nx*Nt);
	dev_alloc((void**) &dev_cc,Nx*Nt);

	dev_alloc((void**) &dev_a2, Nx/BsX*Nt*2);
	dev_alloc((void**) &dev_b2, Nx/BsX*Nt*2);
	dev_alloc((void**) &dev_c2, Nx/BsX*Nt*2);
	dev_alloc((void**) &dev_x2, Nx/BsX*Nt*2);
	dev_alloc((void**) &dev_y2, Nx/BsX*Nt*2);
//
	//FL_DBL gpu_time=step(Nt, Nx);

	FL_DBL* spec_test=new FL_DBL[Nt*Nx];
	dev_d2h(dev_spec,spec_test,Nx*Nt);
//on CPU test
	//goto back from Fourier (using GPU)
		cuFFT_B(planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp);
		cudaThreadSynchronize();
		cutilSafeCall( cudaMemcpy( spec, dev_spec_tmp, Nx*Nt*sizeof(FL_DBL), cudaMemcpyDeviceToHost));
		cudaThreadSynchronize();
	FL_DBL norm=0.0;
	for(i=0;i<Nt*Nx;i++)	norm+=spec[i]*spec[i];
	printf("norm=%g\n",double(norm));
	
//store zero frequency components to restore them on each step later
	store_zero_freq_Kernel<<<-floor(-Nx/Bs), Bs>>>(dev_spec, dev_spec_zero_freq, Nt, Nx);
	cudaThreadSynchronize();

 	printf( "testing\n");

// test with random x matrix but previously set A B C matrixes
	{
		// create test arrays
		FL_DBL* testin=new FL_DBL[Nt*Nx];
		FL_DBL* testout=new FL_DBL[Nt*Nx];
		FL_DBL* testresCPU=new FL_DBL[Nt*Nx];
		for(int i=0;i<Nt;i++)
		{
			for(int j=0;j<Nx;j++)	testin[i+j*Nt]=0.5F+0.5F*FL_DBL(rand())/FL_DBL(RAND_MAX);
		}
		FL_DBL* dev_testin;
		dev_alloc((void**) &dev_testin,	Nx*Nt);

		cudaThreadSynchronize();
		dev_h2d(testin,dev_testin,Nt*Nx);
		cudaThreadSynchronize();
		FL_DBL* dev_testout;
		dev_alloc((void**) &dev_testout,	Nx*Nt);


		// x1(=dev_testout) from A * x1 = x0,  where x0=dev_testin
		SWEEP_1_FIXED_Kernel<<< grid, threads, smem>>>(dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_testin, dev_testout, dev_bb, dev_cc, Nx, Nt/2, BsX, Bs, dev_testin, dev_a2, dev_b2, dev_c2, dev_y2);
		cudaThreadSynchronize();

		simpleSWEEPKernel<<< simple_grid, simple_threads >>>(dev_a2, dev_b2, dev_c2, dev_y2, dev_x2, 2*Nx/BsX, Nt/2, dev_y2);
		cudaThreadSynchronize();

		SWEEP_2_Kernel<<< grid, threads >>>(dev_testout, dev_bb, dev_cc, dev_testin, Nx, Nt/2, BsX, Bs, dev_x2);
		cudaThreadSynchronize();

		// x2 = A * x1, where x1=dev_testout; x2=dev_testin
	  	forward_t_Kernel<<< forward_grid, forward_threads>>>(dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_testout, dev_testin, Nx, Nt/2);
		cudaThreadSynchronize();

		dev_d2h(dev_testin,testout,Nt*Nx);
		cudaThreadSynchronize();

		norm=0.0F;
		FL_DBL diff=0.0F;
		for(int i=0;i<Nt-2;i++)
		{
			for(int j=0;j<Nx;j++)
			{
				diff+=pow(testin[i+j*Nt]-testout[i+j*Nt],2);
				norm+=pow(testout[i+j*Nt],2);
			}
		}

		printf("TEST RESULTS(type=%s): %g (diff=%g, norm=%g) [NOTE: result should be ~e-15 for float and ~e-32 for double]\n",DATATYPE,diff/norm,diff,norm);
		delete[] testin;
		delete[] testout;
		cudaFree(dev_testin);
		cudaFree(dev_testout);
	}

	//printout the A matrix for debuging
	{
		FILE* fA=fopen("Atest","w");
		for(int i=0;i<Nt/2;i++)
		{
			for(int j=0;j<Nx;j++) fprintf(fA,"w=%d: A[%d]=%.10g+i%.10g, B[%d]=%.10g+i%.10g, C[%d]=%.10g+i%.10g\n",i,j,Ax[j]*Aw[i*2]+Ax[j+Nx]*Aw[i*2+Nt],Ax[j]*Aw[i*2+1]+Ax[j+Nx]*Aw[i*2+Nt+1],j,Bx[j]*BCw[i*2],Bx[j]*BCw[i*2+1],j,Cx[j]*BCw[i*2],Cx[j]*BCw[i*2+1]);
		}
		fclose(fA);
	}
}


extern "C" 
FL_DBL step(int Nt, int Nx, FL_DBL z, bool sv)
{

FL_DBL gpu_time;
FL_DBL overal_time=0;

	int i;

	dim3 grid(-floor(-FL_DBL(Nt/2)/FL_DBL(Bs)), Nx/BsX, 1);
	dim3 threads(Bs, 1, 1);
	dim3 simple_grid(-floor(-FL_DBL(Nt/2)/FL_DBL(Bs)), 1, 1);
	dim3 simple_threads(Bs, 1, 1);
	//  dim3 forward_grid( -floor(-FL_DBL((Nt-2)/2)/FL_DBL(256)), Nx, 1);  !!! DOESNOT WORK WITH THIS GEOMETRY
	//  dim3 forward_threads(256, 1, 1);  !!! DOESNOT WORK WITH THIS GEOMETRY

	dim3 forward_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(16)), Nx/16, 1);
	dim3 forward_threads(16, 16, 1);

	dim3 Nl_grid(Nx/64, 1, 1);
	dim3 Nl_threads(64, 1, 1);

	dim3 NlS_grid(Nx/NlS_BSX, 1, 1);
	dim3 NlS_threads(NlS_BSX, 1, 1);



  unsigned int timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));

	SWEEP_1_FIXED_Kernel<<< grid, threads, smem>>>(dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_spec, dev_spec_tmp, dev_bb, dev_cc, Nx, Nt/2, BsX, Bs, dev_spec, dev_a2, dev_b2, dev_c2, dev_y2);
  cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("bd_I:%.2f\t",gpu_time);
		#endif

		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));

  simpleSWEEPKernel<<< simple_grid, simple_threads >>>(dev_a2, dev_b2, dev_c2, dev_y2, dev_x2, 2*Nx/BsX, Nt/2, dev_y2);
  cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("bd_II:%.2f\t",gpu_time);
		#endif
  		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));

  SWEEP_2_Kernel<<< grid, threads >>>(dev_spec_tmp, dev_bb, dev_cc, dev_spec, Nx, Nt/2, BsX, Bs, dev_x2);
  cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("bd_III:%.2f\t",gpu_time);
		#endif

		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));

  	forward_Kernel<<< forward_grid, forward_threads>>>(dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_spec_tmp, dev_spec, Nx, Nt/2);
	  cudaThreadSynchronize();

		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("forward:%.2f\t",gpu_time);
		#endif

//restore zero frequency
	restore_zero_freq_Kernel<<<-floor(-Nx/Bs), Bs>>>(dev_spec_zero_freq, dev_spec, Nt, Nx);
	cudaThreadSynchronize();

//add plasma responce
	//goto back from Fourier (using GPU)
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	cuFFT_B(planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp);
	
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("FT_I:%.2f\t",gpu_time);
		#endif
	//integ nonlinearity
  	dim3 Nl_fast_grid(Nx, 1, 1);
	  dim3 Nl_fast_threads(Nl_fast_NTH, 1, 1);

		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
		
	find_Nl_fast_Kernel<<<Nl_fast_grid, Nl_fast_threads>>>(dev_spec_tmp,Nl_fast_NTH,Nt-2,Nx,_dt,(sv==true)?dev_plasma1:NULL,(sv==true)?dev_THz_src:NULL);

//	find_Nl_SHARED_Kernel<<<NlS_grid, NlS_threads>>>(dev_spec_tmp,Nt-2,Nx,_dt);
//	find_Nl_Kernel<<<Nl_grid, Nl_threads>>>(dev_spec_tmp,Nt-2,Nx,_dt);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("Nl:%.2f\t",gpu_time);
		#endif
	//goto Fourier (need new array for FTT(NL) term)	

	//goto back from Fourier (using GPU)
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	cuFFT_F(planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_Nl_spec);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("FT_F:%.2f\t",gpu_time);	
		#endif
	//add dev_Nl_spec
  	dim3 addNl_fast_grid(Nt/2/Nl_fast_NTH, Nx, 1);
	  dim3 addNl_fast_threads(Nl_fast_NTH, 1, 1);
		#define Nbl_addNl	(240*16)
  	dim3 addNl_ffast_grid(Nbl_addNl, 1, 1);
	  dim3 addNl_ffast_threads(Nl_fast_NTH, 1, 1);

		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	//add_Nl_fast_Kernel<<<addNl_fast_grid, addNl_fast_threads>>>(dev_spec,dev_Nl_spec,Nt/2,Nx,_dz);
	add_Nl_Kernel<<<forward_grid, forward_threads>>>(dev_spec,dev_Nl_spec,Nt/2,Nx,_dz);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("addNl:%.2f\t",gpu_time);
		#endif
// SEC ORDER PLASMA
//add plasma responce
	//goto back from Fourier (using GPU)
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	cuFFT_B(planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("FT_I:%.2f\t",gpu_time);
		#endif
	//integ nonlinearity
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	//find_Nl_Kernel<<<Nl_grid, Nl_threads>>>(dev_spec_tmp,Nt-2,Nx,_dt);
	find_Nl_fast_Kernel<<<Nl_fast_grid, Nl_fast_threads>>>(dev_spec_tmp,Nl_fast_NTH,Nt-2,Nx,_dt,(sv==true)?dev_plasma2:NULL,NULL);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("Nl:%.2f\t",gpu_time);
		#endif
	//goto Fourier (need new array for FTT(NL) term)	

	//goto back from Fourier (using GPU)
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	cuFFT_F(planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_Nl_spec_2);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("FT_F:%.2f\t",gpu_time);	
		#endif
	//add dev_Nl_spec
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	add_Nl_R2_Kernel<<<forward_grid, forward_threads>>>(dev_spec,dev_Nl_spec,dev_Nl_spec_2,Nt/2,Nx,_dz);
	cudaThreadSynchronize();
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		#ifdef PRINT
		printf("addNl:%.2f\t",gpu_time);
		#endif

		#ifdef PRINT
		printf("\n");
		#endif

/*
	int bsx=256;
	int nth=32;
//  dim3 SMPL_grid((Nt-2)*Nx/bsx, 1, 1);
//  dim3 SMPL_threads(nth, 1, 1);
  dim3 SMPL_grid((Nt-2)*Nx/512, 1, 1);
  dim3 SMPL_threads(512, 1, 1);
	//TEST SIMPLE CP
		timer=0;
		cutilCheckError( cutCreateTimer( &timer));
		cutilCheckError( cutStartTimer( timer));
	for(i=0;i<100;i++)
	{
//		simple_cpy_Kernel<<<SMPL_grid, SMPL_threads>>>(dev_Nl_spec,dev_spec,0.1,0.01,bsx/nth,bsx,nth);
		Simple_cpy_Kernel<<<SMPL_grid, SMPL_threads>>>(dev_Nl_spec,dev_spec,0.1,0.01,bsx/nth,bsx,nth);
		cudaThreadSynchronize();
	}
		cutilCheckError(cutStopTimer(timer));
		cutilCheckMsg("Kernel execution failed");
		gpu_time=cutGetTimerValue( timer);
		cutilCheckError( cutDeleteTimer( timer));
		overal_time+=gpu_time;
		printf("SIMPLE : %f ms \n",gpu_time/100);
*/

	return overal_time;
}

extern "C"
void get_spec(FL_DBL* spec ,int Nt, int Nx)	{dev_d2h(dev_spec,spec,Nx*Nt);}
extern "C"
void get_field(FL_DBL* field ,int Nt, int Nx)
{
	//goto back from Fourier (using GPU)
	cuFFT_B(planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp);
	cudaThreadSynchronize();
	cutilSafeCall( cudaMemcpy( field, dev_spec_tmp, Nx*Nt*sizeof(FL_DBL), cudaMemcpyDeviceToHost));
	cudaThreadSynchronize();
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<(Nt-2);j++) field[j+i*(Nt-2)]=field[j+i*(Nt-2)]/double(Nt-2);
	}

}

extern "C"
void get_density(FL_DBL* field ,int Nt, int Nx)
{
	//goto back from Fourier (using GPU)
	cutilSafeCall( cudaMemcpy( field, dev_plasma1, Nx*Nt*sizeof(FL_DBL), cudaMemcpyDeviceToHost));
	cudaThreadSynchronize();
}

__global__ void
filter_Kernel(FL_DBL* spc, int Nx, int Nt, FL_DBL dt)
{
//because the method of invertion of A matrix we need to specify some value for omega when i=0 (zero frequency)
//that is why we need to restore the value of zero frequency on each step
	const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
	const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if((tidx<Nt)&(tidy<Nx))
	{
		const FL_DBL norma=ONE/FL_DBL(Nt);
		FL_DBL freq=FL_DBL(tidx)*((FL_DBL)2.0)*Pi/(dt*FL_DBL(Nt*2));
		freq*=ONE/((FL_DBL)0.25);
		FL_DBL ex=exp(-freq*freq);
		spc[tidx*2+tidy*Nt*2]*=ex*norma;
		spc[tidx*2+1+tidy*Nt*2]*=ex*norma;
	}
}

__global__ void
mask_Kernel(FL_DBL* spc, int Nx, int Nt, FL_DBL dt)
{
//because the method of invertion of A matrix we need to specify some value for omega when i=0 (zero frequency)
//that is why we need to restore the value of zero frequency on each step
	const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
	const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if((tidx<Nt)&(tidy<Nx))
	{
		FL_DBL t=FL_DBL(tidx)*dt;
		FL_DBL mask=tanh(t*(FL_DBL(Nt-1)*dt-t)/((FL_DBL)10000.0));
		//spc[tidx+tidy*Nt]+=0.0000001;
		spc[tidx+tidy*Nt]*=mask;
	}
}

extern "C"
void get_THz_source(FL_DBL* field ,int Nt, int Nx)
{
	//filter source
	//switch to spectrum
	cuFFT_F(planF, (cuFFT_rt*)dev_THz_src, (cuFFT_ct*)dev_Nl_spec);
	//apply filter and integration
	dim3 filter_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(16)), Nx/16, 1);
	dim3 filter_threads(16, 16, 1);
	filter_Kernel<<<filter_grid,filter_threads>>>(dev_Nl_spec,Nx,Nt/2,_dt);
	//switch to real space
	cuFFT_B(planI, (cuFFT_ct*)dev_Nl_spec, (cuFFT_rt*)dev_THz_src);
	//goto back from Fourier (using GPU)
	dim3 mask_grid( -floor(-FL_DBL(Nt-2)/FL_DBL(16)), Nx/16, 1);
	dim3 mask_threads(16, 16, 1);
	mask_Kernel<<<mask_grid,mask_threads>>>(dev_THz_src,Nx,Nt-2,_dt);

//	filter_Kernel<<<forward_grid,forward_threads>>>(dev_THz_src,Nx,Nt/2,_dt);
	cutilSafeCall( cudaMemcpy( field, dev_THz_src, Nx*Nt*sizeof(FL_DBL), cudaMemcpyDeviceToHost));
	cudaThreadSynchronize();
}
