#include <math_constants.h>
#include <cutil_inline.h>
#include "dfs.h"
#include "sweep_kernel.cu"
#include "nl_kernel.cu"
const bool PRINT = false;
const bool KERR_ON = true;
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
FL_DBL* dev_Nl_spec; // for plasma
FL_DBL* dev_Nl_spec_2; // for plasma sec order
FL_DBL* dev_NlKerr_spec = nullptr; // for kerr nonlinearity spectra
FL_DBL* dev_NlKerr_spec_2 = nullptr; // for kerr nonlinearity spectra sec order
FL_DBL* dev_NlKerr = nullptr; // for kerr nonlinearity
FL_DBL* dev_NlKerr_2 = nullptr; // for kerr nonlinearity sec order
// we need different arrays since plasma and kerr nonlinearity cannot be added in real space
// since it is nE + d^2/dtau^2 E^3 there is a need for individual storage for both terms
// JUST to calculate second derivtive with FFT
// it is quite important since E^3 has hight order harmonics
// yep n is integrated in real space, I know
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

FL_DBL* fraction_params_pack;
FL_DBL* dev_fraction_params_pack;

extern "C"
void get_field(FL_DBL*, int, int);

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

#define CALL_KERNEL(kernel, msg, grid, threads, smem, ...)	timer=0; \
															cutilCheckError( cutCreateTimer( &timer)); \
															cutilCheckError( cutStartTimer( timer)); \
															if(!smem) \
																kernel <<< grid, threads >>> (__VA_ARGS__); \
															else \
																kernel <<< grid, threads, smem >>> (__VA_ARGS__); \
															cudaThreadSynchronize(); \
															cutilCheckError(cutStopTimer(timer)); \
															cutilCheckMsg("Kernel execution failed"); \
															gpu_time=cutGetTimerValue( timer); \
															cutilCheckError( cutDeleteTimer( timer)); \
															overal_time+=gpu_time; \
															if(PRINT) printf("%s:%.2f\t",msg,gpu_time);

#define CALL_CUFFT(CUFFT_FUNC, msg, PLAN, ...)	timer=0; \
												cutilCheckError( cutCreateTimer( &timer)); \
												cutilCheckError( cutStartTimer( timer)); \
												CUFFT_FUNC(PLAN, __VA_ARGS__); \
												cudaThreadSynchronize(); \
												cutilCheckError(cutStopTimer(timer)); \
												cutilCheckMsg("Kernel execution failed"); \
												gpu_time=cutGetTimerValue( timer); \
												cutilCheckError( cutDeleteTimer( timer)); \
												overal_time+=gpu_time; \
												if(PRINT) printf("%s:%.2f\t",msg,gpu_time);


extern "C"
double RefIndexAir_um(double lambda) {
	const double lambda_crit = 0.78/5.;
	double lambda_1 = lambda > lambda_crit ? 1./lambda : 1./lambda_crit;
	return 0.05792105/(238.0185-lambda_1*lambda_1) + 0.00167917/(57.362-lambda_1*lambda_1);
}

#define DISPERSION_TEST .0 //.00031

extern "C"
FL_DBL LaplassianA(int i, FL_DBL dr) {
	return -(2. + 1./(FL_DBL(i) + .5))/(dr*dr);
}

extern "C"
FL_DBL LaplassianB(int i, FL_DBL dr) {
	return 1./(dr*dr);
}

extern "C"
FL_DBL LaplassianC(int i, FL_DBL dr) {
	return ((i == 0 ? 2. : 1.) + 1./(FL_DBL(i) + .5))/(dr*dr);
}

extern "C"
FL_DBL LaplassianMetric(int i, FL_DBL dr) {
	return (i==0 ? 3./8. : (FL_DBL(i) + .5)) * dr * 2. * M_PI;
}

extern "C" 
void cuda_load(FL_DBL* spec, int Nt, int Nx, FL_DBL dt, FL_DBL dx, FL_DBL dz, FL_DBL airDensity)
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
	
//	for(i=0;i<Nx;i++) Ax[i+Nx]=(((i==0)||(i==Nx-1))?1.0:2.0)*dz/dx/dx+1.0/(double(i)+0.5)/dx/dx*dz; // Laplassian
	for(i=0;i<Nx;i++) Ax[i+Nx] = -LaplassianA(i,dx) * dz;//2.*dz/dx/dx + 1./(double(i) + .5)/dx/dx*dz; // Laplassian
	if(NOFA>=3)
		for(i=0;i<Nx;i++) Ax[i+2*Nx] = dz; // Dispersion
	
	for(i=0;i<Nx;i++) Bx[i] = -LaplassianB(i,dx)*dz;//-1.0*dz/dx/dx; // Laplassian
	//Bx[Nx-1]=-dz/dx/dx;  //non hermit
	Bx[0]=0.0F;   // BECAUSE FAST SWEEP HAS A BUG THIS SHOULD BE 0.0! ALWAYS!!!
	//for(i=0;i<Nx;i++) Cx[i]=-1.0*dz/dx/dx-1.0/(double(i)+0.5)/dx/dx*dz; // Laplassian
	for(i=0;i<Nx;i++) Cx[i] = -LaplassianC(i,dx)*dz;//-(((i==0)||(i==Nx-1))?2.:1.)*dz/dx/dx - 1./(double(i) + .5)/dx/dx*dz; // Laplassian
	//Cx[0]=-2.0F*dz/dx/dx-1.0F/(0.5F)/dx/dx*dz;  // non hermit
	//Ax[0]=0.5*Ax[Nx-1]; // hermit   USE HERMIT MATRIX TO AVOID PROBLEMS WITH ENERGY CONSERVATION CONCERNED WITH BOUNDARY EFFECTS
	Cx[Nx-1]=0.0;   // BECAUSE FAST SWEEP HAS A BUG THIS SHOULD BE 0.0! ALWAYS!!!
		
	for(i=0;i<Nt/2;i++)
	{
		//NOFA - N of A. here we have 3 different A
		Aw[i*2] = 1.0; // Identity
		Aw[i*2+1] = 0.0; // Identity


		FL_DBL omega = double(i) * 2. * Pi / T;
		FL_DBL lambda_vac = i==0 ? 1.0 : 0.78 / omega; // um
		FL_DBL omega_1 = i==0 ? 1.0 : 1./omega*(.5 - .5 * tanh(-(2. * Pi / 250.*double(i)-.5)/.1)); // be carefull with i=0. FAST SWEEP will not work with omega=1 <=> A=I YES! FAST SWEEP can't invert identity matrix.

		Aw[i*2+Nt] = 0.0; // Laplassian
		Aw[i*2+1+Nt] = omega_1; // Laplassian

		if(NOFA>=3)
		{
			Aw[i*2+2*Nt] = .0; // Dispersion
			Aw[i*2+1+2*Nt] = DISPERSION_TEST * omega_1 - RefIndexAir_um(lambda_vac) * omega; // Dispersion
		}
		
		BCw[i*2] = 0.0; // Laplassian
		BCw[i*2+1] = omega_1; // Laplassian


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
	if(KERR_ON) {
		dev_alloc((void**) &dev_NlKerr_spec,Nx*Nt);
		dev_alloc((void**) &dev_NlKerr_spec_2,Nx*Nt);
		dev_alloc((void**) &dev_NlKerr,Nx*Nt);
		dev_alloc((void**) &dev_NlKerr_2,Nx*Nt);
	}

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

	dev_alloc((void**) &dev_fraction_params_pack, fractions_count * fractions_stride);
	fraction_params_pack = new FL_DBL[fractions_count * fractions_stride];
	//set fractions parameters:
	//O2
	fraction_params_pack[0] = airDensity * .2;
	fraction_params_pack[1] = .557;
	fraction_params_pack[2] = .123;
	fraction_params_pack[3] = 6.02;
  //N2
	fraction_params_pack[0+fractions_stride] = airDensity * .8;
	fraction_params_pack[1+fractions_stride] = .818;
	fraction_params_pack[2+fractions_stride] = .868;
	fraction_params_pack[3+fractions_stride] = 9.69;

	dev_h2d(fraction_params_pack, dev_fraction_params_pack, fractions_count * fractions_stride);

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
	/*
	{
		FILE* fA=fopen("Atest","w");
		for(int i=0;i<Nt/2;i++)
		{
			for(int j=0;j<Nx;j++) fprintf(fA,"w=%d: A[%d]=%.10g+i%.10g, B[%d]=%.10g+i%.10g, C[%d]=%.10g+i%.10g\n",i,j,Ax[j]*Aw[i*2]+Ax[j+Nx]*Aw[i*2+Nt],Ax[j]*Aw[i*2+1]+Ax[j+Nx]*Aw[i*2+Nt+1],j,Bx[j]*BCw[i*2],Bx[j]*BCw[i*2+1],j,Cx[j]*BCw[i*2],Cx[j]*BCw[i*2+1]);
		}
		fclose(fA);
	}*/
	
	#define printKernelData(X) { cudaFuncAttributes attrib; \
															 cudaFuncGetAttributes(&attrib, X); \
															 printf("kernel: %s==========================================\n", #X); \
															 printf("\tbinaryVersion: %d\n", attrib.binaryVersion); \
															 printf("\tcacheModeCA: %d\n", attrib.cacheModeCA); \
															 printf("\tconstSizeBytes: %lu\n", attrib.constSizeBytes); \
															 printf("\tlocalSizeBytes: %lu\n", attrib.localSizeBytes); \
															 printf("\tmaxDynamicSharedSizeBytes: %d\n", attrib.maxDynamicSharedSizeBytes); \
															 printf("\tmaxThreadsPerBlock: %d\n", attrib.maxThreadsPerBlock); \
															 printf("\tnumRegs: %d\n", attrib.numRegs); \
															 printf("\tpreferredShmemCarveout: %d\n", attrib.preferredShmemCarveout); \
															 printf("\tptxVersion: %d\n", attrib.ptxVersion); \
															 printf("\tsharedSizeBytes: %lu\n", attrib.sharedSizeBytes); \
															 printf("====================================================\n");}
	
	cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	printKernelData(find_Nl_fractions_fast_Kernel);
	printKernelData(add_Nl_Kernel);
	printKernelData(add_Nl_R2_Kernel);
	printKernelData(SWEEP_1_FIXED_Kernel);
	printKernelData(simpleSWEEPKernel);
	printKernelData(SWEEP_2_Kernel);
	printKernelData(forward_Kernel);
}


extern "C" 
FL_DBL step(int Nt, int Nx, FL_DBL z, FL_DBL T, bool sv, FL_DBL airDensity, FL_DBL Kerr, FL_DBL tunnelingExponent, FL_DBL tunnelingFactor)
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

  dim3 forward_threads(16, 16, 1);
  dim3 forward_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(forward_threads.x)), Nx/forward_threads.y, 1);

  //dim3 Nl_threads(64, 1, 1);
  //dim3 Nl_grid(Nx/Nl_threads.x, 1, 1);

  //dim3 NlS_threads(NlS_BSX, 1, 1);
  //dim3 NlS_grid(Nx/NlS_threads.x, 1, 1);

  dim3 addNl_threads(64, 1, 1);
  dim3 addNl_grid( -floor(-FL_DBL(Nt/2)/FL_DBL(addNl_threads.x)), Nx/addNl_BSY, 1);
  
  unsigned int timer=0;

  CALL_KERNEL(SWEEP_1_FIXED_Kernel, "bd_I", grid, threads, smem, dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_spec, dev_spec_tmp, dev_bb, dev_cc, Nx, Nt/2, BsX, Bs, dev_spec, dev_a2, dev_b2, dev_c2, dev_y2)
	
  CALL_KERNEL(simpleSWEEPKernel,"bd_II", simple_grid, simple_threads, 0, dev_a2, dev_b2, dev_c2, dev_y2, dev_x2, 2*Nx/BsX, Nt/2, dev_y2)

  CALL_KERNEL(SWEEP_2_Kernel, "bd_III", grid, threads, 0, dev_spec_tmp, dev_bb, dev_cc, dev_spec, Nx, Nt/2, BsX, Bs, dev_x2)

  CALL_KERNEL(forward_Kernel, "fwd", forward_grid, forward_threads, 0, dev_Ax, dev_Bx, dev_Cx, dev_Aw, dev_BCw, dev_BCw, dev_spec_tmp, dev_spec, Nx, Nt/2)

  //restore zero frequency
  //CALL_KERNEL(restore_zero_freq_Kernel, "rz", -floor(-Nx/Bs), Bs, 0, dev_spec_zero_freq, dev_spec, Nt, Nx)

  //add plasma responce
  //goto back from Fourier (using GPU)
  CALL_CUFFT(cuFFT_B, "FT_I", planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp)
  //integ nonlinearity
#ifndef FRACTIONS
  dim3 Nl_fast_grid(Nx, 1, 1);
  dim3 Nl_fast_threads(Nl_fast_NTH, 1, 1);

  CALL_KERNEL(find_Nl_fast_Kernel,"Nl",Nl_fast_grid, Nl_fast_threads, 0, dev_spec_tmp, dev_NlKerr, Nl_fast_NTH, Nt-2, Nx, _dt, (sv==true)?dev_plasma1:NULL,(sv==true)?dev_THz_src:NULL, airDensity , Kerr, tunnelingExponent, tunnelingFactor)
  //goto Fourier (need new array for FTT(NL) term)	
#else
  dim3 Nl_fast_grid(Nx, 1, 1);
  dim3 Nl_fast_threads(Nl_fractions_nth, 1, 1);

	CALL_KERNEL(find_Nl_fractions_fast_Kernel, "Nlf", Nl_fast_grid, Nl_fast_threads, 0, dev_spec_tmp, dev_NlKerr, Nt-2, Nx, _dt, Kerr, dev_fraction_params_pack, (sv==true)?dev_plasma1:NULL)
#endif

  //goto back from Fourier (using GPU)
  CALL_CUFFT(cuFFT_F, "FT_F", planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_Nl_spec)
  if(dev_NlKerr_spec) {
		CALL_CUFFT(cuFFT_F, "FT_F", planF, (cuFFT_rt*)dev_NlKerr, (cuFFT_ct*)dev_NlKerr_spec)
	}
  //add dev_Nl_spec
  dim3 addNl_fast_grid(Nt/2/Nl_fast_NTH, Nx, 1);
  dim3 addNl_fast_threads(Nl_fast_NTH, 1, 1);
  #define Nbl_addNl	(240*16)
  dim3 addNl_ffast_grid(Nbl_addNl, 1, 1);
  dim3 addNl_ffast_threads(Nl_fast_NTH, 1, 1);

  CALL_KERNEL(add_Nl_Kernel, "addNl", addNl_grid, addNl_threads, 0, dev_spec,dev_Nl_spec,dev_NlKerr_spec,Nt/2,Nx,_dz, T)

// SEC ORDER PLASMA
//add plasma responce
//goto back from Fourier (using GPU)
  CALL_CUFFT(cuFFT_B, "FT_I", planI, (cuFFT_ct*)dev_spec, (cuFFT_rt*)dev_spec_tmp)
  //integ nonlinearity
#ifndef FRACTIONS
  CALL_KERNEL(find_Nl_fast_Kernel, "Nl", Nl_fast_grid, Nl_fast_threads, 0, dev_spec_tmp, dev_NlKerr_2, Nl_fast_NTH, Nt-2, Nx, _dt, (sv==true)?dev_plasma2:NULL, NULL, airDensity , Kerr, tunnelingExponent, tunnelingFactor)
  //goto Fourier (need new array for FTT(NL) term)	
#else
	CALL_KERNEL(find_Nl_fractions_fast_Kernel, "Nlf", Nl_fast_grid, Nl_fast_threads, 0, dev_spec_tmp, dev_NlKerr_2, Nt-2, Nx, _dt, Kerr, dev_fraction_params_pack, (sv==true)?dev_plasma1:NULL)
#endif

  //goto back from Fourier (using GPU)
  CALL_CUFFT(cuFFT_F, "FT_F", planF, (cuFFT_rt*)dev_spec_tmp, (cuFFT_ct*)dev_Nl_spec_2)
  if(dev_NlKerr_spec_2) {
		CALL_CUFFT(cuFFT_F, "FT_F", planF, (cuFFT_rt*)dev_NlKerr_2, (cuFFT_ct*)dev_NlKerr_spec_2)
  }
  //add dev_Nl_spec
  CALL_KERNEL(add_Nl_R2_Kernel, "addNl", addNl_grid, addNl_threads, 0, dev_spec,dev_Nl_spec,dev_Nl_spec_2,dev_NlKerr_spec,dev_NlKerr_spec_2,Nt/2,Nx,_dz, T)

  if(PRINT)
    printf("\n");

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
		for(int j=0;j<(Nt-2);j++) field[j+i*(Nt-2)] /= double(Nt-2);
	}

}

extern "C"
void get_density(FL_DBL* field ,int Nt, int Nx)
{
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
