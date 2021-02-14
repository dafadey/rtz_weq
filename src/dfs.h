//#define USEDOUBLE

#ifdef USEDOUBLE
	#define DATATYPE "double"
	#define cuFFT_F cufftExecD2Z	
	#define cuFFT_B	cufftExecZ2D
	#define cuFFT_FP	CUFFT_D2Z
	#define cuFFT_BP	CUFFT_Z2D
	#define cuFFT_rt cufftDoubleReal
	#define cuFFT_ct cufftDoubleComplex
	#define TYP 
	#define ONE 1.0
	#define ZERO 0.0
#endif
#ifndef USEDOUBLE
	#define DATATYPE "float"
	#define cuFFT_F cufftExecR2C	
	#define cuFFT_B	cufftExecC2R
	#define cuFFT_FP	CUFFT_R2C
	#define cuFFT_BP	CUFFT_C2R
	#define cuFFT_rt cufftReal
	#define cuFFT_ct cufftComplex
	#define TYP f
	#define ONE 1.0f
	#define ZERO 0.0f
#endif

#define NNX (512)
#define NNT (512*2+2)
//#define NNX (1024*2)
//#define NNT (512+2)

#define Lx 5000.0F
#define T 250.0F
#define Pi 3.1415926535897932384626433832795F

#define NET
#define	SCR	1

//#define CPU_TRANSCEND

#define NlS_BSX	32

#define Nl_fast_NTH	32
#define Nl_fast_SMyS	16
#define Nl_fast_SMxS	16
#define Nl_fast_SMxSa	17



#define SAVE
//#define DRAW

//---------------------------

#ifdef USEDOUBLE
	#define FL_DBL double
#endif
#ifndef USEDOUBLE
	#define FL_DBL float
#endif

#define NOFA 2

#define REP 400


#define MAX_THREADS_PER_BLOCK 32
#define MIN_BLOCKS_PER_MP 8
