__global__ void
add_Nl_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, int Nt, int Nx, FL_DBL dz)
{
  const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
 	const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if(tidx<Nt)
	{

		#define simple_ADD //simple is faster

		#ifdef simple_ADD
		FL_DBL omega=(((tidx==0)||(tidx==Nt-1))?0.0F:T/2.0F/Pi/FL_DBL(tidx))*dz*(0.5F-0.5F*tanh(-(2.0F*Pi/250.0F*FL_DBL(tidx)-0.5F)/0.1F))*0.5F;
// do not forget to delete this fuckin 0.5! it here because in CPU code there is extra dz. and it is wrog and result depends on dz
// here nothing should not depend in dz but 0.5 not logical here. just for CPU=GPU test.

		spc[tidx*2+tidy*Nt*2]+=spc_Nl[tidx*2+1+tidy*Nt*2]*omega;
		spc[tidx*2+1+tidy*Nt*2]+=-spc_Nl[tidx*2+tidy*Nt*2]*omega;

		#endif
		#ifndef simple_ADD
		const FL_DBL sgn[]={ONE,-ONE};
		const int adr[]={1,-1};
		FL_DBL omega=1.0F/(2.0F*Pi/T)*dz/(((tidx>>1)>0)?FL_DBL((tidx>>1)):sgn[0]);
		spc[tidx+tidy*Nt*2]+=sgn[tidx&1]*spc_Nl[tidx+adr[tidx&1]+tidy*Nt*2]*omega;
		omega=1.0F/(2.0F*Pi/T)*dz/FL_DBL(((tidx+Nt)>>1));
		spc[tidx+Nt+tidy*Nt*2]+=sgn[tidx&1]*spc_Nl[tidx+Nt+adr[tidx&1]+tidy*Nt*2]*omega;



		#endif

	}
}

__global__ void
add_Nl_R2_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, FL_DBL* spc_Nl2, int Nt, int Nx, FL_DBL dz)
{
  const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if(tidx<Nt)
	{
		FL_DBL omega=(((tidx==0)||(tidx==Nt-1))?0.0:T/2.0F/Pi/FL_DBL(tidx))*dz*0.5F*(0.5F-0.5F*tanh(-(2.0F*Pi/250.0F*FL_DBL(tidx)-0.5F)/0.1F));
		spc[tidx*2+tidy*Nt*2]+=(spc_Nl2[tidx*2+1+tidy*Nt*2]-spc_Nl[tidx*2+1+tidy*Nt*2])*omega;
		spc[tidx*2+1+tidy*Nt*2]+=(-spc_Nl2[tidx*2+tidy*Nt*2]+spc_Nl[tidx*2+tidy*Nt*2])*omega;
	}
}


__global__ void
find_Nl_fast_Kernel(FL_DBL* spc, int nth, int Nt, int Nx, FL_DBL dt, FL_DBL* plasma, FL_DBL* THz_src)
{
	#define find_Nl_fast_2way

//	#define find_Nl_fast_nionly
	//consts
	
	#define E0 10.0F
	#define N0 100000.0F
	#define BASE0 0.000001F

	const FL_DBL norma=ONE/FL_DBL(Nt);
	const int PBS=256; // block size under processing
	const int SMxS=Nl_fast_SMxS; // shared X size 16+1 in order to prevent bank conflicts
	const int SMxSa=Nl_fast_SMxSa; // shared X size 16+1 in order to prevent bank conflicts
	const int SMyS=Nl_fast_SMyS; // shared Y size = PBS/(SMxS-1)
	__shared__ FL_DBL sdata[Nl_fast_SMxSa*Nl_fast_SMyS+1]; // contains all data of block //last value @ SMxS*SMyS is curent density (density per block)
	register int i,j;

	if (threadIdx.x==0) sdata[SMxSa*SMyS]=ZERO;
	syncthreads();

	for(j=0;j<Nt/PBS;j++)
	{
#ifndef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=1)//1 way unrolled read to SM
		{
			sdata[threadIdx.x+(threadIdx.x/16)+i*SMxSa*2]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt];
		}
#endif
#ifdef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=2)//2 way unrolled read to SM
		{
			sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt];
			sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2+SMxSa+SMxSa]=spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt];
		}
#endif
		syncthreads();

		/* NOTE threadIdx.x+(threadIdx.x>>SHFT):
		0->0
		1->1
		...
		15->15
		16->17
		18->19
		...
		31->32
		for no bank conflict sequtial reading later
		*/



		//subintegration (16steps each thread)
		register FL_DBL ni;
		register FL_DBL s;		
		if((threadIdx.x&1)==0)
		{
			int tid=(threadIdx.x>>1); // tid=0..15
			ni=ZERO;
			for(i=0;i<SMxS;i++)
			{

				s=sdata[i+tid*SMxSa]*norma;
				//s=s*s*dt;
				//s=N0*exp(-E0/(fabs(s)+BASE0))*dt;
				s=(fabs(s)>0.000001F)?pow(s*s+0.0000000000000001F,-0.3125F)*exp(-1.0F/(fabs(s)+0.000000000000000000001F)):0.0F;
				ni+=s*dt;
				sdata[i+tid*SMxSa]=ni;

			}
		}
		syncthreads();

		//2 collect ni steps in 17-th s elements of each sm rows


		if(threadIdx.x==0)
		{
			ni=sdata[SMxS-1];
			for(i=1;i<SMxS;i++)
			{
				sdata[i*SMxSa-1]=ni;
				ni+=sdata[i*SMxSa+SMxS-1];
			}
			sdata[SMxS*SMxSa-1]=ni;
		}
		syncthreads();


		//3 add step function
		if((threadIdx.x&1)==0)
		{
			int tid=(threadIdx.x>>1); // tid=0..15
			if(tid>0)
			{
				s=sdata[tid*SMxSa-1];
				for(i=0;i<SMxS;i++)	sdata[i+tid*SMxSa]+=s;
			}
		}

		syncthreads();


	//write back to global
		#ifndef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=1)// 1 way write
		{
			#ifndef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=norma*spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]*(sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
			if(THz_src!=NULL) THz_src[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]/norma;
			if(plasma!=NULL) plasma[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=(sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
			#endif
			#ifdef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2];
			#endif
		}
		#endif
		#ifdef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=2) //2 way unrolled r/w to SM
		{
			#ifndef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=norma*spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]*(sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
			spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=norma*spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]*(sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2]+sdata[SMxSa*SMyS]);

			if(THz_src!=NULL) THz_src[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]/norma;
			if(THz_src!=NULL) THz_src[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]/norma;

			if(plasma!=NULL) plasma[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=(sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
			if(plasma!=NULL) plasma[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=(sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2]+sdata[SMxSa*SMyS]);
			#endif
			#ifdef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2];
			spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2];
			#endif
		}
#endif

		syncthreads();

		//dont forget to store curent density
		if (threadIdx.x==0) sdata[SMxSa*SMyS]+=sdata[SMxS*SMxSa-1];
		syncthreads();
	}
}

__global__ void
find_Nl_SHARED_Kernel(FL_DBL* spc, int Nt, int Nx, FL_DBL dt)
{
	//const FL_DBL E0=10.0;
	//const FL_DBL N0=100000.0;
	//const FL_DBL BASE0=0.000001;
	const int SMxS=Nl_fast_SMxS; // shared X size 16+1 in order to prevent bank conflicts
	const int SMxSa=Nl_fast_SMxSa; // shared X size 16+1 in order to prevent bank conflicts
//	const int SMyS=Nl_fast_SMyS; // shared X size 16+1 in order to prevent bank conflicts
	__shared__ FL_DBL sdata[NlS_BSX*NlS_BSX+NlS_BSX];
//	extern __shared__ FL_DBL sdata[];
	register FL_DBL ni=0;
	syncthreads();
	register int i,j;
	for(j=0;j<Nt/NlS_BSX;j++)
	{
		for(i=0;i<NlS_BSX;i++){
//		 sdata[(threadIdx.x&15)+i*SMxSa+(threadIdx.x/16)*SMxSa*SMyS]=spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt];
			sdata[threadIdx.x+(threadIdx.x/16)+i*SMxSa*2]=spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt];
		}
		syncthreads();

/*
		if(threadIdx.x&1==0)
		{
			register int tid=(threadIdx.x>>1);
			register FL_DBL s;
			register int j;
			for(j=0;j<4;j++)
			{
				for(i=0;i<SMxS;i++)
				{
					s=sdata[i+SMxSa*tid+j*SMxSa*SMxS*j]/FL_DBL(Nt);
					//ni+=s*s*s*s*s*s*s*s*s*s*dt;
		//			ni+=N0*exp(-E0/(fabs(s)+BASE0))*dt;
					ni+=s*s*dt;
					sdata[i+SMxSa*tid+j*SMxSa*SMxS*j]=ni*s;

		//					sdata[i+SMxSa*tid+j*SMxSa*SMxS*j]+=1.0;
				}
			}
		}
*/
		syncthreads();
		for(i=0;i<NlS_BSX;i++){
//			spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt]=sdata[(threadIdx.x&15)+i*SMxSa+(threadIdx.x/16)*SMxSa*SMyS];
			spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt]=sdata[threadIdx.x+(threadIdx.x/16)+i*SMxSa*2];
		}

//		syncthreads();
	}
}


__global__ void
find_Nl_Kernel(FL_DBL* spc, int Nt, int Nx, FL_DBL dt)
{
	//const FL_DBL E0=10.0;
	//const FL_DBL N0=100000.0;
	//const FL_DBL BASE0=0.000001;
  	const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
	register int i;
	register FL_DBL ni=0;
	register FL_DBL s;

	for(i=0;i<Nt;i++)
	{
		s=spc[i+tid*Nt]/FL_DBL(Nt);
		//ni+=s*s*s*s*s*s*s*s*s*s*dt;
		ni+=s*s*dt;
//		ni+=N0*exp(-E0/(fabs(s)+BASE0))*dt;
		spc[i+tid*Nt]=ni*s;
	}
}

__global__ void
add_Nl_fast_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, int Nt, int Nx, FL_DBL dz)
{
  const int tidy = blockIdx.y;
	
	__shared__ FL_DBL sdata[(Nl_fast_NTH<<2)];
	sdata[threadIdx.x]=spc[threadIdx.x+blockDim.x*blockIdx.x*2+tidy*Nt*2];
	sdata[threadIdx.x+blockDim.x]=spc[threadIdx.x+blockDim.x+blockDim.x*blockIdx.x*2+tidy*Nt*2];
	sdata[threadIdx.x+blockDim.x*2]=spc_Nl[threadIdx.x+tidy*Nt*2+blockDim.x*blockIdx.x*2];
	sdata[threadIdx.x+blockDim.x*2+blockDim.x]=spc_Nl[threadIdx.x+blockDim.x+tidy*Nt*2+blockDim.x*blockIdx.x*2];
	syncthreads();



	const FL_DBL sgn[]={1.0F,-1.0F};
	const int adr[]={1,-1};
	register int tidx=threadIdx.x+blockDim.x*blockIdx.x*2;
	register FL_DBL omega=1.0F/(2.0F*Pi/T)*dz/(((tidx>>1)>0)?FL_DBL((tidx>>1)):sgn[0]);
	sdata[threadIdx.x]+=sgn[threadIdx.x&1]*sdata[threadIdx.x+adr[threadIdx.x&1]+blockDim.x*2]*omega;
	omega=1.0F/(2.0F*Pi/T)*dz/FL_DBL(((tidx+blockDim.x)>>1));
	sdata[threadIdx.x+blockDim.x]+=sgn[threadIdx.x&1]*sdata[threadIdx.x+blockDim.x+adr[threadIdx.x&1]+blockDim.x*2]*omega;
	syncthreads();

	spc[threadIdx.x+tidy*Nt*2+blockDim.x*blockIdx.x*2]=sdata[threadIdx.x];
	spc[threadIdx.x+blockDim.x+tidy*Nt*2+blockDim.x*blockIdx.x*2]=sdata[threadIdx.x+blockDim.x];
	syncthreads();
}

#ifdef USEDOUBLE
	#define w_O2(x) ((x)*(x)*1.2)
	#define w_N2(x) ((x)*(x)*1.1)
#else
	#define w_O2(x) ((x)*(x)*1.2f)
	#define w_N2(x) ((x)*(x)*1.1f)
#endif

__global__ void
qalc_ionization_part(FL_DBL* A, FL_DBL* B, FL_DBL* C, int Nt, FL_DBL dt, FL_DBL E_O2, FL_DBL E_N2)
{
	// block size
	//set banksize according to major rev. #
	#if ARCH_MR==1
		#define BANK 16 // G80, GT200
		#ifdef USEDOUBLE
			#define BATCH 64 // portion of data (FL_DBL values) for processing
			#define THS 32 // simultaneously running threads
		#else
			#define BATCH 128 // portion of data (FL_DBL values) for processing
			#define THS 64 // simultaneously running threads
		#endif
	#else // Fermi, Kepler
		#define BANK 32
		#ifdef USEDOUBLE
			#define BATCH 256 // portion of data (FL_DBL values) for processing
		#else
			#define BATCH 512 // portion of data (FL_DBL values) for processing
		#endif
		#define THS 128 // simultaneously running threads
	#endif
	// padded address function
	#define PADDED_ADDR(i) ((i)/BANK+(i))
//	#define PADDED_ADDR(i) (i)
	//smem size:
	#ifdef USEDOUBLE
		#define SMALL 0.0000000001
		#define HALF 0.5
		#define UNO 1.0
		//#define DOS 2.0
		//#define TRES 3.0
	#else
		#define SMALL 0.0000000001f
		#define HALF 0.5f
		#define UNO	1.0f
		//#define DOS 2.0f
		//#define TRES 3.0f
	#endif
	__shared__ FL_DBL smem[(BATCH/BANK+BATCH)*3];
	//register int THS=blockDim.x;
	const int base=blockIdx.x*Nt;
	register int i;
	const int tid=threadIdx.x;
	const int BATCH2=BATCH*2;
	for(i=0;i<Nt;i+=BATCH)
	{
		register int j;


		// read data from spc
/*
		for(j=0;j<BATCH;j+=THS){
			register FL_DBL V=A[base+i+j+tid];
			smem[PADDED_ADDR(j+tid)]=V; // should use registers for ARCH_MR>2 (Kepler+)
			register FL_DBL w_O2_c=w_O2(V);
			register FL_DBL w_N2_c=w_N2(V);
			register FL_DBL dt_2=dt*HALF;
			//C[base+i+j+tid]=(E_O2*w_O2_c+E_N2*w_N2_c)/(V*V+SMALL)*V; // out put to global here for ionization damping member
			smem[PADDED_ADDR(j+tid+BATCH)]=w_O2_c*dt_2;
			smem[PADDED_ADDR(j+tid+BATCH*2)]=w_N2_c*dt_2;
		}
		*/
		{
			register FL_DBL dt_2=dt*HALF;
			register FL_DBL v[(BATCH/THS)];
			register FL_DBL wo2[(BATCH/THS)];
			register FL_DBL wn2[(BATCH/THS)];
			for(j=0;j<BATCH;j+=THS) v[j]=A[base+i+j+tid];
			for(j=0;j<BATCH;j+=THS) wo2[j]=w_O2(v[j])*dt_2;
			for(j=0;j<BATCH;j+=THS) wn2[j]=w_N2(v[j])*dt_2;
				//C[base+i+j+tid]=(E_O2*w_O2_c+E_N2*w_N2_c)/(V*V+SMALL)*V; // out put to global here for ionization damping member
			for(j=0;j<BATCH;j+=THS)	smem[PADDED_ADDR(j+tid)]=-v[j]; // should use registers for ARCH_MR>2 (Kepler+)
			for(j=0;j<BATCH;j+=THS)	smem[PADDED_ADDR(j+tid+BATCH)]=wo2[j];
			for(j=0;j<BATCH;j+=THS)	smem[PADDED_ADDR(j+tid+BATCH*2)]=wn2[j];
		}
		__syncthreads();
		//inv data




		register int invBATCH=BATCH/THS;
/*
		for(j=0;j<BATCH/THS;j++)  // should use registers for ARCH_MR>2 (Kepler+)
		{
		
//			smem[PADDED_ADDR(j+tid*invBATCH)]=-smem[PADDED_ADDR(j+tid*invBATCH)];
//			smem[PADDED_ADDR(j+tid*invBATCH+BATCH)]=TRES*smem[PADDED_ADDR(j+tid*invBATCH)];	
//			smem[PADDED_ADDR(j+tid*invBATCH+BATCH*2)]=DOS*smem[PADDED_ADDR(j+tid*invBATCH)];
		}
*/


		{
			register FL_DBL ws[(BATCH/THS+1)*2];
			register FL_DBL wsn[BATCH/THS*2];
			for(j=0;j<BATCH/THS+1;j++) ws[j]=smem[PADDED_ADDR(j-1+tid*invBATCH+BATCH)]; // exclude this smem use for Kepler! use __shfl!!!
			for(j=0;j<BATCH/THS+1;j++) ws[j+BATCH/THS+1]=smem[PADDED_ADDR(j-1+tid*invBATCH+BATCH*2)];
			
			for(j=0;j<BATCH/THS;j++){
				wsn[j]=(UNO+ws[j])/(UNO-ws[j+1]);
			}
			for(j=0;j<BATCH/THS;j++){
				wsn[j+BATCH/THS]=(UNO+ws[j+BATCH/THS+1])/(UNO-ws[j+1+BATCH/THS+1]);
			}
			__syncthreads();
			for(j=0;j<BATCH/THS;j++) smem[PADDED_ADDR(j+tid*invBATCH+BATCH)]=wsn[j+1];
			for(j=0;j<BATCH/THS;j++) smem[PADDED_ADDR(j+tid*invBATCH+BATCH*2)]=wsn[j+1+BATCH/THS+1];
		}

		
		// I-s calculated
		// calculate partial multiplication


		#define THSM (THS) //threads that perform multiplication
		// STAGE I
		// make simple sequential multiplications

		if(tid<THSM)
		{
		  register int i;
		  register int portion=BATCH/THSM/2;
		  if(BATCH>THSM*2){
		    for(i=1;i<portion;i++)	smem[PADDED_ADDR(i+portion*tid+BATCH)]*=smem[PADDED_ADDR(i-1+portion*tid+BATCH)];
		    for(i=1;i<portion;i++)	smem[PADDED_ADDR(i+portion*tid+BATCH*2)]*=smem[PADDED_ADDR(i-1+portion*tid+BATCH*2)];
		  }
		  __syncthreads();

		  // STAGE II
		  // make log(THSM*2) integrations
		  #define ADDRo2(i) (PADDED_ADDR(portion*(i)+portion-1+BATCH))
		  #define ADDRn2(i) (PADDED_ADDR(portion*(i)+portion-1+BATCH*2))

			
		  register int level=0;
		  for(level=0;(1<<level)<THSM*2;level++)
		  {
			register int mask=(1<<level)-1;
		    register int indx0=((tid&(~mask))<<1)+mask;
		    register FL_DBL V0o2=smem[PADDED_ADDR(indx0+BATCH)];
		    register FL_DBL V0n2=smem[PADDED_ADDR(indx0+BATCH2)];
		    register int indx=indx0+1+(tid&mask);
		    smem[PADDED_ADDR(indx+BATCH)]*=V0o2;
		    smem[PADDED_ADDR(indx+BATCH2)]*=V0n2;
/*
			register FL_DBL V0o2=smem[PADDED_ADDR(THS-1+BATCH)];
		    register FL_DBL V0n2=smem[PADDED_ADDR(THS-1+BATCH*2)];
		    smem[PADDED_ADDR(THS+tid+BATCH)]*=V0o2;
		    smem[PADDED_ADDR(THS+tid+BATCH*2)]*=V0n2;
*/
			__syncthreads();
		  }
		  #undef ADDRo2
		  #undef ADDRn2
		  
		  // STAGE III
		  // finish integration
		  if((BATCH>THSM*2)&(tid>0)){
		    for(i=0;i<portion-1;i++)	smem[PADDED_ADDR(i+portion*tid+BATCH)]*=smem[PADDED_ADDR(-1+portion*tid+BATCH)];
		    for(i=0;i<portion-1;i++)	smem[PADDED_ADDR(i+portion*tid+BATCH*2)]*=smem[PADDED_ADDR(-1+portion*tid+BATCH*2)];
		  }
		  __syncthreads();
		}
		
		//put back to output array
		for(j=0;j<BATCH;j+=THS) B[base+i+j+tid]=smem[PADDED_ADDR(j+tid)]; // should use registers for ARCH_MR>2 (Kepler+)
	}
	

	#undef THSM
	#undef BANK
	#undef BATCH
	//keep THS //undef THS
	#undef PADDED_ADDR
	#undef SMALL 0.0000000001
	#undef HALF 0.5
	#undef UNO 1.0
	//#undef DOS 2.0
	//#undef TRES 3.0
}
