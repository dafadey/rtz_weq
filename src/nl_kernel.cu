__global__ void
add_Nl_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, FL_DBL* spc_NlKerr, int Nt, int Nx, FL_DBL dz, FL_DBL T)
{
  const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
 	const unsigned int tidy0 = addNl_BSY*blockIdx.y;
	if(tidx<Nt)
	{
		const FL_DBL omega_1 = tidx == 0 ? 1.0F : T/2.0F/Pi/FL_DBL(tidx)*dz*(0.5F-0.5F*tanh(-(2.0F*Pi/T*FL_DBL(tidx)-0.5F)/0.1F));
		const FL_DBL omega = 3. * tanh(2.0F*Pi*FL_DBL(tidx)/T/3.)*dz;
		const FL_DBL factor = exp(.1 * (1. + tanh((2.0F*Pi*FL_DBL(tidx)/T - 5.)/.3))*dz);
		for(int j=0; j<addNl_BSY;j++)
		{
			const unsigned int tidy = j + tidy0;
			#define simple_ADD //simple is faster

			#ifdef simple_ADD
	//		FL_DBL omega_1=(((tidx==0)||(tidx==Nt-1))?0.0F:T/2.0F/Pi/FL_DBL(tidx))*dz*(0.5F-0.5F*tanh(-(2.0F*Pi/250.0F*FL_DBL(tidx)-0.5F)/0.1F))*0.5F;
	// do not forget to delete this fuckin 0.5! it is here because in CPU code there is extra dz. and it is wrong and result depends on dz
	// here nothing should depend on dz and 0.5 does not make any sense here. just for CPU=GPU test.

      FL_DBL spc_c = spc[tidx*2+tidy*Nt*2];
			spc_c += spc_Nl[tidx*2+1+tidy*Nt*2] * omega_1 + (spc_NlKerr ? spc_NlKerr[tidx*2+1+tidy*Nt*2] * omega : .0F);
      FL_DBL spc_s = spc[tidx*2+1+tidy*Nt*2];
			spc_s += -spc_Nl[tidx*2+tidy*Nt*2] * omega_1 - (spc_NlKerr ? spc_NlKerr[tidx*2+tidy*Nt*2] * omega : .0F);

			spc_c *= factor;
			spc_s *= factor;
			
			spc[tidx*2+tidy*Nt*2] = spc_c;
			spc[tidx*2+1+tidy*Nt*2] = spc_s;

			#endif
			#ifndef simple_ADD
			const FL_DBL sgn[]={ONE,-ONE};
			const int adr[]={1,-1};
			FL_DBL omega_1 = 1.0F/(2.0F*Pi/T)*dz/(((tidx>>1)>0)?FL_DBL((tidx>>1)):sgn[0]);
			FL_DBL omega = 3.* tanh((2.0F*Pi/T) * FL_DBL((tidx>>1))/3.) * dz;
			
			spc[tidx+tidy*Nt*2] += sgn[tidx&1] * (spc_Nl[tidx+adr[tidx&1]+tidy*Nt*2] * omega_1 + (spc_NlKerr ? spc_NlKerr[tidx+adr[tidx&1]+tidy*Nt*2] * omega : .0F));
			omega=1.0F/(2.0F*Pi/T)*dz/FL_DBL(((tidx+Nt)>>1));
			spc[tidx+Nt+tidy*Nt*2] += sgn[tidx&1] * (spc_Nl[tidx+Nt+adr[tidx&1]+tidy*Nt*2] * omega_1 + (spc_NlKerr ? spc_NlKerr[tidx+Nt+adr[tidx&1]+tidy*Nt*2] * omega : .0F));

			#endif
		}
	}
}

__global__ void
add_Nl_R2_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, FL_DBL* spc_Nl2, FL_DBL* spc_NlKerr, FL_DBL* spc_Nl2Kerr, int Nt, int Nx, FL_DBL dz, FL_DBL T)
{
  const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
 	const unsigned int tidy0 = addNl_BSY*blockIdx.y;
	if(tidx<Nt)
	{
		const FL_DBL omega_1 = tidx == 0 ? 1.0 : T/2.0F/Pi/FL_DBL(tidx)*dz*0.5F*(0.5F-0.5F*tanh(-(2.0F*Pi/250.0F*FL_DBL(tidx)-0.5F)/0.1F));
		const FL_DBL omega = 3. * tanh(2.0F*Pi*FL_DBL(tidx)/T/3.)*dz;
		const FL_DBL factor = exp(.1 * (1. + tanh((2.0F*Pi*FL_DBL(tidx)/T - 5.)/.3))*dz);

		for(int j=0; j<addNl_BSY;j++)
		{
			const unsigned int tidy = j + tidy0;
			//FL_DBL omega_1=(((tidx==0)||(tidx==Nt-1))?0.0:T/2.0F/Pi/FL_DBL(tidx))*dz*0.5F*(0.5F-0.5F*tanh(-(2.0F*Pi/250.0F*FL_DBL(tidx)-0.5F)/0.1F));
      FL_DBL spc_c = spc[tidx*2+tidy*Nt*2];
			spc_c += (spc_Nl2[tidx*2+1+tidy*Nt*2] - spc_Nl[tidx*2+1+tidy*Nt*2])*omega_1 + 
																((spc_Nl2Kerr && spc_NlKerr) ? (spc_Nl2Kerr[tidx*2+1+tidy*Nt*2] - spc_NlKerr[tidx*2+1+tidy*Nt*2])*omega : .0F);
      FL_DBL spc_s = spc[tidx*2+1+tidy*Nt*2];
			spc_s += (-spc_Nl2[tidx*2+tidy*Nt*2] + spc_Nl[tidx*2+tidy*Nt*2])*omega_1 +
																	((spc_Nl2Kerr && spc_NlKerr) ? (-spc_Nl2Kerr[tidx*2+tidy*Nt*2] + spc_NlKerr[tidx*2+tidy*Nt*2])*omega : .0F);
			spc_c *= factor;
			spc_s *= factor;
			
			spc[tidx*2+tidy*Nt*2] = spc_c;
			spc[tidx*2+1+tidy*Nt*2] = spc_s;
		}
	}
}

#define NL_DISPERSION_TEST .0  // .00062

__global__ void
find_Nl_fast_Kernel(FL_DBL* spc, FL_DBL* kerr_term, int nth, int Nt, int Nx, FL_DBL dt, FL_DBL* plasma, FL_DBL* THz_src, FL_DBL n0, FL_DBL Kerr, FL_DBL tunnelingExponent, FL_DBL tunnelingFactor)
{
	#define find_Nl_fast_2way

//	#define find_Nl_fast_nionly
	//consts
	
	#define E0 10.0F
	#define N0 100000.0F
	#define BASE0 0.000001F

	const FL_DBL norma=ONE/FL_DBL(Nt);
	const int PBS=256; // block size under processing
	const int SMxS=Nl_fast_SMxS; // shared X size 16 in order to prevent bank conflicts
	const int SMxSa=Nl_fast_SMxSa; // shared X size 16+1 in order to prevent bank conflicts
	const int SMyS=Nl_fast_SMyS; // shared Y size = PBS/(SMxS-1)
	__shared__ FL_DBL sdata[Nl_fast_SMxSa*Nl_fast_SMyS+1]; // contains all data of block //last value @ SMxS*SMyS is curent density (density per block)
	register int i,j;

	if (threadIdx.x==0) sdata[SMxSa*SMyS]=ZERO;
	__syncthreads();

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
		__syncthreads();

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
				register FL_DBL s = sdata[i+tid*SMxSa]*norma;
				//s=n0*s*s;
				//s=N0*exp(-E0/(fabs(s)+BASE0))*dt;
				s = n0 * (fabs(s)>0.000001F)?pow(s*s+0.0000000000000001F,-0.3125F * tunnelingFactor)*exp(-tunnelingExponent/(fabs(s)+0.000000000000000000001F)):0.0F;
				//new way!
				//s = n0 * (fabs(s)>1e-6F)?expf(-tunnelingExponent/(fabs(s)+1e-13F)-0.3125F*tunnelingFactor*logf(s*s+1e-13F)):0.0F;
				ni+=s*dt;
				sdata[i+tid*SMxSa] = ni;

			}
		}
		__syncthreads();

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
		__syncthreads();


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

		__syncthreads();

	//write back to global
		#ifndef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=1) // 1 way write
		{
			#ifndef find_Nl_fast_nionly
			register FL_DBL value = norma * spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt];
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = value * (NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);// - Kerr * value * value);
			if(kerr_term!=NULL)
				kerr_term[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = - Kerr * value * value * value;
			if(THz_src!=NULL) THz_src[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]/norma;
			if(plasma!=NULL) plasma[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=(NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
			#endif
			#ifdef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2];
			#endif
		}
		#endif
		#ifdef find_Nl_fast_2way
		for(i=0;i<PBS/SMyS/2;i+=2) //2 way unrolled r/w to SM
		{
			#ifndef find_Nl_fast_nionly
			register FL_DBL value1 = norma * spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt];
			register FL_DBL value2 = norma * spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt];
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = value1 * (NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);// - Kerr * value1 * value1);
			spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt] = value2 * (NL_DISPERSION_TEST + sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2]+sdata[SMxSa*SMyS]);// - Kerr * value2 * value2);

			if(kerr_term != NULL) {
				kerr_term[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = - Kerr * value1 * value1 * value1;
				kerr_term[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt] = - Kerr * value2 * value2 * value2;
			}

			if(THz_src!=NULL) {
				THz_src[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]/norma;
				THz_src[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]/norma;
			}
			
			if(plasma!=NULL) {
				plasma[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt]=(NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2]);
				plasma[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt]=(NL_DISPERSION_TEST + sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2]+sdata[SMxSa*SMyS]);
			}
			#endif
			#ifdef find_Nl_fast_nionly
			spc[threadIdx.x+i*nth+j*PBS+blockIdx.x*Nt] = NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+i*SMxSa*2];
			spc[threadIdx.x+i*nth+nth+j*PBS+blockIdx.x*Nt] = NL_DISPERSION_TEST + sdata[SMxSa*SMyS]+sdata[threadIdx.x+(threadIdx.x>>4)+(i+1)*SMxSa*2];
			#endif
		}
#endif

		__syncthreads();

		//dont forget to store curent density
		if (threadIdx.x==0) sdata[SMxSa*SMyS]+=sdata[SMxS*SMxSa-1];
		__syncthreads();
	}
}

__global__ void
find_Nl_fractions_fast_Kernel(FL_DBL* spc, FL_DBL* kerr_term, int Nt, int Nx, FL_DBL dt, FL_DBL Kerr, FL_DBL* params_pack, FL_DBL* plasma) {
//params pack format: {n0_0, tunneling_exp_0, tunneling_factor_0}, {n0_1, ...}, {}, ...
	const FL_DBL norma=ONE/FL_DBL(Nt);
	//depends on smem size (POT)
	const int PBS=256; // block size under processing
	//PLEASE NOTE THE ASSUMPTION PBS <= Nl_fractions_nth * SMxS
	
	//externals (used in host code to configure grid,threads):
	//#define Nl_fractions_nth (default is 32 = warp size)
	//externals used to fill up params_pack on host:
	//#define fractions_count
	//#define fractions_stride
	//below are is the fundamental one:
	const int SMxS_exponent = 4; // log2 of smem bank size
	//somewhat derived:
	const int SMxS = 1 << SMxS_exponent; // shared X *row*. basically it is a smem bank size
	const int SMxSa = (SMxS+1); // shared X size 16+1 in order to prevent bank conflicts
	const int SMyS = (PBS/SMxS); // shared 'Y' size
  
	const int halfWarps=(Nl_fractions_nth/SMxS);

	const int fracCount = fractions_count;
	const int fracStride = fractions_stride;
	__shared__ FL_DBL smemraw[SMxSa*SMyS*3+fracCount*(1+fracStride)]; // contains all data of block //last values are curent density (density per block)

	//let us use register file for constants and accumulators

	FL_DBL* sdata = smemraw;
	FL_DBL* sdens = &smemraw[SMxSa*SMyS];
	FL_DBL* sfrac = &smemraw[SMxSa*SMyS*2];
	FL_DBL* sfrac_acc = &smemraw[SMxSa*SMyS*3]; // fractions accumulators, size is Nl_fast_fractions_count
	FL_DBL* sfrac_params = &smemraw[SMxSa*SMyS*3 + fracCount]; // fraction parameters
	
	//first zero fractions accumulators
	for(int i=0; i < fracCount/blockDim.x+1; i++)
	{
		int id = threadIdx.x + blockDim.x * i;
		if (id < fracCount)
			sfrac_acc[id] = ZERO;
	}
	for(int i=0; i < fracCount/blockDim.x+1; i++) //push all parameters to smem
	{
		int id = threadIdx.x + blockDim.x * i;
		if (id < fracCount * fracStride)
			sfrac_params[id] = params_pack[id];
	}
	__syncthreads();

	FL_DBL* spc_row = &spc[blockIdx.x*Nt];
	FL_DBL* kerr_term_row = &kerr_term[blockIdx.x*Nt];
	FL_DBL* plasma_row = plasma ? &plasma[blockIdx.x*Nt] : nullptr;

	for(int j=0;j<Nt/PBS;j++)
	{
		FL_DBL* spc_block = &spc_row[j*PBS];
		FL_DBL* kerr_term_block = &kerr_term_row[j*PBS];
		FL_DBL* plasma_block = plasma_row ? &plasma_row[j*PBS] : nullptr;
		
		for(int i=0;i<PBS/SMxS/halfWarps;i++) {
			const int addr = threadIdx.x+(threadIdx.x>>SMxS_exponent) + i * SMxSa*halfWarps;
			sdata[addr] = spc_block[threadIdx.x+i*blockDim.x] * norma;
			sdens[addr] = ZERO;
		}
		//NOTE the assumption SMyS <= blockDim.x !!!!!!!
		//if you need greater PBS plase add a loop here
//#define PARALLEL_REDUCTION
//#define ALLTHREADS_ADD_STAIRS
//#define NOTRANS_TEST
#define TRANS_TRICK
#define DOKERR
#define DOTHEJOB // undefine it to test pure memory performance
#ifdef DOTHEJOB
//this gives quite a boost for exp and log calculation, why?...

		__syncthreads();
		for(int fracId=0; fracId<fracCount; fracId++) {
			register FL_DBL ni;
			register FL_DBL s;		
			register FL_DBL n0 = sfrac_params[fracId*fracStride];
			register FL_DBL tunneling_exponent = sfrac_params[fracId*fracStride+1];
			register FL_DBL tunneling_power = sfrac_params[fracId*fracStride+2];
			register FL_DBL tunneling_factor = sfrac_params[fracId*fracStride+3];
			
			#ifdef TRANS_TRICK
			//0 prepare ratios
			for(int i=0;i<PBS/SMxS/halfWarps;i++) {
				const int addr = threadIdx.x+(threadIdx.x>>SMxS_exponent) + i * SMxSa*halfWarps;
				const FL_DBL s = sdata[addr];
				#ifdef NOTRANS_TEST
				sfrac[addr] = tunneling_exponent * n0 * s * s;
				#else
				sfrac[addr] = (fabs(s)>1e-6F)?FL_DBL(Adim)*tunneling_factor*exp(-tunneling_exponent/(fabs(s)+1e-13F)-tunneling_power*log(fabs(s)+1e-13F)):0.0F;
				#endif
			}
			__syncthreads();
			#endif
						
			//1 subintegration (16steps each thread)
			if(threadIdx.x < SMyS) {
				ni=ZERO;
				for(int i = 0; i < SMxS; i++) {
					#ifndef TRANS_TRICK
					FL_DBL s = sdata[i + threadIdx.x * SMxSa];
					#ifdef NOTRANS_TEST
					s = tunneling_exponent * n0 * s * s;
					#else
					s = (fabs(s)>1e-6F)?FL_DBL(Adim)*tunneling_factor*exp(-tunneling_exponent/(fabs(s)+1e-13F)-tunneling_power*log(fabs(s)+1e-13F)):0.0F;
					#endif
					#else
					const FL_DBL s = sfrac[i + threadIdx.x * SMxSa];
					#endif
					ni += s * dt;
					sfrac[i + threadIdx.x * SMxSa] = ni;
				}
			}
			__syncthreads();
			//2 collect ni steps in SMxSa-th s elements of each sm rows
#ifdef PARALLEL_REDUCTION
			//collect SMxS-1 value to cell with index SMxS (last cell of SMxSa extended block)
			//since SMyS<blockDim.x we can do it per thread:

			if(threadIdx.x < SMyS)
				sfrac[SMxSa*threadIdx.x+SMxS] = sfrac[SMxSa*threadIdx.x+SMxS-1];
			__syncthreads();
			
			//now do the actual parallel reduction
			for(int level = 0; (1 << level) <= (SMyS >> 1); level++) {
				if(threadIdx.x < (SMyS >> 1)) {
					int id = (threadIdx.x >> level) << (level + 1);
					int subid = 1 << level;
					int maxoff = subid - 1;
					int off = threadIdx.x & maxoff;
					sfrac[(id+subid+off)*SMxSa+SMxS] += sfrac[(id+maxoff)*SMxSa+SMxS];
				}
				//__syncthreads();
			}
			__syncthreads();
#else			
			if(threadIdx.x==0)
			{
				ni = sfrac[SMxS - 1];
				for(int i=1;i<SMyS;i++)
				{
					sfrac[i * SMxSa - 1] = ni;
					ni += sfrac[i * SMxSa + SMxS - 1];
				}
				sfrac[SMyS * SMxSa - 1] = ni;
			}
			__syncthreads();
#endif

			//3 add step function
#ifdef ALLTHREADS_ADD_STAIRS
			int blockid = threadIdx.x >> 1;
			if(blockid != 0) {
				int blocksubsize = SMxS >> 1;
				s = sfrac[blockid * SMxSa-1];
				for(int i = 0; i < blocksubsize; i++)
					sfrac[i + blockid * SMxSa + blocksubsize * (threadIdx.x & 1)] += s;				
			}
#else
			if(threadIdx.x < SMyS) {
				if(threadIdx.x != 0) {
					s = sfrac[threadIdx.x * SMxSa-1];
					for(int i = 0; i < SMxS; i++)
						sfrac[i + threadIdx.x * SMxSa] += s;				
				}
			}
#endif
			__syncthreads();
			
			//4 summ up to sdens (overall plasma density)
/*
			dn / dt = (n0 - n) * f(t)
			d(ln(n0-n)) = - f(t) * dt
			n0-n = n0 exp(-Int(f(t))
			n = n0 (1-exp(-Int(f)))
*/
			FL_DBL base = sfrac_acc[fracId];
			for(int i=0;i<PBS/SMxS/halfWarps;i++) {
				int addr = threadIdx.x+(threadIdx.x>>SMxS_exponent) + i * SMxSa*halfWarps;
				FL_DBL v = sfrac[addr] + base;
#ifdef NOTRANS_TEST
				sdens[addr] += n0 * v;
#else
				sdens[addr] += n0 * (ONE - exp(-v));
#endif
			}
			
			__syncthreads();
			if(threadIdx.x == 0)
				sfrac_acc[fracId] = base + sfrac[SMxSa*SMyS-1]; // store fraction base
			__syncthreads();
		}
#else //skip the real job and just copy
		FL_DBL n0 = sfrac_params[0];
		for(int i=0;i<PBS/SMxS/halfWarps;i++) {
				int addr = threadIdx.x+(threadIdx.x>>SMxS_exponent) + i * SMxSa*halfWarps;
				sdens[addr] = n0;
		}
#endif		
		//apply plasma to the field
		for(int i=0;i<PBS/SMxS/halfWarps;i++) {
			int addr = threadIdx.x+(threadIdx.x>>SMxS_exponent) + i * SMxSa*halfWarps;
			int ram_addr = threadIdx.x+i*blockDim.x;
			FL_DBL v = sdata[addr];
#ifdef DOKERR
			kerr_term_block[ram_addr] = - Kerr * v*v*v;
#endif
			spc_block[ram_addr] = v * sdens[addr];
			if(plasma_block)
				plasma_block[ram_addr] = sdens[addr];
		}
		__syncthreads();
	}
#undef DOTHEJOB
#undef PARALLEL_REDUCTION
#undef ALLTHREADS_ADD_STAIRS
#undef NOTRANS_TEST
#undef TRANS_TRICK
#undef DOKERR

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
	__syncthreads();
	register int i,j;
	for(j=0;j<Nt/NlS_BSX;j++)
	{
		for(i=0;i<NlS_BSX;i++){
//		 sdata[(threadIdx.x&15)+i*SMxSa+(threadIdx.x/16)*SMxSa*SMyS]=spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt];
			sdata[threadIdx.x+(threadIdx.x/16)+i*SMxSa*2]=spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt];
		}
		__syncthreads();

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
		__syncthreads();
		for(i=0;i<NlS_BSX;i++){
//			spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt]=sdata[(threadIdx.x&15)+i*SMxSa+(threadIdx.x/16)*SMxSa*SMyS];
			spc[threadIdx.x+j*NlS_BSX+(i+blockIdx.x*NlS_BSX)*Nt]=sdata[threadIdx.x+(threadIdx.x/16)+i*SMxSa*2];
		}

//		__syncthreads();
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
add_Nl_fast_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, int Nt, int Nx, FL_DBL dz, FL_DBL T)
{
  const int tidy = blockIdx.y;
	
	__shared__ FL_DBL sdata[(Nl_fast_NTH<<2)];
	sdata[threadIdx.x]=spc[threadIdx.x+blockDim.x*blockIdx.x*2+tidy*Nt*2];
	sdata[threadIdx.x+blockDim.x]=spc[threadIdx.x+blockDim.x+blockDim.x*blockIdx.x*2+tidy*Nt*2];
	sdata[threadIdx.x+blockDim.x*2]=spc_Nl[threadIdx.x+tidy*Nt*2+blockDim.x*blockIdx.x*2];
	sdata[threadIdx.x+blockDim.x*2+blockDim.x]=spc_Nl[threadIdx.x+blockDim.x+tidy*Nt*2+blockDim.x*blockIdx.x*2];
	__syncthreads();



	const FL_DBL sgn[]={1.0F,-1.0F};
	const int adr[]={1,-1};
	register int tidx=threadIdx.x+blockDim.x*blockIdx.x*2;
	register FL_DBL omega=1.0F/(2.0F*Pi/T)*dz/(((tidx>>1)>0)?FL_DBL((tidx>>1)):sgn[0]);
	sdata[threadIdx.x]+=sgn[threadIdx.x&1]*sdata[threadIdx.x+adr[threadIdx.x&1]+blockDim.x*2]*omega;
	omega=1.0F/(2.0F*Pi/T)*dz/FL_DBL(((tidx+blockDim.x)>>1));
	sdata[threadIdx.x+blockDim.x]+=sgn[threadIdx.x&1]*sdata[threadIdx.x+blockDim.x+adr[threadIdx.x&1]+blockDim.x*2]*omega;
	__syncthreads();

	spc[threadIdx.x+tidy*Nt*2+blockDim.x*blockIdx.x*2]=sdata[threadIdx.x];
	spc[threadIdx.x+blockDim.x+tidy*Nt*2+blockDim.x*blockIdx.x*2]=sdata[threadIdx.x+blockDim.x];
	__syncthreads();
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
