//#ifndef _TEMPLATE_KERNEL_H_
//#define _TEMPLATE_KERNEL_H_
__device__ FL_DBL mulRe(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return xr*yr-xi*yi;}
__device__ FL_DBL mulIm(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return xi*yr+xr*yi;}

extern "C"
inline FL_DBL mulReCPU(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return xr*yr-xi*yi;}
extern "C"
inline FL_DBL mulImCPU(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return xi*yr+xr*yi;}


//__device__ FL_DBL divRe(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return __fdividef(xr*yr+xi*yi,yr*yr+yi*yi);}
//__device__ FL_DBL divIm(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return __fdividef(xi*yr-xr*yi,yr*yr+yi*yi);}

__device__ FL_DBL divRe(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return (xr*yr+xi*yi)/(yr*yr+yi*yi);}
__device__ FL_DBL divIm(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return (xi*yr-xr*yi)/(yr*yr+yi*yi);}

extern "C"
inline FL_DBL divReCPU(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return (xr*yr+xi*yi)/(yr*yr+yi*yi);}
extern "C"
inline FL_DBL divImCPU(FL_DBL xr, FL_DBL xi, FL_DBL yr, FL_DBL yi) {return (xi*yr-xr*yi)/(yr*yr+yi*yi);}


__device__ FL_DBL getA(FL_DBL* Ay, FL_DBL* Ax, int stp, int nofa)
{
	register int i;
	register FL_DBL res=0.0;
	for(i=0;i<nofa;i++) res+=Ay[i]*Ax[i*stp];
	return res;
}

extern "C"
inline FL_DBL getACPU(FL_DBL* Ay, FL_DBL* Ax, int stp, int nofa)
{
	register int i;
	register FL_DBL res=0.0;
	for(i=0;i<nofa;i++) res+=Ay[i]*Ax[i*stp];
	return res;
}

__device__ FL_DBL get_t_A2(FL_DBL* Ay, FL_DBL* Ax, int stpy, int stpx, int nofa)
{
	register int i;
	register FL_DBL res=0.0;
	for(i=0;i<nofa;i++) res+=Ay[i*stpy]*Ax[i*stpx];
	return res;
}

extern "C"
inline FL_DBL get_t_A2CPU(FL_DBL* Ay, FL_DBL* Ax, int stpy, int stpx, int nofa)
{
	register int i;
	register FL_DBL res=0.0;
	for(i=0;i<nofa;i++) res+=Ay[i*stpy]*Ax[i*stpx];
	return res;
}

__device__ FL_DBL getA2(FL_DBL* Ay, FL_DBL* Ax, int stpy, int stpx, int nofa)
{
	register int i;
	register FL_DBL res=Ay[0]*Ax[0];
	for(i=1;i<nofa;i++) res=res-Ay[i*stpy]*Ax[i*stpx];
	return res;
}

extern "C"
inline FL_DBL getA2CPU(FL_DBL* Ay, FL_DBL* Ax, int stpy, int stpx, int nofa)
{
	register int i;
	register FL_DBL res=Ay[0]*Ax[0];
	for(i=1;i<nofa;i++) res=res-Ay[i*stpy]*Ax[i*stpx];
	return res;
}

__global__ void
simpleSWEEPKernel( FL_DBL* a, FL_DBL* b, FL_DBL* c, FL_DBL* y, FL_DBL* x, int NX, int NY, FL_DBL* yy)
{
//if no need in y array any more please set yy=y
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  register int i;
  if(tid<NY)
  {
		FL_DBL k[2];
		FL_DBL tmp[2];
		x[tid*2]=a[tid*2];
		x[tid*2+1]=a[tid*2+1];
	//forward sweep:
	//use x as temporary variable for new a coeffisients during forward sweep
	
		yy[tid*2]=y[tid*2];
		yy[tid*2+1]=y[tid*2+1];
	
		for(i=0;i<NX-1;i++)
		{
			k[0]=divRe(b[tid*2+NY*(i+1)*2],b[tid*2+1+NY*(i+1)*2],x[tid*2+NY*i*2],x[tid*2+1+NY*i*2]);
			k[1]=divIm(b[tid*2+NY*(i+1)*2],b[tid*2+1+NY*(i+1)*2],x[tid*2+NY*i*2],x[tid*2+1+NY*i*2]);
			x[tid*2+NY*(i+1)*2]=a[tid*2+NY*(i+1)*2]-mulRe(c[tid*2+NY*i*2],c[tid*2+1+NY*i*2],k[0],k[1]);
			x[tid*2+1+NY*(i+1)*2]=a[tid*2+1+NY*(i+1)*2]-mulIm(c[tid*2+NY*i*2],c[tid*2+1+NY*i*2],k[0],k[1]);
		
			yy[tid*2+NY*(i+1)*2]=y[tid*2+NY*(i+1)*2]-mulRe(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
			yy[tid*2+1+NY*(i+1)*2]=y[tid*2+1+NY*(i+1)*2]-mulIm(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
			//__syncthreads();
		}
	//backward sweep:
	//find the first x:
	
		k[0]=x[tid*2+NY*(NX-1)*2];
		k[1]=x[tid*2+1+NY*(NX-1)*2];
		x[tid*2+NY*(NX-1)*2]=divRe(yy[tid*2+NY*(NX-1)*2],yy[tid*2+1+NY*(NX-1)*2],k[0],k[1]);
		x[tid*2+1+NY*(NX-1)*2]=divIm(yy[tid*2+NY*(NX-1)*2],yy[tid*2+1+NY*(NX-1)*2],k[0],k[1]);
	
		for(i=NX-2;i>-1;i--)
		{
			k[0]=mulRe(c[tid*2+NY*i*2], c[tid*2+1+NY*i*2], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);
			k[1]=mulIm(c[tid*2+NY*i*2], c[tid*2+1+NY*i*2], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);

			tmp[0]=x[tid*2+NY*i*2];
			tmp[1]=x[tid*2+1+NY*i*2];

			x[tid*2+NY*i*2]=divRe(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);
			x[tid*2+1+NY*i*2]=divIm(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);
		}
  }
}

extern "C"
void
simpleFIXEDSWEEPKernel_CPU(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, int NX, int NY, FL_DBL* yy)
{
//if no need in y array any more please set yy=y
	for(int tid=0;tid<NY;tid++)
	{
	  register int i;

		register FL_DBL BRe=-BY[tid*2];
		register FL_DBL BIm=-BY[tid*2+1];
		register FL_DBL CRe=-CY[tid*2];
		register FL_DBL CIm=-CY[tid*2+1];
		register FL_DBL ARe=getACPU(&AY[tid*2],&AX[0],NX, NOFA);
		register FL_DBL AIm=getACPU(&AY[(tid*2)+1],&AX[0],NX, NOFA);

		x[tid*2]=ARe;
		x[tid*2+1]=AIm;
	//forward sweep:
	//use x as temporary variable for new a coeffisients during forward sweep
		FL_DBL k[2];
		FL_DBL tmp[2];

		yy[tid*2]=y[tid*2];
		yy[tid*2+1]=y[tid*2+1];
	
		for(i=0;i<NX-1;i++)
		{

			k[0]=divReCPU(BX[i+1]*BRe,BX[i+1]*BIm,x[tid*2+2*NY*i],x[tid*2+1+2*NY*i]);
			k[1]=divImCPU(BX[i+1]*BRe,BX[i+1]*BIm,x[tid*2+2*NY*i],x[tid*2+1+2*NY*i]);
		
			ARe=getACPU(&AY[tid*2],&AX[i+1],NX, NOFA);
			AIm=getACPU(&AY[(tid*2)+1],&AX[i+1],NX, NOFA);

			x[tid*2+2*NY*(i+1)]=ARe-mulReCPU(CRe*CX[i],CIm*CX[i],k[0],k[1]);
			x[tid*2+1+2*NY*(i+1)]=AIm-mulImCPU(CRe*CX[i],CIm*CX[i],k[0],k[1]);
		
		
			yy[tid*2+NY*(i+1)*2]=y[tid*2+NY*(i+1)*2]-mulReCPU(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
			yy[tid*2+1+NY*(i+1)*2]=y[tid*2+1+NY*(i+1)*2]-mulImCPU(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
		}
	//backward sweep:
	//find the first x:
	
		tmp[0]=x[tid*2+NY*(NX-1)*2];
		tmp[1]=x[tid*2+1+NY*(NX-1)*2];
	
		x[tid*2+NY*(NX-1)*2]=divReCPU(yy[tid*2+NY*(NX-1)*2], yy[tid*2+1+NY*(NX-1)*2], tmp[0],  tmp[1]);
		x[tid*2+1+NY*(NX-1)*2]=divImCPU(yy[tid*2+NY*(NX-1)*2], yy[tid*2+1+NY*(NX-1)*2], tmp[0],  tmp[1]);
		for(i=NX-2;i>-1;i--)
		{
			k[0]=mulReCPU(CRe*CX[i], CIm*CX[i], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);
			k[1]=mulImCPU(CRe*CX[i], CIm*CX[i], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);
		
			tmp[0]=x[tid*2+NY*i*2];
			tmp[1]=x[tid*2+1+NY*i*2];
	
			x[tid*2+NY*i*2]=divReCPU(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);
			x[tid*2+1+NY*i*2]=divImCPU(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);

		}
	}
}


__global__ void
simpleFIXEDSWEEPKernel(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, int NX, int NY, FL_DBL* yy)
{
//if no need in y array any more please set yy=y
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  register int i;
  if(tid<NY)
  {
		FL_DBL ARe=AY[tid*2];
		FL_DBL AIm=AY[tid*2+1];
		FL_DBL BRe=BY[tid*2];
		FL_DBL BIm=BY[tid*2+1];
		FL_DBL CRe=CY[tid*2];
		FL_DBL CIm=CY[tid*2+1];

		x[tid*2]=ARe*AX[0];
		x[tid*2+1]=AIm*AX[0];
	//forward sweep:
	//use x as temporary variable for new a coeffisients during forward sweep
		FL_DBL k[2];
		FL_DBL tmp[2];

		yy[tid*2]=y[tid*2];
		yy[tid*2+1]=y[tid*2+1];
	
		for(i=0;i<NX-1;i++)
		{
			k[0]=divRe(BX[i+1]*BRe,BX[i+1]*BIm,x[tid*2+2*NY*i],x[tid*2+1+2*NY*i]);
			k[1]=divIm(BX[i+1]*BRe,BX[i+1]*BIm,x[tid*2+2*NY*i],x[tid*2+1+2*NY*i]);
		
			x[tid*2+2*NY*(i+1)]=AX[i+1]*ARe-mulRe(CRe*CX[i],CIm*CX[i],k[0],k[1]);
			x[tid*2+1+2*NY*(i+1)]=AX[i+1]*AIm-mulIm(CRe*CX[i],CIm*CX[i],k[0],k[1]);
		
		
			yy[tid*2+NY*(i+1)*2]=y[tid*2+NY*(i+1)*2]-mulRe(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
			yy[tid*2+1+NY*(i+1)*2]=y[tid*2+1+NY*(i+1)*2]-mulIm(yy[tid*2+NY*i*2],yy[tid*2+1+NY*i*2],k[0],k[1]);
		}
	//backward sweep:
	//find the first x:
	
		tmp[0]=x[tid*2+NY*(NX-1)*2];
		tmp[1]=x[tid*2+1+NY*(NX-1)*2];
	
		x[tid*2+NY*(NX-1)*2]=divRe(yy[tid*2+NY*(NX-1)*2], yy[tid*2+1+NY*(NX-1)*2], tmp[0],  tmp[1]);
		x[tid*2+1+NY*(NX-1)*2]=divIm(yy[tid*2+NY*(NX-1)*2], yy[tid*2+1+NY*(NX-1)*2], tmp[0],  tmp[1]);
		for(i=NX-2;i>-1;i--)
		{
			k[0]=mulRe(CRe*CX[i], CIm*CX[i], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);
			k[1]=mulIm(CRe*CX[i], CIm*CX[i], x[tid*2+NY*(i+1)*2], x[tid*2+1+NY*(i+1)*2]);
		
			tmp[0]=x[tid*2+NY*i*2];
			tmp[1]=x[tid*2+1+NY*i*2];
	
			x[tid*2+NY*i*2]=divRe(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);
			x[tid*2+1+NY*i*2]=divIm(yy[tid*2+NY*i*2]-k[0],yy[tid*2+1+NY*i*2]-k[1],tmp[0], tmp[1]);

		}
  }
}

__global__ void
SWEEP_1_FIXED_Kernel(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, FL_DBL* bb, FL_DBL* cc, int NX, int NY, int BLOCKX, int BLOCKY, FL_DBL* yy, FL_DBL* a2, FL_DBL* b2, FL_DBL* c2, FL_DBL* y2) 
{

	extern __shared__ FL_DBL sdata[];
		FL_DBL* sAX=sdata;
		FL_DBL* sBX=&sAX[BLOCKX*NOFA];
		FL_DBL* sCX=&sBX[BLOCKX];
		FL_DBL* sAY=&sCX[BLOCKX];
		FL_DBL* sBY=&sAY[2*BLOCKY*NOFA];
		FL_DBL* sCY=&sBY[BLOCKY*2];
	
		const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
		const unsigned int tidy = blockIdx.y*BLOCKX;
		register int i;

  if (tid<NY-1) //caution!!! needs revision.
  {

	for(int i=0;i<NOFA;i++)
	{
		sAY[threadIdx.x+BLOCKY*2*i]=AY[threadIdx.x+blockDim.x*blockIdx.x*2+NY*2*i];
		sAY[threadIdx.x+blockDim.x+BLOCKY*2*i]=AY[threadIdx.x+blockDim.x+blockDim.x*blockIdx.x*2+NY*2*i];
	}

    sBY[threadIdx.x]=BY[threadIdx.x+blockDim.x*blockIdx.x*2];
    sBY[threadIdx.x+blockDim.x]=BY[threadIdx.x+blockDim.x+blockDim.x*blockIdx.x*2];

    sCY[threadIdx.x]=CY[threadIdx.x+blockDim.x*blockIdx.x*2];
    sCY[threadIdx.x+blockDim.x]=CY[threadIdx.x+blockDim.x+blockDim.x*blockIdx.x*2];

	for(i=0;i<-floor(-FL_DBL(BLOCKX)/FL_DBL(BLOCKY));i++)
    {
			if (threadIdx.x<BLOCKX)
			{
					for(int j=0;j<NOFA;j++)
						sAX[threadIdx.x+i*BLOCKY+BLOCKX*j]=AX[threadIdx.x+i*BLOCKY+tidy+NX*j];
					
					sBX[threadIdx.x+i*BLOCKY]=BX[threadIdx.x+i*BLOCKY+tidy];
					sCX[threadIdx.x+i*BLOCKY]=CX[threadIdx.x+i*BLOCKY+tidy];
			}
    }

  }

  __syncthreads();
  if(tid<NY)
  {
		register FL_DBL ARe[NOFA];
		for(int i=0;i<NOFA;i++)
			ARe[i]=sAY[threadIdx.x*2+BLOCKY*2*i];
		register FL_DBL AIm[NOFA];
		for(int i=0;i<NOFA;i++)
			AIm[i]=sAY[threadIdx.x*2+1+BLOCKY*2*i];
		register FL_DBL BRe=sBY[threadIdx.x*2];
		register FL_DBL BIm=sBY[threadIdx.x*2+1];
		register FL_DBL CRe=sCY[threadIdx.x*2];
		register FL_DBL CIm=sCY[threadIdx.x*2+1];

		x[tid*2+2*NY*tidy]=getA(ARe,&sAX[0],BLOCKX,NOFA);
		x[tid*2+1+2*NY*tidy]=getA(AIm,&sAX[0],BLOCKX,NOFA);
	//forward sweep:
	//use x as temporary variable for new a coeffisients during forward sweep

		bb[2*tid+NY*tidy*2]=BRe*sBX[0];
		bb[2*tid+1+NY*tidy*2]=BIm*sBX[0];
		yy[2*tid+NY*tidy*2]=y[2*tid+NY*tidy*2];
		yy[2*tid+1+NY*tidy*2]=y[2*tid+1+NY*tidy*2];
		register FL_DBL k[2];
		register FL_DBL tmp_x[2];
		register FL_DBL tmp_bb[2];
		register FL_DBL tmp_yy[2];
		register FL_DBL tmp_cc[2];
		register int adr=tid*2+NY*tidy*2;
		register int adr2;
		tmp_x[0]=x[adr];
		tmp_x[1]=x[adr+1];
		tmp_yy[0]=yy[adr];
		tmp_yy[1]=yy[adr+1];
		tmp_bb[0]=bb[adr];
		tmp_bb[1]=bb[adr+1];

		for(i=0;i<BLOCKX-1;i++)
		{
			adr2=adr+NY*2;
			register FL_DBL tmp;

	//		k=B*BX[i+1+tidy]/x[adr];

			//k=B*__fdividef(BX[i+1+tidy];tmp_x);
			k[0]=divRe(BRe*sBX[i+1],BIm*sBX[i+1],tmp_x[0],tmp_x[1]);
			k[1]=divIm(BRe*sBX[i+1],BIm*sBX[i+1],tmp_x[0],tmp_x[1]);
		
			//tmp_x=A*sAX[i+1]-C*sCX[i]*k;
			tmp_x[0]=getA(ARe,&sAX[i+1],BLOCKX,NOFA)-mulRe(CRe*sCX[i],CIm*sCX[i],k[0],k[1]);
			tmp_x[1]=getA(AIm,&sAX[i+1],BLOCKX,NOFA)-mulIm(CRe*sCX[i],CIm*sCX[i],k[0],k[1]);

			//tmp_bb=-tmp_bb*k;
			tmp=tmp_bb[0];
			tmp_bb[0]=-mulRe(tmp_bb[0],tmp_bb[1],k[0],k[1]);
			tmp_bb[1]=-mulIm(tmp,tmp_bb[1],k[0],k[1]);

			//tmp_yy=y[adr2]-tmp_yy*k;
			tmp=tmp_yy[0];
			tmp_yy[0]=y[adr2]-mulRe(tmp_yy[0],tmp_yy[1],k[0],k[1]);
			tmp_yy[1]=y[adr2+1]-mulIm(tmp,tmp_yy[1],k[0],k[1]);

			//x[adr2]=tmp_x;
			x[adr2]=tmp_x[0];
			x[adr2+1]=tmp_x[1];

			//yy[adr2]=tmp_yy;
			yy[adr2]=tmp_yy[0];
			yy[adr2+1]=tmp_yy[1];
		
			//bb[adr2]=tmp_bb;
			bb[adr2]=tmp_bb[0];
			bb[adr2+1]=tmp_bb[1];
		
			adr+=NY*2;
		}
	//backward sweep:
	//find the first x:
		//__syncthreads();
		adr=tid*2+NY*(tidy+BLOCKX-2)*2;
		adr2=adr+NY*2;

		cc[adr2]=CRe*sCX[BLOCKX-1];
		cc[adr2+1]=CIm*sCX[BLOCKX-1];

		tmp_yy[0]=yy[adr2];
		tmp_yy[1]=yy[adr2+1];
	
		tmp_bb[0]=bb[adr2];
		tmp_bb[1]=bb[adr2+1];
	
		tmp_cc[0]=cc[adr2];
		tmp_cc[1]=cc[adr2+1];
	
		for(i=BLOCKX-2;i>-1;i--)
		{
			adr2=adr+NY*2;
			register FL_DBL tmp;
			//k=C*sCX[i]/x[adr2];
			//k=C*__fdividef(sCX[i];x[adr2]);
			k[0]=divRe(CRe*sCX[i],CIm*sCX[i],x[adr2],x[adr2+1]);
			k[1]=divIm(CRe*sCX[i],CIm*sCX[i],x[adr2],x[adr2+1]);
		
			//tmp_yy=yy[adr]-tmp_yy*k;
			tmp=tmp_yy[0];
			tmp_yy[0]=yy[adr]-mulRe(tmp_yy[0],tmp_yy[1],k[0],k[1]);
			tmp_yy[1]=yy[adr+1]-mulIm(tmp,tmp_yy[1],k[0],k[1]);
				
			//tmp_cc=-tmp_cc*k;
			tmp=tmp_cc[0];
			tmp_cc[0]=-mulRe(tmp_cc[0],tmp_cc[1],k[0],k[1]);
			tmp_cc[1]=-mulIm(tmp,tmp_cc[1],k[0],k[1]);
		
			//tmp_bb=bb[adr]-tmp_bb*k;
			tmp=tmp_bb[0];
			tmp_bb[0]=bb[adr]-mulRe(tmp_bb[0],tmp_bb[1],k[0],k[1]);
			tmp_bb[1]=bb[adr+1]-mulIm(tmp,tmp_bb[1],k[0],k[1]);
		
			cc[adr]=tmp_cc[0];
			cc[adr+1]=tmp_cc[1];
		
			bb[adr]=tmp_bb[0];
			bb[adr+1]=tmp_bb[1];

			yy[adr]=tmp_yy[0];
			yy[adr+1]=tmp_yy[1];

			adr-=NY*2;
		}


		//fill small secondary matrix:
		y2[tid*2+(blockIdx.y*2)*NY*2]=yy[tid*2+NY*tidy*2];
		y2[tid*2+1+(blockIdx.y*2)*NY*2]=yy[tid*2+1+NY*tidy*2];

		y2[tid*2+(blockIdx.y*2+1)*NY*2]=yy[tid*2+NY*(tidy+BLOCKX-1)*2];
		y2[tid*2+1+(blockIdx.y*2+1)*NY*2]=yy[tid*2+1+NY*(tidy+BLOCKX-1)*2];
	
		if(blockIdx.y!=0)
		{
			a2[tid*2+(blockIdx.y*2)*NY*2]=bb[tid*2+NY*tidy*2];
			a2[tid*2+1+(blockIdx.y*2)*NY*2]=bb[tid*2+1+NY*tidy*2];
		}
		else
		{
			a2[tid*2+(blockIdx.y*2)*NY*2]=x[tid*2+NY*tidy*2];
			a2[tid*2+1+(blockIdx.y*2)*NY*2]=x[tid*2+1+NY*tidy*2];
		}
	
		if(blockIdx.y!=NX/BLOCKX-1)
		{
			a2[tid*2+(blockIdx.y*2+1)*NY*2]=cc[tid*2+NY*(tidy+BLOCKX-1)*2];
			a2[tid*2+1+(blockIdx.y*2+1)*NY*2]=cc[tid*2+1+NY*(tidy+BLOCKX-1)*2];
		}
		else
		{
			a2[tid*2+(blockIdx.y*2+1)*NY*2]=x[tid*2+NY*(tidy+BLOCKX-1)*2];
			a2[tid*2+1+(blockIdx.y*2+1)*NY*2]=x[tid*2+1+NY*(tidy+BLOCKX-1)*2];
		}
	
		b2[tid*2+(blockIdx.y*2)*NY*2]=x[tid*2+NY*tidy*2];
		b2[tid*2+1+(blockIdx.y*2)*NY*2]=x[tid*2+1+NY*tidy*2];

		b2[tid*2+(blockIdx.y*2+1)*NY*2]=bb[tid*2+NY*(tidy+BLOCKX-1)*2];
		b2[tid*2+1+(blockIdx.y*2+1)*NY*2]=bb[tid*2+1+NY*(tidy+BLOCKX-1)*2];
	
		c2[tid*2+(blockIdx.y*2)*NY*2]=cc[tid*2+NY*tidy*2];
		c2[tid*2+1+(blockIdx.y*2)*NY*2]=cc[tid*2+1+NY*tidy*2];

		c2[tid*2+(blockIdx.y*2+1)*NY*2]=x[tid*2+NY*(tidy+BLOCKX-1)*2];
		c2[tid*2+1+(blockIdx.y*2+1)*NY*2]=x[tid*2+1+NY*(tidy+BLOCKX-1)*2];
  }
}

__global__ void
SWEEP_2_Kernel( FL_DBL* a, FL_DBL* bb, FL_DBL* cc, FL_DBL* y, int NX, int NY, int BLOCKX, int BLOCKY, FL_DBL* x2)
{
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tidy = blockIdx.y*BLOCKX;
  register int i;
  
  if(tid<NY)
  {
		register FL_DBL k[2];
		register FL_DBL l[2];
		register FL_DBL r[2];
		register FL_DBL m1[2];
		register FL_DBL m2[2];
		l[0]=x2[tid*2+(blockIdx.y*2)*NY*2];
		l[1]=x2[tid*2+1+(blockIdx.y*2)*NY*2];
		r[0]=x2[tid*2+(blockIdx.y*2+1)*NY*2];
		r[1]=x2[tid*2+1+(blockIdx.y*2+1)*NY*2];
///last sweep:
		for(i=0;i<BLOCKX;i++)
		{
			m1[0]=mulRe(l[0],l[1],bb[tid*2+NY*(tidy+i)*2],bb[tid*2+1+NY*(tidy+i)*2]);
			m1[1]=mulIm(l[0],l[1],bb[tid*2+NY*(tidy+i)*2],bb[tid*2+1+NY*(tidy+i)*2]);
		
			m2[0]=mulRe(r[0],r[1],cc[tid*2+NY*(tidy+i)*2],cc[tid*2+1+NY*(tidy+i)*2]);
			m2[1]=mulIm(r[0],r[1],cc[tid*2+NY*(tidy+i)*2],cc[tid*2+1+NY*(tidy+i)*2]);
		
			k[0]=a[tid*2+NY*(tidy+i)*2];
			k[1]=a[tid*2+1+NY*(tidy+i)*2];
		
			a[tid*2+NY*(tidy+i)*2]=divRe(y[tid*2+NY*(tidy+i)*2]-m1[0]-m2[0],y[tid*2+1+NY*(tidy+i)*2]-m1[1]-m2[1],k[0],k[1]);
			a[tid*2+1+NY*(tidy+i)*2]=divIm(y[tid*2+NY*(tidy+i)*2]-m1[0]-m2[0],y[tid*2+1+NY*(tidy+i)*2]-m1[1]-m2[1],k[0],k[1]);
		}
  }
}

__global__ void
forward_t_Kernel(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, int NX, int NY) 
{
	//extern __shared__ FL_DBL sdata[];
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	register FL_DBL ARe=get_t_A2(&AY[tid*2],&AX[tidy],NY*2,NX, NOFA);
	register FL_DBL AIm=get_t_A2(&AY[tid*2+1],&AX[tidy],NY*2,NX, NOFA);
	register FL_DBL BRe=BY[tid*2]*BX[tidy];
	register FL_DBL BIm=BY[tid*2+1]*BX[tidy];
	register FL_DBL CRe=CY[tid*2]*CX[tidy];
	register FL_DBL CIm=CY[tid*2+1]*CX[tidy];
	if(tid<NY)
	{
		if((tidy<NX-1)&(tidy>0))
		{
			x[tid*2+tidy*NY*2]=mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulRe(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
			x[tid*2+1+tidy*NY*2]=mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulIm(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
		}
		if(tidy==0)
		{
			x[tid*2+tidy*NY*2]=mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
			x[tid*2+1+tidy*NY*2]=mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
		}
		if(tidy==NX-1)
		{
			x[tid*2+tidy*NY*2]=mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]);
			x[tid*2+1+tidy*NY*2]=mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]);
		}
	}
}

extern "C"
void
forward_t_Kernel_CPU(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, int NX, int NY) 
{
	for(int tid=0;tid<NY;tid++)
	{
		for(int tidy=0;tidy<NX;tidy++)
		{

			register FL_DBL ARe=get_t_A2CPU(&AY[tid*2],&AX[tidy],NY*2,NX, NOFA);
			register FL_DBL AIm=get_t_A2CPU(&AY[tid*2+1],&AX[tidy],NY*2,NX, NOFA);
			register FL_DBL BRe=BY[tid*2]*BX[tidy];
			register FL_DBL BIm=BY[tid*2+1]*BX[tidy];
			register FL_DBL CRe=CY[tid*2]*CX[tidy];
			register FL_DBL CIm=CY[tid*2+1]*CX[tidy];
			if((tidy<NX-1)&(tidy>0))
			{
				x[tid*2+tidy*NY*2]=mulReCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulReCPU(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulReCPU(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
				x[tid*2+1+tidy*NY*2]=mulImCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulImCPU(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulImCPU(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
			}
			if(tidy==0)
			{
				x[tid*2+tidy*NY*2]=mulReCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulReCPU(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
				x[tid*2+1+tidy*NY*2]=mulImCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulImCPU(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]);
			}
			if(tidy==NX-1)
			{
				x[tid*2+tidy*NY*2]=mulReCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulReCPU(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]);
				x[tid*2+1+tidy*NY*2]=mulImCPU(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulImCPU(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]);
			}
		}
	}
}

__global__ void
forward_Kernel(FL_DBL* AX, FL_DBL* BX, FL_DBL* CX, FL_DBL* AY, FL_DBL* BY, FL_DBL* CY, FL_DBL* y, FL_DBL* x, int NX, int NY) 
{
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if(tid<NY)
	{
		register FL_DBL ARe=getA2(&AY[tid*2],&AX[tidy],NY*2,NX, NOFA);
		register FL_DBL AIm=getA2(&AY[(tid*2)+1],&AX[tidy],NY*2,NX, NOFA);
		register FL_DBL BRe=-BY[tid*2]*BX[tidy];
		register FL_DBL BIm=-BY[tid*2+1]*BX[tidy];
		register FL_DBL CRe=-CY[tid*2]*CX[tidy];
		register FL_DBL CIm=-CY[tid*2+1]*CX[tidy];
		if((tidy<NX-1)&(tidy>0))
		{
			x[tid*2+tidy*NY*2]=(mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulRe(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]));
			x[tid*2+1+tidy*NY*2]=(mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2])+mulIm(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]));
		}

		if(tidy==0)
		{
			x[tid*2+tidy*NY*2]=(mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]));
			x[tid*2+1+tidy*NY*2]=(mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(CRe,CIm,y[tid*2+(tidy+1)*NY*2],y[tid*2+1+(tidy+1)*NY*2]));
		}
		if(tidy==NX-1)
		{
			x[tid*2+tidy*NY*2]=(mulRe(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulRe(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]));
			x[tid*2+1+tidy*NY*2]=(mulIm(ARe,AIm,y[tid*2+tidy*NY*2],y[tid*2+1+tidy*NY*2])+mulIm(BRe,BIm,y[tid*2+(tidy-1)*NY*2],y[tid*2+1+(tidy-1)*NY*2]));
		}

		x[tid*2+tidy*NY*2]=(tid==0)?0:x[tid*2+tidy*NY*2];
		x[tid*2+tidy*NY*2+1]=(tid==0)?0:x[tid*2+tidy*NY*2+1];
		x[tid*2+tidy*NY*2]=(tid==NY-1)?0:x[tid*2+tidy*NY*2];
		x[tid*2+tidy*NY*2+1]=(tid==NY-1)?0:x[tid*2+tidy*NY*2+1];
	}

}


__global__ void
simple_cp_Kernel(FL_DBL* spc, FL_DBL* spc_Nl, int Nt)
{
	const unsigned int tidx = threadIdx.x+blockDim.x*blockIdx.x;
	const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if(tidx<Nt)	
	{
		spc[tidx+tidy*(Nt-2)]=spc_Nl[tidx+tidy*(Nt-2)]/((tidx>0)?FL_DBL(tidx):1);
//		spc[tidx*2+tidy*(Nt)*2]=spc_Nl[tidx*2+tidy*(Nt)*2]/((tidx>0)?FL_DBL(tidx):1);
//		spc[tidx*2+1+tidy*(Nt)*2]=spc_Nl[tidx*2+1+tidy*(Nt)*2]/((tidx>0)?FL_DBL(tidx):1);
	}
}



__global__ void
cpy_Kernel(FL_DBL* y, FL_DBL* x, int NX, int NY) 
{
  const unsigned int tid = threadIdx.x+blockDim.x*blockIdx.x;
  const unsigned int tidy = threadIdx.y+blockDim.y*blockIdx.y;
	if(tid<NY)
	{
			x[tid*2+tidy*NY*2]=y[tid*2+tidy*NY*2];
			x[tid*2+1+tidy*NY*2]=y[tid*2+1+tidy*NY*2];
	}
}

__global__ void __launch_bounds__( MAX_THREADS_PER_BLOCK,MIN_BLOCKS_PER_MP )
simple_cpy_Kernel(FL_DBL* x, FL_DBL* y, FL_DBL CC, FL_DBL BB,int r, int bsx, int nth) 
{
//	__shared__ FL_DBL sdata[4000];
	__shared__ FL_DBL sdata[272+1];
//	x[threadIdx.x+blockDim.x*blockIdx.x]=B+C*y[threadIdx.x+blockDim.x*blockIdx.x]+y[threadIdx.x+blockDim.x*blockIdx.x]*y[threadIdx.x+blockDim.x*blockIdx.x];//+sdata[0];
	int i,adr;
	FL_DBL B=1.0;
	FL_DBL C=0.1;
	FL_DBL v=0.01;

/*
	for(i=0;i<r;i++)
	{
		adr=threadIdx.x+i*nth+bsx*blockIdx.x;
		x[adr]=B+C*y[adr]+(B+C)/(y[adr]*y[adr]+B)+(B+C*B)/(y[adr]+B*B*B)+sdata[0];

//		x[adr]=B+C*y[adr]+(B+C)/(y[adr]*y[adr]+B)+sdata[0];
//		x[adr]=B+C*y[adr]+sdata[0];
//		x[adr]=y[adr];
	}
*/
	//read to shared


//1way
//	for(i=0;i<r;i+=1)	sdata[threadIdx.x+i*nth]=y[threadIdx.x+i*nth+bsx*blockIdx.x];

//2 way - OPTIMAL

	for(i=0;i<r;i+=2)
	{
		sdata[threadIdx.x+i*nth]=y[threadIdx.x+i*nth+bsx*blockIdx.x];
		sdata[threadIdx.x+i*nth+nth]=y[threadIdx.x+i*nth+nth+bsx*blockIdx.x];
	}

/*
//4 way
	for(i=0;i<r;i+=4)
	{
		sdata[threadIdx.x+i*nth]=y[threadIdx.x+i*nth+bsx*blockIdx.x];
		sdata[threadIdx.x+i*nth+nth]=y[threadIdx.x+i*nth+nth+bsx*blockIdx.x];
		sdata[threadIdx.x+i*nth+nth+nth]=y[threadIdx.x+i*nth+nth+nth+bsx*blockIdx.x];
		sdata[threadIdx.x+i*nth+nth+nth+nth]=y[threadIdx.x+i*nth+nth+nth+nth+bsx*blockIdx.x];
	}
*/


	__syncthreads();

/*
	for(i=0;i<r;i+=two)
	{
		adr=threadIdx.x+i*nth+bsx*blockIdx.x;
		x[adr]=B+C*y[adr]+(B+C)/(y[adr]*y[adr]+B)+(B+C*B)/(y[adr]+B*B*B)+sdata[0];
		adr+=nth;
		x[adr]=B+C*y[adr]+(B+C)/(y[adr]*y[adr]+B)+(B+C*B)/(y[adr]+B*B*B)+sdata[0];
	}
*/


	FL_DBL n=0.0;
//1way - OPTIMAL (nodiff)

	for(i=0;i<r;i+=1)
	{
		adr=threadIdx.x+i*nth;
		sdata[adr]=C*sdata[adr];
//		n+=sdata[adr]*sdata[adr];
//		sdata[adr]=v+C*sdata[adr]+(B+C)/(sdata[adr]*sdata[adr]+B)+(B+C*B)/(sdata[adr]+B*B*B)+n*sdata[adr];
	}


/*
//2way
	for(i=0;i<r;i+=2)
	{
		adr=threadIdx.x+i*nth;
		n+=sdata[adr]*sdata[adr];
		sdata[adr]=v+C*sdata[adr]+(B+C)/(sdata[adr]*sdata[adr]+B)+(B+C*B)/(sdata[adr]+B*B*B)+n*sdata[adr];
		adr+=nth;
		n+=sdata[adr]*sdata[adr];
		sdata[adr]=v+C*sdata[adr]+(B+C)/(sdata[adr]*sdata[adr]+B)+(B+C*B)/(sdata[adr]+B*B*B)+n*sdata[adr];
	}
*/





//	for(i=0;i<r;i+=1)	x[threadIdx.x+i*nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth];
	for(i=0;i<r;i+=2)
	{
		x[threadIdx.x+i*nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth];
		x[threadIdx.x+i*nth+nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth+nth];
	}

/*
	for(i=0;i<r;i+=4)
	{
		x[threadIdx.x+i*nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth];
		x[threadIdx.x+i*nth+nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth+nth];
		x[threadIdx.x+i*nth+nth+nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth+nth+nth];
		x[threadIdx.x+i*nth+nth+nth+nth+bsx*blockIdx.x]=sdata[threadIdx.x+i*nth+nth+nth+nth];
	}

*/
}



__global__ void
Simple_cpy_Kernel(FL_DBL* x, FL_DBL* y, FL_DBL CC, FL_DBL BB,int r, int bsx, int nth) 
{
	int adr=threadIdx.x+blockDim.x*blockIdx.x;
	x[adr]=y[adr];
}

//#endif
