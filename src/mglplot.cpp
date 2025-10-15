#include "mglplot.h"

#include <omp.h>
#include <mgl2/mgl.h>

#include "fft_real_3.cpp"

void postprocess::prepare_spectrum() {
  if(!tmp)
  tmp=new FL_DBL[Nt*2];
  FL_DBL* sig = tmp + Nt;
	for(int j=0;j<Nx;j++)
	{
		for(int i=0;i<Nt;i++)
			sig[i+nx] = field[i+j*Nt];
		qqFFT_freal(Nt, sig, tmp);
		for(int i=0;i<Nt;i++)
			spec[i+j*Nt] = sig[i+nx];
  }
}

fil


FL_DBL* filter(FL_DBL* in, int nx, int ny, int freq, int gap)
{
	FL_DBL* out=new FL_DBL[nx*ny];
	FL_DBL* tmp=new FL_DBL[nx*omp_get_num_threads()];
	FL_DBL* sig=new FL_DBL[nx*omp_get_num_threads()];
  #pragma omg prallel for
	for(int j=0;j<ny;j++)
	{
		for(int i=0;i<nx;i++)
			sig[i+omp_get_thread_num()*nx]=in[i+j*nx];
		qqFFT_freal(nx, &sig[omp_get_thread_num()*nx], &tmp[omp_get_thread_num()*nx]);

		sig[omp_get_thread_num()*nx+0]=0;		
		sig[omp_get_thread_num()*nx+1]=0;		
		for(int i=0;i<nx/2;i++)
		{
			FL_DBL f=FL_DBL(i-freq);
			FL_DBL w=FL_DBL(gap);
			sig[omp_get_thread_num()*nx+i*2]*=exp(-pow(f/w,2));
			sig[omp_get_thread_num()*nx+i*2+1]*=exp(-pow(f/w,2));
		}

		qqFFT_freal_1(nx, &sig[omp_get_thread_num()*nx], &tmp[omp_get_thread_num()*nx]);
		for(int i=0;i<nx;i++)
			out[i+j*nx]=sig[i+omp_get_thread_num()*nx];
	}
	delete[] tmp;
	delete[] sig;
	return out;
}

postprocess::postprocess(int _Nt, int _Nx, double _DT, double _DX) : Nt(_Nt), Nx(_Nx), DT(_DT), DX(_DX) {
  field = new FL_DBL[Nt * Nx];
}

postprocess::~postprocess() {
  delete[] field;
}

void postprocess::w2w3w_plot()
{
  HMDT hmw, hm2w, hm3w;
  HMGL gr;
  int stepx=1;
  int stepy=1;
  double MAXw=.0;
  double MAX2w=.0;
  double MAX3w=.0;
  hmw = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  hm2w = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  hm3w = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  FL_DBL* value_w1 = filter(field,Nt,Nx,40,20);
  FL_DBL* value_w2 = filter(field,Nt,Nx,80,20);
  FL_DBL* value_w3 = filter(field,Nt,Nx,120,20);

  for(int i=0;i<Nx/stepy;i++)
  {
    for(int j=0;j<Nt/stepx;j++)
    {
      double vw=fabs(value_w1[i*Nt*stepy+j*stepx]);
      double v2w=fabs(value_w2[i*Nt*stepy+j*stepx]);
      double v3w=fabs(value_w3[i*Nt*stepy+j*stepx]);
      mgl_data_put_val(hmw, vw,j,i,0);
      mgl_data_put_val(hm2w, v2w,j,i,0);
      mgl_data_put_val(hm3w, v3w,j,i,0);
      MAXw = std::max(MAXw, vw);
      MAX2w = std::max(MAX2w, v2w);
      MAX3w = std::max(MAX3w, v3w);
    }
  }
  mgl_set_color('W',1,1,1);

  gr = mgl_create_graph(768*2, 1024*2);

  mgl_subplot(gr,1,3,0,"E_\\omega");
  mgl_title(gr,"E_\\omega","",6);
  mgl_rotate(gr,0,0,0);
  mgl_set_range_val(gr,'z',0,MAXw);
  mgl_set_range_val(gr,'c',0,MAXw);
  mgl_set_range_val(gr,'y',0,400);
  mgl_set_range_val(gr,'x',0,100);
  mgl_surf(gr,hmw,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry<");
  mgl_box(gr);		

  mgl_subplot(gr,1,3,1,"E_{2\\omega}");
  mgl_title(gr,"E_{2\\omega}","",6);
  mgl_rotate(gr,0,0,0);
  mgl_set_range_val(gr,'z',0,MAX2w);
  mgl_set_range_val(gr,'c',0,MAX2w);
  mgl_set_range_val(gr,'y',0,400);
  mgl_set_range_val(gr,'x',0,100);
  mgl_surf(gr,hm2w,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry<");
  mgl_box(gr);		


  mgl_subplot(gr,1,3,2,"E_{3\\omega}");
  mgl_title(gr,"E_{3\\omega}","",6);
  mgl_rotate(gr,0,0,0);
  mgl_set_range_val(gr,'z',0,MAX3w);
  mgl_set_range_val(gr,'c',0,MAX3w);
  mgl_set_range_val(gr,'y',0,400);
  mgl_set_range_val(gr,'x',0,100);
  mgl_surf(gr,hm3w,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry<");
  mgl_box(gr);

  mgl_write_jpg(gr,filename.c_str(),0);
  mgl_delete_graph(gr);

  mgl_delete_data(hmw);
  mgl_delete_data(hm2w);
  mgl_delete_data(hm3w);

}
