#include "postprocess.h"

#include <omp.h>
#include <iostream>
#include <mgl2/mgl.h>

#include "fft_real_3.cpp"
#include "host_methods.h"

void postprocess::prepare_spectrum() {
  
  if(status & postprocess::SPEC_READY)
    return;

  if(!spec)
    spec = new FL_DBL[Nt*Nx];

  if(tmp.size() < Nt*2)
    tmp.resize(Nt*2);
  FL_DBL* sig = tmp.data() + Nt;
	for(int j=0;j<Nx;j++)
	{
		for(int i=0;i<Nt;i++)
			sig[i] = field[i+j*Nt];
		qqFFT_freal(Nt, sig, tmp.data());
		for(int i=0;i<Nt;i++)
			spec[i+j*Nt] = sig[i];
  }
  
  if(!spec_on_aperture)
    spec_on_aperture = new FL_DBL[Nt/2];
  FL_DBL dr = DX / (FL_DBL) Nx;
  for(int j=0;j<Nt/2;j++)
  {
    spec_on_aperture[j] = .0;
    for(int i=0;i<Nx;i++) {
      FL_DBL r = ((FL_DBL)i+.5) * DX / (FL_DBL) Nx;
      spec_on_aperture[j] += (std::pow(spec[j*2+i*Nt],2.) + std::pow(spec[j*2+1+i*Nt],2.)) * 2. * Pi * r * dr; 
    }
  }    

  status |= (int) postprocess::SPEC_READY;
  status &= (int) ~postprocess::CLEAN;
}

void postprocess::filter(FL_DBL** out, double freq, double gap) {
  if(!*out)
    *out = new FL_DBL[Nt*Nx];

  for(int j=0;j<Nx;j++)
	{
    if(tmp.size() < Nt*2)
      tmp.resize(Nt*2);
    FL_DBL* sig = tmp.data() + Nt;
		for(int i=0;i<Nt;i++)
			sig[i] = spec[i+j*Nt];
		sig[0] = .0;		
		sig[1] = .0;
		for(int i=0;i<Nt/2;i++)
		{
			const FL_DBL f = FL_DBL(i)-freq;
			sig[i*2]*=exp(-pow(f/gap,4));
			sig[i*2+1]*=exp(-pow(f/gap,4));
		}

		qqFFT_freal_1(Nt, sig, tmp.data());
		for(int i=0;i<Nt;i++)
			(*out)[i+j*Nt]=sig[i];
	}
}

void postprocess::w2w3w_filter()
{
  if(status & postprocess::FILTERS_READY)
    return;
  prepare_spectrum();
  double freq1 = round(w0/(2.*M_PI/DT));
  double freq2 = round(2.*w0/(2.*M_PI/DT));
  double freq3 = round(3.*w0/(2.*M_PI/DT));
  double gap = 3*freq1/4;
  std::cout << "freq=" << freq1 << '\n';
  filter(&filtered_w,freq1,gap);
  filter(&filtered_2w,freq2,gap);
  filter(&filtered_3w,freq3,gap);
    
  status |= postprocess::FILTERS_READY;
  status &= ~postprocess::CLEAN;
}

void postprocess::reset() {
  status = postprocess::CLEAN;
}

postprocess::postprocess(int _Nt, int _Nx, double _DT, double _DX) : Nt(_Nt), Nx(_Nx), DT(_DT), DX(_DX) {
  field = new FL_DBL[Nt * Nx];
	plasma = new FL_DBL[Nt * Nx];
}

postprocess::~postprocess() {
  #define SAFE_DEL(X) if(X) delete[] X;
  SAFE_DEL(field)
  SAFE_DEL(spec_on_aperture)
  SAFE_DEL(spec)
  SAFE_DEL(filtered_w)
  SAFE_DEL(filtered_2w)
  SAFE_DEL(filtered_3w)
}

void postprocess::saveField() {
  FILE* fp = fopen(signalFilename.c_str(), "w");
  fprintf(fp,"head\n1d\t1\n%d\ndata\n", Nt);
  for(int j = 0; j < Nt; j++)
    fprintf(fp, "%g\n", field[j]);
  fclose(fp);
}

void postprocess::savePlasma() {
  FILE* fp = fopen(plasmaProfileFilename.c_str(), "w");
  fprintf(fp,"head\n1d\t1\n%d\ndata\n", Nx);
  for(int i = 0; i < Nx; i++)
    fprintf(fp, "%g\n", plasma[Nt-1+Nt*i]);
  fclose(fp);
}

void postprocess::saveEnergies() {
  w2w3w_filter();
  std::vector<double> ewr(Nx);
  std::vector<double> e2wr(Nx);
  std::vector<double> e3wr(Nx);
  double ew{};
  double e2w{};
  double e3w{};
  double dr = DX/(double)Nx;
  double dt = DT/(double)Nt;
  for(int i = 0; i < Nx; i++)
  {
    ewr[i] = .0;
    e2wr[i] = .0;
    e3wr[i] = .0;
    double r_2PI = LaplassianMetric(i, dr);
    for(int j = 0; j < Nt; j++)
    {
      double w = pow(filtered_w[j+i*Nt], 2);
      double w2 = pow(filtered_2w[j+i*Nt], 2);
      double w3 = pow(filtered_3w[j+i*Nt], 2);
      
      ewr[i] += w;
      e2wr[i] += w2;
      e3wr[i] += w3;

      ew += w * r_2PI * dr * dt;
      e2w += w2 * r_2PI * dr * dt;
      e3w += w3 * r_2PI * dr * dt;
    }
  }
  
  #define GETR(S,R) double R; \
          for(int i=0;i<Nx;i++) { \
            if(S[i] < S[0] * .5) { \
              R = (double) i * dr; \
              break; \
            } \
          }
  GETR(ewr, rw)
  GETR(e2wr, r2w)
  GETR(e3wr, r3w)
  #undef GETR

  #define GETWIDTH(WIDTH, HARM) double WIDTH = .0; \
          { double norm = .0; \
          for(int i=0;i<Nt/2;i++) { \
            double omega = 2. * M_PI / DT * (double) i; \
            if(std::abs(omega - w0 * HARM) < w0 * .5) { \
              WIDTH += spec_on_aperture[i] * std::pow(omega - w0 * HARM, 2); \
              norm += spec_on_aperture[i]; \
            } \
          } \
          WIDTH = std::sqrt(WIDTH/norm); }
  GETWIDTH(w, 1.);
  GETWIDTH(w2, 2.);
  GETWIDTH(w3, 3.);
  #undef GETWIDTH

  double total = .0;
  for(int i=0;i<Nx;i++)
  {
    double r_2PI = LaplassianMetric(i,dr);
    for(int j=0;j<Nt;j++)
      total += std::pow(field[j+i*Nt], 2.) * r_2PI * dr * dt;
  }

/*
  for(int i=0;i<Nt/2;i++)
    total += spec_on_aperture[i] * 2. * M_PI / DT;
*/
  
  std::cout << "=======" << "w=" << w << " w2=" << w2 << " w3=" << w3 << '\n';
  FILE* fp = fopen(energyFilename.c_str(), "w");
  fprintf(fp,"head\n0d\t10\n1\ndata\n");
  fprintf(fp, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", ew, e2w, e3w, rw, r2w, r3w, w, w2, w3, total);
  fclose(fp);
}

void postprocess::saveSources()
{
	FILE* ofwww{nullptr};
	FILE* of2w2w_w{nullptr};
	FILE* ofpw{nullptr};
	FILE* ofW3env{nullptr};
	FILE* ofw_env{nullptr};
	FILE* of2w_env{nullptr};
	
  if(saveSource) {
	  ofwww = fopen((sourcePrefix+"_src_www.dat").c_str(), "w");
	  of2w2w_w = fopen((sourcePrefix+"src_2w2w_w.dat").c_str(), "w");
	  ofpw = fopen((sourcePrefix+"src_pw.dat").c_str(), "w");
  }

	if(saveW3envelope)
	  ofW3env = fopen((sourcePrefix+"env_w3.dat").c_str(), "w");
  
  if(saveFieldEnvelope) {
	  ofw_env = fopen((sourcePrefix+"env_w.dat").c_str(), "w");
	  of2w_env = fopen((sourcePrefix+"env_2w.dat").c_str(), "w");
  }

	double freq1 = round(w0/(2.*M_PI/DT));
  double freq2 = round(2.*w0/(2.*M_PI/DT));
  double freq3 = round(3.*w0/(2.*M_PI/DT));
  int gap = freq1/2;
	
	if(tmp.size() < Nt*3)
		tmp.resize(Nt*3);
	FL_DBL* sigS = tmp.data() + Nt;
	FL_DBL* sigC = tmp.data() + Nt*2;
	
	auto shiftSpec = [&](FL_DBL central_frequency) {
		qqFFT_freal(Nt, sigS, tmp.data());
		//
		int icf = central_frequency;
		sigS[0]=2.*sigS[icf*2+1];
		sigS[1]=.0;
		sigC[0]=2.*sigS[icf*2];
		sigC[1]=.0;
		for(int f=1;f<gap;f++)
		{
			sigS[f*2] = sigS[(icf+f)*2+1] + sigS[(icf-f)*2+1];
			sigS[f*2+1] = -sigS[(icf+f)*2] + sigS[(icf-f)*2];
			sigC[f*2] = sigS[(icf+f)*2] + sigS[(icf-f)*2];
			sigC[f*2+1] = sigS[(icf+f)*2+1] - sigS[(icf-f)*2+1];
		}
		for(int i=gap*2;i<Nt;i++)
		{
			sigS[i]=.0;
			sigC[i]=.0;
		}
		qqFFT_freal_1(Nt, sigS, tmp.data());
		qqFFT_freal_1(Nt, sigC, tmp.data());
	};
	
  for(int j=0;j<Nx/8;j++)
	{
		const int step = 8;
		
		auto dump=[&](FILE* of) {
			fprintf(of, "%d\t%d\n", j, Nt/step);
			for(int i=0;i<Nt;i+=step)
				fprintf(of, "%g\t%g\n", sigC[i], sigS[i]);
		};
		
		//www channel
		// sin^3 wt = sin wt * sin^2 wt = sin wt  * (cos 2w - 1)/2 = .5 sin wt * cos 2wt =.25 * ( sin (w + 2w) + sin (w - 2w))
    if(ofwww) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = 9. * Kerr * std::pow(filtered_w[i+j*Nt], 3); // 9 is for d^2 / dt^2 -> omega_0^2 = 9
		  shiftSpec(freq3);
		  dump(ofwww);
    }
    
		//2w2w_w channel
		// sin2w*sin2w*sinw=.5(cos 4w-1) * sin w = .25 * sin(3w) + ...
		// sin2w*sin2w*cosw=.5(cos 4w-1) * cos w = .25 * cos(3w) + ...
		if(of2w2w_w) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = 9. * Kerr * std::pow(filtered_2w[i+j*Nt], 2) * filtered_w[i+j*Nt];
		  shiftSpec(freq3);
		  dump(of2w2w_w);
    }
    
		//plasma n * w
		//sin wt * sin 2wt = .5 (cos(3wt)-cos(wt)), but I do not need this 0.5 coefficient since field, field_w field_2w and plasma are not just envelopes but full fields
		if(ofpw) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = plasma[i+j*Nt] * field[i+j*Nt];
		  shiftSpec(freq3);
		  dump(ofpw);
    }
    
		//field w3
		if(ofW3env) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = filtered_3w[i+j*Nt];
		  shiftSpec(freq3);
		  dump(ofW3env);
		}
		
		//field w
		if(ofw_env) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = filtered_w[i+j*Nt];
		  shiftSpec(freq1);
		  dump(ofw_env);			
		}

		//field 2w
		if(of2w_env) {
		  for(int i=0;i<Nt;i++)
			  sigS[i] = filtered_2w[i+j*Nt];
		  shiftSpec(freq2);
		  dump(of2w_env);
		}
	}
	
	auto safe_close = [](FILE* f) {
		if(f)
			fclose(f);
	};
	
	safe_close(ofwww);
	safe_close(of2w2w_w);
	safe_close(ofpw);
	safe_close(ofW3env);
	safe_close(ofw_env);
	safe_close(of2w_env);

}

void postprocess::plot()
{
	//return;
  HMDT hmw, hm2w, hm3w, hmspec2d;
  HMGL gr;
  int stepx=1;
  int stepy=1;
  double MAXw=.0;
  double MAX2w=.0;
  double MAX3w=.0;
  hmw = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  hm2w = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  hm3w = mgl_create_data_size(Nt/stepx,Nx/stepy,1);
  hmspec2d = mgl_create_data_size(Nt/2/3,1,1);
  
  w2w3w_filter();

  double MAXspec = .0;
  for(int j=0;j<Nt/2;j++)
    MAXspec = std::max(MAXspec, (double) spec_on_aperture[j]);

  for(int j=0;j<Nt/2;j++)
    mgl_data_put_val(hmspec2d, spec_on_aperture[j] + MAXspec * 1e-7,j,0,0);

  for(int i=0;i<Nx/stepy;i++)
  {
    for(int j=0;j<Nt/stepx;j++)
    {
      double vw=fabs(filtered_w[i*Nt*stepy+j*stepx]);
      double v2w=fabs(filtered_2w[i*Nt*stepy+j*stepx]);
      double v3w=fabs(filtered_3w[i*Nt*stepy+j*stepx]);
      mgl_data_put_val(hmw, vw,j,i,0);
      mgl_data_put_val(hm2w, v2w,j,i,0);
      mgl_data_put_val(hm3w, v3w,j,i,0);
      MAXw = std::max(MAXw, vw);
      MAX2w = std::max(MAX2w, v2w);
      MAX3w = std::max(MAX3w, v3w);
    }
  }
  
  //std::cout << "MAXw=" << MAXw << " MAX2w=" << MAX2w << " MAX3w=" << MAX3w << '\n' << std::flush;
  mgl_set_color('W',1,1,1);

  gr = mgl_create_graph(1024*2, 768*2);

  //mgl_set_plotfactor(gr, 1.3);

  mgl_subplot(gr,2,2,0,"E_\\omega");
  mgl_title(gr,"E_\\omega","",6);
  mgl_rotate(gr,0,0,0);
  //mgl_set_range_val(gr,'z',0,MAXw);
  mgl_set_range_val(gr,'c',0,MAXw);
  mgl_set_range_val(gr,'y',0,DX*dimx*std::pow(10.,4.));
  mgl_set_range_val(gr,'x',0,DT*dimt*std::pow(10.,15.));
  mgl_surf(gr,hmw,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  //mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry>");
  mgl_box(gr);		

  mgl_subplot(gr,2,2,1,"E_{2\\omega}");
  mgl_title(gr,"E_{2\\omega}","",6);
  mgl_rotate(gr,0,0,0);
  //mgl_set_range_val(gr,'z',0,MAX2w);
  mgl_set_range_val(gr,'c',0,MAX2w);
  mgl_set_range_val(gr,'y',0,DX*dimx*std::pow(10.,4.));
  mgl_set_range_val(gr,'x',0,DT*dimt*std::pow(10.,15.));
  mgl_surf(gr,hm2w,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  //mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry>");
  mgl_box(gr);		

  mgl_subplot(gr,2,2,2,"E_{3\\omega}");
  mgl_title(gr,"E_{3\\omega}","",6);
  mgl_rotate(gr,0,0,0);
  //mgl_set_range_val(gr,'z',0,MAX3w);
  mgl_set_range_val(gr,'c',0,MAX3w);
  mgl_set_range_val(gr,'y',0,DX*dimx*std::pow(10.,4.));
  mgl_set_range_val(gr,'x',0,DT*dimt*std::pow(10.,15.));
  mgl_surf(gr,hm3w,"WRry","");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\tau \\ [fs]",0,"");
  mgl_label(gr,'y',"r \\ [\\mu m]",0,"");
  //mgl_label(gr,'z',"E_\\omega",0,"");
  mgl_colorbar(gr,"WRry>");
  mgl_box(gr);

  mgl_subplot(gr,2,2,3,"|E(\\omega)|^2");
  mgl_title(gr,"|E(\\omega)|^2","",6);
  mgl_rotate(gr,0,0,0);
  mgl_set_func(gr, "", "lg(y)", nullptr, nullptr);
  //std::cout << "MAXspec: " << MAXspec*1e-7 << " : " << MAXspec+MAXspec*1e-7 << '\n';
  mgl_set_range_val(gr,'y',MAXspec*1e-7,MAXspec+MAXspec*1e-7);
  mgl_set_range_val(gr,'x',0,2.*Pi/(DT*dimt*std::pow(10.,15.))*(Nt/2)/3);
  mgl_plot(gr,hmspec2d, "k", "");
  mgl_axis(gr,"","","");
  mgl_label(gr,'x',"\\omega \\ [fs^{-1}]",0,"");
  mgl_label(gr,'y',"|E(\\omega)|^2",0,"");
  mgl_box(gr);


  std::cout << "saving to file " << pictureFilename << '\n' << std::flush;

  mgl_write_jpg(gr,pictureFilename.c_str(),0);
  mgl_delete_graph(gr);

  mgl_delete_data(hmw);
  mgl_delete_data(hm2w);
  mgl_delete_data(hm3w);

}
