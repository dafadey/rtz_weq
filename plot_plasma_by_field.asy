import graph;
import descr;
import fft;

real Adim = 17.1;

real w_o2(real E) {
	real F = abs(E)+1e-7;
	return 6.018834097463142 * exp( - 0.12307858699746754 * log(F) - 0.5573147246191912 / F);
}

real w_n2(real E) {
	real F = abs(E)+1e-7;
	return 9.691581554890364 * exp( - 0.8679102160158225 * log(F) - 0.8183342643785264 / F);
}

real[] filter(real[] a, real freq, real dt) {
	int n = a.length;
	real[] out = new real[n];
	for(int i = 0; i < n; ++i)
		out[i] = a[i];
	qqFFT_freal(n, out);
	real dw = 2 * pi / (n * dt);
	for(int i = 0; i < floor(n/2); ++i) {
		real f = i * dw;
		real fac = exp(-((f-freq)/(.5))^4);
		out[i*2] *= fac;
		out[i*2+1] *= fac;
	}
	qqFFT_freal_1(n, out);
	return out;
}

struct data {
	real[] plasma;
	real[] field;
	real[] src_www_env_sim;
	real[] src_ne_env_sim;
	int n;
	int nsrc;
	descr d;
};

data get_plasma(string path, int file_id, real dens=0) {
	data dat;
	dat.d = read_descr(path+"/descr.txt");
	if(dens==0)
		dens=dat.d.airDensity;
	file inf = input(path+"/pics/res/field_axis"+number(file_id)+".dat").line();
	string[] raw = inf;
	dat.n = raw.length-4;
	int n = dat.n;
	real[] rawfield = new real[n];
	for(int i = 0; i < n; ++i)
		rawfield[i] = (real) raw[i+4];
	dat.field = filter(rawfield, 1., dat.d.dt);
	dat.field[0]=0;
	dat.field[dat.field.length-1]=0;
	
	real[] o2 = new real[n];
	real[] n2 = new real[n];
	o2[0] = 0;
	n2[0] = 0;
	for(int i = 1; i < n; ++i) {
		o2[i] = o2[i-1] + Adim * w_o2(dat.field[i-1]) * dat.d.dt;
		n2[i] = n2[i-1] + Adim * w_n2(dat.field[i-1]) * dat.d.dt;
	}
	for(int i = 0; i < n; ++i) {
		o2[i] = dens*.2 * (1. - exp(-o2[i]));
		n2[i] = dens*.8 * (1. - exp(-n2[i]));
	}

	dat.plasma = new real[n];
	for(int i = 0; i < n; ++i)
		dat.plasma[i] = o2[i]+n2[i];
	write(max(dat.plasma));
	write("-------------------");
	real[] read_sim_src(string srcfilepostfix) {
		real[] res;
		string fullname = path+"/pics/res/src"+number(file_id)+srcfilepostfix+".dat";
		file infs = input(fullname).line();
		string[] raw = infs;
		dat.nsrc = (int) split(raw[0],'\t')[1];
		int nsrc = dat.nsrc;
		for(int i=0;i<nsrc;++i) {
			real v = (((real) split(raw[i+1],'\t')[0])^2 + ((real) split(raw[i+1],'\t')[1])^2)^.5;
			res.push(v);
		}
		return res;
	}
	dat.src_www_env_sim = read_sim_src("_src_www"); 
	dat.src_ne_env_sim = read_sim_src("src_pw"); 
	return dat;
}

string[] amps = {"0.00005", "0.000065", "0.000075", "0.000085", "0.0001", "0.00011", "0.00013", "0.00015", "0.00017", "0.000185", "0.0002", "0.000225", "0.00025", "0.0003"};
//int[] file_ids = {15543, 15543, 15543, 15510, 15411, 15345, 15246, 15246, 15213, 15147, 15114, 14949, 14850, 14685};

int[] file_ids;
for(int i=0;i<amps.length;++i)
	file_ids.push(15708);
	
picture pic, pics;
picture pict;
picture picpt;
picture picet;
picture picst;

pen[] colors = {blue, .5 * green, red};

struct record {
	string name;
	int id;
	void operator init(string _n, int _i) {
		this.name = _n;
		this.id = _i;
	}
};

record[][] records;
{
	record[] r;
	for(int id = 0; id < amps.length; ++id)
		r.push(record("NEW_harm_ampW0="+amps[id]+"_Kerr=0.00754_airDensity=0.0122_", file_ids[id]));
	records.push(r);
}

{
	record[] r;
	for(int id = 0; id < amps.length; ++id)
		r.push(record("NEW_harm_ampW0="+amps[id]+"_Kerr=0.00754_", file_ids[id]));
		//r.push(record("SINGLE_STEP_harm_ampW0="+amps[id]+"_Kerr=0.00754_airDensity=0.0122_", 0));
	records.push(r);
}

real defDens = 0.0122; // 0 for auto
real pows = .5;
for(int plotid=0;plotid<2;++plotid) {
	real color_factor = (2+plotid)/3;
	int color_id = 0;
	guide gpl, gwww, gpw;
	for(int id = 0; id < amps.length; ++id) {
		string amp = amps[id];
		
		//real[] plasma = get_plasma("SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_ampW0="+amp+"_airDensity=0.0122_", 0);
		//data dat = get_plasma("NEW_harm_ampW0="+amp+"_Kerr=0.00754_airDensity=0.0122_", file_ids[id]);
		//data dat = get_plasma("SINGLE_STEP_harm_ampW0="+amp+"_Kerr=0.00754_airDensity=0.0122_", 0);
		data dat = get_plasma(records[plotid][id].name, records[plotid][id].id, defDens);
		int n = dat.plasma.length;
		real res_p = dat.plasma[n-1];
		write(amp+ " : " + string(res_p));
		
		//df/dt = -nu*f + src
		
		real[] f=new real[n];
		f[0]=0;
		real nu = 0.5;
		for(int i=1;i<n;++i)
			f[i] = f[i-1] - nu * f[i-1] + nu * dat.plasma[i];
		for(int i=0;i<n;++i)
			f[i] = dat.plasma[i] - f[i];
		
		real[] src_ne = new real[n];
		for(int i=0; i < n; ++i) src_ne[i] = 0;
		for(int i=0; i<n; ++i)
			src_ne[i] = dat.plasma[i] * dat.field[i];
		real[] src_ne_3 = filter(src_ne, 3., dat.d.dt);
		real src_ne_res = 0;
		for(int i = 0; i < n; ++i)
			src_ne_res += src_ne_3[i]^2 * dat.d.dt;

		real[] src_www = new real[n];
		for(int i = 0; i < n; ++i)
			src_www[i] = 0.00754 * 9. * dat.field[i]^3;
		real[] src_www_3 = filter(src_www, 3., dat.d.dt);
		real src_www_res = 0;
		for(int i = 0; i < n; ++i)
			src_www_res += src_www_3[i]^2 * dat.d.dt;
		
		if(dat.d.airDensity!=0) {
		  dot(pic, (2*log10((real)amp), (res_p/dat.d.airDensity*2.2*1e19)));
			gpl = gpl -- (2*log10((real)amp), (res_p/dat.d.airDensity*2.2*1e19));
		} else {
			dot(pic, (2*log10((real)amp), (res_p/defDens*2.2*1e19)));
			gpl = gpl -- (2*log10((real)amp), (res_p/defDens*2.2*1e19));
		}
		
		dot(pics, (2*log10((real)amp), log10(src_ne_res)), color_factor * purple);
		dot(pics, (2*log10((real)amp), log10(src_www_res)), color_factor * .5*green);
		
		gpw = gpw -- (2*log10((real)amp), log10(src_ne_res));
		gwww = gwww -- (2*log10((real)amp), log10(src_www_res));
		
		if(amp == "0.00013" || amp == "0.00017" || amp == "0.0003") {
		//if(amp == "0.0003") {
			guide g;
			for(int i=floor(n*.15); i<floor(n*.85); ++i)
				g = g -- (i*dat.d.dt, f[i]);
			draw(pict, g, color_factor * colors[color_id % colors.length]);

			guide gpl;
			for(int i=0; i<n; ++i)
				gpl = gpl -- (i*dat.d.dt, dat.plasma[i]);
			draw(picpt, gpl, color_factor * colors[color_id % colors.length]);
			

			guide ge;

			for(int i=floor(n*.15); i<floor(n*.85); ++i) {
				if(dat.field[i] > dat.field[i+1] && dat.field[i] > dat.field[i-1])
					ge = ge -- (i*dat.d.dt, dat.field[i]);
			}
			draw(picet, ge, color_factor * colors[color_id % colors.length]);
			
			{
				guide gs;
				for(int i=floor(n*.15); i<floor(n*.85); ++i) {
					if(src_ne_3[i] > src_ne_3[i-1] && src_ne_3[i] > src_ne_3[i+1])
						gs = gs -- (i*dat.d.dt, abs(src_ne_3[i])^pows);
				}
				draw(picst, gs, color_factor * colors[color_id % colors.length]);
			}
			/*
			{
				guide gs;
				int nsrc = dat.nsrc;
				for(int i=floor(nsrc*.15); i<floor(nsrc*.85); ++i)
					gs = gs -- (i*dat.d.dt*n/nsrc, abs(dat.src_ne_env_sim[i])^pows);
				draw(picst, gs, color_factor * colors[color_id % colors.length]);
			}
			*/
			{
				guide gs;
				for(int i=floor(n*.15); i<floor(n*.85); ++i) {
					if(src_www_3[i] > src_www_3[i-1] && src_www_3[i] > src_www_3[i+1])
						gs = gs -- (i*dat.d.dt, abs(src_www_3[i])^pows);
				}
				draw(picst, gs, color_factor * colors[color_id % colors.length] + dashed);
			}
			/*
			{
				guide gs;
				int nsrc = dat.nsrc;
				for(int i=floor(nsrc*.15); i<floor(nsrc*.85); ++i)
					gs = gs -- (i*dat.d.dt*n/nsrc, abs(dat.src_www_env_sim[i])^pows);
				draw(picst, gs, color_factor * colors[color_id % colors.length] + dashed);
			}
			*/
			
			color_id+=1;
		}
	}
	
	if(gpl != nullpath)
		draw(pic, gpl, color_factor * defaultpen);
	
	draw(pics, gwww, color_factor * .5*green);
	draw(pics, gpw, color_factor * purple);
}

scale(pic, Log, Linear);
scale(pics, Log, Log);
void addpic(picture pic, pair pos = (0,0), pair sz = (7cm, 5cm), string xname = "$I_\omega$", string yname) {
	xaxis(pic, xname, BottomTop, LeftTicks);
	yaxis(pic, yname, LeftRight, RightTicks);
	size(pic, sz.x, sz.y, point(pic, SW), point(pic, NE));
	pic = shift(-point(pic, SW)) * pic;
	add(pic.fit(), pos);

}

addpic(pic, xname = "$I_\omega$", yname = "$n_\Sigma$");
addpic(pics, (10cm, 0), xname = "$I_\omega$", yname = "src");
addpic(pict, (0cm, -6.3cm), (17cm,5cm), xname = "$\tau$", yname = "$n_\Sigma$");
addpic(picet, (0cm, 6cm), (17cm,5cm), xname = "$\tau$", yname = "$E$");
addpic(picst, (0cm, 18.7cm), (17cm,5cm), xname = "$\tau$", yname = "$\mathrm{src}^{0.5}$");
addpic(picpt, (0cm, 12.3cm), (17cm,5cm), xname = "$\tau$", yname = "$n$");

