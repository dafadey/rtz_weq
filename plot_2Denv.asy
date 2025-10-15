import palette;
import contour;
import graph;
import descr;

real c0 = 3.*10^10; //cm/s
real k0 = 80553.65779; // cm^-1
real x0 = 1. / (2.^.5 * k0);
real t0 = 1 / k0 / c0;
real z0 = 2/k0;
real dz = .5*z0;
real n0 = 1.83485883*10.^21;
real Ea_Vcm = 5.14e9;

defaultpen(fontsize(8pt));

real tmax_fs;
real rmax_um;

real[][] getData(string fn) {
	write(fn);
	file inf = input(fn).word();
	real[] raw = inf;
	pair[][] data;
	int id=0;
	int rowSize=0;
	int colCount=0;
	while(true) {
		int rowId = floor(raw[id]);
		id+=1;
		rowSize = floor(raw[id]);
		//write("rowId=", rowId);
		//write("rowSize=", rowSize);
		id+=1;
		pair[] row = new pair[rowSize];
		for(int i=0; i < rowSize; ++i) {
			row[i] = (raw[id], raw[id+1]);
			id+=2;
		}
		data.push(row);
		colCount+=1;
		if(id==raw.length)
			break;
	}
	
	real[][] env = new real[rowSize][colCount*2-1];
	for(int i=0;i<rowSize;i+=1) {
		pair v = data[0][i];
		real rv = sqrt(v.x^2+v.y^2);
		env[i][colCount-1] = rv;
	}
	for(int j=1;j<colCount;++j) {
		for(int i=0;i<rowSize;++i) {
			pair v = data[j][i];
			real rv = sqrt(v.x^2+v.y^2);
			env[i][colCount-1+j] = rv;
			env[i][colCount-1-j] = rv;
		}
	}
	for(int i=0;i<rowSize;i+=1) {
		env[i][0] = 0;
		env[i][2*colCount-1] = 0;
	}

	return env;
}

triple plot(picture pic, real dr, real T, real[][] env, pen pn = defaultpen) {
	
	//image(pic, env, (0,0), (rowSize, 2*colCount-1), BWRainbow());
	int nlevels=4;
	real[] levels = new real[nlevels];
	real bot = min(env);
	real top = max(env);
	real step = (top-bot)/(nlevels);
	step = (top-bot-1e-2*step)/(nlevels);
	for(int i=0; i<nlevels; ++i)
		levels[i] = bot + step*(i+1);
	
	guide[][] cs = contour(env, (0, -rmax_um), (tmax_fs, rmax_um), levels);
	for(guide[] g : cs)
		filldraw(pic, g, pn+opacity(1/nlevels), pn);
	return (point(cs[nlevels-1][0],0).x, point(cs[nlevels-1][0],0).y, top);
}


void make(string prefix, picture p, string szi) {
	descr d = read_descr(prefix+"/descr.txt");
	real z0_mm = d.dz * z0 * 10;
	real zshift_mm = z0 * 10 * d.z_shift;
	real z_mm = (int) szi * z0_mm + zshift_mm;

	real[][] wdata = getData(prefix + "/pics/res/src"+szi+"env_w.dat");
	real[][] w2data = getData(prefix + "/pics/res/src"+szi+"env_2w.dat");

	rmax_um = d.dr * wdata.length * x0 * 10^4; // use d.dr because each r is dumped
	tmax_fs = d.T * t0 * 10^15; // use T because tau points are dumped with stride

	draw(p, shift(tmax_fs/5,0)*scale(tmax_fs*3/5, 2*rmax_um)*shift(0,-.5)*unitsquare);
	//write(tmax_fs);
	//write(rmax_um);
	triple maxw = plot(p, d.dr, d.T, wdata, red);
	triple max2w = plot(p, d.dr, d.T, w2data, blue);
	dot(p, (maxw.x, maxw.y));
	label(p, rotate(-45)*format("$%.3f$", maxw.z), (maxw.x, maxw.y), SE);//, red, UnFill());
	dot(p, (max2w.x, max2w.y));
	label(p, rotate(45)*format("$%.3f$", max2w.z), (max2w.x, max2w.y), NE);//, blue, UnFill());
	xaxis(p, shift(0,.1cm)*"$\tau$ [fs]", BottomTop, LeftTicks(shift(0, .1cm)*"$%g$"));
	yaxis(p, shift(.4cm,0)*rotate(90)*"$r$ [$\mu$m]", LeftRight, RightTicks(shift(.05cm, 0)*"$%g$"));
	size(p, 5cm, 4cm, point(p, SW), point(p, NE));
	label(p, "z=F" + (z_mm > 0 ? "+":"") + format("$%.2f$", z_mm) + " mm", point(p, NW), 3*SE, UnFill());
}

string[] amps={"0.00013", "0.00015","0.00017", "0.000185", "0.0002", "0.000225", "0.00025", "0.0003"};

for(int i=0; i<amps.length; i+=1) {
	string a = amps[i];
	string prefix = "NEW_harm_ampW0="+a+"_ampW2=0.00005_phase2w=max_Kerr=0.00754_airDensity=0.0122_";

	label("$A_\omega=" + format("$%.3g$", (real) a) + "$ $E_a$", (5.7cm*i+2.5cm, 10cm), N);
	
	picture p0;
	make(prefix, p0, "0014322");
	p0 = shift(-point(p0, SW)) * p0;
	add(p0.fit(), (5.7cm*i,-5cm));

	picture p1;
	make(prefix, p1, "0015576");
	p1 = shift(-point(p1, SW)) * p1;
	add(p1.fit(), (5.7cm*i, 0));

	picture p2;
	make(prefix, p2, "0016995");
	p2 = shift(-point(p2, SW)) * p2;
	add(p2.fit(), (5.7cm*i,5cm));
	
}

//size(5.7cm*8, 15cm);
