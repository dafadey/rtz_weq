import graph;
import patterns;
import descr;

add("crosshatch", crosshatch(.7mm, .8*cyan+.1*blue));
add("crosshatchdense", crosshatch(.3mm, .4*cyan+.1*blue+.5*red));

defaultpen(linewidth(0.3pt));
defaultpen(fontsize(8pt));

real k0=80553.65779; // cm^-1
real x0 = 1. / (2.^.5 * k0);
real t0 = 1. / k0;
real z0 = 2/k0;
real dz = .5*z0;
real n0 = 1.83485883*10.^21;
real Ea_Vcm = 5.14e9;

int[] iterations;

//updates iterations and tells if data path has iterations.txt
bool read_iterations(string data_path) {
	iterations.delete();
	//try to open iterations.txt
	string filename = data_path+"/iterations.txt";
  file inf = input(filename, check = false).word();
  if(!error(inf)) { // file exists read iterations from there
		iterations = inf;
		close(inf);
		iterations.cyclic = true;
		return true;
	} else { // file does not exist build iterations from pics/res/energy*.dat
		int i=0;
		while(true) {
			filename = data_path+"/pics/res/"+"field_energy"+number(i)+".dat";
			file infe = input(filename, check = false);
			if(error(infe))
				break;
			else
				close(infe);
			iterations.push(i);
			++i;
		}
		iterations.cyclic = true;
		return false;
	}
}

void read_point(real[][] data, string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  string[] sinput = split(raw[4], '\t');
  real[] input;
  for(string v : sinput)
    input.push((real) v);
  data.push(input);
}

real get_w3_energy(string data_path, int current = -1) // -1 means whole energy
{
	write("get_w3_energy current=", current);
  write("data_path=\"" + data_path + "\"");
  write("itearations[current]=", iterations[current]);
	string filename = data_path+"/pics/res/"+"field_energy"+number(iterations[current])+".dat";
  write("reading file "+filename);
  real[][] data;
  read_point(data, filename);
  filename = data_path+"/pics/res/"+"field_energy"+number(iterations[0])+".dat";
  //write("reading file "+filename);
  real[][] data0;
  read_point(data0, filename);

  return data[0][2] - data0[0][2];
}

real[][] read_1ddata(string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  int n = (int) raw[2];
  int items_count = (int) split(raw[1], '\t')[1];
  real[][] data;
  for(int i = 0; i < n; ++i) {
    string[] sitems = split(raw[i+4], '\t');
    real[] items;
    for(int j = 0; j < sitems.length; ++j)
      items.push((real) sitems[j]);
    data.push(items);
  }
  return data;
}

pair[][] read_2dcomplexdata(string filename) {
  file inf = input(filename).line();
  string[] raw = inf;
  close(inf);
  //write(raw);
  pair[][] res;
  int lineid=0;
	while(true) {
		int n = (int) split(raw[lineid], '\t')[1];
		lineid += 1;
		pair[] data;
		for(int i = 0; i < n; ++i) {
			string[] sitems = split(raw[lineid], '\t');
			lineid += 1;
			data.push(((real) sitems[0], (real) sitems[1]));
		}
		res.push(data);
		if(lineid == raw.length)
			break;
	}
  return res;
}


real[][] get_plasma(string path) {
  real[][] res;
  for(int i=0; i<iterations.length; ++i) {
    real[][] data = read_1ddata(path+"/pics/res/plasma"+number(iterations[i])+".dat");
    real[] slice = new real[data.length];
    for(int j=0; j<slice.length; ++j)
      slice[j] = data[j][0];
    res.push(slice);
  }
  return res;
}

real[] get_src(string path, string postfix) {
  real[] res;
  int istart = floor(iterations[-1]*.45);
  int iend = floor(iterations[-1]*.55);
  int i0=0; while(iterations[i0] < istart) i0+=1;
  int i1=0; while(iterations[i1] < iend) i1+=1;
  
	descr desc = read_descr(path+"/descr.txt");
  for(int i=i0; i<=i1; ++i) {
		pair[][] data = read_2dcomplexdata(path+"/pics/res/"+"src"+number(iterations[i])+postfix+".dat");
		//write("read cplx data ",i);

		real energy = .0;
		for(int xx=0; xx<data.length; ++xx) {
			real r = xx;
			real DT = desc.T / data[xx].length;
			real au_factor = 0.41 / 2.42e-2;
			real dt_au = DT * au_factor / desc.airDensity^2; // dt in au
			for(int tt=0; tt<data[xx].length; ++tt)
				energy += .5 * abs(data[xx][tt])^2 * dt_au;// * 2. * pi * (3/8+xx) * (desc.dr * x0)^2;
			break;
		}
    res.push(energy);
  }
  return res;
}


int type=0; // 0 -- W3(W,W2), 1 -- W(phase2w)

if(type == 0)
  defaultfilename = "W3W_W2";
if(type == 1)
  defaultfilename = "W3_phase";

string[] ampsW0 = {"0.0001","0.0002","0.0003","0.0004","0.0005","0.0006","0.00065","0.0007"};


string[] phases2w = {"0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9",
                     "1.0", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9",
                     "2.0", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "2.9",
                     "3.0", "3.1"};

picture pic, pic_plasma_r, pic_plasma_z, pic_plasma_axis, pic_plasma_acc, pic_src;
picture pic_w3dynamic0, pic_w3dynamic0_;
picture pic_w3dynamic1, pic_w3dynamic1_;
picture pic_w3dynamic2, pic_w3dynamic2_;
picture pic_w3dynamic3, pic_w3dynamic3_;
picture pic_w3dynamic4, pic_w3dynamic4_;

void plotw2w(string ampW2, string[] phases2w, string[] ampsW0, string postfix = "", string prefix="", pen color = defaultpen) {
  guide gw2;
  if(phases2w.length != ampsW0.length)
    return;
  int n = phases2w.length;
  for(int i=0; i< n; ++i) {
    string ampW0 =  ampsW0[i];
    string phase2w =  phases2w[i];
    string pth = "ampW0="+ampW0+"_ampW2="+ampW2+"_phase2w="+phase2w+postfix+"_";
		read_iterations(pth);
    real v = get_w3_energy(pth);
    descr desc = read_descr(pth+"/descr.txt");
    read_iterations(pth);
    real ampW0_focal = (real) ampW0 * desc.r_0/(desc.dw*.5);
    gw2 = gw2 -- (log10((ampW0_focal/27.49*Ea_Vcm)^2), log10(v+1e-17));
  }
  draw(pic, gw2, color);
  dot(pic, gw2, color);
}

real mylog10(real x)
{
  return x > 1e-5 ? log10(x) : -5;
}

void plotw2wrange(string ampW2, string config, string postfix="", string prefix="", pen color = defaultpen) { //config="0.0005:0.3:2.2 ..." eg ampW:phaseMin:phaseMax
  string[] items = split(config, ' ');
  real [][] ampW0MinAvgMax;
  for(string item : items) {
    string[] tokens = split(item, ':');
    string ampW0 = tokens[0];
    string phase2wMin = tokens[1];
    string phase2wMax = tokens[2];
    
    string pth_www = prefix+"ampW0="+ampW0+postfix+"_";
    read_iterations(pth_www);
    real vwww = get_w3_energy(pth_www);
    string pth_2w_min = prefix+"ampW0="+ampW0+"_ampW2="+ampW2+"_phase2w="+phase2wMin+postfix+"_";
    read_iterations(pth_2w_min);
    real vmin = get_w3_energy(pth_2w_min);
    string pth_2w_max = prefix+"ampW0="+ampW0+"_ampW2="+ampW2+"_phase2w="+phase2wMax+postfix+"_";
    read_iterations(pth_2w_max);
    real vmax = get_w3_energy(pth_2w_max);
    
    descr desc = read_descr(prefix+"ampW0="+ampW0+"_ampW2="+ampW2+"_phase2w="+phase2wMax+postfix+"_/descr.txt");
    real ampW0_focal = desc.ampW * desc.r_0/(desc.dw*.5);
    real [] data = {ampW0_focal/27.49*Ea_Vcm, vmin - vwww, .5 * (vmin + vmax) - vwww, vmax - vwww};
    ampW0MinAvgMax.push(data);
  }
  int n = ampW0MinAvgMax.length;
  guide gw2;
  for(int i=0; i < n; ++i) {
    gw2 = gw2 -- (log10(((real)ampW0MinAvgMax[i][0])^2), mylog10(ampW0MinAvgMax[i][1]));
  }
  for(int i=n-1; i >= 0; --i) {
    gw2 = gw2 -- (log10(((real)ampW0MinAvgMax[i][0])^2), mylog10(ampW0MinAvgMax[i][3]));
  }
  gw2 = gw2 -- cycle;
  pair[] points;
  guide gavg;
  for(int i=0; i < n; ++i) {
    pair pt = (log10(((real)ampW0MinAvgMax[i][0])^2), mylog10(ampW0MinAvgMax[i][2]));
    points.push(pt);
    gavg = gavg -- pt;
  }
  fill(pic, gw2, color+opacity(.17));
  dot(pic, gw2, color);
  draw(pic, gavg, color, marker(scale(0.5mm)*rotate(45)*shift(-.5,-.5)*unitsquare, color, FillDraw(color)));
  dot(pic, gavg, color);
  for(int i=0;i<points.length-1;++i) {
    pair d = points[i+1] - points[i];
    pair avg = .5 * (points[i+1] + points[i]);
    label(pic, string(d.y/d.x,3), avg, color+fontsize(7pt));//, UnFill);
  }
  
  //make a label
  string lbl = "$A_{2\omega}$="+ampW2;
  string[] par = split(postfix, '_');
  for(string p : par)
  {
		write("------>"+p);
		string[] tokens = split(p, '=');
		write(tokens);
		if(tokens[0] == "Kerr")
			lbl+=" $K$="+tokens[1];
		if(tokens[0] == "airDensity")
			lbl+=" $n_0$="+tokens[1];
	}
	label(pic, lbl, points[0], W);
}

void plotw(string ampW2="", string ampsW0s, string postfix="", string prefix="", pen color = defaultpen)
{
  string[] ampsW0 = split(ampsW0s, ' ');
  guide gw0;
  pair[] points;
  for(int i=0; i< ampsW0.length; ++i) {
    string item =  ampsW0[i];
    string[] tokens = split(item, ':');
    real v;
    if(tokens.length==1) { // w
      string pth = prefix+"ampW0="+tokens[0]+postfix+"_";
			read_iterations(pth);
      v = get_w3_energy(pth);
    } else { // 2w
      string pth = prefix+"ampW0="+tokens[0]+"_ampW2="+ampW2+"_phase2w="+tokens[1]+postfix+"_";
      read_iterations(pth);
      real e2w = get_w3_energy(pth);
      
      string pthwww = prefix+"ampW0="+tokens[0]+postfix+"_";
      read_iterations(pthwww);
      real ewww = get_w3_energy(pthwww);
      
      v = e2w - ewww;
    }

    descr desc = read_descr(prefix+"ampW0="+tokens[0]+postfix+"_/descr.txt");
    
    real ampW0_focal = desc.ampW * desc.r_0/(desc.dw*.5);
		if(v>0) {
			points.push((log10((ampW0_focal/27.49*Ea_Vcm)^2), log10(v+1e-17)));
			gw0 = gw0 -- (log10((ampW0_focal/27.49*Ea_Vcm)^2), log10(v+1e-17));
		}
  }
  draw(pic, gw0, color, marker(scale(0.5mm)*rotate(45)*shift(-.5,-.5)*unitsquare, color, FillDraw(color)));
  dot(pic, gw0, color);
  for(int i=0;i<points.length-1;++i) {
    pair d = points[i+1] - points[i];
    pair avg = .5 * (points[i+1] + points[i]);
    label(pic, string(d.y/d.x,3), avg, color+fontsize(7pt));//, UnFill);
  }

  //make a label  
  string lbl;
  string[] par = split(postfix, '_');
  for(string p : par)
  {
		string[] tokens = split(p, '=');
		if(tokens[0] == "Kerr")
			lbl+=" $K$="+tokens[1];
		if(tokens[0] == "airDensity")
			lbl+=" $n_0$="+tokens[1];
	}
	label(pic, lbl, points[0], W);
}

void plot_w3dynamic(picture pic_, string pth, pen p, bool lbl = false) {
	write("plot_w3dynamic at " + pth);
	bool hasIter = read_iterations(pth);
	descr desc = read_descr(pth+"/descr.txt");
	real dz = desc.dz * z0 * (hasIter ? 1 : desc.save_interval);
	real z_shift = desc.z_shift * z0;
	guide g;
	for(int iterId=0; iterId < iterations.length; ++iterId) {
		real v = get_w3_energy(pth, iterId);
		real z = z_shift + dz * iterations[iterId];
		//write((z, v));
		g = g--(z, v);
	}

  real aw_focal = desc.ampW * desc.r_0/(desc.dw*.5);// lens: (real) tokens[0];
  real I_Wcm2 = (aw_focal/27.49*Ea_Vcm)^2; 

	draw(pic_, g, p);
	if(lbl)
		label(pic_, "$"+string(I_Wcm2, 3)+"$ $\mathrm{W}/\mathrm{cm}^2$", (point(pic_, NW).x, point(pic_, NW).y*1.3), SE);
	write("plot w3 dynamic is done");
}


void plotPlasma(string ampW2, string config, string postfix="", string prefix="", pen color = defaultpen) {
  string[] items = split(config, ' ');

  guide g_plasma_axis;
  guide g_plasma_amount;

  for(string item : items) {
    string[] tokens = split(item, ':');
    string pth="";
    if(tokens.length==1) // w
      pth = prefix+"ampW0="+tokens[0]+ampW2+postfix+"_";
    else // 2w
      pth = prefix+"ampW0="+tokens[0]+"_ampW2="+ampW2+"_phase2w="+tokens[1]+postfix+"_";
    read_iterations(pth);
    descr desc = read_descr(pth+"/descr.txt");
    write("plotting plasma from "+pth);
    real[][] p = get_plasma(pth);
    real max_amount = .0;
    int zimax=0;
    for(int zi = 0; zi < p.length; ++zi) {
      real amount = .0;
      for(int ri = 0; ri < p[zi].length; ++ri)
        amount += p[zi][ri] * 2. * pi * (3/8+ri) * (desc.dr * x0)^2;
      if(amount > max_amount) {
        max_amount = amount;
        zimax = zi;
      }
    }
    
    guide g_plasma_z;
    guide g_plasma_r;

    real zmin=10^10;
    real zmax=-zmin;
    for(int zi=0;zi<p.length;++zi)
    {
      real z = desc.z_shift + desc.dz * zi * desc.save_interval;
      zmin = min(zmin, z);
      zmax = max(zmax, z);
      g_plasma_z = g_plasma_z -- (z*z0, p[zi][0]*n0);
    }
    draw(pic_plasma_z, g_plasma_z, color);
    draw(pic_plasma_z, shift(0, desc.airDensity * n0) * ((zmin*z0, .0)--(zmax*z0, .0)), color+dashed);

    real rmax=0;
    real p_amount=0;
    for(int ri=0;ri<p[zimax].length/8;++ri)
    {
      real r = ri * desc.dr;
      rmax = max(rmax, r);
      g_plasma_r = g_plasma_r -- (r * x0, p[zimax][ri]*n0);
      p_amount += r * p[zimax][ri] * desc.dr * 2. * pi;
    }
    draw(pic_plasma_r, g_plasma_r, color);
    draw(pic_plasma_r, shift(0, desc.airDensity * n0) * ((.0, .0)--(rmax*x0, .0)), color+dashed);
    
    
    real aw_focal = desc.ampW * desc.r_0/(desc.dw*.5);// lens: (real) tokens[0];
    g_plasma_axis = g_plasma_axis -- (log10((aw_focal/27.49*Ea_Vcm)^2), log10(p[zimax][0]*n0)); 
    g_plasma_amount = g_plasma_amount -- (log10((aw_focal/27.49*Ea_Vcm)^2), log10(p_amount*n0/*/2.2e4*/)); 
  }
  draw(pic_plasma_axis, g_plasma_axis, color);
  draw(pic_plasma_acc, g_plasma_amount, color+dashed);
}


void plotSRC(string ampW2, string config, string postfix="", string prefix="", string src_type, pen color = defaultpen) {
  string[] items = split(config, ' ');

  guide g_src;
  
  for(string item : items) {
    string[] tokens = split(item, ':');
    string pth="";
    if(tokens.length==1) // w
      pth = prefix+"ampW0="+tokens[0]+ampW2+postfix+"_";
    else // 2w
      pth = prefix+"ampW0="+tokens[0]+"_ampW2="+ampW2+"_phase2w="+tokens[1]+postfix+"_";
    read_iterations(pth);
    descr desc = read_descr(pth+"/descr.txt");
    write("plotting src from "+pth);
    real[] p = get_src(pth, src_type);
    real max_p = .0;
    int zimax=0;
    for(int zi = 0; zi < p.length; ++zi) {
      real amount = .0;
      if(p[zi] > max_p) {
        max_p = p[zi];
        zimax = zi;
      }
    }
    
    real aw_focal = desc.ampW * desc.r_0/(desc.dw*.5);// lens: (real) tokens[0];
    g_src = g_src -- (log10((aw_focal/27.49*Ea_Vcm)^2), log10(p[zimax])); 
  }
  draw(pic_src, g_src, color);
}



int[][] getItemsFromRangesString(string ranges_string) {
  int[][] res;
  
  int[] irange;
  if(ranges_string=="") {
    for(int i = 0; i < phases2w.length; ++i)
      irange.push(i);
    res.push(irange);
    return res;
  }
    
  string[] ranges = split(ranges_string, ' ');
  for(string range : ranges) {
    string[] items = split(range, '-');
    if(items.length==2) {
      if(irange.length != 0) {
        res.push(copy(irange));
        irange.delete();
      }
      for(int i = (int) items[0]; i <= (int) items[1]; ++i) {
        irange.push(i);
      }
      res.push(copy(irange));
      irange.delete();
    }
    else {
      int i = (int)items[0];
      if(irange.length==0 || i == irange[irange.length-1]+1)
        irange.push(i);
      else {
        res.push(copy(irange));
        irange.delete();
        irange.push(i);
      }
    }
  }
  if(irange.length!=0)
    res.push(irange);
  write("RANGES:");
  write(res);
  return res;
}

void plotw2wphase(string ampW0, string ampW2, string ranges="", string postfix="", string prefix="", pen color = defaultpen) {
  int[][] items_collections = getItemsFromRangesString(ranges);
  real maxW3=0;
  int imaxW3=0;
  real minW3=1e11;
  int iminW3=0;
  for(int[] items : items_collections) {
    guide g;
    for(int i : items) {
      real phase = (real) phases2w[i];
      string pth = prefix+"ampW0="+ampW0+"_ampW2="+ampW2+"_phase2w="+phases2w[i]+postfix+"_";
      read_iterations(pth);
      real v = get_w3_energy(pth);
      if(v>maxW3) {
        maxW3=v;
        imaxW3=i;
      }
      if(v<minW3) {
        minW3=v;
        iminW3=i;
      }
      g = g -- (phase, log10(v));
    }
    draw(pic, g, color, marker(scale(0.5mm)*rotate(45)*shift(-.5,-.5)*unitsquare, color, FillDraw(color)));
    //dot(pic, g, color);    
  
  }
  picture picmax;
  draw(picmax, scale(.7mm)*unitcircle, color);
  label(picmax, rotate(45) * shift(.7mm+8pt,0) * phases2w[imaxW3], color);
  add(pic, picmax.fit(), ((real)phases2w[imaxW3], log10(maxW3)));
  picture picmin;
  draw(picmin, scale(.7mm)*unitcircle, color);
  label(picmin, rotate(-45) * shift(.7mm+8pt,0) * phases2w[iminW3], color);
  add(pic, picmin.fit(), ((real)phases2w[iminW3], log10(minW3)));
}

if(type==0) {
	/*
  
  //plotw2wrange("0.0002", "0.0001:1.2:2.8 0.0002:1.2:2.7 0.0003:1.1:2.7 0.0004:1.0:2.5 0.0005:0.8:2.4 0.0006:0.7:2.2 0.00065:0.6:2.1", .8*cyan+.1*blue);
  
  //plotw2wrange("0.0003", "0.0001:0.8:2.4 0.0002:0.8:2.4 0.0003:0.7:2.3 0.0004:0.6:2.2 0.0005:0.5:2.0 0.0006:0.3:1.9 0.00065:0.2:1.9", blue);

  //plotw2wrange("0.0004", "0.0001:0.3:1.9 0.0002:0.2:1.8 0.0003:0.2:1.8 0.0004:0.1:1.7 0.0005:3.1:1.5 0.0006:2.9:1.4", .5 * green);
  
  //plotw2wrange("0.0005", "0.0001:1.3:2.8 0.0002:1.2:2.7 0.0003:1.2:2.8 0.0004:1.1:2.6 0.0005:0.9:2.3", red);

  //plotw("", "0.00005 0.000075 0.0001 0.000125 0.00015 0.0002 0.0003 0.0004 0.0005 0.0006 0.00065 0.0007");
  
  //plotw("", "0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.0011", "_airDensity=0.0137_tunExp=3", dashed);
  
  
  //plotw("0.0002", "0.0001:2.8 0.0002:2.8 0.0003:2.7 0.0004:2.1 0.0005:1.4 0.0006:1.4 0.0007:1.5 0.0008:1.5", "_airDensity=0.0137_tunExp=3", dashed+cyan);
  
  

  //plotPlasma("", "0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.0011" , "_airDensity=0.0137_tunExp=3");
  
  //plotPlasma("0.0002", "0.0001:2.8 0.0002:2.8 0.0003:2.7 0.0004:2.1 0.0005:1.4 0.0006:1.4 0.0007:1.5 0.0008:1.5" , "_airDensity=0.0137_tunExp=3", cyan);
  
  //plotPlasma("0.0002", "0.0001:1.2:2.8 0.0002:1.2:2.8 0.0003:1.2:2.7 0.0004:1.0:2.6 0.00045:0.3:2.1 0.0005:0.0:1.4 0.0006:3.0:1.4 0.00065:3.0:1.4 0.0007:3.0:1.5 0.0008:3.0:1.5 0.0009:3.0:1.5 0.001:3.0:1.5" , "_airDensity=0.0137_tunExp=3", .8*cyan+.1*blue+dashed);


  
  plotw2wrange("0.00005", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);

//NEW_harm_Nt=2048_TAU=70_T=500_


  //plotw2wrange("0.0002", "0.0001:1.2:2.8 0.0002:1.2:2.8 0.0003:1.2:2.7 0.0004:1.0:2.6 0.00045:0.3:2.1 0.0005:0.0:1.4 0.0006:3.0:1.4 0.00065:3.0:1.4 0.0007:3.0:1.5 0.0008:3.0:1.5 0.0009:3.0:1.5 0.001:3.0:1.5" , "_airDensity=0.0137_tunExp=3", .8*cyan+.1*blue);

  //plotPlasma("", "0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001" , "_airDensity=0.0137_tunExp=3", orange);

  //plotw("" ,"0.0001 0.0002 0.0003 0.0004 0.00045 0.0005 0.0006 0.00065 0.0007 0.0008 0.0009 0.001", "_airDensity=0.0137_tunExp=3", dashed+orange);
 
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_airDensity=0.0122", "NEW_harm_", dashed+purple);
  
  
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", dashed+cyan);
  
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "", "NEW_harm_", black+dotted);

  write("----------------------------");

  //plotPlasma("_ampW2=0.00005_phase2w=max", "0.00005 0.000075 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);
  //plotPlasma("_ampW2=0.00005_phase2w=min", "0.00005 0.000075 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", orange);
  
  //plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", "src_2w2w_w", .5*green+dashed);


  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", "src_pw", .7*magenta);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", "src_pw", .7*magenta+dashed);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", "_src_www", .5*green);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", "_src_www", .5*green+dashed);

  plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", purple);
  
  plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", cyan);
  */



/*
	lastFile=92;
  
  plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", cyan);

  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_airDensity=0.0122", "NEW_harm_Nt=2048_TAU=70_T=500_", dashed+cyan);

	lastFile=0;
  plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_", purple);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_", "src_pw", .7*magenta+dashed);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_", "_src_www", .5*green+dashed);
*/








/*
	plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_", "src_pw", .7*magenta);
  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_", "_src_www", .5*green);

  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_Kerr=0.00754_", "src_pw", .7*magenta+dashed);
  plotSRC("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "SINGLE_STEP_harm_Nt=2048_TAU=70_T=500_Kerr=0.00754_", "_src_www", .5*green+dashed);
*/

//  plotw2wrange("0.000025", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.00754_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);

//  plotw2wrange("0.0000125", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.00754_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);

//  plotw2wrange("0.00005", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.017_airDensity=0.0122", "NEW_harm_", .4*green+.7*red);

  plotw2wrange("0.00005", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.00754_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);

  plotw2wrange("0.000025", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.00754_airDensity=0.0122", "NEW_harm_", .8*blue+.7*red);

  //plotw2wrange("0.00005", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_airDensity=0.0122", "NEW_harm_", .8*blue+.7*green);

//  plotw2wrange("0.00005", "0.00005:min:max 0.000075:min:max 0.0001:min:max 0.00013:min:max 0.00015:min:max 0.00017:min:max 0.000185:min:max 0.0002:min:max 0.000225:min:max 0.00025:min:max 0.0003:min:max" , "_Kerr=0.0377_airDensity=0.0122", "NEW_check_harm_", .1*blue+.8*green);

  //plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_airDensity=0.0122", "NEW_harm_", .7*green);

	//plotPlasma("", "0.00005 0.000065 0.000075 0.000085 0.0001 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003" , "_Kerr=0_airDensity=0.0122", "NEW_harm_", blue);

  //plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_airDensity=0.0122", "NEW_harm_", .7*green+opacity(0.5));
  //plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_Kerr=0.0377_airDensity=0.0122", "NEW_check_harm_", .7*blue+opacity(0.5));
  
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_Kerr=0.00754_airDensity=0.0122", "NEW_harm_", .7*green+dashed);
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_Kerr=0_airDensity=0.0122", "NEW_harm_", blue);
  //plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "", "NEW_harm_", red);
  plotw("" ,"0.00005 0.000065 0.000075 0.000085 0.0001 0.00011 0.00013 0.00015 0.00017 0.000185 0.0002 0.000225 0.00025 0.0003", "_Kerr=0.00754", "NEW_harm_", red+dashed);

/*
	plot_w3dynamic(pic_w3dynamic0, "NEW_harm_ampW0=0.000065_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic0, "NEW_harm_ampW0=0.000065_", red, true);

	plot_w3dynamic(pic_w3dynamic1, "NEW_harm_ampW0=0.0001_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic1, "NEW_harm_ampW0=0.0001_", red, true);

	plot_w3dynamic(pic_w3dynamic2, "NEW_harm_ampW0=0.00013_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic2, "NEW_harm_ampW0=0.00013_", red, true);

	plot_w3dynamic(pic_w3dynamic3, "NEW_harm_ampW0=0.000185_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic3, "NEW_harm_ampW0=0.000185_", red, true);

	plot_w3dynamic(pic_w3dynamic4, "NEW_harm_ampW0=0.00025_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic4, "NEW_harm_ampW0=0.00025_", red, true);
*/
	plot_w3dynamic(pic_w3dynamic0_, "NEW_harm_ampW0=0.000065_Kerr=0.00754_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic0_, "NEW_harm_ampW0=0.000065_Kerr=0.00754_", red, true);

	plot_w3dynamic(pic_w3dynamic1_, "NEW_harm_ampW0=0.0001_Kerr=0.00754_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic1_, "NEW_harm_ampW0=0.0001_Kerr=0.00754_", red, true);

	plot_w3dynamic(pic_w3dynamic2_, "NEW_harm_ampW0=0.00013_Kerr=0.00754_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic2_, "NEW_harm_ampW0=0.00013_Kerr=0.00754_", red, true);

	plot_w3dynamic(pic_w3dynamic3_, "NEW_harm_ampW0=0.000185_Kerr=0.00754_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic3_, "NEW_harm_ampW0=0.000185_Kerr=0.00754_", red, true);

	plot_w3dynamic(pic_w3dynamic4_, "NEW_harm_ampW0=0.00025_Kerr=0.00754_airDensity=0.0122_", .7*green);
	plot_w3dynamic(pic_w3dynamic4_, "NEW_harm_ampW0=0.00025_Kerr=0.00754_", red, true);

  //string[] amps = {"0.00005", "0.000075", "0.0001", "0.00013", "0.00015", "0.00017", "0.000185", "0.0002", "0.000225", "0.00025", "0.0003"};
  string[] amps = {"0.00013", "0.00015", "0.00017", "0.000185", "0.0002", "0.000225", "0.00025", "0.0003"};
  real bot = point(pic, SW).y;
  real top = point(pic, NE).y;
  for(string a : amps) {
    descr desc = read_descr("NEW_harm_"+"ampW0="+a+"_Kerr=0.00754_airDensity=0.0122"+"_/descr.txt");
    real ampW0_focal = desc.ampW * desc.r_0/(desc.dw*.5);
		real v = log10((ampW0_focal/27.49*Ea_Vcm)^2);
		draw(pic,(v, bot)--(v, top), gray+dashed+opacity(0.333));
		label(pic, rotate(90)*("$A_\omega=$"+a), (v,bot), N, darkgrey, UnFill());
	}
	
	//add margins
	real margin = 0.025;
	pair pic_dim = point(pic,NE)-point(pic,SW);
	pair pic_origin = point(pic,SW);
	draw(pic, shift(pic_origin - pic_dim*margin)*scale(pic_dim.x * (1+2*margin), pic_dim.y * (1+2*margin))*unitsquare);
}
if(type==1) {
  string postfix="";
  
  string ampW2="0.0002";
  postfix = "_airDensity=0.0137_tunExp=3";
  plotw2wphase("0.0001", ampW2, "11-13 27-29", postfix, cyan);
  plotw2wphase("0.0002", ampW2, "11-13 26-29", postfix, blue);
  plotw2wphase("0.0003", ampW2, "10-14 26-28", postfix, 0.7*blue+.3* green);
  plotw2wphase("0.0004", ampW2, "9-11 24-27", postfix, .5*green);
  plotw2wphase("0.00045", ampW2, "0-9 11 12 16-28", postfix, .5*green+0.3*red);
  plotw2wphase("0.0005", ampW2, "0-2 7-19 24-31", postfix, red*0.7+0.3*green);
  plotw2wphase("0.0006", ampW2, "0-2 5-15 21-31", postfix, red);
  plotw2wphase("0.00065", ampW2, "0 1 6 12-16 18 21 24 28-31", postfix, orange);
  plotw2wphase("0.0007", ampW2, "0 4 8 12 14-16 20 24 28-31", postfix, black);
  plotw2wphase("0.0008", ampW2, "0 4 8 14-16 24 29-31", postfix, gray);
  plotw2wphase("0.0009", ampW2, "0 8 14-16 24 29-31", postfix, blue);
  plotw2wphase("0.001", ampW2, "0 8 14-16 24 29-31", postfix, red);
  
/*  
  string ampW2="0.0002";
  plotw2wphase("0.0001", ampW2, "7-13 23-29", cyan);
  plotw2wphase("0.0002", ampW2, "7-13 23-28", blue);
  plotw2wphase("0.0003", ampW2, "6-12 22-28", 0.7*blue+.3* green);
  plotw2wphase("0.0004", ampW2, "5-12 21-26", .5*green);
  plotw2wphase("0.0005", ampW2, "4-9 19-25", red*0.7+0.3*green);
  plotw2wphase("0.0006", ampW2, "2-8 18-23", red);
  plotw2wphase("0.00065", ampW2, "1-7 18-23", orange);
*/
/*
  string ampW2="0.0003";
  plotw2wphase("0.0001", ampW2, "2-4 7-9 17-20 22-25", cyan);
  plotw2wphase("0.0002", ampW2, "2-4 7-9 16-19 22-25", blue);
  plotw2wphase("0.0003", ampW2, "1-3 5-8 16-19 22-25", 0.7*blue+.3* green);
  plotw2wphase("0.0004", ampW2, "0-1 5-7 15-18 21-23", .5*green);
  plotw2wphase("0.0005", ampW2, "0-31", red*0.7+0.3*green);
  plotw2wphase("0.0006", ampW2, "0-6 17-23", red);
  plotw2wphase("0.00065", ampW2, "0-31", orange);
*/
/*
  string ampW2="0.0004";
  plotw2wphase("0.0001", ampW2, "2-4 16-23", cyan);
  plotw2wphase("0.0002", ampW2, "0-5 14-21", blue);
  plotw2wphase("0.0003", ampW2, "0-4 13-20 30-31", 0.7*blue+.3* green);
  plotw2wphase("0.0004", ampW2, "0-2 10 14-21 30-31", .5*green);
  plotw2wphase("0.0005", ampW2, "0-1 6-21 29-31", red*0.7+0.3*green);
  plotw2wphase("0.0006", ampW2, "0-18 23-31", red);
*/
/*  
  string ampW2="0.0005";
  plotw2wphase("0.0001", ampW2, "11-14 27-30", cyan);
  plotw2wphase("0.0002", ampW2, "11-14 26-29", blue);
  plotw2wphase("0.0003", ampW2, "9-15 26-30", 0.7*blue+.3* green);
  plotw2wphase("0.0004", ampW2, "7-14 23-30", .5*green);
  plotw2wphase("0.0005", ampW2, "0-31", red*0.7+0.3*green);
*/
  defaultfilename+="_W2="+ampW2+postfix;
}



if(type==0) {
  scale(pic_plasma_axis, Log, Log);
  xaxis(pic_plasma_axis, "$I$ [$\mathrm{W}/\mathrm{cm}^2$]", BottomTop, LeftTicks);
  yaxis(pic_plasma_axis, "$n$ [$\mathrm{cm}^{-3}$]", LeftRight, RightTicks);
  size(pic_plasma_axis, 4.5cm, 4.9cm, point(pic_plasma_axis, SW), point(pic_plasma_axis, NE));

  scale(pic_plasma_acc, Log, Log);
  xaxis(pic_plasma_acc, "$I$ [$\mathrm{W}/\mathrm{cm}^2$]", BottomTop, LeftTicks);
  yaxis(pic_plasma_acc, "$n$ [$\mathrm{cm}^{-1}$]", LeftRight, RightTicks);
  size(pic_plasma_acc, 4.5cm, 4.9cm, point(pic_plasma_acc, SW), point(pic_plasma_acc, NE));

  scale(pic_src, Log, Log);
  xaxis(pic_src, "$I$ [$\mathrm{W}/\mathrm{cm}^2$]", BottomTop, LeftTicks);
  yaxis(pic_src, "$\mathrm{SRC}$", LeftRight, RightTicks);
  size(pic_src, 4.5cm, 4.9cm, point(pic_src, SW), point(pic_src, NE));

  scale(pic, Log, Log);
  xaxis(pic, "$I$ [$\mathrm{W}/\mathrm{cm}^2$]", BottomTop, LeftTicks);
  yaxis(pic, "$W_{3\omega}$", LeftRight, RightTicks);
  size(pic, 17cm, 13cm, point(pic, SW), point(pic, NE));

  xaxis(pic_plasma_z, "$z$ [cm]", BottomTop, LeftTicks);
  yaxis(pic_plasma_z, "$n$ [$\mathrm{cm}^{-3}$]", LeftRight, RightTicks);
  size(pic_plasma_z, 7.5cm, 5cm, point(pic_plasma_z, SW), point(pic_plasma_z, NE));

  xaxis(pic_plasma_r, "$r$ [cm]", BottomTop, LeftTicks);
  yaxis(pic_plasma_r, "$n$ [$\mathrm{cm}^{-3}$]", LeftRight, RightTicks);
  size(pic_plasma_r, 7.5cm, 5cm, point(pic_plasma_r, SW), point(pic_plasma_r, NE));

	void scale_dynamic_pic(picture pic_dyn, bool tcks = false) {
		ylimits(pic_dyn, 0, point(pic_dyn, NE).y*1.03);
		xaxis(pic_dyn, tcks ? "$z$ [$cm$]" : "", BottomTop, tcks ? LeftTicks : LeftTicks("$$"));
		yaxis(pic_dyn, "$W_{3\omega}$", LeftRight, RightTicks);
		size(pic_dyn, 4.7cm, 2.5cm, point(pic_dyn, SW), point(pic_dyn, NE));
	}	

	scale_dynamic_pic(pic_w3dynamic0, true);
	scale_dynamic_pic(pic_w3dynamic1);
	scale_dynamic_pic(pic_w3dynamic2);
	scale_dynamic_pic(pic_w3dynamic3);
	scale_dynamic_pic(pic_w3dynamic4);
	
	
	scale_dynamic_pic(pic_w3dynamic0_, true);
	scale_dynamic_pic(pic_w3dynamic1_);
	scale_dynamic_pic(pic_w3dynamic2_);
	scale_dynamic_pic(pic_w3dynamic3_);
	scale_dynamic_pic(pic_w3dynamic4_);

}
if(type==1) {
  //ylimits(pic, 0);
  xlimits(pic, 0, 3.1);
  scale(pic, Linear, Log);
  xaxis(pic, "$\varphi$ [rad]", BottomTop, LeftTicks);
  yaxis(pic, "$W_{3\omega}$", LeftRight, RightTicks);
  size(pic, 17cm, 13cm, point(pic, SW), point(pic, NE));
}

void myaddpic(picture p, pair pos=(.0, .0)) {
  p = shift(-point(p,SW))*p;
  add(p.fit(), pos);
}

myaddpic(pic_plasma_axis, (0cm,14.2cm));

myaddpic(pic_plasma_acc, (6.1cm,14.2cm));

myaddpic(pic_src, (12.2cm,14.2cm));

myaddpic(pic);

myaddpic(pic_w3dynamic0, (18.4cm, 0cm));
myaddpic(pic_w3dynamic1, (18.4cm, 2.5cm));
myaddpic(pic_w3dynamic2, (18.4cm, 5cm));
myaddpic(pic_w3dynamic3, (18.4cm, 7.5cm));
myaddpic(pic_w3dynamic4, (18.4cm, 10cm));

myaddpic(pic_w3dynamic0_, (24.2cm, 0cm));
myaddpic(pic_w3dynamic1_, (24.2cm, 2.5cm));
myaddpic(pic_w3dynamic2_, (24.2cm, 5cm));
myaddpic(pic_w3dynamic3_, (24.2cm, 7.5cm));
myaddpic(pic_w3dynamic4_, (24.2cm, 10cm));


myaddpic(pic_plasma_z, (0cm,-6cm));

myaddpic(pic_plasma_r, (9.5cm,-6cm));

write("writing output to "+defaultfilename);
