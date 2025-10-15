#include "rtz_harm_reader.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mutex>

static std::string snum(int num) {
  std::stringstream ss;
  ss << std::setw(7) << std::setfill('0') << num;
  return ss.str();
}

static void get2ints(std::string& s, int& i0, int& i1)
{
  char* buff = const_cast<char*>(s.c_str()); // yep this is aggressive and dirty optimization, harrr!
  int i=0;
  for(;i<s.size();i++) {
    if(buff[i] == '\t') {
      buff[i] = char(0);
      break;
    }
  }
  i++;
  i0 = atoi(buff);
  i1 = atoi(buff+i);
}

//this is code duplication. oh shi...
static void get2floats(std::string& s, float& f0, float& f1)
{
  char* buff = const_cast<char*>(s.c_str());
  int i=0;
  for(;i<s.size();i++) {
    if(buff[i] == '\t') {
      buff[i] = char(0);
      break;
    }
  }
  i++;
  f0 = atof(buff);
  f1 = atof(buff+i);
}

static std::vector<float> get_floats(std::string& s)
{
  std::vector<float> res;
  char* buff = const_cast<char*>(s.c_str());
  int ip=0;
  for(int i=0;i<s.size();i++) {
    if(buff[i] == '\t') {
      buff[i] = char(0);
      res.push_back(atof(buff+ip));
      ip=i+1;
    }
  }
  res.push_back(atof(buff+ip));
  return res;
}

float maxval(const std::vector<float>& vec) {
  float res = .0f;
  for(auto& v : vec)
    res = std::max(v, res);
  return res;
}

float src_record::maxval() const {
  float res = .0;
  for(auto& v : (*this)[2])
    res = std::max(res, v);
  return res;
}

std::array<std::string, 2> split(const std::string& s)
{
  std::array<std::string, 2> res;
  int id=0;
  for(int i=0;i<s.size();i++) {
    if(s.c_str()[i] == '=') {
      id=1;
      continue;
    }
    res[id] += s.c_str()[i];
  }
  return res;
}

float RefIndexAir_nm(float lam) {
  return 0.05792105f/(238.0185f-lam*lam)+0.00167917f/(57.362f-lam*lam);
}

void rtz_harm_reader::init() {
  std::ifstream ifs((path+"/descr.txt").c_str());
  std::string item;
  double dz_dimless;
  while(getline(ifs, item))
  {
    auto items = split(item);
    if(items[0]=="save_interval")
      save_interval = atoi(items[1].c_str());
    if(items[0]=="k0")
      k0 = atof(items[1].c_str());
    if(items[0]=="dz")
      dz_dimless = std::fabs(atof(items[1].c_str()));
  }
  x0 = 1.f / (std::sqrt(2.f) * k0);
  t0 = 1.f / k0;
  z0 = 2.f / k0;
  dz = dz_dimless * z0;
  ifs.close();
  std::cout << "k0=" << k0 << '\n' << "dz=" << dz << '\n';
  records_count=0;

	if(save_interval == -1) { // dynamic
		std::ifstream ifs((path+"/iterations.txt").c_str());
		std::string item;
		while(getline(ifs, item))
			iterations.push_back(atoi(item.c_str()));
		ifs.close();
	  records_count = iterations.size();
	} else {
		//read energies
		const std::filesystem::path pth{(path+"/pics/res/").c_str()};
		// directory_iterator can be iterated using a range-for loop
		for (auto const& dir_entry : std::filesystem::directory_iterator{pth}) {
			std::string filename = dir_entry.path().string();
			if(filename.find("field_energy") == std::string::npos)
				continue;
			iterations.push_back(records_count);
			records_count++;
		}
		records_count--;
	}
	
	records_count = records_count ? records_count : 1;
		
	std::cout << "save interval=" << save_interval << " records count=" << records_count << '\n';

  for(int r=0;r<records_count;r++) {
    std::string filename = path+"/pics/res/"+"field_energy"+snum(iterations[r])+".dat";
    std::ifstream ifs(filename.c_str());
    std::string entry;
    for(int i=0;i<5;i++)
      getline(ifs,entry);
    auto items = get_floats(entry);
    if(items.size()!=10)
      std::cerr << "ERROR reding energy data!\n";
    std::cout << "file:" << filename << '\n';
    e_rec.ew.push_back(items[0]);
    e_rec.ew2.push_back(items[1]);
    e_rec.ew3.push_back(items[2]);
    e_rec.rw.push_back(items[3]);
    e_rec.rw2.push_back(items[4]);
    e_rec.rw3.push_back(items[5]);
    e_rec.sw.push_back(items[6]);
    e_rec.sw2.push_back(items[7]);
    e_rec.sw3.push_back(items[8]);
    e_rec.total.push_back(items[9]);
    ifs.close();
  }
  float _max_src_www = .0f;
  float _max_src_2w2w_w = .0f;
  float _max_src_pw = .0f;
  float _max_w3 = .0f;
  
  www_src_integral.resize(records_count);
  www_src_instant.resize(records_count);
  w2w2w_src_integral.resize(records_count);
  w2w2w_src_instant.resize(records_count);
  pw_src_integral.resize(records_count);
  pw_src_instant.resize(records_count);
  
  std::mutex maxLock;
  
  #pragma omp parallel for
  for(int r=0;r<records_count;r++) {
    auto www_src = get_src(iterations[r], rtz_harm_reader::src_rec_type::WWW);
    auto w2w2w_src = get_src(iterations[r], rtz_harm_reader::src_rec_type::W2W2_W);
    auto pw_src = get_src(iterations[r], rtz_harm_reader::src_rec_type::PW);
    auto w3 = get_src(iterations[r], rtz_harm_reader::src_rec_type::W3);
    int nx = www_src.nx;
    int ny = www_src.ny;
    
    auto integrate = [&nx, &ny](src_record& src_integral, src_record& src_instant) {
      if(src_integral[2].size()==0)
      {
        for(int i=0;i<3;i++)
          src_integral[i].resize(nx * ny, .0f);
        src_integral.nx = nx;
        src_integral.ny = ny;
      }
      for(int i=0;i<2;i++)
      {
        for(int j=0;j<nx*ny;j++)
          src_integral[i][j] += src_instant[i][j];
      }
      for(int j=0;j<nx*ny;j++)
        src_integral[2][j] = std::sqrt(std::pow(src_integral[0][j],2)+std::pow(src_integral[1][j],2));
    };

		src_record www_intgral;
		src_record w2w2w_intgral;
		src_record pw_intgral;

    integrate(www_intgral, www_src);
    integrate(w2w2w_intgral, w2w2w_src);
    integrate(pw_intgral, pw_src);
    
    auto integrate_RT=[](src_record& rec) {
      float res = .0f;
      for(int i=0;i<rec.ny;i++) {
        for(int j=0;j<rec.nx;j++)
          res += std::pow(rec[2][j+i*rec.nx],2) * (float) i;    
      }
      return res;
    };
    
    www_src_integral[r] = integrate_RT(www_intgral);
    w2w2w_src_integral[r] = integrate_RT(w2w2w_intgral);
    pw_src_integral[r] = integrate_RT(pw_intgral);
    www_src_instant[r] = integrate_RT(www_src);
    w2w2w_src_instant[r] = integrate_RT(w2w2w_src);
    pw_src_instant[r] = integrate_RT(pw_src);
    
    maxLock.lock();
    _max_src_www = std::max(_max_src_www, www_src.maxval());
    _max_src_2w2w_w = std::max(_max_src_2w2w_w, w2w2w_src.maxval()) ;
    _max_src_pw = std::max(_max_src_pw, pw_src.maxval()) ;
		_max_w3 = std::max(_max_w3, w3.maxval());
    maxLock.unlock();
  }
  
  max_src_www = _max_src_www;
  max_src_2w2w_w = _max_src_2w2w_w;
  max_src_pw = _max_src_pw;
  max_w3 = _max_w3;
  
  //global max:
  float global_max = std::max(std::max(max_src_www, max_src_2w2w_w), max_src_pw);
  std::cout << "number of files is " << records_count << " max_src_www=" << max_src_www << " max_src_2w2w_w=" << max_src_2w2w_w << " max_src_pw=" << max_src_pw << ", max_w3=" << max_w3 << '\n';
  max_src_www = global_max;
  max_src_2w2w_w = global_max;
  max_src_pw = global_max;
  std::cout << "use global max, so:\n";
  std::cout << "number of files is " << records_count << " max_src_www=" << max_src_www << " max_src_2w2w_w=" << max_src_2w2w_w << " max_src_pw=" << max_src_pw << ", max_w3=" << max_w3 << '\n';

  float timeline_src_instant_max = std::max(std::max(maxval(www_src_instant), maxval(w2w2w_src_instant)), maxval(pw_src_instant));
  float timeline_src_integral_max = std::max(std::max(maxval(www_src_integral), maxval(w2w2w_src_integral)), maxval(pw_src_integral));
  float timeline_ew3_max = maxval(e_rec.ew3);
  
  std::cout << "timeline_ew3_max=" << timeline_ew3_max << " timeline_src_integral_max=" << timeline_src_integral_max << " timeline_src_instant_max=" << timeline_src_instant_max << ", max_w3=" << max_w3 << '\n';
  
  for(int i=0;i<records_count;i++) {
    www_src_instant[i] *= timeline_ew3_max/timeline_src_instant_max;
    w2w2w_src_instant[i] *= timeline_ew3_max/timeline_src_instant_max;
    pw_src_instant[i] *= timeline_ew3_max/timeline_src_instant_max;
    www_src_integral[i] *= timeline_ew3_max/timeline_src_integral_max;
    w2w2w_src_integral[i] *= timeline_ew3_max/timeline_src_integral_max;
    pw_src_integral[i] *= timeline_ew3_max/timeline_src_integral_max;
  }
  
}

src_record rtz_harm_reader::get_src(int record_id, src_rec_type t, double fac) const {
  float norma = .0f;

  std::string filename = path + "/pics/res/src" + snum(record_id);

  if(t == rtz_harm_reader::src_rec_type::WWW) {
    filename += "_src_www";
    norma = -max_src_www;
  }
  if(t == rtz_harm_reader::src_rec_type::W2W2_W) {
    filename += "src_2w2w_w";
    norma = -max_src_2w2w_w;
  }
  if(t == rtz_harm_reader::src_rec_type::PW) {
    filename += "src_pw";
    norma = -max_src_pw;
  }
  if(t == rtz_harm_reader::src_rec_type::W3) {
    filename += "env_w3";
    //norma = 0.001;//max_w3;
    norma = max_w3;
  }
  filename += ".dat";
  std::cout << "filename is " << filename << '\n';

  auto read_src=[&](std::ifstream& ifs, src_record& r)
  {
    std::string entry;
    if(norma == .0)
      norma = 1.;
    else
      norma = 1./norma;
    int ri;
    int tcount;
    
    float z = float(iterations[record_id]) * dz;
    float airPhase = 0.f;//3. * k0 * z * RefIndexAir_nm(.78*3);
    float csap = 1.f;//std::cos(airPhase);
    float snap = 0.f;//std::sin(airPhase);
        
    while(true) {
      getline(ifs, entry);
      get2ints(entry, ri, tcount);
      
      for(int i=0;i<tcount;i++) {
        getline(ifs, entry);
        float s,c;
        get2floats(entry, c, s);

        float cs = c * csap + s * snap;
        float sn = s * csap - c * snap; 
				cs *= fac;
				sn *= fac;
        r[0].push_back(cs*norma);
        r[1].push_back(sn*norma);
        r[2].push_back(std::sqrt(sn*sn+cs*cs)*norma);
      }
      if(ifs.peek() == EOF)
        break;
    }
    
    r.nx = tcount;
    r.ny = ri+1;
  };

  src_record rec;
  std::ifstream ifs(filename.c_str());
  read_src(ifs, rec);
  ifs.close();

  if(!refpath.empty()) {
    src_record ref_rec;
    std::string reffilename = refpath+ "/pics/res/src" + snum(record_id)+"env_w3.dat";
    std::cout << "OPENING REFERENCE " << reffilename << '\n';
    ifs.open(reffilename);
    read_src(ifs, ref_rec);
    ifs.close();
    
    if(ref_rec.nx == rec.nx && ref_rec.ny == rec.ny ) {
      int n = rec[0].size();
      for(int i=0;i<n;i++) {
        double cs = ref_rec[0][i]/ref_rec[2][i];
        double sn = ref_rec[1][i]/ref_rec[2][i];
        double x = rec[0][i];
        double y = rec[1][i];
        rec[0][i] = x * cs + y * sn;
        rec[1][i] = x * sn - y * cs;
      }
    } else
      std::cerr << "ERROR! ref_rec.nx != rec.nx || ref_rec.ny != rec.ny\n";
  }
  
  std::cout << "got src " << record_id << " of " << rec.nx << 'x' << rec.ny << " size, max_src_www=" << max_src_www << ", max_src_2w2w_w=" << max_src_2w2w_w << ", max_src_pw=" << max_src_pw << ", max_w3=" << max_w3 << '\n'; 
  return rec;
}
