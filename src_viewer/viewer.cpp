#include <iostream>
#include <vector>
#include <cmath>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include "picture_widget.h"
#include "graph_widget.h"
#include "3Dgraph_widget.h"
#include "tex.h"
#include "rtz_harm_reader.h"

//std::vector<widget*> widgets;
int winNx=1700;
int winNy=768;

struct dashboard {
  std::vector<widget*> widgets;
  picture_widget* pic_src_www{};
  picture_widget* pic_src_w2w2_w{};
  picture_widget* pic_src_pw{};
	picture_widget* pic_w3{};
  picture_widget* pic_src_www_rel{};
  picture_widget* pic_src_w2w2_w_rel{};
  picture_widget* pic_src_pw_rel{};
	
  graph_widget* timeline_view{};
  graph3D_widget* phase_view{};
  rtz_harm_reader data;
  src_record src_rec, src2w_rec, srcpw_rec, w3_rec;
  src_record src_rec_rel, src2w_rec_rel, srcpw_rec_rel, w3_rec_rel;
  
  double getPercentage(int i) const {
		return (double) data.iterations[i<data.iterations.size()-1 ? i+1 : i] / (double) data.iterations.back(); 
	}
  
  int getIdbyPercentage(double p) const {
		int i=0;
		for(;i<data.iterations.size()-1;i++) {
			double p0 = getPercentage(i);
			double p1 = getPercentage(i+1);
			if(p>=p0 && p<=p1)
				return i;
		}
		return i;
	}
  
  std::vector<float> xs;
  std::vector<RGB> colors_www;
  std::vector<RGB> colors_w2w2_w;
  std::vector<RGB> colors_pw;
  std::vector<RGB> colors_w3;

  int last_record_id{7};
  void set_colors_and_xs(int n) {
    if(xs.size() != n) {
      xs.resize(n);
      for(int i=0;i<n;i++)
        xs[i] = 3.*(float(i)/float(n) - .5f);
    }

    auto set_colors=[&n](std::vector<RGB>& colors, RGB (*f) (float)) {
      if(colors.size() != n) {
        colors.resize(n);
        for(int i=0;i<n;i++) {
          float x=(float) i / (float) n;
          colors[i] = f(x);
        }
      }
    };
    set_colors(colors_www, [](float v) {return RGB{v,1.f-v,.0f};});
    set_colors(colors_w2w2_w, [](float v) {return RGB{.0f,v,1.f-v};});
    set_colors(colors_pw, [](float v) {return RGB{v,.0f,1.f-v};});
    set_colors(colors_w3, [](float v) {return RGB{.3f+v*.4f,.3f+v*.4f,.3f+v*.4f};});
  }
  
  void upate_section() {
    int wsec = pic_src_www_rel->hline * src_rec_rel.ny;
    int w2sec = pic_src_w2w2_w_rel->hline * src2w_rec_rel.ny;
    int wpsec = pic_src_pw_rel->hline * srcpw_rec_rel.ny;
    set_colors_and_xs(src_rec_rel.nx);
    phase_view->add_graph("www", src_rec_rel.nx, xs.data(), src_rec[0].data() + wsec * src_rec_rel.nx, src_rec[1].data() + wsec * src_rec_rel.nx, colors_www.data());
    phase_view->add_graph("w2w2_w", src_rec_rel.nx, xs.data(), src2w_rec_rel[0].data() + w2sec * src2w_rec_rel.nx, src2w_rec_rel[1].data() + w2sec * src2w_rec_rel.nx, colors_w2w2_w.data());
    phase_view->add_graph("pw", src_rec_rel.nx, xs.data(), srcpw_rec_rel[0].data() + wpsec * srcpw_rec_rel.nx, srcpw_rec_rel[1].data() + wpsec * srcpw_rec_rel.nx, colors_pw.data());
    if(w3_rec_rel.nx > 1)
			phase_view->add_graph("w3", src_rec_rel.nx, xs.data(), w3_rec[0].data() + wpsec * w3_rec_rel.nx, w3_rec[1].data() + wpsec * w3_rec_rel.nx, colors_w3.data());    
  }
  
  
  void subtract(float& ax, float& ay, float refx, float refy) {
		const float ax0 = ax;
		const float ay0 = ay;
		const float norm = std::sqrt(refx * refx + refy * refy);
		ax = (ax0 * refx + ay0 * refy) / norm;
		ay = (ay0 * refx - ax0 * refy) / norm;
		//a -> (a * ref*) / |ref|^2
	}
  
  void subtract(src_record& rec1, const src_record& ref, float c=1, float s=0) {
		if(rec1.nx != ref.nx || rec1.ny != ref.ny) {
			std::cout << "records have non matching geometries\n";
			return;
		}
		int nx = rec1.nx;
		int ny = rec1.ny;
		for(int j=0;j<ny;j++) {
			for(int i=0;i<nx;i++) {
				int addr = i+j*nx;
				subtract(rec1[0][addr], rec1[1][addr], ref[0][addr], ref[1][addr]);
				subtract(rec1[0][addr], rec1[1][addr], c, s);
			}
		}
	}
  
  void switch_record(int record_id)
  {
    if(record_id != last_record_id)
    {
      last_record_id = record_id;
      src_rec = data.get_src(last_record_id, rtz_harm_reader::src_rec_type::WWW);
      src2w_rec = data.get_src(last_record_id, rtz_harm_reader::src_rec_type::W2W2_W);
      srcpw_rec = data.get_src(last_record_id, rtz_harm_reader::src_rec_type::PW);
      w3_rec = data.get_src(last_record_id, rtz_harm_reader::src_rec_type::W3);
			
			src_rec_rel = src_rec;
			src2w_rec_rel = src2w_rec;
			srcpw_rec_rel = srcpw_rec;
			w3_rec_rel = w3_rec;
			if(w3_rec.nx == src_rec.nx && w3_rec.ny == src_rec.ny)
				subtract(src_rec_rel, w3_rec,0,1);
			if(w3_rec.nx == src2w_rec.nx && w3_rec.ny == src2w_rec.ny)
				subtract(src2w_rec_rel, w3_rec,0,1);
			if(w3_rec.nx == srcpw_rec.nx && w3_rec.ny == srcpw_rec.ny)
				subtract(srcpw_rec_rel, w3_rec,0,1);
			subtract(w3_rec_rel, w3_rec);
			
			
      pic_src_www->add_pic(src_rec.nx, src_rec.ny, src_rec[0].data(), src_rec[1].data());
      pic_src_w2w2_w->add_pic(src2w_rec.nx, src2w_rec.ny, src2w_rec[0].data(), src2w_rec[1].data());
      pic_src_pw->add_pic(srcpw_rec.nx, srcpw_rec.ny, srcpw_rec[0].data(), srcpw_rec[1].data());

      pic_src_www_rel->add_pic(src_rec_rel.nx, src_rec_rel.ny, src_rec_rel[0].data(), src_rec_rel[1].data());
      pic_src_w2w2_w_rel->add_pic(src2w_rec_rel.nx, src2w_rec_rel.ny, src2w_rec_rel[0].data(), src2w_rec_rel[1].data());
      pic_src_pw_rel->add_pic(srcpw_rec_rel.nx, srcpw_rec_rel.ny, srcpw_rec_rel[0].data(), srcpw_rec_rel[1].data());

      pic_w3->add_pic(w3_rec.nx, w3_rec.ny, w3_rec[0].data(), w3_rec[1].data());

      upate_section();
    }
  }
  
  widget* get_by_name(const std::string& name) {
    for(auto w : widgets)
    {
      if(w->name == name)
        return w;
    }
    return nullptr;
  }
  
  void toggle_integ() {
    auto gw=dynamic_cast<graph_widget*>(get_by_name("timeline"));
    gw->graphs["www_integ"].visible = !gw->graphs["www_integ"].visible;
    gw->graphs["w2w2w_integ"].visible = !gw->graphs["w2w2w_integ"].visible;
    gw->graphs["pw_integ"].visible = !gw->graphs["pw_integ"].visible;
    gw->update_graphs();
  }

  void toggle_inst() {
    auto gw=dynamic_cast<graph_widget*>(get_by_name("timeline"));
    gw->graphs["www_inst"].visible = !gw->graphs["www_inst"].visible;
    gw->graphs["w2w2w_inst"].visible = !gw->graphs["w2w2w_inst"].visible;
    gw->graphs["pw_inst"].visible = !gw->graphs["pw_inst"].visible;
    gw->update_graphs();
  }

};

void display(dashboard* db)
{
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glClearColor(.7,.7,.7,1.);
  glClear(GL_COLOR_BUFFER_BIT);
  glClear(GL_DEPTH_BUFFER_BIT);
  for(auto& w : db->widgets)
    w->display();
}

std::array<float, 4> v2c(float v) {
  std::array<float, 4> res{0,0,0,1};
  res[v>.0? 0 : 2] = std::pow(std::abs(v),.25f);
  res[1] = v*v;
  res[1] *= res[1];
  return res;
}

std::array<float, 4> cv2c(float v_re, float v_im) {
  std::array<float, 4> res{0,0,0,1};
  res[v_re > .0 ? 0 : 2] = std::pow(std::abs(v_re),.5);
  res[1] = std::pow(std::abs(v_im), 1.);
  return res;
}

void kbd_callback(GLFWwindow* win, int key, int scancode, int action, int mods)
{
  auto db = (dashboard*) glfwGetWindowUserPointer(win);

	auto update_timeline = [&](int i) {
		for(auto w : db->widgets) {
			if(w->name == "timeline") {
				auto gw = dynamic_cast<graph_widget*>(w);
				if(!gw)
					break;
				gw->vline = db->getPercentage(i);
				gw->update_lines();
			}
		}
	};

  if(key == GLFW_KEY_1 && action == GLFW_PRESS)
    db->toggle_integ();
  if(key == GLFW_KEY_2 && action == GLFW_PRESS)
    db->toggle_inst();
  if(key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
    int i=0;
    for(i;i<db->data.iterations.size();i++) {
			if(db->data.iterations[i] == db->last_record_id)
				break;
		}
		i = i < db->data.iterations.size() - 1 ? i+1 : i;
    db->switch_record(db->data.iterations[i]);
		update_timeline(i);
  }
  if(key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
    int i=0;
    for(i;i<db->data.iterations.size();i++) {
			if(db->data.iterations[i] == db->last_record_id)
				break;
		}
		i = i > 0 ? i-1 : i;
    db->switch_record(db->data.iterations[i]);
		update_timeline(i);
  }
  
  //find 3d widget
  graph3D_widget* w3d{};
  for(auto w : db->widgets) {
		w3d = dynamic_cast<graph3D_widget*>(w);
		if(w3d)
			break;
	}
  w3d->alt_mode = false;
	if((key == GLFW_KEY_LEFT_CONTROL || key == GLFW_KEY_RIGHT_CONTROL) && action == GLFW_PRESS)
		w3d->alt_mode = true;
}

void mouse_button_clbk(GLFWwindow* win, int button, int PressRelease, int mod) {
  auto db = (dashboard*) glfwGetWindowUserPointer(win);
  double xi, yi;
	glfwGetCursorPos(win, &xi, &yi);
  double x = xi / (double) winNx * 2. - 1.;
  double y = 1. - yi / (double) winNy * 2.;
  for(auto w : db->widgets) {
    if(w->in(x, y))
      w->mouse_button_clbk(win, button, PressRelease, mod, x, y);
  }
}

void mouse_move_clbk(GLFWwindow* win, double xi, double yi) {
  auto db = (dashboard*) glfwGetWindowUserPointer(win);
  static widget* last_widget_in = nullptr;
  double x = xi / (double) winNx * 2. - 1.;
  double y = 1. - yi / (double) winNy * 2.;
  for(auto w : db->widgets) {
    if(w->in(x, y)) {
      w->mouse_move_clbk(win, x, y);
      if(last_widget_in != w && last_widget_in)
        last_widget_in->cursor_got_away();
      last_widget_in = w;
    }
  }
}

void mouse_wheel_clbk(GLFWwindow* win, double xoffset, double yoffset) {
  auto db = (dashboard*) glfwGetWindowUserPointer(win);
  double xi, yi;
	glfwGetCursorPos(win, &xi, &yi);
  double x = xi / (double) winNx * 2. - 1.;
  double y = 1. - yi / (double) winNy * 2.;
  for(auto w : db->widgets) {
    if(w->in(x, y))
      w->mouse_wheel_clbk(win, xoffset, yoffset, x, y);
  }
}

static void win_sz_clbk(GLFWwindow* win, int x, int y)
{
		winNx = x;
		winNy = y;
		glfwMakeContextCurrent(win);
		glViewport(0,0,x,y);
}

void timeline_clbk(widget* w, void* owners_data) {
  auto timeline = dynamic_cast<graph_widget*>(w);
  auto db = (dashboard*) owners_data;
  int record_id = db->getIdbyPercentage(timeline->vline);//timeline->vline * db->data.records_count;
  if(timeline && timeline->LMB == false)
    std::cout << "timeline was modified! new vline=" << timeline->vline << ", new record id is " << record_id  << '\n';
  db->switch_record(db->data.iterations[record_id]);
}

void sec_clbk(widget* w, void* owners_data) {
  auto db = (dashboard*) owners_data;

  auto pic = dynamic_cast<picture_widget*>(w);
  
  if(pic) {
    db->pic_src_www->hline = pic->hline;
    db->pic_src_www_rel->hline = pic->hline;
    db->pic_src_w2w2_w->hline = pic->hline;
    db->pic_src_w2w2_w_rel->hline = pic->hline;
    db->pic_src_pw->hline = pic->hline;
    db->pic_src_pw_rel->hline = pic->hline;
    
    db->pic_src_www->vline = pic->vline;
    db->pic_src_www_rel->vline = pic->vline;
    db->pic_src_w2w2_w->vline = pic->vline;
    db->pic_src_w2w2_w_rel->vline = pic->vline;
    db->pic_src_pw->vline = pic->vline;
    db->pic_src_pw_rel->vline = pic->vline;
    
    db->pic_src_www->update_lines();
    db->pic_src_www_rel->update_lines();
    db->pic_src_w2w2_w->update_lines();
    db->pic_src_w2w2_w_rel->update_lines();
    db->pic_src_pw->update_lines();
    db->pic_src_pw_rel->update_lines();
  }
  
  db->upate_section();

}

int main(int argc, char* argv[]) {
  
  dashboard db;
  if(argc==1)
    db.data.path = "/home/dan/rtz_weq_harm/harm_ampW0=0.0007_ampW2=0.0002_phase2w=1.5_airDensity=0.0137_tunExp=3_";
  else {
    db.data.path = argv[1]; 
    if(argc==3)
      db.data.refpath = argv[2];
  }
  db.data.init();
  db.src_rec = db.data.get_src(7, rtz_harm_reader::src_rec_type::WWW);
  db.src2w_rec = db.data.get_src(7, rtz_harm_reader::src_rec_type::W2W2_W);
  db.srcpw_rec = db.data.get_src(7, rtz_harm_reader::src_rec_type::PW);
  
  if (!glfwInit()) {
    std::cerr << "ERROR: failed init glfw\n";
    return -1;
  }
  
  GLFWwindow* window = glfwCreateWindow(winNx, winNy, "3D geometry viewer", NULL, NULL);
  if (window == NULL) {
    std::cerr << "ERROR: failed to create glfw window\n";
    return -1;
  }
  
  glfwMakeContextCurrent(window);
  auto glew_err = glewInit();
  if (glew_err != GLEW_OK) {
    std::cerr << "renderer::init: ERROR: failed init glew with error " << glewGetErrorString(glew_err) << '\n';
    return -1;
  }
  
  glfwSetWindowUserPointer(window, (void*) &db);
  glfwSetMouseButtonCallback(window, mouse_button_clbk);
	glfwSetCursorPosCallback(window, mouse_move_clbk);
	glfwSetScrollCallback(window, mouse_wheel_clbk);
  glfwSetWindowSizeCallback(window, win_sz_clbk);
  glfwSetKeyCallback(window, kbd_callback);
  
  /*
  shader = loadShadersFromFile("sha.vs", "sha.fs");
  */

  int nx=113;
  int ny=117;
  
  std::vector<float> data1(nx*ny), data2(nx*ny), data3(nx*ny);
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) {
      float y = (float) j / (float) ny - .5f;
      float x = (float) i / (float) nx - .5f;
      data1[i+j*nx] = std::exp(-(x*x+y*y)*33.f);
      data2[i+j*nx] = std::exp(-(x*x+y*y)*133.f);
      data3[i+j*nx] = std::exp(-(x*x+y*y)*13.f);
    }
  }

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_w3";
    pw->x0 = -1.;
    pw->y0 = -1.;
    pw->x1 = -.6;
    pw->y1 = -0.33;
    pw->vline = .22;
    pw->vlineVisible = false;
    pw->hlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    db.widgets.push_back(pw);
    db.pic_w3 = pw;
    pw->callback = nullptr;
    pw->callback_owner_data = (void*) &db;
  }

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_srcwww";
    pw->x0 = .6;
    pw->y0 = -1.;
    pw->x1 = 0.97;
    pw->y1 = -0.33;
    pw->vline = .22;
    pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    pw->add_pic(db.src_rec.nx, db.src_rec.ny, db.src_rec[2].data());
    db.widgets.push_back(pw);
    db.pic_src_www = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }  

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_srcwww_rel";
    pw->x0 = .2;
    pw->y0 = -1.;
    pw->x1 = .57;
    pw->y1 = -0.33;
    pw->vline = .22;
    pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    db.widgets.push_back(pw);
    db.pic_src_www_rel = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }  

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_srcw2w2_w";
    pw->x0 = .6;
    pw->y0 = -.3;
    pw->x1 = .97;
    pw->y1 = .3;
    pw->vline = .22;
    pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    pw->add_pic(db.src2w_rec.nx, db.src2w_rec.ny, db.src2w_rec[2].data());
    db.widgets.push_back(pw);
    db.pic_src_w2w2_w = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }  

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_srcw2w2_w_rel";
    pw->x0 = .2;
    pw->y0 = -.3;
    pw->x1 = .57;
    pw->y1 = .3;
    pw->vline = .22;
    pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    db.widgets.push_back(pw);
    db.pic_src_w2w2_w_rel = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }  

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_src_pw";
    pw->x0 = .6;
    pw->y0 = .33;
    pw->x1 = .97;
    pw->y1 = 1.;
    pw->vline = .22;
		pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    pw->add_pic(db.srcpw_rec.nx, db.srcpw_rec.ny, db.srcpw_rec[2].data());
    db.widgets.push_back(pw);
    db.pic_src_pw = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }

  {
    picture_widget* pw = new picture_widget;
    pw->name = "pic_src_pw_rel";
    pw->x0 = .2;
    pw->y0 = .33;
    pw->x1 = .57;
    pw->y1 = 1.;
    pw->vline = .22;
		pw->vlineVisible = false;
    pw->hline = 0.03;
    pw->tm.value_to_color = &v2c;
    pw->tm.value_to_color_cplx = &cv2c;
    db.widgets.push_back(pw);
    db.pic_src_pw_rel = pw;
    pw->callback = &sec_clbk;
    pw->callback_owner_data = (void*) &db;
  }
    
  //planar plot:
  {
    graph_widget* gw = new graph_widget;
    gw->name = "timeline";
    gw->x0 = -1.;
    gw->y0 = .3;
    gw->x1 = .2;
    gw->y1 = 1.;
		gw->hlineVisible = false;
    gw->vline = .22;
    std::vector<float> xxx(db.data.records_count);
    for(int i=0;i<xxx.size();i++)
      xxx[i] = db.getPercentage(i);
    gw->add_graph("ew3", db.data.records_count, xxx.data(), db.data.e_rec.ew3.data(), RGB{.3, .22, .11});
    gw->add_graph("www_integ", db.data.records_count, xxx.data(), db.data.www_src_integral.data(), RGB{.8, .11, .11});
    gw->add_graph("w2w2w_integ", db.data.records_count, xxx.data(), db.data.w2w2w_src_integral.data(), RGB{.11, .11, .7});
    gw->add_graph("pw_integ", db.data.records_count, xxx.data(), db.data.pw_src_integral.data(), RGB{.6, .11, .75});
    gw->add_graph("www_inst", db.data.records_count, xxx.data(), db.data.www_src_instant.data(), RGB{.8, .11, .11});
    gw->add_graph("w2w2w_inst", db.data.records_count, xxx.data(), db.data.w2w2w_src_instant.data(), RGB{.11, .11, .7});
    gw->add_graph("pw_inst", db.data.records_count, xxx.data(), db.data.pw_src_instant.data(), RGB{.6, .11, .75});

    gw->callback = &timeline_clbk;
    gw->callback_owner_data = (void*) &db;
    db.widgets.push_back(gw);
    db.timeline_view = gw;
  }  
  
  //3d plot
  {
    graph3D_widget* gw = new graph3D_widget;
    gw->name = "3d";
    gw->x0 = -1.;
    gw->y0 = -1.;
    gw->x1 = .2;
    gw->y1 = .33;
    //gw->add_graph("abc", nx, xs.data(), ys.data(), ysc.data(), colors.data());
    //gw->add_graph("abc1", nx, xs.data(), ysc.data(), ys.data(), colors1.data());
    db.widgets.push_back(gw);
    db.phase_view = gw;
  }
  
  for(auto w : db.widgets) {
    w->init();
    w->update();
  }
  
  while (!glfwWindowShouldClose(window)) {
    glfwWaitEvents();
    //glfwPollEvents();
    //std::cout << "shady loop, cur=(" << xcursor << ", " << ycursor << ") view mode is " << viewMode << '\n';
    //std::cout << "posCurrent=" << posCurrent << " posSteady=" << posSteady << " fpCurrent=" << fpCurrent << " fpSteady=" << fpSteady << '\n';
    display(&db);
    glFlush();
    glfwSwapBuffers(window);
    
    //glfwPollEvents();
  }
  
  
  return 0;
}
