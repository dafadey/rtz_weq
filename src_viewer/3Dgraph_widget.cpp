#include <GL/glew.h>
#include "3Dgraph_widget.h"
#include <GLFW/glfw3.h>

#include <iostream>
#include <limits>
#include "shaders.h"

unsigned int graph3D_widget::sha3D = 0;

void graph3D_widget::init() {
  if(!sha3D) {
    std::string vs = 
                    "#version 150\n"
                    "in vec3 vertex_pos;\n"
                    "in vec3 normal_pos;\n"
                    "in vec3 color_pos;\n"
                    "uniform vec3 light_pos;\n"
                    "uniform mat4 view_matrix;\n"
                    "uniform mat4 proj_matrix;\n"
                    "out vec3 color;\n"
                    "void main() {\n"
                    " float f = .97 * abs(dot(light_pos, normal_pos));\n"
                    " float s = f*f;\n"
                    " s = s*s;\n"
                    " s = s*s;\n"
                    " s = s*s;\n"
                    " s = s*s;\n"
                    " color = color_pos * (.3+.7*f) * (1.-s) +  s*vec3(1,1,1);\n"
                    " gl_Position = proj_matrix * (view_matrix * vec4(vertex_pos,1));\n"
                    "}\n";
    std::string fs =
                    "#version 150\n"
                    "in vec3 color;"
                    "out vec4 outCol;\n"
                    "void main() {\n"
                    "  outCol = vec4(color,0.8);\n"
                    "}\n";
  
    sha3D = loadShaders(vs, fs);
  }
  
  if(!vao)
    glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  if(!vbo)
    glGenBuffers(1, &vbo);
  
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  verts_loc = VAPlocation(glGetAttribLocation(sha3D, "vertex_pos"), 3, GL_FLOAT, 9 * sizeof(float), (void*) 0);
  normals_loc = VAPlocation(glGetAttribLocation(sha3D, "normal_pos"), 3, GL_FLOAT, 9 * sizeof(float), (void*) (3*sizeof(float)));
  colors_loc = VAPlocation(glGetAttribLocation(sha3D, "color_pos"), 3, GL_FLOAT, 9 * sizeof(float), (void*) (6*sizeof(float)));
  proj_matrix_location = glGetUniformLocation(sha3D, "proj_matrix");
  view_matrix_location = glGetUniformLocation(sha3D, "view_matrix");
  light_pos_location = glGetUniformLocation(sha3D, "light_pos");
  
  std::cout << "graph_widget::verts_location=" << verts_loc.loc << " graph_widget::normals_location=" << normals_loc.loc << " graph_widget::colors_location=" << colors_loc.loc << '\n';

}

void graph3D_widget::add_graph(const std::string& name, int n, float* buffxyz, const RGB* col) {
  graph3D g(n);
  for(int i=0; i<n; i++) {
    g[i] = vec3f{buffxyz[i*3], buffxyz[i*3+1], buffxyz[i*3+2]};
    g.colors[i] = col[i];
  }
  add_graph(name, g);
}

void graph3D_widget::add_graph(const std::string& name, int n, float* buffx, float* buffy, float* buffz, const RGB* col) {
  graph3D g(n);
  for(int i=0; i<n; i++) {
    g[i]=vec3f{buffx[i], buffy[i], buffz[i]};
    g.colors[i]=col[i];
  }
  add_graph(name, g);
}

void b_box3D::update(const vec3f& v) {
  for(int i=0;i<3;i++) {
    vmin[i] = std::min(vmin[i], v[i]);
    vmax[i] = std::max(vmax[i], v[i]);
  }
}

vec3f b_box3D::get(int id) const {
  bool is2 = (id >> 2) == 1;
  bool is1 = ((id & 2) >> 1) == 1;
  bool is0 = (id & 1) == 1;
  return vec3f{is0 ? vmax[0] : vmin[0],
               is1 ? vmax[1] : vmin[1],
               is2 ? vmax[2] : vmin[2]};
}


void b_box3D::reset() {
  vmin = vec3f{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  vmax = vec3f{-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
}

float sf(float x) {
	return x * x * (x > .0f ? 1.f : -1.f);
}

void graph3D_widget::update_graphs() {
  bb.reset();
  nquads=0;
  for(auto& it : graphs) {
    nquads += it.second.size();
    nquads--;
  }

	if(nquads<=0)
		return;
  //n lines i.e. 2 n verts each vertex has 2 coords + 3 color components
  int divs=16;
  nquads *= divs;
  std::vector<float> buff((nquads+1)*4*9);
  
  int pos=0;
  int cell=0;
	auto addvec=[&](const vec3f& v) {
		for(int i=0;i<3;i++,pos++)
			buff[pos] = v[i];
	};
  for(auto& it : graphs) {
    auto& g = it.second;
    
    for(int i=0;i<g.size()-1;i++) {
      vec3f p0_{g[i][0], tran_factor*sf(g[i][1]), tran_factor*sf(g[i][2])};
      vec3f p1_{g[i+1][0], tran_factor*sf(g[i+1][1]), tran_factor*sf(g[i+1][2])};
      vec3f a0_{g[i][0], .0, .0};
      vec3f a1_{g[i+1][0], .0, .0};
      vec3f c0{g.colors[i]};
      vec3f c1{g.colors[i+1]};
      cell++;
      for(int j=0;j<divs;j++)
      {
        cell++;
        const float shadow_fac = .7f;
        if((cell % 2)==0) {
          c0 = RGB{g.colors[i][0]*shadow_fac,g.colors[i][1]*shadow_fac,g.colors[i][2]*shadow_fac};
          c1 = RGB{g.colors[i+1][0]*shadow_fac,g.colors[i+1][1]*shadow_fac,g.colors[i+1][2]*shadow_fac};
        } else {
          c0 = g.colors[i];
          c1 = g.colors[i+1];
        }

        float a = float(j) / float(divs);
        float d = 1. / float(divs);
        vec3f p0 = p0_*(a+d)+a0_*(1.-a-d);
        vec3f a0 = p0_*a+a0_*(1-a);
        vec3f p1 = p1_*(a+d)+a1_*(1.-a-d);
        vec3f a1 = p1_*a+a1_*(1-a);
        
        vec3f n0=cross_prod(a0-p0, p1-p0);
        normalize(n0);
        vec3f n1=cross_prod(p0-p1, a1-p1);
        normalize(n1);
        vec3f na1=cross_prod(a0-a1, p1-a1);
        normalize(na1);
        vec3f na0=cross_prod(p1-a0, a1-a0);
        normalize(na0);

        addvec(a0);
        addvec(na0);
        addvec(c0);
        addvec(p0);
        addvec(n0);
        addvec(c0);
        addvec(p1);
        addvec(n1);
        addvec(c1);
        addvec(a1);
        addvec(na1);
        addvec(c1);
        bb.update(a0);
        bb.update(p0);
        bb.update(p1);
        bb.update(a1);
      }
    }
  }
  
  double x0 = .5*std::sqrt(std::pow(bb.vmax[1] - bb.vmin[1], 2.) + std::pow(bb.vmax[2] - bb.vmin[2], 2.));
  vec3f vx00 = ((bb.vmax[0]*.13 + bb.vmin[0]*.87)+.01*(bb.vmax[0] - bb.vmin[0]))*vec3f{1,0,0};
  vec3f vx01 = ((bb.vmax[0]*.13 + bb.vmin[0]*.87)+.01*(bb.vmax[0] - bb.vmin[0]))*vec3f{1,0,0} + x0 * vec3f{0,1,0};
  vec3f vx11 = ((bb.vmax[0]*.13 + bb.vmin[0]*.87)-.01*(bb.vmax[0] - bb.vmin[0]))*vec3f{1,0,0} + x0 * vec3f{0,1,0};
  vec3f vx10 = ((bb.vmax[0]*.13 + bb.vmin[0]*.87)-.01*(bb.vmax[0] - bb.vmin[0]))*vec3f{1,0,0};
  addvec(vx00);
  addvec(vec3f{0.f,0.f,1.f});
  addvec(vec3f{1.f,0.f,0.f});
  addvec(vx01);
  addvec(vec3f{0.f,0.f,1.f});
  addvec(vec3f{1.f,0.f,0.f});
  addvec(vx11);
  addvec(vec3f{0.f,0.f,1.f});
  addvec(vec3f{1.f,0.f,0.f});
  addvec(vx10);
  addvec(vec3f{0.f,0.f,1.f});
  addvec(vec3f{1.f,0.f,0.f});
  
  //update vbo
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, buff.size() * sizeof(float), (void*) buff.data(), GL_DYNAMIC_DRAW);
}

void graph3D_widget::add_graph(const std::string& name, graph3D& g) {
  graphs[name] = g;
  update_graphs();
}

void graph3D_widget::del_graph(const std::string& name) {
  graphs.erase(name);
  update_graphs();
}

void graph3D_widget::update() {
  update_graphs();
}

vec4f graph3D_widget::project(const vec3f& pt) const {
  vec4f in{pt[0], pt[1], pt[2], 1.f};
  vec4f v{.0f, .0f, .0f, .0f};
  for(int j = 0; j < 4; j++) {
    v[j] = .0f;
    for(int i = 0; i < 4; i++)
      v[j] += in[i] * view_matrix[i * 4 + j];
  }

  vec4f o=v;
  for(int j = 0; j < 4; j++) {
    o[j] = .0f;
    for(int i = 0; i < 4; i++) {
      o[j] += v[i] * proj_matrix[i * 4 + j];
    }
  }
  o[0] /= o[3];
  o[1] /= o[3];
  o[2] /= o[3];
  o[3] /= o[3];
  return o;
}

void graph3D_widget::update_view() {
  vec3f cam_y = cam_up;
  normalize(cam_y);
  vec3f cam_z = cam_pos - focal_point;
  normalize(cam_z);
  vec3f cam_x = cross_prod(cam_y, cam_z);

  vec3f shift = vec3f{cam_pos * cam_x, cam_pos * cam_y, cam_pos * cam_z};
  
  for(int i=0;i<16;i++)
    view_matrix[i]=0;
  view_matrix[15] = 1.;
  
  for(int i=0;i<3;i++) {
    view_matrix[i*4] = cam_x[i];
    view_matrix[i*4+1] = cam_y[i];
    view_matrix[i*4+2] = cam_z[i];
  }
  for(int i=0;i<3;i++)
    view_matrix[3*4+i] = - shift[i];

  float f = std::numeric_limits<float>::max();
  float n = -std::numeric_limits<float>::max();
  
  for(int i=0; i<8; i++) {
    float proj = static_cast<GLfloat>((bb.get(i) - cam_pos) * cam_z);
    if(!std::isnan(proj) && !std::isinf(proj)) {
      f = f < proj ? f : proj;
      n = n > proj ? n : proj;
    }
  }
  bool parallel = true;
  float scale;
  if(parallel)
    scale = 4.f * std::tan(56.f / 2.f * 3.14f / 180.f) / std::sqrt((cam_pos - focal_point)*(cam_pos - focal_point));
  else
    scale = 1.f/std::tan(56.f / 2.f * 3.14f / 180.f); 

  //std::cout << "near=" << n << " far=" << f << " scale=" << scale << '\n';

  for(int i=0;i<16;i++)
    proj_matrix[i]=0;
  
  int dims[4] = {0};
  glGetIntegerv(GL_VIEWPORT, dims);
  float fbw = dims[2];
  float fbh = dims[3];
  
  if(parallel) {
    proj_matrix[0] = scale;
    proj_matrix[5] = scale / fbh * fbw;
    proj_matrix[10] = 2./(f-n);
    proj_matrix[14] = -(f+n)/(f-n);
    proj_matrix[11] = .0f;
    proj_matrix[15] = 1.f;
  } else { // perspective
    proj_matrix[0] = scale;
    proj_matrix[5] = scale / fbh * fbw;
    proj_matrix[10] = -f/(f-n);
    proj_matrix[11] = -1.f;
    proj_matrix[14] = f*n/(f-n);
    proj_matrix[15] = 0.f;
  }
  proj_matrix[12] = .5*(x0+x1);
  proj_matrix[13] = .5*(y0+y1);
  
}

void graph3D_widget::display() {
  glEnable(GL_DEPTH_TEST);

  glBindVertexArray(vao);
  
  glUseProgram(sha3D);

  update_view();
  
  /*
  std::cout << "bb:=\n";
  for(int i=0;i<8;i++)
    std::cout << bb.get(i) << '\n';
  
  std::cout << "nquads=" << nquads << '\n';
  std::cout << "vmin : " << bb.vmin << " vmax : " << bb.vmax << '\n';
  */
  
  vec3f light = cam_pos + cam_up*3.;
  normalize(light);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  verts_loc.enable();
  normals_loc.enable();
  colors_loc.enable();
  glUniformMatrix4fv(view_matrix_location, 1, false, view_matrix);
  glUniformMatrix4fv(proj_matrix_location, 1, false, proj_matrix);
  glUniform3fv(light_pos_location, 1, light.data());
  
  glDrawArrays(GL_QUADS, 0, 4 * (nquads+1));

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glUseProgram(0);
  
}

void graph3D_widget::mouse_button_clbk(GLFWwindow* win, int button, int action, int mod, float xi, float yi) {
  widget::mouse_button_clbk(win, button, action, mod, xi, yi);

  WIDGET_CALLBACK
}

void graph3D_widget::mouse_move_clbk(GLFWwindow* win, float xi, float yi) {
  if(LMB) {
    float phi = static_cast<float>(-(xi-x_click)/.3);
    float theta = static_cast<float>(-(yi-y_click)/.3);
    //std::cout << "phi=" << phi << " theta=" << theta << '\n';
    
    vec3f cam_y = cam_up;
    normalize(cam_y);
    vec3f cam_z = cam_pos - focal_point;
    normalize(cam_z);
    vec3f cam_x = cross_prod(cam_y, cam_z);
    
    
    vec3f cam_up_rel{-std::sin(phi) * std::sin(theta), std::cos(theta), -std::cos(phi) * std::sin(theta)};
    cam_up = cam_up_rel[0] * cam_x + cam_up_rel[1] * cam_y + cam_up_rel[2] * cam_z;

    normalize(cam_up);

    vec3f cam_pos_rel = cam_pos - focal_point;
    cam_pos_rel = vec3f{cam_pos_rel * cam_x, cam_pos_rel * cam_y, cam_pos_rel * cam_z};
    cam_pos_rel = vec3f{cam_pos_rel[2] * std::sin(phi) * std::cos(theta), cam_pos_rel[2] * std::sin(theta), cam_pos_rel[2] * std::cos(phi) * std::cos(theta)};
    
    cam_pos = cam_pos_rel[0] * cam_x + cam_pos_rel[1] * cam_y + cam_pos_rel[2] * cam_z;
    cam_pos = cross_prod(cam_up, cross_prod(cam_pos, cam_up));
    cam_pos = cam_pos + focal_point;
    //std::cout << "cam_pos : " << cam_pos << '\n';
    //std::cout << "cam_up : " << cam_up << '\n';
  }
  x_click = xi;
  y_click = yi;

  WIDGET_CALLBACK
}

void graph3D_widget::mouse_wheel_clbk(GLFWwindow* win, double xoffset, double yoffset, float xi, float yi) {
  if(alt_mode) {
		tran_factor *= std::exp(.1*yoffset);
		update_graphs();
  } else
		cam_pos = (cam_pos - focal_point) * (1-.1*yoffset) + focal_point;

  WIDGET_CALLBACK
}
