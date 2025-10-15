#pragma once

#include "widget.h"
#include <map>
#include <vector>
#include <string>

#include "vectors.h"

typedef std::array<float,3> RGB;

struct b_box3D {
  vec3f vmin{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
  vec3f vmax{-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
  
  void update(const vec3f& v);
  
  void reset();
  
  vec3f get(int id) const;
};


struct graph3D : public std::vector<vec3f> {
  graph3D() = default;
  graph3D(int n) : std::vector<vec3f>(n) {colors.resize(n);}
  std::vector<RGB> colors;
};

struct graph3D_widget : public widget {
  static unsigned int sha3D;
  VAPlocation verts_loc, normals_loc, colors_loc;
  float view_matrix[16];
  float proj_matrix[16];
  unsigned int proj_matrix_location, view_matrix_location, light_pos_location;
  unsigned int vbo = 0;
  int nquads;
  b_box3D bb;
  double tran_factor = 1.;
  bool alt_mode = false;

  vec3f focal_point{0,0,0};
  vec3f cam_pos{3., 3.3, .0};
  vec3f cam_up{0,0,1};
  std::map<std::string, graph3D> graphs;
  
  void add_graph(const std::string& name, int n, float*, const RGB*);
  void add_graph(const std::string& name, int n, float*, float*, float*, const RGB*);
  void add_graph(const std::string& name, graph3D&);
  void del_graph(const std::string& name);
  void update_graphs();
  
  vec4f project(const vec3f& pt) const;
  
  void update_view();
  void display() override;
  void update() override;
  void init() override;
  void mouse_button_clbk(GLFWwindow* win, int button, int PressRelease, int mod, float xi, float yi) override;
  void mouse_move_clbk(GLFWwindow* win, float xi, float yi) override;
  void mouse_wheel_clbk(GLFWwindow* win, double xoffset, double yoffset, float xi, float yi) override;
};
