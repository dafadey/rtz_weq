#pragma once

#include "widget.h"
#include <map>
#include <vector>
#include <string>

#include "vectors.h"

typedef std::array<float,3> RGB;

struct graph : public std::vector<vec2f> {
  graph() = default;
  graph(int n) : std::vector<vec2f>(n) {}
  RGB color;
  bool visible = true;
};

struct graph_widget : public widget {
  float vline;
  float hline;
  bool vlineVisible = true;
  bool hlineVisible = true;  
  static unsigned int sha_lines;
  VAPlocation verts_loc, colors_loc;
  unsigned int lines_vbo = 0;
  unsigned int vline_vbo = 0;
  int nlines;
  void update_line();
  
  std::map<std::string, graph> graphs;
  
  void add_graph(const std::string& name, int n, float*, const RGB&);
  void add_graph(const std::string& name, int n, float*, float*, const RGB&);
  void add_graph(const std::string& name, graph&);
  void del_graph(const std::string& name);
  void update_graphs();
  void update_lines();
  
  void display() override;
  void update() override;
  void init() override;
  void mouse_button_clbk(GLFWwindow* win, int button, int PressRelease, int mod, float xi, float yi) override;
  void mouse_move_clbk(GLFWwindow* win, float xi, float yi) override;
};
