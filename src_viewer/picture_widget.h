#pragma once

#include "widget.h"
#include "tex.h"

struct picture_widget : public widget {
  tex_manager tm;
  float vline;
  float hline;
  bool vlineVisible = true;
  bool hlineVisible = true;
  static unsigned int sha_pic;
  static unsigned int sha_lines;
  VAPlocation pic_verts_loc, pic_tex_coords_loc, lines_verts_loc;
  unsigned int tex_loc;
  unsigned int pic_vbo = 0;
  unsigned int lines_vbo = 0;
  void add_pic(int nx, int ny, float*);
  void add_pic(int nx, int ny, float*, float*);
  void update_pic();
  void update_lines();
  void display() override;
  void update() override;
  void init() override;
  void mouse_button_clbk(GLFWwindow* win, int button, int PressRelease, int mod, float xi, float yi) override;
  void mouse_move_clbk(GLFWwindow* win, float xi, float yi) override;
};
