#pragma once

#include <array>
#include <string>

#include <GL/gl.h>

class GLFWwindow;

struct VAPlocation {
  VAPlocation() = default;
  VAPlocation(GLuint, GLint, GLenum, GLsizei, void*);
  GLuint loc;
 	GLint count;
 	GLenum type;
 	GLsizei stride;
 	void* offset;
 	void enable();
};


struct widget {
  float x0,y0,x1,y1;
  std::string name;
  static unsigned int vao;
  virtual void display() = 0;
  virtual void update() = 0;
  virtual void init() = 0;
  virtual void cursor_got_away();
  bool in(float x, float y) const;
  virtual void mouse_button_clbk(GLFWwindow* win, int button, int PressRelease, int mod, float xi, float yi);
  virtual void mouse_move_clbk(GLFWwindow* win, float xi, float yi) {}
  virtual void mouse_wheel_clbk(GLFWwindow* win, double xoffset, double yoffset, float xi, float yi) {}
  bool LMB{};
  bool RMB{};
  float x_click;
  float y_click;
  std::array<float, 2> get_local_xy(float x, float y) const;
  void (*callback)(widget* caller, void* owner_data) = nullptr;
  void* callback_owner_data{};
};

#define WIDGET_CALLBACK if(callback) (*callback)(this, callback_owner_data);