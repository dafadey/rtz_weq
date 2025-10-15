#include <GL/glew.h>
#include "widget.h"
#include <GLFW/glfw3.h>
#include <iostream>

unsigned int widget::vao = 0;

VAPlocation::VAPlocation(GLuint loc_, GLint cnt, GLenum t, GLsizei s, void* ptr) : loc(loc_), count(cnt), type(t), stride(s), offset(ptr) {}

void VAPlocation::enable() {
  glVertexAttribPointer(loc, count, type, GL_FALSE, stride, offset);
  glEnableVertexAttribArray(loc);
}

bool widget::in(float x, float y) const {
  return (x-x0) * (x1-x) > .0 && (y-y0) * (y1-y) > .0;
}

std::array<float, 2> widget::get_local_xy(float x, float y) const {
  return std::array<float, 2> {(x-x0)/(x1-x0), (y-y0)/(y1-y0)};
}

void widget::mouse_button_clbk(GLFWwindow* win, int button, int action, int mod, float xi, float yi) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if(action == GLFW_PRESS)
      LMB = true;
    if(action == GLFW_RELEASE)
      LMB = false;
  }
  if (button == GLFW_MOUSE_BUTTON_RIGHT) {
    if(action == GLFW_PRESS)
      RMB = true;
    if(action == GLFW_RELEASE)
      RMB = false;
  }
  if(action == GLFW_PRESS) {
    x_click = xi;
    y_click = yi;
  }
}

void widget::cursor_got_away() {
  LMB = false;
  RMB = false;
}