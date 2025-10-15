#include <GL/glew.h>
#include "picture_widget.h"
#include <GLFW/glfw3.h>

#include <iostream>
#include "shaders.h"

unsigned int picture_widget::sha_pic = 0;

unsigned int picture_widget::sha_lines = 0;

void picture_widget::init() {
  if(!sha_pic) {
    std::string vs = 
                    "#version 150\n"
                    "in vec2 vertex_pos;\n"
                    "in vec2 tex_coord;\n"
                    "out vec2 tex_coord_fs;\n"
                    "void main() {\n"
                    "  gl_Position = vec4(vertex_pos,0,1.);\n"
                    "  tex_coord_fs = tex_coord;\n"
                    "}\n";
    std::string fs =
                    "#version 150\n"
                    "in vec2 tex_coord_fs;\n"
                    "out vec4 outCol;\n"
                    "uniform sampler2D tex;\n"
                    "void main() {\n"
                    "  outCol = texture(tex, tex_coord_fs);\n"
                    "}\n";
    sha_pic = loadShaders(vs, fs);
  }
  if(!sha_lines) {
    std::string vs = 
                    "#version 150\n"
                    "in vec2 vertex_pos;\n"
                    "void main() {\n"
                    " gl_Position = vec4(vertex_pos,0,1);\n"
                    "}\n";
    std::string fs =
                    "#version 150\n"
                    "out vec4 outCol;\n"
                    "void main() {\n"
                    "  outCol = vec4(0,1,0,.5);\n"
                    "}\n";
  
    sha_lines = loadShaders(vs, fs);
  }
  
  if(!vao)
    glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  if(!pic_vbo)
    glGenBuffers(1, &pic_vbo);
  if(!lines_vbo)
    glGenBuffers(1, &lines_vbo);

  glBindBuffer(GL_ARRAY_BUFFER, pic_vbo);
  pic_verts_loc = VAPlocation(glGetAttribLocation(sha_pic, "vertex_pos"), 2, GL_FLOAT, 4 * sizeof(float), (void*) 0);
  pic_tex_coords_loc = VAPlocation(glGetAttribLocation(sha_pic, "tex_coord"), 2, GL_FLOAT, 4 * sizeof(float), (void*) (2 * sizeof(float)));
  tex_loc = glGetUniformLocation(sha_pic, "tex");
  
  std::cout << "picture_widget::pic_verts_location=" << pic_verts_loc.loc << " picture_widget::pic_tex_coords_location=" << pic_tex_coords_loc.loc << " picture_widget::tex_loc=" << tex_loc << '\n';
  
  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  lines_verts_loc = VAPlocation(glGetAttribLocation(sha_lines, "vertex_pos"), 2, GL_FLOAT, 2 * sizeof(float), (void*) 0);
  std::cout << "picture_widget::lines_verts_location=" << lines_verts_loc.loc << '\n';
}


void picture_widget::update_pic() {
  glBindVertexArray(vao);
  
  glBindBuffer(GL_ARRAY_BUFFER, pic_vbo);
  float verts[] = {
                   x0, y0, .0, .0,
                   x1, y0, 1., .0,
                   x1, y1, 1., 1.,
                   x0, y1, .0, 1.,
                  };  

  glBufferData(GL_ARRAY_BUFFER, 16 * sizeof(float), (void*) verts, GL_DYNAMIC_DRAW);
}

void picture_widget::update_lines() {
  glBindVertexArray(vao);
   
  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  float lines[] = {x0+vline*(x1-x0), y0, x0+vline*(x1-x0), y1,
                   x0, y0+hline*(y1-y0), x1, y0+hline*(y1-y0)};
  glBufferData(GL_ARRAY_BUFFER, 8 * sizeof(float), (void*) lines, GL_DYNAMIC_DRAW);
  
}

void picture_widget::update() {
  update_pic();
  update_lines();
}

void picture_widget::add_pic(int nx, int ny, float* data) {
  tm.update_on_host(nx, ny, data);
  tm.update_on_device();
}

//cplx version
void picture_widget::add_pic(int nx, int ny, float* data_re, float* data_im) {
  tm.update_on_host(nx, ny, data_re, data_im);
  tm.update_on_device();
}

void picture_widget::display() {
  glDisable(GL_DEPTH_TEST);
  glUseProgram(sha_pic);
  glBindVertexArray(vao);
  
  glBindTexture(GL_TEXTURE_2D, tm.tex);
  glActiveTexture(GL_TEXTURE0);
  glUniform1i(tex_loc, 0);

  glBindBuffer(GL_ARRAY_BUFFER, pic_vbo);
  pic_verts_loc.enable();
  pic_tex_coords_loc.enable();

	glDrawArrays(GL_QUADS, 0, 4);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glUseProgram(sha_lines);

  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  lines_verts_loc.enable();
  
	if(vlineVisible)
		glDrawArrays(GL_LINES, 0, 2);
	if(hlineVisible)
		glDrawArrays(GL_LINES, 2, 2);
		
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glBindVertexArray(0);
  glUseProgram(0);
  
}

void picture_widget::mouse_button_clbk(GLFWwindow* win, int button, int action, int mod, float xi, float yi) {
  //std::cout << "***xi=" << xi << " yi=" << yi << '\n';
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE && LMB)
  {
    vline = get_local_xy(xi, yi)[0];
    hline = get_local_xy(xi, yi)[1];
    update_lines();
  }
  //std::cout << "vline=" << vline << " LMB=" << LMB << " RMB=" << RMB << '\n';
  widget::mouse_button_clbk(win, button, action, mod, xi, yi);
  
  WIDGET_CALLBACK
}

void picture_widget::mouse_move_clbk(GLFWwindow* win, float xi, float yi) {
  if(LMB) {
    vline = get_local_xy(xi, yi)[0];
    hline = get_local_xy(xi, yi)[1];
    update_lines();
  }

  WIDGET_CALLBACK
}

