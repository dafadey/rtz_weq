#include <GL/glew.h>
#include "graph_widget.h"
#include <GLFW/glfw3.h>

#include <iostream>
#include <limits>
#include "shaders.h"

unsigned int graph_widget::sha_lines = 0;

void graph_widget::init() {
  if(!sha_lines) {
    std::string vs = 
                    "#version 150\n"
                    "in vec2 vertex_pos;\n"
                    "in vec3 color_pos;\n"
                    "out vec3 color;\n"
                    "void main() {\n"
                    " color = color_pos;"
                    " gl_Position = vec4(vertex_pos,0,1);\n"
                    "}\n";
    std::string fs =
                    "#version 150\n"
                    "in vec3 color;"
                    "out vec4 outCol;\n"
                    "void main() {\n"
                    "  outCol = vec4(color,.5);\n"
                    "}\n";
  
    sha_lines = loadShaders(vs, fs);
  }
  
  if(!vao)
    glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  if(!lines_vbo)
    glGenBuffers(1, &lines_vbo);
  
  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  verts_loc = VAPlocation(glGetAttribLocation(sha_lines, "vertex_pos"), 2, GL_FLOAT, 5 * sizeof(float), (void*) 0);
  colors_loc = VAPlocation(glGetAttribLocation(sha_lines, "color_pos"), 3, GL_FLOAT, 5 * sizeof(float), (void*) (2*sizeof(float)));
  std::cout << "graph_widget::verts_location=" << verts_loc.loc << " graph_widget::colors_location=" << colors_loc.loc << '\n';

  if(!vline_vbo)
    glGenBuffers(1, &vline_vbo);
}


void graph_widget::add_graph(const std::string& name, int n, float* buffxy, const RGB& col) {
  graph g(n);
  g.color = col;
  for(int i=0; i<n; i++) {
    g[i][0]=buffxy[i*2];
    g[i][1]=buffxy[i*2+1];
  }
  add_graph(name, g);
}

void graph_widget::add_graph(const std::string& name, int n, float* buffx, float* buffy, const RGB& col) {
  graph g(n);
  g.color = col;
  for(int i=0; i<n; i++) {
    g[i][0]=buffx[i];
    g[i][1]=buffy[i];
  }
  add_graph(name, g);
}

struct b_box {
  float xmin{std::numeric_limits<float>::max()};
  float xmax{-std::numeric_limits<float>::max()};
  float ymin{std::numeric_limits<float>::max()};
  float ymax{-std::numeric_limits<float>::max()};
  void update(float x, float y) {
    xmin = std::min(x, xmin);
    xmax = std::max(x, xmax);
    ymin = std::min(y, ymin);
    ymax = std::max(y, ymax);
  }
};

void graph_widget::update_lines() {
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, vline_vbo);
  float lines[] = {x0+vline*(x1-x0), y0, 1., .0, .0, x0+vline*(x1-x0), y1, 1., 0., .0,
                   x0, y0+hline*(y1-y0), 1., .0, .0, x1, y0+hline*(y1-y0), 1., .0, .0};
  glBufferData(GL_ARRAY_BUFFER, 20 * sizeof(float), (void*) lines, GL_DYNAMIC_DRAW);
}

void graph_widget::update_graphs() {
  b_box bb;
  nlines=0;
  for(auto& it : graphs) {
    for(vec2f& pt : it.second) {
      bb.update(pt[0], pt[1]);
      nlines++;
    }
    nlines--;
  }
  //n lines i.e. 2 n verts each vertex has 2 coords + 3 color components
  if(nlines <= 0)
		return;
  std::vector<float> buff(nlines*2*5);
  
  int pos=0;
  for(auto& it : graphs) {
    auto& g = it.second;
    if(!g.visible)
      continue;
    for(int i=0;i<g.size()-1;i++) {
      for(int j=0;j<2;j++) {
        buff[pos+0] = (g[i+j][0] - bb.xmin) / (bb.xmax - bb.xmin) * (x1-x0)+x0;
        buff[pos+1] = (g[i+j][1] - bb.ymin) / (bb.ymax - bb.ymin) * (y1-y0)+y0;
        buff[pos+2] = g.color[0];
        buff[pos+3] = g.color[1];
        buff[pos+4] = g.color[2];
        pos+=5;
      }
    }
  }
  //update vbo
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  glBufferData(GL_ARRAY_BUFFER, buff.size() * sizeof(float), (void*) buff.data(), GL_DYNAMIC_DRAW);
}

void graph_widget::add_graph(const std::string& name, graph& g) {
  graphs[name] = g;
  update_graphs();
}

void graph_widget::del_graph(const std::string& name) {
  graphs.erase(name);
  update_graphs();
}

void graph_widget::update() {
  update_graphs();
  update_lines();
}

void graph_widget::display() {
  glDisable(GL_DEPTH_TEST);

  glBindVertexArray(vao);
  
  glUseProgram(sha_lines);

  glBindBuffer(GL_ARRAY_BUFFER, lines_vbo);
  verts_loc.enable();
  colors_loc.enable();
  
  glDrawArrays(GL_LINES, 0, 2 * nlines);

  glBindBuffer(GL_ARRAY_BUFFER, vline_vbo);
  verts_loc.enable();
  colors_loc.enable();
	if(vlineVisible)
		glDrawArrays(GL_LINES, 0, 2);
	if(hlineVisible)
		glDrawArrays(GL_LINES, 2, 2);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glUseProgram(0);
  
}

void graph_widget::mouse_button_clbk(GLFWwindow* win, int button, int action, int mod, float xi, float yi) {
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


void graph_widget::mouse_move_clbk(GLFWwindow* win, float xi, float yi) {
  if(LMB) {
    vline = get_local_xy(xi, yi)[0];
    hline = get_local_xy(xi, yi)[1];
    update_lines();
  }

  WIDGET_CALLBACK
}
