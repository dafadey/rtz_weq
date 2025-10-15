#include "tex.h"

#include <GL/glew.h>
#include <GL/gl.h>

void tex_manager::init(int nx_, int ny_) {
  if(nx==nx_ && ny==ny_)
    return;
  nx=nx_;
  ny=ny_;
  if(tex_buff)
    delete[] tex_buff;
  tex_buff = new float[nx*ny*4];
  if(tex)
    glDeleteTextures(1, &tex);
  glGenTextures(1, &tex);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, nx, ny, 0, GL_RGBA, GL_FLOAT, tex_buff);
  glGenerateMipmap(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, 0);
}

void tex_manager::update_on_host(int Dx, int Dy, float* data) {
  init(Dx,Dy);
    
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++) {
      std::array<float, 4> col;
      
      if(value_to_color)
        col = value_to_color(data[j*nx+i]);
      else {
        for(int c=0;c<4;c++)
          col[i] = data[(j*nx+i)*4+c];
      }
      
      for(int c=0;c<4;c++)
        tex_buff[j*nx*4+i*4+c] = col[c];
    }
  }
}

//cplx version
void tex_manager::update_on_host(int Dx, int Dy, float* data_re, float* data_im) {
  init(Dx,Dy);
    
  for(int j=0; j<ny; j++) {
    for(int i=0; i<nx; i++) {
      std::array<float, 4> col;
      
      if(value_to_color_cplx)
        col = value_to_color_cplx(data_re[j*nx+i], data_im[j*nx+i]);
      else {
        for(int c=0;c<4;c++)
          col[i] = data_re[(j*nx+i)*4+c];
      }
      
      for(int c=0;c<4;c++)
        tex_buff[j*nx*4+i*4+c] = col[c];
    }
  }
}

void tex_manager::update_on_device() {
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, nx, ny, 0, GL_RGBA, GL_FLOAT, tex_buff);
  glGenerateMipmap(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, 0);
}
