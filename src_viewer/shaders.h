#pragma once
#include <string>

#include <GL/gl.h>

GLuint loadShaders(const std::string& VertexShaderCode, const std::string& FragmentShaderCode, const std::string& GeometryShaderCode = "");

GLuint loadShadersFromFile(const char* vsfile, const char* fsfile, const char* gsfile = nullptr);
