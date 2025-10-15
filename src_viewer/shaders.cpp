#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <GL/gl.h>

#include "shaders.h"

GLuint loadShadersFromFile(const char* vsFileName, const char* fsFileName, const char* gsFileName) {
  auto getContents = [](const char* fileName) -> std::string {
    std::string out;
    if(!fileName)
      return out;
    std::string s;
    std::ifstream file(fileName);
    while(getline(file,s))
      out += s + '\n';
    file.close();
    return out;
  };
  auto vsContents = getContents(vsFileName);
  auto fsContents = getContents(fsFileName);
  auto gsContents = getContents(gsFileName);
  return loadShaders(vsContents, fsContents, gsContents);
}

GLuint loadShaders(const std::string& VertexShaderCode, const std::string& FragmentShaderCode, const std::string& GeometryShaderCode) {

    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint GeometryShaderID = 0;
    if(!GeometryShaderCode.empty())
      GeometryShaderID = glCreateShader(GL_GEOMETRY_SHADER);

    GLint Result = GL_FALSE;
    int InfoLogLength;

    #ifndef QUIET
    std::cout << "building vertex shader...\n";
    #endif
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);

    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> VertexShaderErrorMessage(InfoLogLength+1);
      glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
      std::cerr << "building vertex shader failed with " << VertexShaderErrorMessage.data() << '\n';
    }

    if(GeometryShaderID) {
      std::cout << "building geometry shader...\n";
      char const* GeometrySourcePointer = GeometryShaderCode.c_str();
      glShaderSource(GeometryShaderID, 1, &GeometrySourcePointer, NULL);
      glCompileShader(GeometryShaderID);

      glGetShaderiv(GeometryShaderID, GL_COMPILE_STATUS, &Result);
      glGetShaderiv(GeometryShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
      if (InfoLogLength > 0) {
        std::vector<char> GeometryShaderErrorMessage(InfoLogLength + 1);
        glGetShaderInfoLog(GeometryShaderID, InfoLogLength, NULL, &GeometryShaderErrorMessage[0]);
        std::cerr << "building geometry shader failed with " << GeometryShaderErrorMessage.data() << '\n';
      }
    }

    #ifndef QUIET
    std::cout << "building fragment shader...\n";
    #endif
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);

    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> FragmentShaderErrorMessage(InfoLogLength+1);
      glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
      std::cerr << "building fragment shader failed with " << FragmentShaderErrorMessage.data() << '\n';
    }

    #ifndef QUIET
    std::cout << "building shader program...\n";
    #endif
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    if (GeometryShaderID)
      glAttachShader(ProgramID, GeometryShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);

    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    if ( InfoLogLength > 0 ){
      std::vector<char> ProgramErrorMessage(InfoLogLength+1);
      glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
      std::cerr << "building shader program failed with " << ProgramErrorMessage.data() << '\n';
    }

    glDeleteShader(VertexShaderID);
    glDeleteShader(GeometryShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}
