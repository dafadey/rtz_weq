all: viewer

f := -fopenmp -g -O3 --std=c++17 -I/home/dan/glew-2.1.0/include/ -I/home/dan/glu-master/include -I/home/dan/mesa-24.2.3/include -I/home/dan/distr/glfw3/glfw-3.2.1/include

glew := /home/dan/glew-2.1.0/lib
glfw := /home/dan/distr/glfw3/build/src
gl := /home/dan/mesa-24.2.3/install/lib/x86_64-linux-gnu/

libs := -L${glew} -L${glfw} -L${gl} -Wl,-rpath=${glew} -Wl,-rpath=${glfw} -Wl,-rpath=${gl}

viewer: rtz_harm_reader.obj viewer.obj tex.obj vectors.obj widget.obj picture_widget.obj graph_widget.obj 3Dgraph_widget.obj shaders.obj
	g++ $f $^ ${libs} -lGL -lGLEW -lglfw -lpng -o $@

3Dgraph_widget.obj: 3Dgraph_widget.cpp 3Dgraph_widget.h vectors.h widget.h
	g++ $f $< -c -o $@

graph_widget.obj: graph_widget.cpp graph_widget.h vectors.h widget.h
	g++ $f $< -c -o $@

picture_widget.obj: picture_widget.cpp picture_widget.h tex.h
	g++ $f $< -c -o $@

widget.obj: widget.cpp widget.h
	g++ $f $< -c -o $@

vectors.obj: vectors.cpp vectors.h
	g++ $f $< -c -o $@

tex.obj: tex.cpp tex.h
	g++ $f $< -c -o $@

shaders.obj: shaders.cpp shaders.h
	g++ $f $< -c -o $@

viewer.obj: viewer.cpp tex.h widget.h picture_widget.h graph_widget.h 3Dgraph_widget.h vectors.h rtz_harm_reader.h
	g++ $f $< -c -o $@
	
rtz_harm_reader.obj: rtz_harm_reader.cpp rtz_harm_reader.h
	g++ $f $< -c -o $@
