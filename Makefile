SRCPATH := src

arch := CC30
mrch := sm_30
vmrch := compute_30
CUDA := /usr/local/cuda-9.0
CUTIL:= /home/dan/cuda3cutil_fPIE
#CUPTI:= /usr/local/cuda-9.0/extras/CUPTI
GLFW3:=  /home/dan/distr/glfw3/build/
TOOLS:= /home/dan/tools/
MGL := /home/dan/distr/mathgl-2.5/install
INC := -I$(GLFW3)/include -I$(CUTIL)/common/inc -I$(CUDA)/include -I$(CUDA)/samples/common/inc -I$(MGL)/include -I$(TOOLS)/include
LIB := $(CUDA)/lib64
CUTIL_LIBS := $(CUTIL)/lib/libcutil.a

fm := 

harm : thz_kernels.obj harm.obj queue.obj postprocess.obj
	g++ -fPIE -L$(LIB) -L$(GLFW3)/src -L$(TOOLS)/lib $^ $(CUTIL_LIBS) -L$(MGL)/lib -lgetTemps -lgetGPUtemp -fopenmp -lsimpledraw_glfw3 -lGL -lglfw -lpthread -lcudart -lcufft -lmgl -Wl,-rpath=$(LIB) -Wl,-rpath=$(TOOLS)/lib -Wl,-rpath=$(MGL)/lib -o $@

thz_kernels.obj : $(SRCPATH)/thz_kernels.cu $(SRCPATH)/dfs.h $(SRCPATH)/nl_kernel.cu $(SRCPATH)/sweep_kernel.cu
	$(CUDA)/bin/nvcc -g -Wno-deprecated-gpu-targets -w $(INC) -Xcompiler -fPIC -O3 -c --gpu-architecture=$(vmrch) --gpu-code=$(mrch) $(fm) -lineinfo $< -o $@

thz_kernels.ptx : $(SRCPATH)/thz_kernels.cu $(SRCPATH)/dfs.h $(SRCPATH)/nl_kernel.cu $(SRCPATH)/sweep_kernel.cu
	$(CUDA)/bin/nvcc -ptx -g -Wno-deprecated-gpu-targets -w $(INC) -Xcompiler -fPIC -O3 -c --gpu-architecture=$(vmrch) --gpu-code=$(mrch) $(fm) -lineinfo $(SRCPATH)/thz_kernels.cu

font.o : $(SRCPATH)/font.dat
	objcopy -I binary -O elf64-x86-64 --binary-architecture i386 $^ $@

harm.obj : $(SRCPATH)/harm.cpp $(SRCPATH)/dfs.h $(SRCPATH)/fft_real_3.cpp $(SRCPATH)/thz_funcs.cpp $(SRCPATH)/queue.h $(SRCPATH)/postprocess.h $(SRCPATH)/simpledraw2D.h
	g++ --std=c++17 $(INC) -g -Wno-deprecated-declarations $< -c -O3 -o $@

queue.obj : $(SRCPATH)/queue.cpp $(SRCPATH)/queue.h
	g++ $(INC) -g -Wno-deprecated-declarations -fopenmp $< -c -O3 -o $@

postprocess.obj : $(SRCPATH)/postprocess.cpp $(SRCPATH)/postprocess.h $(SRCPATH)/dfs.h
	g++ $(INC) -g -Wno-deprecated-declarations -fopenmp $< -c -O3 -o $@

simpledraw2d_glfw3.obj : $(SRCPATH)/simpledraw2d_glfw3.cpp $(SRCPATH)/simpledraw2D.h
	g++ $(INC) -g -Wno-deprecated-declarations $< -c -O3 -o $@
	
.PHONY: clean
	rm -f *.obj *.ptx
