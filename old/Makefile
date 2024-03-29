############### Makefile #######################################################

############### Source, object & binary files path #############################
srcdir := src
basedir := .
objdir := $(basedir)/obj
bindir := $(basedir)/bin
outdir := $(basedir)/out
magmadir := /usr/local/magma/lib
cudadir := /usr/local/cuda/lib64
#atlasdir := /usr/lib64/atlas
atlasdir := /usr/local/atlas/lib

############### Names  #########################################################
program := heat
compiler := gfortran

############### Libraries ######################################################
magmalib := -L$(magmadir) -lmagma
atlaslib := $(atlasdir)/libcblas.a $(atlasdir)/libatlas.a
cudalib :=  -L$(cudadir) -lcublas -lcudart 
#lib := -lblas -llapack $(magmalib) $(atlaslib) $(cudalib) -lstdc++ -lpthread  
lib := $(magmalib) $(cudalib) -llapack -lblas -lstdc++

############### Object files names #############################################
object := $(objdir)/$(program).o
source := $(srcdir)/main.f90

############### Flags ##########################################################
#flags := -Wall -cpp -O2 -pg -g -pthread -fopenmp
flags := -Wall -cpp -O3 -pthread -fopenmp

############### Targets ########################################################
world: $(program)

$(program): $(source) $(object)
	$(compiler) $(flags) -o $(bindir)/$(program) $(object) $(lib)

$(object): $(source)
	$(compiler) $(flags) -c $(source) -o $(object)

run: 
	$(bindir)/$(program)
	mv *.out $(outdir)

plot:
	./plot.py $(outdir)/*.out

video:
	cd $(outdir); ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p heat.mp4; cd -

clean:
	rm -f *.o *.mod *.2 $(program) 

