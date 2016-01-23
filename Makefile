############### Makefile #######################################################

############### Source, object & binary files path #############################
basedir := /tmp/heat
srcdir := ./src
objdir := $(basedir)/obj
bindir := $(basedir)/bin
outdir := $(basedir)/out

#cudadir := /usr/local/cuda/lib64

############### Names  #########################################################
program := heat
compiler := g++

############### Libraries ######################################################
#cudalib :=  -L$(cudadir) -lcublas -lcudart
lapacklib := -llapacke -lopenblas

#-lstdc++ -lpthread
#lib :=  $(cudalib) -llapack -lblas -lstdc++
lib := $(lapacklib)

############### Object files names #############################################
object := $(objdir)/$(program).o
source := $(srcdir)/main.cpp

############### Flags ##########################################################
debugflags := -Wall -g
profflags := -pg
optflags := -O3 -fopenmp
flags := -O2
############### Targets ########################################################
world: flags := $(optflags)
world: $(program)

debug: flags += $(debugflags)
debug: $(program)

prof: flags += $(profflags)
prof: $(program)

$(program): $(source) $(object)
	mkdir -p $(bindir)
	$(compiler) $(flags) -o $(bindir)/$(program) $(object) $(lib)

$(object): $(source)
	mkdir -p $(objdir)
	$(compiler) $(flags) -c $(source) -o $(object)

run:
	$(bindir)/$(program)

plot:
	mkdir -p $(outdir)
	./plot.py $(outdir)/*.out

video: $(outdir)
	cd $(outdir); ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p heat.mp4; cd -

clean:
	rm -f *.o *.mod *.2 $(program)
