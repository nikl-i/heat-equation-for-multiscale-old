############### Makefile #######################################################

############### Source, object & binary files path #############################
basedir := /tmp/heat
srcdir := /home/cf/Develop/heat/src
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
objects := $(objdir)/main.o $(objdir)/solution.o $(objdir)/problem.o $(objdir)/solver.o
sources := $(srcdir)/main.cpp $(srcdir)/problem.cpp $(srcdir)/solution.cpp $(srcdir)/solver.cpp

############### Flags ##########################################################
debugflags := -Wall -g -D DEBUG
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

$(program): $(sources) $(objects)
	mkdir -p $(bindir)
	$(compiler) $(flags) -o $(bindir)/$(program) $(objects) $(lib)

$(objdir)/main.o: $(srcdir)/main.cpp
	mkdir -p $(objdir)
	$(compiler) $(flags) -c $(srcdir)/main.cpp -o $(objdir)/main.o

$(objdir)/solution.o: $(srcdir)/solution.cpp
	$(compiler) $(flags) -c $(srcdir)/solution.cpp -o $(objdir)/solution.o

$(objdir)/problem.o: $(srcdir)/problem.cpp
	$(compiler) $(flags) -c $(srcdir)/problem.cpp -o $(objdir)/problem.o

$(objdir)/solver.o: $(srcdir)/solver.cpp
	$(compiler) $(flags) -c $(srcdir)/solver.cpp -o $(objdir)/solver.o

run:
	$(bindir)/$(program)

plot:
	./plot.py
video: $(outdir)
	cd $(outdir); ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p heat.mp4; cd -

clean:
	rm -f *.o *.mod *.2 $(program)
