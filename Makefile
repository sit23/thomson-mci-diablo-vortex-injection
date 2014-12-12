#
# This is the makefile for qgbt_diablo.
# To compile the code, just type make.  Such an approach
# makes recompilation of the code easy, recompiling only 
# as necessary to account for recent changes to the code.
#

#COMPILER = /Applications/gccfolder/usr/local/bin/gfortran
COMPILER = /usr/bin/gfortran
c COMPILER = gfortran
COMPOPTS = -O3
LINKOPTS = -L/usr/lib -lrfftw -lfftw -lm
c LINKOPTS = -L/usr/lib64 -lrfftw -lfftw -lfftw -lm
c LINKOPTS = -L/usr/lib64 -lrfftw -lfftw -lm

qgbt_diablo: qgbt_diablo.f qgbt_periodic.o qg_fft.o qgbt_diablo.h grid_def
	$(COMPILER) $(COMPOPTS) qgbt_diablo.f -o qgbt_test \
        qgbt_periodic.o qg_fft.o $(LINKOPTS)

qgbt_periodic.o: qgbt_periodic.f qgbt_diablo.h grid_def
	$(COMPILER) $(COMPOPTS) -c qgbt_periodic.f 

qg_fft.o: qg_fft.f qgbt_diablo.h grid_def
	$(COMPILER) $(COMPOPTS) -c qg_fft.f

clean:
	m -f *.o  *~ qgbt_diablo core
