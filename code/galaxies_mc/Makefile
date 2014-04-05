#CFLAGS = -Wall -Wno-deprecated -I${GSLINCROOT} -O -arch i386
CFLAGS = -Wall -Wno-deprecated -I${GSLINCROOT} -O
SRC = driver.cc atomic.cc dnumrecipes.cc dcomplex.cc dcosmology.cc
OBJ = dnumrecipes.o dcomplex.o dcosmology.o astrophysics.o twentyonecm.o reionization.o haloDensity.o spline.o radio.o observation.o fisherCMB.o fisher.o fisherGAL.o fisher21CM.o neutrinos.o foregrounds.o inflation.o lymanforest.o ionization.o spline2D.o galaxies.o fisherLF.o astronomy.o
#CC = g++
CC = g++-4.2

#GSLLIBDIR = /sw/lib/
#GSLINCROOT = /sw/include/
GSLLIBDIR = /opt/local/lib/
GSLINCROOT = /opt/local/include/
LIB = -L${GSLLIBDIR} -lgsl -lgslcblas -lm

#GSLLIBDIR = /newman/user1/jp/code/gsl-1.6/.libs/
#GSLCBLASLIBDIR = /newman/user1/jp/code/gsl-1.6/cblas/.libs/
#GSLINCROOT = /newman/user1/jp/code/gsl-1.6/

#LIB = -L${GSLLIBDIR} -L$(GSLCBLASLIBDIR) -lgsl -lgslcblas -lm
#MISC = -Wl,-rpath=$(GSLLIBDIR) -Wl,-rpath=$(GSLCBLASLIBDIR)



.SUFFIXES : .o .cc
.cc.o :
		$(CC) $(CFLAGS) -c $<

all: driver.x test.x full.x like_driver.x like_driver_mc.x like_driver_vegas.x pdf_driver.x zeta_driver.x postprocess.x contourprep.x lf_driver.x lfz_driver.x

driver.x : $(OBJ) driver.o
	$(CC) $(CFLAGS) -o $@ driver.o $(OBJ) $(LIB) $(MISC)

lf_driver.x : $(OBJ) lf_driver.o
	$(CC) $(CFLAGS) -o $@ lf_driver.o $(OBJ) $(LIB) $(MISC)

lfz_driver.x : $(OBJ) lfz_driver.o
	$(CC) $(CFLAGS) -o $@ lfz_driver.o $(OBJ) $(LIB) $(MISC)

zeta_driver.x : $(OBJ) zeta_driver.o
	$(CC) $(CFLAGS) -o $@ zeta_driver.o $(OBJ) $(LIB) $(MISC)

like_driver.x : $(OBJ) like_driver.o
	$(CC) $(CFLAGS) -o $@ like_driver.o $(OBJ) $(LIB) $(MISC)

pdf_driver.x : $(OBJ) pdf_driver.o
	$(CC) $(CFLAGS) -o $@ pdf_driver.o $(OBJ) $(LIB) $(MISC)

like_driver_mc.x : $(OBJ) like_driver_mc.o
	$(CC) $(CFLAGS) -o $@ like_driver_mc.o $(OBJ) $(LIB) $(MISC)

like_driver_vegas.x : $(OBJ) like_driver_vegas.o
	$(CC) $(CFLAGS) -o $@ like_driver_vegas.o $(OBJ) $(LIB) $(MISC)

test.x : $(OBJ) test.o
	$(CC) $(CFLAGS) -o $@ test.o $(OBJ) $(LIB) $(MISC)

full.x : $(OBJ) full.o
	$(CC) $(CFLAGS) -o $@ full.o $(OBJ) $(LIB) $(MISC)

postprocess.x : $(OBJ) postprocess.o
	$(CC) $(CFLAGS) -o $@ postprocess.o $(OBJ) $(LIB) $(MISC)

contourprep.x : $(OBJ) contourprep.o
	$(CC) $(CFLAGS) -o $@ contourprep.o $(OBJ) $(LIB) $(MISC)

clean : 
	rm -f *.o

new :
	rm -f *.o
	rm -f *spline*.dat

nion :
	./like_driver_mc.x -n0 -c0 -p4 -t5 -g0<standard.ini &
	./like_driver_mc.x -n1 -c0 -p4 -t6 -g0<standard.ini &
	./like_driver_mc.x -n0 -c1 -p4 -t7 -g0<standard.ini &
	./like_driver_mc.x -n0 -c2 -p4 -t8 -g0<standard.ini &

ftwostep :
	./like_driver_mc.x -n0 -c0 -p2 -t1 -g0<standard.ini &
	./like_driver_mc.x -n1 -c0 -p2 -t2 -g0<standard.ini &
	./like_driver_mc.x -n0 -c1 -p2 -t3 -g0<standard.ini &
	./like_driver_mc.x -n0 -c2 -p2 -t4 -g0<standard.ini &

fpolyz :
	./like_driver_mc.x -n0 -c0 -p5 -t5 -g0<standard.ini &
	./like_driver_mc.x -n1 -c0 -p5 -t6 -g0<standard.ini &
	./like_driver_mc.x -n0 -c1 -p5 -t7 -g0<standard.ini &
	./like_driver_mc.x -n0 -c2 -p5 -t8 -g0<standard.ini &

full :
	#make nion
	make twostep

full_d:
	make nion_d
	make twostep_d

nion_d :
	./like_driver_mc.x -n0 -c0 -p4 -d0 -t11<standard.ini &
	./like_driver_mc.x -n1 -c0 -p4 -d0 -t12<standard.ini &
	./like_driver_mc.x -n0 -c1 -p4 -d0 -t13<standard.ini &
	./like_driver_mc.x -n0 -c2 -p4 -d0 -t14<standard.ini &

twostep_d :
	./like_driver_mc.x -n0 -c0 -p2 -d0 -t15<standard.ini &
	./like_driver_mc.x -n1 -c0 -p2 -d0 -t16<standard.ini &
	./like_driver_mc.x -n0 -c1 -p2 -d0 -t17<standard.ini &
	./like_driver_mc.x -n0 -c2 -p2 -d0 -t18<standard.ini &


polyz_d :
	./like_driver_mc.x -n0 -c0 -p5 -d0 -t15<standard.ini
	./like_driver_mc.x -n1 -c0 -p5 -d0 -t16<standard.ini
	./like_driver_mc.x -n0 -c1 -p5 -d0 -t17<standard.ini
	./like_driver_mc.x -n0 -c2 -p5 -d0 -t18<standard.ini

nion_nolls :
	#./like_driver_mc.x -n0 -c0 -p3 -l0 -t1 <standard.ini
	#./like_driver_mc.x -n1 -c0 -p3 -l0 -t2 <standard.ini
	#./like_driver_mc.x -n0 -c1 -p3 -l0 -t3 <standard.ini
	#./like_driver_mc.x -n0 -c2 -p3 -l0 -t4 <standard.ini

	./like_driver_mc.x -n0 -c0 -p2 -l0 -t5 <standard.ini &
	./like_driver_mc.x -n1 -c0 -p2 -l0 -t6 <standard.ini &
	./like_driver_mc.x -n0 -c1 -p2 -l0 -t7 <standard.ini &
	./like_driver_mc.x -n0 -c2 -p2 -l0 -t8 <standard.ini &

pdffull:
	make pdfA
	make pdfB

pdfA :
	./pdf_driver.x -z4 -a1 <standard.ini &
	./pdf_driver.x -z5 -a1 <standard.ini &
	./pdf_driver.x -z6 -a1 <standard.ini &

pdfB :
	./pdf_driver.x -z4 -a2<standard.ini &
	./pdf_driver.x -z5 -a2<standard.ini &
	./pdf_driver.x -z6 -a2<standard.ini &

pdfC :
	./pdf_driver.x -z4 -a3<standard.ini &
	./pdf_driver.x -z5 -a3<standard.ini &
	./pdf_driver.x -z6 -a3<standard.ini &

pdf :
	./pdf_driver.x -z4 -e4 <standard.ini &
	./pdf_driver.x -z5 -e4 <standard.ini &
	./pdf_driver.x -z6 -e4 <standard.ini &

tocm :
	./like_driver_mc.x -n0 -c0 -p4 -m1 -t21<standard.ini 
	#./like_driver_mc.x -n0 -c0 -p2 -m1 -t22<standard.ini &

tocm_d :
	./like_driver_mc.x -n0 -c0 -p4 -m1 -t21 -d0<standard.ini 
	#./like_driver_mc.x -n0 -c0 -p2 -m1 -t22 -d0<standard.ini 

tocm_p :
	./postprocess.x -n0 -c0 -p4 -m1 -t31 -d0<standard.ini &
	./postprocess.x -n0 -c0 -p2 -m1 -t32 -d0<standard.ini &

nion_p :
	./postprocess.x -n0 -c0 -p4 -d0 -t11<standard.ini &
	./postprocess.x -n1 -c0 -p4 -d0 -t12<standard.ini &
	./postprocess.x -n0 -c1 -p4 -d0 -t13<standard.ini &
	./postprocess.x -n0 -c2 -p4 -d0 -t14<standard.ini &

twostep_p :
	./postprocess.x -n0 -c0 -p2 -d0 -t15<standard.ini
	./postprocess.x -n1 -c0 -p2 -d0 -t16<standard.ini
	./postprocess.x -n0 -c1 -p2 -d0 -t17<standard.ini
	./postprocess.x -n0 -c2 -p2 -d0 -t18<standard.ini

post_p :
	./postprocess.x -n0 -c0 -p4 -d0 -t11<standard.ini &
	./postprocess.x -n0 -c0 -p2 -d0 -t15<standard.ini &

contourprep :
	./contourprep.x -s0
	./contourprep.x -s1


	./contourprep.x -s0 -p2
	./contourprep.x -s1 -p2

	./contourprep.x -s2
	./contourprep.x -s2 -p2

contourprep2 :
	./contourprep.x -s0 -c2
	./contourprep.x -s1 -c2


	./contourprep.x -s0 -p2 -c2
	./contourprep.x -s1 -p2 -c2

	./contourprep.x -s2 -c2
	./contourprep.x -s2 -p2 -c2

contourprep_tocm :
	./contourprep.x -s0 -m1
	./contourprep.x -s1 -m1
	./contourprep.x -s2 -m1

	./contourprep.x -s0 -p2 -m1
	./contourprep.x -s1 -p2 -m1
	./contourprep.x -s2 -p2 -m1

nolya :
	./like_driver_mc.x -n2 -c0 -p4 -t41<standard.ini & 
	./like_driver_mc.x -n2 -c0 -p2 -t42<standard.ini &

nolya_d :
	./like_driver_mc.x -n2 -c0 -p4 -t43 -d0<standard.ini 
	./like_driver_mc.x -n2 -c0 -p2 -t44 -d0<standard.ini 

nolya_p :
	./postprocess.x -n2 -c0 -p4 -t51 -d0<standard.ini 
	./postprocess.x -n2 -c0 -p2 -t52 -d0<standard.ini 

contourprep_nolya :
	./contourprep.x -s0 -n2
	./contourprep.x -s1 -n2
	./contourprep.x -s2 -n2

	./contourprep.x -s0 -p2 -n2
	./contourprep.x -s1 -p2 -n2
	./contourprep.x -s2 -p2 -n2

betterg :
	./like_driver_mc.x -n3 -c0 -p4 -t41<standard.ini & 
	./like_driver_mc.x -n3 -c0 -p2 -t42<standard.ini &

betterg_d :
	./like_driver_mc.x -n3 -c0 -p4 -t43 -d0<standard.ini 
	./like_driver_mc.x -n3 -c0 -p2 -t44 -d0<standard.ini 

betterg_p :
	./postprocess.x -n3 -c0 -p4 -t51 -d0<standard.ini 
	./postprocess.x -n3 -c0 -p2 -t52 -d0<standard.ini 

contourprep_betterg :
	./contourprep.x -s0 -n3
	./contourprep.x -s1 -n3
	./contourprep.x -s2 -n3

	./contourprep.x -s0 -p2 -n3
	./contourprep.x -s1 -p2 -n3
	./contourprep.x -s2 -p2 -n3

nocut :
	./like_driver_mc.x -n0 -k1 -c0 -p4 -t41<standard.ini & 
	./like_driver_mc.x -n0 -k1 -c0 -p2 -t42<standard.ini &

nocut_d :
	./like_driver_mc.x -n0 -k1 -c0 -p4 -t43 -d0<standard.ini &
	./like_driver_mc.x -n0 -k1 -c0 -p2 -t44 -d0<standard.ini &

nocut_p :
	./postprocess.x -n0 -k1 -c0 -p4 -t51 -d0<standard.ini 
	./postprocess.x -n0 -k1 -c0 -p2 -t52 -d0<standard.ini 

contourprep_nocut :
	./contourprep.x -s0 -n0 -k1
	./contourprep.x -s1 -n0 -k1
	./contourprep.x -s2 -n0 -k1

	./contourprep.x -s0 -p2 -n0 -k1
	./contourprep.x -s1 -p2 -n0 -k1
	./contourprep.x -s2 -p2 -n0 -k1

nocmb :
	./like_driver_mc.x -n0 -q1 -c0 -p4 -t41<standard.ini & 
	./like_driver_mc.x -n0 -q1 -c0 -p2 -t42<standard.ini &

nocmb_d :
	./like_driver_mc.x -n0 -q1 -c0 -p4 -t43 -d0<standard.ini &
	./like_driver_mc.x -n0 -q1 -c0 -p2 -t44 -d0<standard.ini &

nocmb_p :
	./postprocess.x -n0 -q1 -c0 -p4 -t51 -d0<standard.ini 
	./postprocess.x -n0 -q1 -c0 -p2 -t52 -d0<standard.ini 

pdfgal :
	./pdf_driver.x -z3.8 -e0<standard.ini 
	./pdf_driver.x -z5 -e0<standard.ini
	./pdf_driver.x -z5.9 -e0<standard.ini 
	./pdf_driver.x -z7.4 -e0<standard.ini 
	./pdf_driver.x -z10.2  -e0<standard.ini 

pdfgalA :
	./pdf_driver.x -z3.8 -a1 <standard.ini 
	./pdf_driver.x -z5 -a1  <standard.ini
	./pdf_driver.x -z5.9 -a1  <standard.ini 
	./pdf_driver.x -z7.4 -a1  <standard.ini 
	./pdf_driver.x -z10.2 -a1  <standard.ini 

pdfgalB :
	./pdf_driver.x -z3.8 -a2 <standard.ini 
	./pdf_driver.x -z5 -a2  <standard.ini
	./pdf_driver.x -z5.9 -a2  <standard.ini 
	./pdf_driver.x -z7.4 -a2  <standard.ini 
	./pdf_driver.x -z10.2 -a2  <standard.ini 

pdfgalC :
	./pdf_driver.x -z3.8 -a3 <standard.ini 
	./pdf_driver.x -z5 -a3  <standard.ini
	./pdf_driver.x -z5.9 -a3  <standard.ini 
	./pdf_driver.x -z7.4 -a3  <standard.ini 
	./pdf_driver.x -z10.2 -a3  <standard.ini 

pdfc :
	./pdf_driver.x -z10.2 -a1  <standard.ini 
	./pdf_driver.x -z10.2 -a2  <standard.ini 
	./pdf_driver.x -z10.2 -a3  <standard.ini 
