

FC = gfortran
FFLAGS = -O3 #-fcheck=bounds # -Wunused 

OBJECTS = mcmc.o \
	rdparam.o \
	BDsignal.o \
	random_flat.o \
	r8_normal_ab.o \
	r8_uniform_01.o\
	brack.o \
	ftst.o \
	ftstb.o \
	init_random_seed.o \
	gmh_on.o \
	randnum.o \
	control.o \
	rdinput.o \
	wrout.o

mcmcexe: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o 0348d4

rdparam.o: rdparam.f const_incl
BDsignal.o: BDsignal.f const_incl
r8_normal_ab.o: r8_normal_ab.f
r8_uniform_01.o: r8_uniform_01.f
brack.o: brack.f
ftst.o: ftst.f
ftstb.o: ftstb.f
init_random_seed.o: init_random_seed.f
gmh_on.o: gmh_on.f const_incl
randnum.o: randnum.f
control.o: control.f const_incl
rdinput.o: rdinput.f const_incl
wrout.o: wrout.f const_incl
mcmc.o: mcmc.f const_incl
