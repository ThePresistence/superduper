noinst_LIBRARIES = libguga.a

libguga_a_SOURCES =      \
        basissets.f90    \
	davidson.f90     \
        gaussoei.f90     \
        gtointerface.f90 \
	gugaci.f90	 \
	gugadirect.f90	 \
	gugaglobal.f90	 \
	gugaincore.f90	 \
	gugaloops.f90	 \
	gugalower.f90	 \
	gugaplots.f90	 \
	gugasetup.f90	 \
	gugashape.f90	 \
	gugashape1.f90	 \
	gugashape2.f90	 \
	gugashape3.f90	 \
	gugashape4.f90	 \
	gugashape5.f90	 \
	gugashape6.f90	 \
	gugashape7.f90	 \
	gugashape8.f90	 \
	gugashape9.f90	 \
	gugashape10.f90  \
	gugashape11.f90  \
	gugasrtgen.f90	 \
	gugatime.c	 \
	gugaupper.f90	 \
	gugautils.f90	 \
	matutils.f90	 \
	property.f90	 \
        sphrharm.f90

if COND_CSRMM
libguga_a_SOURCES += gugamkl.f90
GUGASPARSE_O       = gugamkl.o
else
libguga_a_SOURCES += gugasparse.f90
GUGASPARSE_O       = gugasparse.o
endif

if COND_CXX
libguga_a_SOURCES += cpudavliu.cxx
else
libguga_a_SOURCES += guganocxx.f90
endif

if COND_CUDA
libguga_a_SOURCES += gpudavliu.cu
else
libguga_a_SOURCES += guganocuda.f90
endif

.cu.o:
	$(CXX) $(NVFLAGS_GUGACI) -c $<

CONST.o: ../CONST.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c ../CONST.f90

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c ../LIMIT.f90

basissets.o: basissets.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c basissets.f90

cpudavliu.o: cpudavliu.cxx
	$(CXX) $(CXXFLAGS_GUGACI) -c cpudavliu.cxx

davidson.o: davidson.f90 gugaglobal.o gugaloops.o gugashape3.o gugashape4.o \
	 gugashape5.o gugashape6.o $(GUGASPARSE_O)
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c davidson.f90

gaussoei.o: gaussoei.f90 basissets.o sphrharm.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gaussoei.f90

gtointerface.o: gtointerface.f90 gaussoei.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gtointerface.f90

gugaci.o: gugaci.f90 gugadirect.o gugaglobal.o gugaincore.o gugasetup.o \
	gugashape.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaci.f90

gugadirect.o: gugadirect.f90 davidson.o gugaglobal.o gugaloops.o gugashape.o \
	gugashape7.o gugashape8.o gugashape9.o gugashape10.o gugashape11.o \
	gugautils.o matutils.o property.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugadirect.f90

gugaglobal.o: gugaglobal.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaglobal.f90

gugaincore.o: gugaincore.f90 davidson.o gugaglobal.o gugaloops.o gugashape.o \
	gugasrtgen.o gugautils.o matutils.o property.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaincore.f90

gugaloops.o: gugaloops.f90 gugaglobal.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaloops.f90

gugalower.o: gugalower.f90 gugaglobal.o gugaloops.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugalower.f90

gugamkl.o: gugamkl.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugamkl.f90

guganocuda.o: guganocuda.f90 gugaglobal.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c guganocuda.f90

guganocxx.o: guganocxx.f90 gugaglobal.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c guganocxx.f90

gugaplots.o: gugaplots.f90 gugaglobal.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaplots.f90

gugasetup.o: gugasetup.f90 gugaglobal.o gugaplots.o gugashape.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugasetup.f90

gugashape.o: gugashape.f90 gugaglobal.o gugalower.o gugashape1.o gugashape2.o \
	gugasrtgen.o gugaupper.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape.f90

gugashape1.o: gugashape1.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape1.f90

gugashape2.o: gugashape2.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape2.f90

gugashape3.o: gugashape3.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape3.f90

gugashape4.o: gugashape4.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape4.f90

gugashape5.o: gugashape5.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape5.f90

gugashape6.o: gugashape6.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape6.f90

gugashape7.o: gugashape7.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape7.f90

gugashape8.o: gugashape8.f90 gugaglobal.o gugaloops.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape8.f90

gugashape9.o: gugashape9.f90 gugaglobal.o gugaloops.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape9.f90

gugashape10.o: gugashape10.f90 gugaglobal.o gugaloops.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape10.f90

gugashape11.o: gugashape11.f90 gugaglobal.o gugaloops.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugashape11.f90

gugasrtgen.o: gugasrtgen.f90 gugaglobal.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugasrtgen.f90

gugasparse.o: gugasparse.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugasparse.f90

gugatime.o: gugatime.c
	$(CC) $(CFLAGS_GUGACI) -c gugatime.c

gugaupper.o: gugaupper.f90 gugaglobal.o gugaloops.o gugautils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugaupper.f90

gugautils.o: gugautils.f90 gugaglobal.o matutils.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c gugautils.f90

matutils.o: matutils.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c matutils.f90

property.o: property.f90 CONST.o matutils.o
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c property.f90

sphrharm.o: sphrharm.f90
	$(F90) $(FREE)$(FFLAGS2_GUGACI) -c sphrharm.f90
