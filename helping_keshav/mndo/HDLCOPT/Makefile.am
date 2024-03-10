noinst_LIBRARIES = libhdlc.a

libhdlc_a_SOURCES =				\
constraint.f90	global.f90	hdlclib.f90	\
matrixlib.f90	mndolib.f90	primitive.f90	\
search.f90	searchlib.f90

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c ../LIMIT.f90

global.o: global.f90
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c global.f90

matrixlib.o: matrixlib.f90 global.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c matrixlib.f90

primitive.o: primitive.f90 global.o matrixlib.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c primitive.f90

constraint.o: constraint.f90 global.o matrixlib.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c constraint.f90

hdlclib.o: hdlclib.f90 global.o matrixlib.o primitive.o constraint.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c hdlclib.f90

searchlib.o: searchlib.f90 global.o matrixlib.o hdlclib.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c searchlib.f90

search.o: search.f90 global.o matrixlib.o hdlclib.o searchlib.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c search.f90

mndolib.o: mndolib.f90 LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_HDLCOPT) -c mndolib.f90
