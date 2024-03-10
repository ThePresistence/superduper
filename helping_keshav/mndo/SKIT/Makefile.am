noinst_LIBRARIES = libskit.a

libskit_a_SOURCES = formats.f aplbs.f blassm.f unary.f

formats.o: formats.f
	$(F90) $(FFLAGS2_SKIT) -c formats.f

aplbs.o: aplbs.f
	$(F90) $(FFLAGS2_SKIT) -c aplbs.f

blassm.o: blassm.f
	$(F90) $(FFLAGS2_SKIT) -c blassm.f

unary.o: unary.f
	$(F90) $(FFLAGS2_SKIT) -c unary.f
