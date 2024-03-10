noinst_LIBRARIES = libbl2011.a

if COND_BLAS3
MAYBE_BLAS3 = dblas3.f
endif

libbl2011_a_SOURCES = dblas1.f dblas2.f dgemm.f $(MAYBE_BLAS3) lsame.f xerbla.f zblas.f

dblas1.o: dblas1.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c dblas1.f

dblas2.o: dblas2.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c dblas2.f

dblas3.o: dblas3.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c dblas3.f

dgemm.o: dgemm.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c dgemm.f

lsame.o: lsame.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c lsame.f

xerbla.o: xerbla.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c xerbla.f

zblas.o: zblas.f
	$(F90) $(FIXED)$(FFLAGS2_BLAS) -c zblas.f
