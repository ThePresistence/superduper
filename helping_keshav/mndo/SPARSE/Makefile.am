noinst_LIBRARIES = libsparse.a

libsparse_a_SOURCES =					      \
      module1.f90 module2.f90 module3.f90 traceap.f aplbp.f   \
      amubp.f     amubdgp.f   aplbdgp.f   filterp.f transpp.f \
      tracenp.f   csrcscp.f   IDAMAXP.f   coicsrp.f copmatp.f \
      BORDS1.f    CGDMSS.f    CGROOTS.f   DCARTS.f  DDOTS1.f  \
      DDOTS2.f    DIPOLS.f    DIRCGS.f    DMSCOFS.f DMSGRDS.f \
      HFOCKS.f    INDEXN.f    MATDEVS.f   PONES.f   PRTSCFS.f \
      PSORTS.f    PURIFYN.f   PURIFYP.f   SPAPRT.f  SPARSP.f  \
      SPAUPM.f    XERALL.f    XERSPA.f    DIAPRC.f  DFPUPD.f  \
      CGUPDS.f    UTIL.f      SCFDMS.f    FRAGMS.f  FILFRG.f  \
      SRTFRG.f

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_SPARSE) -c ../LIMIT.f90

module1.o: module1.f90 module3.o
	$(F90) $(FREE)$(FFLAGS2_SPARSE) -c module1.f90

module2.o: module2.f90 module3.o LIMIT.o
	$(F90) $(FREE)$(FFLAGS2_SPARSE) -c module2.f90

module3.o: module3.f90
	$(F90) $(FREE)$(FFLAGS2_SPARSE) -c module3.f90

traceap.o: traceap.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c traceap.f

aplbp.o: aplbp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c aplbp.f

amubp.o: amubp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c amubp.f

amubdgp.o: amubdgp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c amubdgp.f

aplbdgp.o: aplbdgp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c aplbdgp.f

filterp.o: filterp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c filterp.f

transpp.o: transpp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c transpp.f

tracenp.o: tracenp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c tracenp.f

csrcscp.o: csrcscp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c csrcscp.f

IDAMAXP.o: IDAMAXP.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c IDAMAXP.f

coicsrp.o: coicsrp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c coicsrp.f

copmatp.o: copmatp.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c copmatp.f

BORDS1.o: BORDS1.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c BORDS1.f

CGDMSS.o: CGDMSS.f module1.o module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c CGDMSS.f

CGROOTS.o: CGROOTS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c CGROOTS.f

DCARTS.o: DCARTS.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DCARTS.f

DDOTS1.o: DDOTS1.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DDOTS1.f

DDOTS2.o: DDOTS2.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DDOTS2.f

DIPOLS.o: DIPOLS.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DIPOLS.f

DIRCGS.o: DIRCGS.f LIMIT.o module2.o module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DIRCGS.f

DMSCOFS.o: DMSCOFS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DMSCOFS.f

DMSGRDS.o: DMSGRDS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DMSGRDS.f

HFOCKS.o: HFOCKS.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c HFOCKS.f

INDEXN.o: INDEXN.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c INDEXN.f

MATDEVS.o: MATDEVS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c MATDEVS.f

PONES.o: PONES.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c PONES.f

PRTSCFS.o: PRTSCFS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c PRTSCFS.f

PSORTS.o: PSORTS.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c PSORTS.f

PURIFYN.o: PURIFYN.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c PURIFYN.f

PURIFYP.o: PURIFYP.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c PURIFYP.f

SPAPRT.o: SPAPRT.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c SPAPRT.f

SPARSP.o: SPARSP.f LIMIT.o module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c SPARSP.f

SPAUPM.o: SPAUPM.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c SPAUPM.f

XERALL.o: XERALL.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c XERALL.f

XERSPA.o: XERSPA.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c XERSPA.f

DIAPRC.o: DIAPRC.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DIAPRC.f

DFPUPD.o: DFPUPD.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c DFPUPD.f

CGUPDS.o: CGUPDS.f module3.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c CGUPDS.f

UTIL.o: UTIL.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c UTIL.f

SCFDMS.o: SCFDMS.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c SCFDMS.f

FRAGMS.o: FRAGMS.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c FRAGMS.f

FILFRG.o: FILFRG.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c FILFRG.f

SRTFRG.o: SRTFRG.f
	$(F90) $(FIXED)$(FFLAGS2_SPARSE) -c SRTFRG.f
