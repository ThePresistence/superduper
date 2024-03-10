noinst_LIBRARIES = libanalytic.a

libanalytic_a_SOURCES = 			 \
d1mch.f   dexint.f  psbbla.f  psbpsz.f  pscfca.f \
psci.f    psdadm.f  psdext.f  psdipd.f  psdmpr.f \
psdprt.f  psdrv.f   psds1a.f  psdst1.f  psdst2.f \
psdst3.f  psdst4.f  psdxtr.f  psfock.f  psgvml.f \
pshalf.f  pshini.f  psidr.f   psidrd.f  psidre.f \
psidrp.f  psidrq.f  psidrs.f  psidrt.f  psil.f   \
psildd.f  psildp.f  psilds.f  psildx.f  psilgr.f \
psilpd.f  psilpp.f  psilps.f  psilpx.f  psilsd.f \
psilsp.f  psilss.f  psilsx.f  psint1.f  psints.f \
psir1.f   psir1d.f  psir1p.f  psir1s.f  psir2.f  \
psistb.f  psistd.f  psistf.f  psistp.f  psists.f \
pskrow.f  pslins.f  psmmgr.f  psmutl.f  psnmr.f  \
psnmrc.f  psnmrs.f  psov3b.f  psov3c.f  psov3d.f \
psov3e.f  psover.f  psovgs.f  psovx.f   psrm.f   \
pstran.f  psxutl.f  psomx.f   psgsbp.f  COMPHF.f

LIMIT.o: ../LIMIT.f90
	$(F90) $(FREE)$(FFLAGS2_PS) -c ../LIMIT.f90

d1mch.o: d1mch.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c d1mch.f

# For the following file the optimization level may cause
# errors with some compilers. Reduce level if necessary.
dexint.o: dexint.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c dexint.f

psbbla.o: psbbla.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psbbla.f

psbpsz.o: psbpsz.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psbpsz.f

pscfca.o: pscfca.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pscfca.f

psci.o: psci.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psci.f

psdadm.o: psdadm.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdadm.f

psdext.o: psdext.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdext.f

psdipd.o: psdipd.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdipd.f

psdmpr.o: psdmpr.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdmpr.f

psdprt.o: psdprt.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdprt.f

psdrv.o: psdrv.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdrv.f

psomx.o: psomx.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psomx.f

psds1a.o: psds1a.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psds1a.f

psdst1.o: psdst1.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdst1.f

psdst2.o: psdst2.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdst2.f

psdst3.o: psdst3.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdst3.f

psdst4.o: psdst4.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdst4.f

psdxtr.o: psdxtr.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psdxtr.f

psfock.o: psfock.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psfock.f

psgvml.o: psgvml.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psgvml.f

pshalf.o: pshalf.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pshalf.f

pshini.o: pshini.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pshini.f

psidr.o: psidr.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidr.f

psidrd.o: psidrd.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidrd.f

psidre.o: psidre.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidre.f

psidrp.o: psidrp.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidrp.f

psidrq.o: psidrq.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidrq.f

psidrs.o: psidrs.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidrs.f

psidrt.o: psidrt.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psidrt.f

psil.o: psil.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psil.f

psildd.o: psildd.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psildd.f

psildp.o: psildp.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psildp.f

psilds.o: psilds.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilds.f

psildx.o: psildx.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psildx.f

psilgr.o: psilgr.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilgr.f

psilpd.o: psilpd.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilpd.f

psilpp.o: psilpp.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilpp.f

psilps.o: psilps.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilps.f

psilpx.o: psilpx.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilpx.f

psilsd.o: psilsd.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilsd.f

psilsp.o: psilsp.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilsp.f

psilss.o: psilss.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilss.f

psilsx.o: psilsx.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psilsx.f

psint1.o: psint1.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psint1.f

psints.o: psints.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psints.f

psir1.o: psir1.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psir1.f

psir1d.o: psir1d.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psir1d.f

psir1p.o: psir1p.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psir1p.f

psir1s.o: psir1s.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psir1s.f

psir2.o: psir2.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psir2.f

psistb.o: psistb.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psistb.f

psistd.o: psistd.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psistd.f

psistf.o: psistf.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psistf.f

psistp.o: psistp.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psistp.f

psists.o: psists.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psists.f

pskrow.o: pskrow.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pskrow.f

pslins.o: pslins.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pslins.f

psmmgr.o: psmmgr.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psmmgr.f

psmutl.o: psmutl.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psmutl.f

psnmr.o: psnmr.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psnmr.f

psnmrc.o: psnmrc.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psnmrc.f

psnmrs.o: psnmrs.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psnmrs.f

psov3b.o: psov3b.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psov3b.f

psov3c.o: psov3c.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psov3c.f

psov3d.o: psov3d.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psov3d.f

psov3e.o: psov3e.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psov3e.f

psover.o: psover.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psover.f

psovgs.o: psovgs.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psovgs.f

psovx.o: psovx.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psovx.f

# For the following file the optimization level may cause
# errors with some compilers. Reduce level if necessary.
psrm.o: psrm.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psrm.f

pstran.o: pstran.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c pstran.f

psxutl.o: psxutl.f
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psxutl.f

psgsbp.o: psgsbp.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c psgsbp.f

COMPHF.o: COMPHF.f LIMIT.o
	$(F90) $(FIXED)$(FFLAGS2_PS) -c COMPHF.f
