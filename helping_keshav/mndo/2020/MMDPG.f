      SUBROUTINE MMDPG (PRT)
C     *
C     CALCULATION OF DISPERSION CORRECTIONS TO THE GRADIENT.
C     *
C     INPUT AND OUTPUT VIA PARAMETER LIST AND COMMON BLOCKS.
C     PRT        PRINTING FLAG DEFINED IN SCFCAL (I)
C     VS6,VCD    GLOBAL DISPERSION CORRECTION PARAMETERS (I)
C     CG(3,N)    DISPERSION CORRECTION TO THE CARTESIAN GRADIENT (O)
C     *
      USE LIMIT, ONLY: LM1, LM1M
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PRT
      LOGICAL noabc
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT
     ./CGRAD / CG(3,LM1+LM1M)
     ./INOPT2/ IN2(300)
     ./OMXDC / VS6,VCD
C *** INITIALIZATION.
      IOP    = IN2(2)
      IMMDP  = IN2(30)
C *** COMPUTE CARTESIAN DERIVATIVES CG OF DISPERSION CORRECTION
C     AM1-D AND PM3-D. SEE PCCP 9, 2362 (2007).
      IF((IOP.EQ.-2.OR.IOP.EQ.-7).AND.IMMDP.EQ.1) THEN
        DO I=1,NUMAT
          DO J=1,NUMAT
            IF (J.EQ.I) CYCLE
            TEMP  = EXP(-VCD*(RIJNM(I,J)/R0IJ(I,J)-1.0D0))
            CGTMP =
     1      ((VS6*C6IJ(I,J))/((1.0D0+TEMP)*RIJNM(I,J)**7))
     2      *((VCD*TEMP)/((1.0+TEMP)*R0IJ(I,J))-6.0D0/RIJNM(I,J))
C           CONVERT FROM J MOL-1 NM-1 TO KCAL MOL-1 ANGSTROM-1
            CGTMP = CGTMP/418404.3D0
            DO K=1,3
              CG(K,I)=CG(K,I)+CGTMP*(COORD(K,J)-COORD(K,I))
            ENDDO
          ENDDO
        ENDDO
C     ELSTNER DAMPING FUNCTION IN JCP 114, 5149-5155 (2001).
C     IMPLEMENTED FOR OM2 AND OM3.
      ELSE IF((IOP.EQ.-6.OR.IOP.EQ.-8).AND.IMMDP.EQ.1) THEN
        DO I=1,NUMAT
          DO J=1,NUMAT
            IF (J.EQ.I) CYCLE
            TEMP  = EXP(-VCD*(RIJNM(I,J)/R0IJ(I,J))**7)
            CGTMP  =
     1      ((VS6*C6IJ(I,J)*(1.0D0-TEMP)**3)/(RIJNM(I,J)))
     2      *((28.0D0*VCD*TEMP)/(R0IJ(I,J)**7)
     3      -6.0D0*(1.0D0-TEMP)/(RIJNM(I,J)**7))
C           CONVERT FROM J MOL-1 NM-1 TO KCAL MOL-1 ANGSTROM-1
            CGTMP = CGTMP/418404.3D0
            DO K=1, 3
              CG(K,I)=CG(K,I)+CGTMP*(COORD(K,J)-COORD(K,I))
            ENDDO
          ENDDO
        ENDDO
C     D2 CORRECTION WITH YANG DAMPING FUNCTION. JCC 27, 1787 (2006).
C     IMPLEMENTED FOR OM2 AND OM3.
      ELSE IF((IOP.EQ.-6.OR.IOP.EQ.-8).AND.IMMDP.EQ.2) THEN
        DO I=1,NUMAT
          DO J=1,NUMAT
            IF (J.EQ.I) CYCLE
            TEMP  = EXP(-VCD*(RIJNM(I,J)/R0IJ(I,J)-1.0D0))
            CGTMP =
     1      ((VS6*C6IJ(I,J))/((1.0D0+TEMP)*RIJNM(I,J)**7))
     2      *((VCD*TEMP)/((1.0+TEMP)*R0IJ(I,J))-6.0D0/RIJNM(I,J))
C           CONVERT FROM J MOL-1 NM-1 TO KCAL MOL-1 ANGSTROM-1
            CGTMP = CGTMP/418404.3D0
            DO K=1, 3
              CG(K,I)=CG(K,I)+CGTMP*(COORD(K,J)-COORD(K,I))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
C     D3 correction from Grimme, JCP 132, 154104 (2010)
C     Implemented by Lasse Spoerkel; Copied and adapted from
C     http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=DFT-D3&lang=english
C     Used for OM2, OM3, ODM2 and ODM3.
      IF((IOP.EQ.-6.OR.IOP.EQ.-8.OR.IOP.EQ.-22.OR.IOP.EQ.-23)
     1        .AND.(IMMDP.EQ.3.OR.IMMDP.EQ.-3)) THEN
        IF (IMMDP.EQ.3) THEN
          noabc  = .TRUE.
        ELSE
          noabc  = .FALSE.
        ENDIF
        CALL DFTD3(.FALSE.,.TRUE.,noabc,PRT)
      ENDIF
      RETURN
      END SUBROUTINE MMDPG
C     ******************************************************************
C     The following routines were downloaded from the Grimme website.
C     Changes in DFTD3: Input of printing flag PRT via parameter list.
C     Changes in DFTD3: Input and output data via mndo common blocks.
C     Changes in DFTD3: Printing governed by PRT, IN2(41), and IN2(71).
C     Further changes are merely cosmetic - modified comments.
C     ******************************************************************
      SUBROUTINE DFTD3(ene,grad,noabc,PRT)
C     *
C     ene    flag for energy calculation (I)
c     grad   flag for gradient calculation (I)
c     noabc  flag for inclusion of three-center contributions (I)
c     PRT    printing flag (I)
C     *
      USE LIMIT, ONLY: LM1, LM1M, LMZ
      IMPLICIT NONE
      COMMON
     ./ATOMC / COORD(3,LM1)
     ./ATOMS / NUMAT,NAT(LM1)
     ./CGRAD / CG(3,LM1+LM1M)
     ./INOPT2/ IN2(300)
     ./MMDP  / EMMDP
     ./MMDPP / C6(LMZ),R0(LMZ),D3S6,D3S8,D3A1,D3A2
      INTEGER NUMAT,NAT,IN2
      REAL*8 COORD,CG,EMMDP,C6,R0,D3S6,D3S8,D3A1,D3A2
      REAL*8,DIMENSION(:,:), ALLOCATABLE :: g
c
c     input variables
      LOGICAL ene, grad, noabc,PRT
      INTEGER max_elem,maxc
      PARAMETER (max_elem=94)
c     maximum coordination number references per element
      PARAMETER (maxc=5)
c     conversion factors
      REAL*8 autoang,autokcal,autoev
      PARAMETER (autoang = 0.52917726d0)
      PARAMETER (autokcal= 627.509541d0)
      PARAMETER (autoev  = 27.21138505d0)
c
c     number of atoms
      INTEGER n
c     coordinates in au
      REAL*8,DIMENSION(:,:), ALLOCATABLE :: xyz
c     cardinal numbers of elements
      INTEGER,DIMENSION(:), ALLOCATABLE :: iz
c     cut-off radii for all element pairs
      REAL*8 r0ab(max_elem,max_elem)
c     c6 for all element pairs
      REAL*8 c6ab(max_elem,max_elem,maxc,maxc,3)
c     how many dIFferent c6 for one element
      INTEGER mxc(max_elem)
c     coordination numbers of the atoms
      REAL*8,DIMENSION(:), ALLOCATABLE :: cn
c     energies
      REAL*8 e6, e8, e6abc, disp
c     the PARAMETERs of the method (not all a "free")
      REAL*8 rs6, rs8, s6, s18, rs18, alp, alp6, alp8
c     r^2 distance neglect threshold (important for speed in CASE of large systems)
      REAL*8 rthr,rthr2
c     r^2 distance to cutoff for cn_calculation
      REAL*8 cn_thr
c     local and dummy variables
      INTEGER i,j,k,z,nn,iat,jat,i1,i2
      REAL*8  x,gdsp
      INTEGER iadr,jadr,kk,nlines
c
      INCLUDE 'MMDPPAR.f'
c     cutoff r^2 thresholds for the gradient in bohr^2.
c     rthr influences n^2 part of the gradient.
c     rthr2 influences the n^3 part of the gradient. when using
c     When using dftd3 in combination with semi-empirical methods or ffs,
c     and large (>1000 atoms) systems, rthr2 is crucial for speed:
c     recommended values are 20^2 to 25^2 bohr.
      rthr   = 9000.0d0   ! ur, se
      rthr2  = 9000.0d0
      cn_thr = 9000.0d0
c     set radii
      CALL setr0ab(max_elem,autoang,r0ab)
c
      c6ab   =-1
      mxc    = 0
      kk     = 1
      DO i=1,nlines
        iat=INT(pars(kk+1))
        jat=INT(pars(kk+2))
        CALL limit_adr(iat,jat,iadr,jadr)
        mxc(iat)=MAX(mxc(iat),iadr)
        mxc(jat)=MAX(mxc(jat),jadr)
c
        c6ab(iat,jat,iadr,jadr,1)=pars(kk)
        c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
        c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)
c
        c6ab(jat,iat,jadr,iadr,1)=pars(kk)
        c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
        c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
        kk=(i*5)+1
      ENDDO
c
      n = numat
      ALLOCATE(iz(n))
      DO i=1,n
        iz(i)=nat(i)
      ENDDO
      ALLOCATE(xyz(3,n))
      DO i=1,n
        DO k=1,3
          xyz(k,i)=coord(k,i)/autoang
        ENDDO
      ENDDO
      ALLOCATE(g(3,n))
      ALLOCATE(cn(n))
c     set constants
      s6   = D3S6
      s18  = D3S8
      rs6  = D3A1
      rs18 = D3A2
      alp  = 14.0D0
c
      rs8  = rs18
      alp6 = alp
      alp8 = alp+2.d0
c     *
c *** compute d3 dispersion energy correction.
c     *
      IF (ene .eq. .TRUE.) THEN
        CALL edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     1     rs6,rs8,alp6,alp8,noabc,rthr,cn_thr,e6,e8,e6abc)
        e6   = e6*s6
        e6abc= e6abc*s6
        e8   = e8*s18
        disp = -e6-e8-e6abc
c       printing section.
        IF (PRT .and. IN2(72).ge.0) THEN
          WRITE(*,*) ""
          WRITE(*,'(A,f11.4)') ' D3-disp E (kcal/mol)= ',disp*autokcal
          WRITE(*,'(A,f11.4)') '   E6      (kcal/mol)= ',-e6*autokcal
          WRITE(*,'(A,f11.4)') '   E8      (kcal/mol)= ',-e8*autokcal
          IF (.NOT.noabc) THEN
           WRITE(*,'(A,f11.4)')'   E6(ABC) (kcal/mol)= ',-e6abc*autokcal
          END IF
        END IF
        emmdp=disp*autokcal
      END IF
c     *
c *** compute d3 dispersion energy gradient correction.
c     *
      IF (grad .eq. .TRUE.) THEN
        g=0.0d0
        CALL gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     1           s6,s18,rs6,rs8,alp6,alp8,noabc,rthr,rthr2,
     2           g,gdsp,x,cn_thr)
        DO i=1,n
          DO k=1,3
            cg(k,i)=cg(k,i) + g(k,i)*autokcal/autoang
          ENDDO
        ENDDO
c       printing section.
        IF (PRT .and. IN2(41).ge.5) THEN
          WRITE(*,*) ""
          WRITE(*,'(A)') " D3 dispersion correction to gradient"
          WRITE(*,'(a5,3(a15))') "NI","X","Y","Z"
          DO i=1,n
            WRITE(*,'(i5,3(f15.8))') i,g(1:3,i)*autokcal/autoang
          ENDDO
        END IF
        IF (PRT .and. IN2(41).ge.2) THEN
          WRITE(*,*)
          WRITE(*,'(A,f11.4)')' D3-disp |g 2-norm| (kcal/(mol*ang))= ',
     .          x*autokcal/autoang
        END IF
      END IF
c
      DEALLOCATE(xyz,g,iz,cn)
      END SUBROUTINE dftd3
C     ******************************************************************
c     Compute dispersion energy
C     ******************************************************************
      SUBROUTINE edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     1           rs6,rs8,alp6,alp8,noabc,rthr,cn_thr,e6,e8,e63)
      IMPLICIT NONE
      INTEGER n,iz(*),max_elem,maxc,mxc(max_elem)
      REAL*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 rs6,rs8,alp8,alp6,rcov(max_elem)
      REAL*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      REAL*8 e6, e8 ,e63
      logical noabc
c
      integer iat,jat,kat
      real*8 r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,ang,rav
      real*8 rr,c9,rthr,cn_thr
      real*8 cn(n)
      real*8 r2ab(n*n),cc6ab(n*n),dmp(n*n),t1,t2,t3,a1,a2,tmp2
      real*8 abcthr
      integer*2 icomp(n*n)
      integer ij,ik,jk
      INTEGER,EXTERNAL :: lin
c
      e6 =0.0d0
      e8 =0.0d0
      abcthr=cn_thr ! the threebody term uses the same threshold as the CN-calculation
c     becke-johnson PARAMETERs
      a1=rs6
      a2=rs8
c dft-d3
      CALL ncoord(n,rcov,iz,xyz,cn,cn_thr)

      icomp=0
      DO iat=1,n-1
         DO jat=iat+1,n
            dx=xyz(1,iat)-xyz(1,jat)
            dy=xyz(2,iat)-xyz(2,jat)
            dz=xyz(3,iat)-xyz(3,jat)
            r2=dx*dx+dy*dy+dz*dz
cthr
            IF(r2.GT.rthr) CYCLE
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r
c get c6
            CALL getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                    cn(iat),cn(jat),c6)

            r6=r2**3
            r8=r6*r2
c r2r4 stored in main as SQRT
            c8 =3.0d0*c6*r2r4(iz(iat))*r2r4(iz(jat))

            tmp=SQRT(c8/c6)
            e6=e6+  c6/(r6+(a1*tmp+a2)**6)
            e8=e8+  c8/(r8+(a1*tmp+a2)**8)

            if((.not.noabc).and.(r2.lt.abcthr)) then
              ij=lin(jat,iat)
              icomp(ij)=1
c store C6 for C9, calc as sqrt
              cc6ab(ij)=sqrt(c6)
c store R^2 for abc
              r2ab(ij)=r2
c store for abc damping
              dmp(ij)=(1.d0/rr)**(1.d0/3.d0)
            endif !noabc
         ENDDO
      ENDDO
      if(noabc)return

C compute non-additive third-order energy using averaged C6
      do iat=1,n
         do jat=1,iat-1
            ij=lin(jat,iat)
            if(icomp(ij).eq.1)then
            do kat=1,jat-1

              ik=lin(kat,iat)
              jk=lin(kat,jat)
              if((icomp(ik).eq.0).or.(icomp(jk).eq.0)) cycle
c damping func product
              rav=(4.d0/3.d0)/(dmp(ik)*dmp(jk)*dmp(ij))
              tmp=1.d0/( 1.d0+6.d0*rav**alp8 )
c triple C6 coefficient (stored as sqrt)
              c9=cc6ab(ij)*cc6ab(ik)*cc6ab(jk)
c angular terms with law of cosines

              t1 = (r2ab(ij)+r2ab(jk)-r2ab(ik))
              t2 = (r2ab(ij)+r2ab(ik)-r2ab(jk))
              t3 = (r2ab(ik)+r2ab(jk)-r2ab(ij))
              tmp2=r2ab(ij)*r2ab(jk)*r2ab(ik)
              ang=(0.375d0*t1*t2*t3/tmp2+1.0d0)/tmp2**1.5d0

              e63=e63-tmp*c9*ang

            enddo
            endif
         enddo
      enddo
      RETURN
      END SUBROUTINE edisp
C     ******************************************************************
c     Compute gradient of dispersion energy
C     ******************************************************************
      SUBROUTINE gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .                 s6,s18,rs6,rs8,alp6,alp8,noabc,rthr,
     .                 rthr2,g,disp,gnorm,cn_thr)
      IMPLICIT NONE
c k1-k3
      REAL*8 k1
c global ad hoc PARAMETERs
      PARAMETER (k1=16.0d0)

      integer n,iz(*),max_elem,maxc,mxc(max_elem)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 g(3,*),s6,s18,rcov(max_elem)
      real*8 rs6,rs8,alp8,alp6,a1,a2
      logical noabc

      integer iat,jat,i,j,kat,k
      real*8 R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc
      real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
      real*8 damp6,damp8,damp10,e6,e8,e10,e12,gnorm,tmp1
      real*8 s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
      real*8 cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,c8,cn_thr
      real*8 rthr2,rthr3

      real*8 rij(3),rik(3),rjk(3),r7,r9
      real*8 rik_dist,rjk_dist
      real*8 drij(n*(n+1)/2)  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
      real*8 drik,drjk
      real*8 rcovij
      real*8 dc6,c6chk !d(C6ij)/d(r_ij)
      real*8 expterm,dcni
      real*8 dcn               !dCN(iat)/d(r_ij) is equal to
                               !dCN(jat)/d(r_ij)
      real*8 dc6_rest ! saves (1/r^6*f_dmp + 3*r4r2/r^8*f_dmp) for kat loop
      integer,external :: lin
      integer linij,linik,linjk
      real*8 vec(3),vec2(3)
      real*8 dc6i(n)       ! dE_disp/dCN(iat) in dc6i(iat)
      real*8 dc6ij(n,n)    ! saves dC6(ij)/dCN(iat)
      real*8 dc6iji,dc6ijj
      logical abccalc(n*(n+1)/2)
      real*8 abcthr
      real*8 labc,rabc
      real*8 c6abc(n*(n+1)/2)
      real*8 r2abc(n*(n+1)/2)
      real*8 r3abc(n*(n+1)/2)
      real*8 c9,rav,rav3,fdmp,ang,angr9,eabc,dc9,dfdmp,dang
      real*8 r2ij,r2jk,r2ik,mijk,imjk,ijmk,rijk3
      integer mat,linim,linjm,linkm,kk

      dc6i=0.0d0
      abccalc=.FALSE.
      abcthr=rthr2

c Becke-Johnson finite damping
      a1 =rs6
      a2 =rs8
      s8 =s18

      disp=0.d0
      call ncoord(n,rcov,iz,xyz,cn,cn_thr)

      drij=0.0d0
      dc6_rest=0.0d0
      kat=0

      do iat=1,n
        do jat=1,iat-1
          rij=xyz(:,jat)-xyz(:,iat)
          r2=sum(rij*rij)
          if (r2.gt.rthr) cycle

          linij=lin(iat,jat)
          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
!
!      get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
!      dC6(iat,jat)/dCN(jat).
!
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),
     .          mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),
     .          c6,dc6iji,dc6ijj)

          r=dsqrt(r2)
          r4=r2*r2
          r6=r4*r2
          r8=r6*r2

          if ((.not.noabc).and.(r2.lt.abcthr)) then
            abccalc(linij)=.TRUE.
            dc6ij(iat,jat)=dc6iji
            dc6ij(jat,iat)=dc6ijj
            c6abc(linij)=c6
            r2abc(linij)=r2
            r3abc(linij)=(r/r0)**(1.0d0/3.0d0)
          endif !noabc
c use BJ radius
          R0=a1*DSQRT(3.0d0*r42)+a2

          t6=(r6+R0**6)
          t8=(r8+R0**8)

          drij(linij)=drij(linij)
     .        -s6*C6*6.0d0*r4*r/(t6*t6)
     .        -s8*C6*24.0d0*r42*r6*r/(t8*t8)

          dc6_rest=s6/t6+3.d0*s8*r42/t8
          disp=disp-dc6_rest*c6  ! calculate E_disp for sanity check

!     saving all (1/r^6...)* dC6/dCN(i) for each atom
          dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
          dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj

        enddo !jat
      enddo !iat

      if(.not.noabc)then

        do iat=1,n
          do jat=1,iat-1
            linij=lin(iat,jat)
            if(.NOT.abccalc(linij))cycle
              r2ij=r2abc(linij)
            do kat=1,jat-1

              linik=lin(iat,kat)
              linjk=lin(jat,kat)
              if(.NOT.(abccalc(linjk).AND.abccalc(linik)))cycle    !cutoff crit
! calculating the 3body energy:
              r2jk=r2abc(linjk)
              r2ik=r2abc(linik)
              c9=c6abc(linij)*c6abc(linjk)*c6abc(linik)
              c9=dsqrt(c9)
              rav=r3abc(linij)*r3abc(linjk)*r3abc(linik)
              fdmp=1.0d0/(1+6.0d0*(0.75d0*rav)**(-alp8))
              mijk=-r2ij+r2jk+r2ik
              imjk= r2ij-r2jk+r2ik
              ijmk= r2ij+r2jk-r2ik
              rijk3=r2ij*r2jk*r2ik
              rav3=rijk3**1.5d0
              ang=0.375d0*ijmk*imjk*mijk/rijk3
              angr9=(ang +1.0d0)/rav3

              eabc=eabc+c9*angr9*fdmp
!end of 3body energy calculation

!start calculating the derivatives of each part w.r.t. r_ij
              r=dsqrt(r2ij)
              dfdmp=-2.d0*alp8*(0.75d0*rav)**(-alp8)*fdmp*fdmp

              dang=-0.375d0*(r2ij**3+r2ij**2*(r2jk+r2ik)
     .             +r2ij*(3.0d0*r2jk**2+2.0d0*r2jk*r2ik+3.0d0*r2ik**2)
     .             -5.0d0*(r2jk-r2ik)**2*(r2jk+r2ik))
     .             /(r*rijk3*rav3)

              tmp1=dfdmp/r*c9*angr9-dang*c9*fdmp
              drij(linij)=drij(linij)+tmp1

!start calculating the derivatives of each part w.r.t. r_jk
              r=dsqrt(r2jk)

              dang=-0.375d0*(r2jk**3+r2jk**2*(r2ik+r2ij)
     .             +r2jk*(3.0d0*r2ik**2+2.0d0*r2ik*r2ij+3.0d0*r2ij**2)
     .             -5.0d0*(r2ik-r2ij)**2*(r2ik+r2ij))
     .             /(r*rijk3*rav3)

              drij(linjk)=drij(linjk)
     .                       +dfdmp/r*c9*angr9-dang*c9*fdmp

!start calculating the derivatives of each part w.r.t. r_ik
              r=dsqrt(r2abc(linik))

              dang=-0.375d0*(r2ik**3+r2ik**2*(r2jk+r2ij)
     .             +r2ik*(3.0d0*r2jk**2+2.0d0*r2jk*r2ij+3.0d0*r2ij**2)
     .             -5.0d0*(r2jk-r2ij)**2*(r2jk+r2ij))
     .             /(r*rijk3*rav3)

              drij(linik)=drij(linik)
     .                       +dfdmp/r*c9*angr9-dang*c9*fdmp

! calculate rest* dc9/dcn(iat)  and sum it up for every atom ijk
              dc6_rest=angr9*fdmp

             dc9=dc6ij(iat,jat)/c6abc(linij)+dc6ij(iat,kat)/c6abc(linik)
              dc9=-0.5d0*c9*dc9
              dc6i(iat)=dc6i(iat)+dc6_rest*dc9

             dc9=dc6ij(jat,iat)/c6abc(linij)+dc6ij(jat,kat)/c6abc(linjk)
              dc9=-0.5d0*c9*dc9
              dc6i(jat)=dc6i(jat)+dc6_rest*dc9

             dc9=dc6ij(kat,iat)/c6abc(linik)+dc6ij(kat,jat)/c6abc(linjk)
              dc9=-0.5d0*c9*dc9
              dc6i(kat)=dc6i(kat)+dc6_rest*dc9

            ENDDO !kat
          ENDDO !jat
        ENDDO !iat

        disp=disp+eabc
      endif !noabc

! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
      do iat=2,n
        do jat=1,iat-1
          linij=lin(iat,jat)
          rij=xyz(:,jat)-xyz(:,iat)

          r2=sum(rij*rij)
          r=dsqrt(r2)
          if (r2.lt.cn_thr) then
            rcovij=rcov(iz(iat))+rcov(iz(jat))
            expterm=exp(-k1*(rcovij/r-1.d0))
            dcn=-k1*rcovij*expterm/
     .                   (r*r*(expterm+1.d0)*(expterm+1.d0))
          else
            dcn=0.d0
          endif
          x1=drij(linij)+dcn*(dc6i(iat)+dc6i(jat))

          g(:,iat)=g(:,iat)+x1*rij/r
          g(:,jat)=g(:,jat)-x1*rij/r

        enddo !iat
      enddo !jat

      gnorm=sum(abs(g(1:3,1:n)))

      END SUBROUTINE gdisp

cccccccccccccccccccccccccccccccccccccccc
c      the   n e w   gradc6 routine    c
cccccccccccccccccccccccccccccccccccccccc
!
      SUBROUTINE get_dc6_dcnij(maxc,max_elem,c6ab,mxci,mxcj,cni,cnj,
     .           izi,izj,c6check,dc6i,dc6j)

      IMPLICIT NONE
c k1-k3
      REAL*8 k3
c global ad hoc PARAMETERs
      PARAMETER (k3=-4.0d0)

      INTEGER maxc,max_elem
      REAL*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      INTEGER mxci,mxcj
      REAL*8 cni,cnj,term
      INTEGER izi,izj
      REAL*8  dc6i,dc6j,c6check

      INTEGER a,b
      REAL*8 zaehler,nenner,dzaehler_i,dnenner_i,dzaehler_j,dnenner_j
      REAL*8 expterm,cn_refi,cn_refj,c6ref,r
      REAL*8 c6mem,r_save

      c6mem=-1.d99
      r_save=9999.0d0
      zaehler=0.0d0
      nenner=0.0d0

      dzaehler_i=0.d0
      dnenner_i=0.d0
      dzaehler_j=0.d0
      dnenner_j=0.d0

      DO a=1,mxci
        DO b=1,mxcj
          c6ref=c6ab(izi,izj,a,b,1)
          IF (c6ref.GT.0) THEN
            cn_refi=c6ab(izi,izj,a,b,2)
            cn_refj=c6ab(izi,izj,a,b,3)
            r=(cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj)
            IF (r.LT.r_save) THEN
               r_save=r
               c6mem=c6ref
            ENDIF
            expterm=exp(k3*r)
            zaehler=zaehler+c6ref*expterm
            nenner=nenner+expterm
            expterm=expterm*2.d0*k3
            term=expterm*(cni-cn_refi)
            dzaehler_i=dzaehler_i+c6ref*term
            dnenner_i =dnenner_i +      term

            term=expterm*(cnj-cn_refj)
            dzaehler_j=dzaehler_j+c6ref*term
            dnenner_j =dnenner_j +      term
          ENDIF
        ENDDO !b
      ENDDO !a

      IF (nenner.GT.1.0d-99) THEN
        c6check=zaehler/nenner
        dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler))
     .    /(nenner*nenner)
        dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler))
     .    /(nenner*nenner)
      ELSE
        c6check=c6mem
        dc6i=0.0d0
        dc6j=0.0d0
      ENDIF
      END SUBROUTINE get_dc6_dcnij

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INTerpolate c6
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
      IMPLICIT NONE
      INTEGER maxc,max_elem
      INTEGER iat,jat,i,j,mxc(max_elem)
      REAL*8  nci,ncj,c6,c6mem
      REAL*8  c6ab(max_elem,max_elem,maxc,maxc,3)
c the exponential is sensitive to numerics
c when nci or ncj is much larger than cn1/cn2
      REAL*8  cn1,cn2,r,rsum,csum,tmp1
      REAL*8  r_save
c k1-k3
      REAL*8 k3
c global ad hoc PARAMETERs
      PARAMETER (k3=-4.0d0)

      c6mem=-1.d+99
      rsum=0.0d0
      csum=0.0d0
      c6  =0.0d0
      r_save=1.0d99
      DO i=1,mxc(iat)
      DO j=1,mxc(jat)
         c6=c6ab(iat,jat,i,j,1)
         IF(c6.GT.0)THEN
            cn1=c6ab(iat,jat,i,j,2)
            cn2=c6ab(iat,jat,i,j,3)
c distance
            r=(cn1-nci)**2+(cn2-ncj)**2
            IF (r.LT.r_save) THEN
               r_save=r
               c6mem=c6
            ENDIF
            tmp1=EXP(k3*r)
            rsum=rsum+tmp1
            csum=csum+tmp1*c6
         ENDIF
      ENDDO
      ENDDO

      IF(rsum.GT.1.0d-99)THEN
         c6=csum/rsum
      ELSE
         c6=c6mem
      ENDIF

      END SUBROUTINE getc6

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c compute coordination numbers by adding an inverse damping FUNCTION
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
      IMPLICIT NONE
c global ad hoc PARAMETERs
      REAL*8 k1
      PARAMETER (k1=16.0d0)

      INTEGER iz(*),natoms,i
      REAL*8 xyz(3,*),cn(*),rcov(94)
      REAL*8 cn_thr

      INTEGER iat
      REAL*8 dx,dy,dz,r,damp,xn,rr,rco,r2

      DO i=1,natoms
        xn=0.0d0
        DO iat=1,natoms
          IF(iat.ne.i)THEN
              dx=xyz(1,iat)-xyz(1,i)
              dy=xyz(2,iat)-xyz(2,i)
              dz=xyz(3,iat)-xyz(3,i)
              r2=dx*dx+dy*dy+dz*dz
              IF (r2.GT.cn_thr) CYCLE
              r=SQRT(r2)
c covalent distance in bohr
              rco=rcov(iz(i))+rcov(iz(iat))
              rr=rco/r
c counting FUNCTION exponential has a better long-range behavior than mhgs inverse damping
              damp=1.d0/(1.d0+EXP(-k1*(rr-1.0d0)))
              xn=xn+damp
          ENDIF
        ENDDO
        cn(i)=xn
      ENDDO
      END SUBROUTINE ncoord


      INTEGER FUNCTION lin(i1,i2)
        INTEGER i1,i2,idum1,idum2
        idum1=MAX(i1,i2)
        idum2=MIN(i1,i2)
        lin=idum2+idum1*(idum1-1)/2
        RETURN
      END FUNCTION lin


      SUBROUTINE limit_adr(iat,jat,iadr,jadr)
        IMPLICIT NONE
        INTEGER iat,jat,iadr,jadr,i
        iadr=1
        jadr=1
        i=100
  10   IF(iat.GT.100) THEN
          iat=iat-100
          iadr=iadr+1
          GOTO 10
        ENDIF

        i=100
  20   IF(jat.GT.100) THEN
          jat=jat-100
          jadr=jadr+1
          GOTO 20
        ENDIF
      END SUBROUTINE limit_adr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set cut-off radii
c in parts due to INTel compiler bug
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE setr0ab(max_elem,autoang,r)
      IMPLICIT NONE
      INTEGER max_elem,i,j,k
      REAL*8 r(max_elem,max_elem),autoang
      REAL*8 r0ab(4465)
      r0ab(   1:  70)=(/
     . 2.1823d0,1.8547d0,1.7347d0,2.9086d0,2.5732d0,3.4956d0,2.3550d0
     .,2.5095d0,2.9802d0,3.0982d0,2.5141d0,2.3917d0,2.9977d0,2.9484d0
     .,3.2160d0,2.4492d0,2.2527d0,3.1933d0,3.0214d0,2.9531d0,2.9103d0
     .,2.3667d0,2.1328d0,2.8784d0,2.7660d0,2.7776d0,2.7063d0,2.6225d0
     .,2.1768d0,2.0625d0,2.6395d0,2.6648d0,2.6482d0,2.5697d0,2.4846d0
     .,2.4817d0,2.0646d0,1.9891d0,2.5086d0,2.6908d0,2.6233d0,2.4770d0
     .,2.3885d0,2.3511d0,2.2996d0,1.9892d0,1.9251d0,2.4190d0,2.5473d0
     .,2.4994d0,2.4091d0,2.3176d0,2.2571d0,2.1946d0,2.1374d0,2.9898d0
     .,2.6397d0,3.6031d0,3.1219d0,3.7620d0,3.2485d0,2.9357d0,2.7093d0
     .,2.5781d0,2.4839d0,3.7082d0,2.5129d0,2.7321d0,3.1052d0,3.2962d0
     ./)
      r0ab(  71: 140)=(/
     . 3.1331d0,3.2000d0,2.9586d0,3.0822d0,2.8582d0,2.7120d0,3.2570d0
     .,3.4839d0,2.8766d0,2.7427d0,3.2776d0,3.2363d0,3.5929d0,3.2826d0
     .,3.0911d0,2.9369d0,2.9030d0,2.7789d0,3.3921d0,3.3970d0,4.0106d0
     .,2.8884d0,2.6605d0,3.7513d0,3.1613d0,3.3605d0,3.3325d0,3.0991d0
     .,2.9297d0,2.8674d0,2.7571d0,3.8129d0,3.3266d0,3.7105d0,3.7917d0
     .,2.8304d0,2.5538d0,3.3932d0,3.1193d0,3.1866d0,3.1245d0,3.0465d0
     .,2.8727d0,2.7664d0,2.6926d0,3.4608d0,3.2984d0,3.5142d0,3.5418d0
     .,3.5017d0,2.6190d0,2.4797d0,3.1331d0,3.0540d0,3.0651d0,2.9879d0
     .,2.9054d0,2.8805d0,2.7330d0,2.6331d0,3.2096d0,3.5668d0,3.3684d0
     .,3.3686d0,3.3180d0,3.3107d0,2.4757d0,2.4019d0,2.9789d0,3.1468d0
     ./)
      r0ab( 141: 210)=(/
     . 2.9768d0,2.8848d0,2.7952d0,2.7457d0,2.6881d0,2.5728d0,3.0574d0
     .,3.3264d0,3.3562d0,3.2529d0,3.1916d0,3.1523d0,3.1046d0,2.3725d0
     .,2.3289d0,2.8760d0,2.9804d0,2.9093d0,2.8040d0,2.7071d0,2.6386d0
     .,2.5720d0,2.5139d0,2.9517d0,3.1606d0,3.2085d0,3.1692d0,3.0982d0
     .,3.0352d0,2.9730d0,2.9148d0,3.2147d0,2.8315d0,3.8724d0,3.4621d0
     .,3.8823d0,3.3760d0,3.0746d0,2.8817d0,2.7552d0,2.6605d0,3.9740d0
     .,3.6192d0,3.6569d0,3.9586d0,3.6188d0,3.3917d0,3.2479d0,3.1434d0
     .,4.2411d0,2.7597d0,3.0588d0,3.3474d0,3.6214d0,3.4353d0,3.4729d0
     .,3.2487d0,3.3200d0,3.0914d0,2.9403d0,3.4972d0,3.7993d0,3.6773d0
     .,3.8678d0,3.5808d0,3.8243d0,3.5826d0,3.4156d0,3.8765d0,4.1035d0
     ./)
      r0ab( 211: 280)=(/
     . 2.7361d0,2.9765d0,3.2475d0,3.5004d0,3.4185d0,3.4378d0,3.2084d0
     .,3.2787d0,3.0604d0,2.9187d0,3.4037d0,3.6759d0,3.6586d0,3.8327d0
     .,3.5372d0,3.7665d0,3.5310d0,3.3700d0,3.7788d0,3.9804d0,3.8903d0
     .,2.6832d0,2.9060d0,3.2613d0,3.4359d0,3.3538d0,3.3860d0,3.1550d0
     .,3.2300d0,3.0133d0,2.8736d0,3.4024d0,3.6142d0,3.5979d0,3.5295d0
     .,3.4834d0,3.7140d0,3.4782d0,3.3170d0,3.7434d0,3.9623d0,3.8181d0
     .,3.7642d0,2.6379d0,2.8494d0,3.1840d0,3.4225d0,3.2771d0,3.3401d0
     .,3.1072d0,3.1885d0,2.9714d0,2.8319d0,3.3315d0,3.5979d0,3.5256d0
     .,3.4980d0,3.4376d0,3.6714d0,3.4346d0,3.2723d0,3.6859d0,3.8985d0
     .,3.7918d0,3.7372d0,3.7211d0,2.9230d0,2.6223d0,3.4161d0,2.8999d0
     ./)
      r0ab( 281: 350)=(/
     . 3.0557d0,3.3308d0,3.0555d0,2.8508d0,2.7385d0,2.6640d0,3.5263d0
     .,3.0277d0,3.2990d0,3.7721d0,3.5017d0,3.2751d0,3.1368d0,3.0435d0
     .,3.7873d0,3.2858d0,3.2140d0,3.1727d0,3.2178d0,3.4414d0,2.5490d0
     .,2.7623d0,3.0991d0,3.3252d0,3.1836d0,3.2428d0,3.0259d0,3.1225d0
     .,2.9032d0,2.7621d0,3.2490d0,3.5110d0,3.4429d0,3.3845d0,3.3574d0
     .,3.6045d0,3.3658d0,3.2013d0,3.6110d0,3.8241d0,3.7090d0,3.6496d0
     .,3.6333d0,3.0896d0,3.5462d0,2.4926d0,2.7136d0,3.0693d0,3.2699d0
     .,3.1272d0,3.1893d0,2.9658d0,3.0972d0,2.8778d0,2.7358d0,3.2206d0
     .,3.4566d0,3.3896d0,3.3257d0,3.2946d0,3.5693d0,3.3312d0,3.1670d0
     .,3.5805d0,3.7711d0,3.6536d0,3.5927d0,3.5775d0,3.0411d0,3.4885d0
     ./)
      r0ab( 351: 420)=(/
     . 3.4421d0,2.4667d0,2.6709d0,3.0575d0,3.2357d0,3.0908d0,3.1537d0
     .,2.9235d0,3.0669d0,2.8476d0,2.7054d0,3.2064d0,3.4519d0,3.3593d0
     .,3.2921d0,3.2577d0,3.2161d0,3.2982d0,3.1339d0,3.5606d0,3.7582d0
     .,3.6432d0,3.5833d0,3.5691d0,3.0161d0,3.4812d0,3.4339d0,3.4327d0
     .,2.4515d0,2.6338d0,3.0511d0,3.2229d0,3.0630d0,3.1265d0,2.8909d0
     .,3.0253d0,2.8184d0,2.6764d0,3.1968d0,3.4114d0,3.3492d0,3.2691d0
     .,3.2320d0,3.1786d0,3.2680d0,3.1036d0,3.5453d0,3.7259d0,3.6090d0
     .,3.5473d0,3.5327d0,3.0018d0,3.4413d0,3.3907d0,3.3593d0,3.3462d0
     .,2.4413d0,2.6006d0,3.0540d0,3.1987d0,3.0490d0,3.1058d0,2.8643d0
     .,2.9948d0,2.7908d0,2.6491d0,3.1950d0,3.3922d0,3.3316d0,3.2585d0
     ./)
      r0ab( 421: 490)=(/
     . 3.2136d0,3.1516d0,3.2364d0,3.0752d0,3.5368d0,3.7117d0,3.5941d0
     .,3.5313d0,3.5164d0,2.9962d0,3.4225d0,3.3699d0,3.3370d0,3.3234d0
     .,3.3008d0,2.4318d0,2.5729d0,3.0416d0,3.1639d0,3.0196d0,3.0843d0
     .,2.8413d0,2.7436d0,2.7608d0,2.6271d0,3.1811d0,3.3591d0,3.3045d0
     .,3.2349d0,3.1942d0,3.1291d0,3.2111d0,3.0534d0,3.5189d0,3.6809d0
     .,3.5635d0,3.5001d0,3.4854d0,2.9857d0,3.3897d0,3.3363d0,3.3027d0
     .,3.2890d0,3.2655d0,3.2309d0,2.8502d0,2.6934d0,3.2467d0,3.1921d0
     .,3.5663d0,3.2541d0,3.0571d0,2.9048d0,2.8657d0,2.7438d0,3.3547d0
     .,3.3510d0,3.9837d0,3.6871d0,3.4862d0,3.3389d0,3.2413d0,3.1708d0
     .,3.6096d0,3.6280d0,3.6860d0,3.5568d0,3.4836d0,3.2868d0,3.3994d0
     ./)
      r0ab( 491: 560)=(/
     . 3.3476d0,3.3170d0,3.2950d0,3.2874d0,3.2606d0,3.9579d0,2.9226d0
     .,2.6838d0,3.7867d0,3.1732d0,3.3872d0,3.3643d0,3.1267d0,2.9541d0
     .,2.8505d0,2.7781d0,3.8475d0,3.3336d0,3.7359d0,3.8266d0,3.5733d0
     .,3.3959d0,3.2775d0,3.1915d0,3.9878d0,3.8816d0,3.5810d0,3.5364d0
     .,3.5060d0,3.8097d0,3.3925d0,3.3348d0,3.3019d0,3.2796d0,3.2662d0
     .,3.2464d0,3.7136d0,3.8619d0,2.9140d0,2.6271d0,3.4771d0,3.1774d0
     .,3.2560d0,3.1970d0,3.1207d0,2.9406d0,2.8322d0,2.7571d0,3.5455d0
     .,3.3514d0,3.5837d0,3.6177d0,3.5816d0,3.3902d0,3.2604d0,3.1652d0
     .,3.7037d0,3.6283d0,3.5858d0,3.5330d0,3.4884d0,3.5789d0,3.4094d0
     .,3.3473d0,3.3118d0,3.2876d0,3.2707d0,3.2521d0,3.5570d0,3.6496d0
     ./)
      r0ab( 561: 630)=(/
     . 3.6625d0,2.7300d0,2.5870d0,3.2471d0,3.1487d0,3.1667d0,3.0914d0
     .,3.0107d0,2.9812d0,2.8300d0,2.7284d0,3.3259d0,3.3182d0,3.4707d0
     .,3.4748d0,3.4279d0,3.4182d0,3.2547d0,3.1353d0,3.5116d0,3.9432d0
     .,3.8828d0,3.8303d0,3.7880d0,3.3760d0,3.7218d0,3.3408d0,3.3059d0
     .,3.2698d0,3.2446d0,3.2229d0,3.4422d0,3.5023d0,3.5009d0,3.5268d0
     .,2.6026d0,2.5355d0,3.1129d0,3.2863d0,3.1029d0,3.0108d0,2.9227d0
     .,2.8694d0,2.8109d0,2.6929d0,3.1958d0,3.4670d0,3.4018d0,3.3805d0
     .,3.3218d0,3.2815d0,3.2346d0,3.0994d0,3.3937d0,3.7266d0,3.6697d0
     .,3.6164d0,3.5730d0,3.2522d0,3.5051d0,3.4686d0,3.4355d0,3.4084d0
     .,3.3748d0,3.3496d0,3.3692d0,3.4052d0,3.3910d0,3.3849d0,3.3662d0
     ./)
      r0ab( 631: 700)=(/
     . 2.5087d0,2.4814d0,3.0239d0,3.1312d0,3.0535d0,2.9457d0,2.8496d0
     .,2.7780d0,2.7828d0,2.6532d0,3.1063d0,3.3143d0,3.3549d0,3.3120d0
     .,3.2421d0,3.1787d0,3.1176d0,3.0613d0,3.3082d0,3.5755d0,3.5222d0
     .,3.4678d0,3.4231d0,3.1684d0,3.3528d0,3.3162d0,3.2827d0,3.2527d0
     .,3.2308d0,3.2029d0,3.3173d0,3.3343d0,3.3092d0,3.2795d0,3.2452d0
     .,3.2096d0,3.2893d0,2.8991d0,4.0388d0,3.6100d0,3.9388d0,3.4475d0
     .,3.1590d0,2.9812d0,2.8586d0,2.7683d0,4.1428d0,3.7911d0,3.8225d0
     .,4.0372d0,3.7059d0,3.4935d0,3.3529d0,3.2492d0,4.4352d0,4.0826d0
     .,3.9733d0,3.9254d0,3.8646d0,3.9315d0,3.7837d0,3.7465d0,3.7211d0
     .,3.7012d0,3.6893d0,3.6676d0,3.7736d0,4.0660d0,3.7926d0,3.6158d0
     ./)
      r0ab( 701: 770)=(/
     . 3.5017d0,3.4166d0,4.6176d0,2.8786d0,3.1658d0,3.5823d0,3.7689d0
     .,3.5762d0,3.5789d0,3.3552d0,3.4004d0,3.1722d0,3.0212d0,3.7241d0
     .,3.9604d0,3.8500d0,3.9844d0,3.7035d0,3.9161d0,3.6751d0,3.5075d0
     .,4.1151d0,4.2877d0,4.1579d0,4.1247d0,4.0617d0,3.4874d0,3.9848d0
     .,3.9280d0,3.9079d0,3.8751d0,3.8604d0,3.8277d0,3.8002d0,3.9981d0
     .,3.7544d0,4.0371d0,3.8225d0,3.6718d0,4.3092d0,4.4764d0,2.8997d0
     .,3.0953d0,3.4524d0,3.6107d0,3.6062d0,3.5783d0,3.3463d0,3.3855d0
     .,3.1746d0,3.0381d0,3.6019d0,3.7938d0,3.8697d0,3.9781d0,3.6877d0
     .,3.8736d0,3.6451d0,3.4890d0,3.9858d0,4.1179d0,4.0430d0,3.9563d0
     .,3.9182d0,3.4002d0,3.8310d0,3.7716d0,3.7543d0,3.7203d0,3.7053d0
     ./)
      r0ab( 771: 840)=(/
     . 3.6742d0,3.8318d0,3.7631d0,3.7392d0,3.9892d0,3.7832d0,3.6406d0
     .,4.1701d0,4.3016d0,4.2196d0,2.8535d0,3.0167d0,3.3978d0,3.5363d0
     .,3.5393d0,3.5301d0,3.2960d0,3.3352d0,3.1287d0,2.9967d0,3.6659d0
     .,3.7239d0,3.8070d0,3.7165d0,3.6368d0,3.8162d0,3.5885d0,3.4336d0
     .,3.9829d0,4.0529d0,3.9584d0,3.9025d0,3.8607d0,3.3673d0,3.7658d0
     .,3.7035d0,3.6866d0,3.6504d0,3.6339d0,3.6024d0,3.7708d0,3.7283d0
     .,3.6896d0,3.9315d0,3.7250d0,3.5819d0,4.1457d0,4.2280d0,4.1130d0
     .,4.0597d0,3.0905d0,2.7998d0,3.6448d0,3.0739d0,3.2996d0,3.5262d0
     .,3.2559d0,3.0518d0,2.9394d0,2.8658d0,3.7514d0,3.2295d0,3.5643d0
     .,3.7808d0,3.6931d0,3.4723d0,3.3357d0,3.2429d0,4.0280d0,3.5589d0
     ./)
      r0ab( 841: 910)=(/
     . 3.4636d0,3.4994d0,3.4309d0,3.6177d0,3.2946d0,3.2376d0,3.2050d0
     .,3.1847d0,3.1715d0,3.1599d0,3.5555d0,3.8111d0,3.7693d0,3.5718d0
     .,3.4498d0,3.3662d0,4.1608d0,3.7417d0,3.6536d0,3.6154d0,3.8596d0
     .,3.0301d0,2.7312d0,3.5821d0,3.0473d0,3.2137d0,3.4679d0,3.1975d0
     .,2.9969d0,2.8847d0,2.8110d0,3.6931d0,3.2076d0,3.4943d0,3.5956d0
     .,3.6379d0,3.4190d0,3.2808d0,3.1860d0,3.9850d0,3.5105d0,3.4330d0
     .,3.3797d0,3.4155d0,3.6033d0,3.2737d0,3.2145d0,3.1807d0,3.1596d0
     .,3.1461d0,3.1337d0,3.4812d0,3.6251d0,3.7152d0,3.5201d0,3.3966d0
     .,3.3107d0,4.1128d0,3.6899d0,3.6082d0,3.5604d0,3.7834d0,3.7543d0
     .,2.9189d0,2.6777d0,3.4925d0,2.9648d0,3.1216d0,3.2940d0,3.0975d0
     ./)
      r0ab( 911: 980)=(/
     . 2.9757d0,2.8493d0,2.7638d0,3.6085d0,3.1214d0,3.4006d0,3.4793d0
     .,3.5147d0,3.3806d0,3.2356d0,3.1335d0,3.9144d0,3.4183d0,3.3369d0
     .,3.2803d0,3.2679d0,3.4871d0,3.1714d0,3.1521d0,3.1101d0,3.0843d0
     .,3.0670d0,3.0539d0,3.3890d0,3.5086d0,3.5895d0,3.4783d0,3.3484d0
     .,3.2559d0,4.0422d0,3.5967d0,3.5113d0,3.4576d0,3.6594d0,3.6313d0
     .,3.5690d0,2.8578d0,2.6334d0,3.4673d0,2.9245d0,3.0732d0,3.2435d0
     .,3.0338d0,2.9462d0,2.8143d0,2.7240d0,3.5832d0,3.0789d0,3.3617d0
     .,3.4246d0,3.4505d0,3.3443d0,3.1964d0,3.0913d0,3.8921d0,3.3713d0
     .,3.2873d0,3.2281d0,3.2165d0,3.4386d0,3.1164d0,3.1220d0,3.0761d0
     .,3.0480d0,3.0295d0,3.0155d0,3.3495d0,3.4543d0,3.5260d0,3.4413d0
     ./)
      r0ab( 981:1050)=(/
     . 3.3085d0,3.2134d0,4.0170d0,3.5464d0,3.4587d0,3.4006d0,3.6027d0
     .,3.5730d0,3.4945d0,3.4623d0,2.8240d0,2.5960d0,3.4635d0,2.9032d0
     .,3.0431d0,3.2115d0,2.9892d0,2.9148d0,2.7801d0,2.6873d0,3.5776d0
     .,3.0568d0,3.3433d0,3.3949d0,3.4132d0,3.3116d0,3.1616d0,3.0548d0
     .,3.8859d0,3.3719d0,3.2917d0,3.2345d0,3.2274d0,3.4171d0,3.1293d0
     .,3.0567d0,3.0565d0,3.0274d0,3.0087d0,2.9939d0,3.3293d0,3.4249d0
     .,3.4902d0,3.4091d0,3.2744d0,3.1776d0,4.0078d0,3.5374d0,3.4537d0
     .,3.3956d0,3.5747d0,3.5430d0,3.4522d0,3.4160d0,3.3975d0,2.8004d0
     .,2.5621d0,3.4617d0,2.9154d0,3.0203d0,3.1875d0,2.9548d0,2.8038d0
     .,2.7472d0,2.6530d0,3.5736d0,3.0584d0,3.3304d0,3.3748d0,3.3871d0
     ./)
      r0ab(1051:1120)=(/
     . 3.2028d0,3.1296d0,3.0214d0,3.8796d0,3.3337d0,3.2492d0,3.1883d0
     .,3.1802d0,3.4050d0,3.0756d0,3.0478d0,3.0322d0,3.0323d0,3.0163d0
     .,3.0019d0,3.3145d0,3.4050d0,3.4656d0,3.3021d0,3.2433d0,3.1453d0
     .,3.9991d0,3.5017d0,3.4141d0,3.3520d0,3.5583d0,3.5251d0,3.4243d0
     .,3.3851d0,3.3662d0,3.3525d0,2.7846d0,2.5324d0,3.4652d0,2.8759d0
     .,3.0051d0,3.1692d0,2.9273d0,2.7615d0,2.7164d0,2.6212d0,3.5744d0
     .,3.0275d0,3.3249d0,3.3627d0,3.3686d0,3.1669d0,3.0584d0,2.9915d0
     .,3.8773d0,3.3099d0,3.2231d0,3.1600d0,3.1520d0,3.4023d0,3.0426d0
     .,3.0099d0,2.9920d0,2.9809d0,2.9800d0,2.9646d0,3.3068d0,3.3930d0
     .,3.4486d0,3.2682d0,3.1729d0,3.1168d0,3.9952d0,3.4796d0,3.3901d0
     ./)
      r0ab(1121:1190)=(/
     . 3.3255d0,3.5530d0,3.5183d0,3.4097d0,3.3683d0,3.3492d0,3.3360d0
     .,3.3308d0,2.5424d0,2.6601d0,3.2555d0,3.2807d0,3.1384d0,3.1737d0
     .,2.9397d0,2.8429d0,2.8492d0,2.7225d0,3.3875d0,3.4910d0,3.4520d0
     .,3.3608d0,3.3036d0,3.2345d0,3.2999d0,3.1487d0,3.7409d0,3.8392d0
     .,3.7148d0,3.6439d0,3.6182d0,3.1753d0,3.5210d0,3.4639d0,3.4265d0
     .,3.4075d0,3.3828d0,3.3474d0,3.4071d0,3.3754d0,3.3646d0,3.3308d0
     .,3.4393d0,3.2993d0,3.8768d0,3.9891d0,3.8310d0,3.7483d0,3.3417d0
     .,3.3019d0,3.2250d0,3.1832d0,3.1578d0,3.1564d0,3.1224d0,3.4620d0
     .,2.9743d0,2.8058d0,3.4830d0,3.3474d0,3.6863d0,3.3617d0,3.1608d0
     .,3.0069d0,2.9640d0,2.8427d0,3.5885d0,3.5219d0,4.1314d0,3.8120d0
     ./)
      r0ab(1191:1260)=(/
     . 3.6015d0,3.4502d0,3.3498d0,3.2777d0,3.8635d0,3.8232d0,3.8486d0
     .,3.7215d0,3.6487d0,3.4724d0,3.5627d0,3.5087d0,3.4757d0,3.4517d0
     .,3.4423d0,3.4139d0,4.1028d0,3.8388d0,3.6745d0,3.5562d0,3.4806d0
     .,3.4272d0,4.0182d0,3.9991d0,4.0007d0,3.9282d0,3.7238d0,3.6498d0
     .,3.5605d0,3.5211d0,3.5009d0,3.4859d0,3.4785d0,3.5621d0,4.2623d0
     .,3.0775d0,2.8275d0,4.0181d0,3.3385d0,3.5379d0,3.5036d0,3.2589d0
     .,3.0804d0,3.0094d0,2.9003d0,4.0869d0,3.5088d0,3.9105d0,3.9833d0
     .,3.7176d0,3.5323d0,3.4102d0,3.3227d0,4.2702d0,4.0888d0,3.7560d0
     .,3.7687d0,3.6681d0,3.6405d0,3.5569d0,3.4990d0,3.4659d0,3.4433d0
     .,3.4330d0,3.4092d0,3.8867d0,4.0190d0,3.7961d0,3.6412d0,3.5405d0
     ./)
      r0ab(1261:1330)=(/
     . 3.4681d0,4.3538d0,4.2136d0,3.9381d0,3.8912d0,3.9681d0,3.7909d0
     .,3.6774d0,3.6262d0,3.5999d0,3.5823d0,3.5727d0,3.5419d0,4.0245d0
     .,4.1874d0,3.0893d0,2.7917d0,3.7262d0,3.3518d0,3.4241d0,3.5433d0
     .,3.2773d0,3.0890d0,2.9775d0,2.9010d0,3.8048d0,3.5362d0,3.7746d0
     .,3.7911d0,3.7511d0,3.5495d0,3.4149d0,3.3177d0,4.0129d0,3.8370d0
     .,3.7739d0,3.7125d0,3.7152d0,3.7701d0,3.5813d0,3.5187d0,3.4835d0
     .,3.4595d0,3.4439d0,3.4242d0,3.7476d0,3.8239d0,3.8346d0,3.6627d0
     .,3.5479d0,3.4639d0,4.1026d0,3.9733d0,3.9292d0,3.8667d0,3.9513d0
     .,3.8959d0,3.7698d0,3.7089d0,3.6765d0,3.6548d0,3.6409d0,3.5398d0
     .,3.8759d0,3.9804d0,4.0150d0,2.9091d0,2.7638d0,3.5066d0,3.3377d0
     ./)
      r0ab(1331:1400)=(/
     . 3.3481d0,3.2633d0,3.1810d0,3.1428d0,2.9872d0,2.8837d0,3.5929d0
     .,3.5183d0,3.6729d0,3.6596d0,3.6082d0,3.5927d0,3.4224d0,3.2997d0
     .,3.8190d0,4.1865d0,4.1114d0,4.0540d0,3.6325d0,3.5697d0,3.5561d0
     .,3.5259d0,3.4901d0,3.4552d0,3.4315d0,3.4091d0,3.6438d0,3.6879d0
     .,3.6832d0,3.7043d0,3.5557d0,3.4466d0,3.9203d0,4.2919d0,4.2196d0
     .,4.1542d0,3.7573d0,3.7039d0,3.6546d0,3.6151d0,3.5293d0,3.4849d0
     .,3.4552d0,3.5192d0,3.7673d0,3.8359d0,3.8525d0,3.8901d0,2.7806d0
     .,2.7209d0,3.3812d0,3.4958d0,3.2913d0,3.1888d0,3.0990d0,3.0394d0
     .,2.9789d0,2.8582d0,3.4716d0,3.6883d0,3.6105d0,3.5704d0,3.5059d0
     .,3.4619d0,3.4138d0,3.2742d0,3.7080d0,3.9773d0,3.9010d0,3.8409d0
     ./)
      r0ab(1401:1470)=(/
     . 3.7944d0,3.4465d0,3.7235d0,3.6808d0,3.6453d0,3.6168d0,3.5844d0
     .,3.5576d0,3.5772d0,3.5959d0,3.5768d0,3.5678d0,3.5486d0,3.4228d0
     .,3.8107d0,4.0866d0,4.0169d0,3.9476d0,3.6358d0,3.5800d0,3.5260d0
     .,3.4838d0,3.4501d0,3.4204d0,3.3553d0,3.6487d0,3.6973d0,3.7398d0
     .,3.7405d0,3.7459d0,3.7380d0,2.6848d0,2.6740d0,3.2925d0,3.3386d0
     .,3.2473d0,3.1284d0,3.0301d0,2.9531d0,2.9602d0,2.8272d0,3.3830d0
     .,3.5358d0,3.5672d0,3.5049d0,3.4284d0,3.3621d0,3.3001d0,3.2451d0
     .,3.6209d0,3.8299d0,3.7543d0,3.6920d0,3.6436d0,3.3598d0,3.5701d0
     .,3.5266d0,3.4904d0,3.4590d0,3.4364d0,3.4077d0,3.5287d0,3.5280d0
     .,3.4969d0,3.4650d0,3.4304d0,3.3963d0,3.7229d0,3.9402d0,3.8753d0
     ./)
      r0ab(1471:1540)=(/
     . 3.8035d0,3.5499d0,3.4913d0,3.4319d0,3.3873d0,3.3520d0,3.3209d0
     .,3.2948d0,3.5052d0,3.6465d0,3.6696d0,3.6577d0,3.6388d0,3.6142d0
     .,3.5889d0,3.3968d0,3.0122d0,4.2241d0,3.7887d0,4.0049d0,3.5384d0
     .,3.2698d0,3.1083d0,2.9917d0,2.9057d0,4.3340d0,3.9900d0,4.6588d0
     .,4.1278d0,3.8125d0,3.6189d0,3.4851d0,3.3859d0,4.6531d0,4.3134d0
     .,4.2258d0,4.1309d0,4.0692d0,4.0944d0,3.9850d0,3.9416d0,3.9112d0
     .,3.8873d0,3.8736d0,3.8473d0,4.6027d0,4.1538d0,3.8994d0,3.7419d0
     .,3.6356d0,3.5548d0,4.8353d0,4.5413d0,4.3891d0,4.3416d0,4.3243d0
     .,4.2753d0,4.2053d0,4.1790d0,4.1685d0,4.1585d0,4.1536d0,4.0579d0
     .,4.1980d0,4.4564d0,4.2192d0,4.0528d0,3.9489d0,3.8642d0,5.0567d0
     ./)
      r0ab(1541:1610)=(/
     . 3.0630d0,3.3271d0,4.0432d0,4.0046d0,4.1555d0,3.7426d0,3.5130d0
     .,3.5174d0,3.2884d0,3.1378d0,4.1894d0,4.2321d0,4.1725d0,4.1833d0
     .,3.8929d0,4.0544d0,3.8118d0,3.6414d0,4.6373d0,4.6268d0,4.4750d0
     .,4.4134d0,4.3458d0,3.8582d0,4.2583d0,4.1898d0,4.1562d0,4.1191d0
     .,4.1069d0,4.0639d0,4.1257d0,4.1974d0,3.9532d0,4.1794d0,3.9660d0
     .,3.8130d0,4.8160d0,4.8272d0,4.6294d0,4.5840d0,4.0770d0,4.0088d0
     .,3.9103d0,3.8536d0,3.8324d0,3.7995d0,3.7826d0,4.2294d0,4.3380d0
     .,4.4352d0,4.1933d0,4.4580d0,4.2554d0,4.1072d0,5.0454d0,5.1814d0
     .,3.0632d0,3.2662d0,3.6432d0,3.8088d0,3.7910d0,3.7381d0,3.5093d0
     .,3.5155d0,3.3047d0,3.1681d0,3.7871d0,3.9924d0,4.0637d0,4.1382d0
     ./)
      r0ab(1611:1680)=(/
     . 3.8591d0,4.0164d0,3.7878d0,3.6316d0,4.1741d0,4.3166d0,4.2395d0
     .,4.1831d0,4.1107d0,3.5857d0,4.0270d0,3.9676d0,3.9463d0,3.9150d0
     .,3.9021d0,3.8708d0,4.0240d0,4.1551d0,3.9108d0,4.1337d0,3.9289d0
     .,3.7873d0,4.3666d0,4.5080d0,4.4232d0,4.3155d0,3.8461d0,3.8007d0
     .,3.6991d0,3.6447d0,3.6308d0,3.5959d0,3.5749d0,4.0359d0,4.3124d0
     .,4.3539d0,4.1122d0,4.3772d0,4.1785d0,4.0386d0,4.7004d0,4.8604d0
     .,4.6261d0,2.9455d0,3.2470d0,3.6108d0,3.8522d0,3.6625d0,3.6598d0
     .,3.4411d0,3.4660d0,3.2415d0,3.0944d0,3.7514d0,4.0397d0,3.9231d0
     .,4.0561d0,3.7860d0,3.9845d0,3.7454d0,3.5802d0,4.1366d0,4.3581d0
     .,4.2351d0,4.2011d0,4.1402d0,3.5381d0,4.0653d0,4.0093d0,3.9883d0
     ./)
      r0ab(1681:1750)=(/
     . 3.9570d0,3.9429d0,3.9112d0,3.8728d0,4.0682d0,3.8351d0,4.1054d0
     .,3.8928d0,3.7445d0,4.3415d0,4.5497d0,4.3833d0,4.3122d0,3.8051d0
     .,3.7583d0,3.6622d0,3.6108d0,3.5971d0,3.5628d0,3.5408d0,4.0780d0
     .,4.0727d0,4.2836d0,4.0553d0,4.3647d0,4.1622d0,4.0178d0,4.5802d0
     .,4.9125d0,4.5861d0,4.6201d0,2.9244d0,3.2241d0,3.5848d0,3.8293d0
     .,3.6395d0,3.6400d0,3.4204d0,3.4499d0,3.2253d0,3.0779d0,3.7257d0
     .,4.0170d0,3.9003d0,4.0372d0,3.7653d0,3.9672d0,3.7283d0,3.5630d0
     .,4.1092d0,4.3347d0,4.2117d0,4.1793d0,4.1179d0,3.5139d0,4.0426d0
     .,3.9867d0,3.9661d0,3.9345d0,3.9200d0,3.8883d0,3.8498d0,4.0496d0
     .,3.8145d0,4.0881d0,3.8756d0,3.7271d0,4.3128d0,4.5242d0,4.3578d0
     ./)
      r0ab(1751:1820)=(/
     . 4.2870d0,3.7796d0,3.7318d0,3.6364d0,3.5854d0,3.5726d0,3.5378d0
     .,3.5155d0,4.0527d0,4.0478d0,4.2630d0,4.0322d0,4.3449d0,4.1421d0
     .,3.9975d0,4.5499d0,4.8825d0,4.5601d0,4.5950d0,4.5702d0,2.9046d0
     .,3.2044d0,3.5621d0,3.8078d0,3.6185d0,3.6220d0,3.4019d0,3.4359d0
     .,3.2110d0,3.0635d0,3.7037d0,3.9958d0,3.8792d0,4.0194d0,3.7460d0
     .,3.9517d0,3.7128d0,3.5474d0,4.0872d0,4.3138d0,4.1906d0,4.1593d0
     .,4.0973d0,3.4919d0,4.0216d0,3.9657d0,3.9454d0,3.9134d0,3.8986d0
     .,3.8669d0,3.8289d0,4.0323d0,3.7954d0,4.0725d0,3.8598d0,3.7113d0
     .,4.2896d0,4.5021d0,4.3325d0,4.2645d0,3.7571d0,3.7083d0,3.6136d0
     .,3.5628d0,3.5507d0,3.5155d0,3.4929d0,4.0297d0,4.0234d0,4.2442d0
     ./)
      r0ab(1821:1890)=(/
     . 4.0112d0,4.3274d0,4.1240d0,3.9793d0,4.5257d0,4.8568d0,4.5353d0
     .,4.5733d0,4.5485d0,4.5271d0,2.8878d0,3.1890d0,3.5412d0,3.7908d0
     .,3.5974d0,3.6078d0,3.3871d0,3.4243d0,3.1992d0,3.0513d0,3.6831d0
     .,3.9784d0,3.8579d0,4.0049d0,3.7304d0,3.9392d0,3.7002d0,3.5347d0
     .,4.0657d0,4.2955d0,4.1705d0,4.1424d0,4.0800d0,3.4717d0,4.0043d0
     .,3.9485d0,3.9286d0,3.8965d0,3.8815d0,3.8500d0,3.8073d0,4.0180d0
     .,3.7796d0,4.0598d0,3.8470d0,3.6983d0,4.2678d0,4.4830d0,4.3132d0
     .,4.2444d0,3.7370d0,3.6876d0,3.5935d0,3.5428d0,3.5314d0,3.4958d0
     .,3.4730d0,4.0117d0,4.0043d0,4.2287d0,3.9939d0,4.3134d0,4.1096d0
     .,3.9646d0,4.5032d0,4.8356d0,4.5156d0,4.5544d0,4.5297d0,4.5083d0
     ./)
      r0ab(1891:1960)=(/
     . 4.4896d0,2.8709d0,3.1737d0,3.5199d0,3.7734d0,3.5802d0,3.5934d0
     .,3.3724d0,3.4128d0,3.1877d0,3.0396d0,3.6624d0,3.9608d0,3.8397d0
     .,3.9893d0,3.7145d0,3.9266d0,3.6877d0,3.5222d0,4.0448d0,4.2771d0
     .,4.1523d0,4.1247d0,4.0626d0,3.4530d0,3.9866d0,3.9310d0,3.9115d0
     .,3.8792d0,3.8641d0,3.8326d0,3.7892d0,4.0025d0,3.7636d0,4.0471d0
     .,3.8343d0,3.6854d0,4.2464d0,4.4635d0,4.2939d0,4.2252d0,3.7169d0
     .,3.6675d0,3.5739d0,3.5235d0,3.5126d0,3.4768d0,3.4537d0,3.9932d0
     .,3.9854d0,4.2123d0,3.9765d0,4.2992d0,4.0951d0,3.9500d0,4.4811d0
     .,4.8135d0,4.4959d0,4.5351d0,4.5105d0,4.4891d0,4.4705d0,4.4515d0
     .,2.8568d0,3.1608d0,3.5050d0,3.7598d0,3.5665d0,3.5803d0,3.3601d0
     ./)
      r0ab(1961:2030)=(/
     . 3.4031d0,3.1779d0,3.0296d0,3.6479d0,3.9471d0,3.8262d0,3.9773d0
     .,3.7015d0,3.9162d0,3.6771d0,3.5115d0,4.0306d0,4.2634d0,4.1385d0
     .,4.1116d0,4.0489d0,3.4366d0,3.9732d0,3.9176d0,3.8983d0,3.8659d0
     .,3.8507d0,3.8191d0,3.7757d0,3.9907d0,3.7506d0,4.0365d0,3.8235d0
     .,3.6745d0,4.2314d0,4.4490d0,4.2792d0,4.2105d0,3.7003d0,3.6510d0
     .,3.5578d0,3.5075d0,3.4971d0,3.4609d0,3.4377d0,3.9788d0,3.9712d0
     .,4.1997d0,3.9624d0,4.2877d0,4.0831d0,3.9378d0,4.4655d0,4.7974d0
     .,4.4813d0,4.5209d0,4.4964d0,4.4750d0,4.4565d0,4.4375d0,4.4234d0
     .,2.6798d0,3.0151d0,3.2586d0,3.5292d0,3.5391d0,3.4902d0,3.2887d0
     .,3.3322d0,3.1228d0,2.9888d0,3.4012d0,3.7145d0,3.7830d0,3.6665d0
     ./)
      r0ab(2031:2100)=(/
     . 3.5898d0,3.8077d0,3.5810d0,3.4265d0,3.7726d0,4.0307d0,3.9763d0
     .,3.8890d0,3.8489d0,3.2706d0,3.7595d0,3.6984d0,3.6772d0,3.6428d0
     .,3.6243d0,3.5951d0,3.7497d0,3.6775d0,3.6364d0,3.9203d0,3.7157d0
     .,3.5746d0,3.9494d0,4.2076d0,4.1563d0,4.0508d0,3.5329d0,3.4780d0
     .,3.3731d0,3.3126d0,3.2846d0,3.2426d0,3.2135d0,3.7491d0,3.9006d0
     .,3.8332d0,3.8029d0,4.1436d0,3.9407d0,3.7998d0,4.1663d0,4.5309d0
     .,4.3481d0,4.2911d0,4.2671d0,4.2415d0,4.2230d0,4.2047d0,4.1908d0
     .,4.1243d0,2.5189d0,2.9703d0,3.3063d0,3.6235d0,3.4517d0,3.3989d0
     .,3.2107d0,3.2434d0,3.0094d0,2.8580d0,3.4253d0,3.8157d0,3.7258d0
     .,3.6132d0,3.5297d0,3.7566d0,3.5095d0,3.3368d0,3.7890d0,4.1298d0
     ./)
      r0ab(2101:2170)=(/
     . 4.0190d0,3.9573d0,3.9237d0,3.2677d0,3.8480d0,3.8157d0,3.7656d0
     .,3.7317d0,3.7126d0,3.6814d0,3.6793d0,3.6218d0,3.5788d0,3.8763d0
     .,3.6572d0,3.5022d0,3.9737d0,4.3255d0,4.1828d0,4.1158d0,3.5078d0
     .,3.4595d0,3.3600d0,3.3088d0,3.2575d0,3.2164d0,3.1856d0,3.8522d0
     .,3.8665d0,3.8075d0,3.7772d0,4.1391d0,3.9296d0,3.7772d0,4.2134d0
     .,4.7308d0,4.3787d0,4.3894d0,4.3649d0,4.3441d0,4.3257d0,4.3073d0
     .,4.2941d0,4.1252d0,4.2427d0,3.0481d0,2.9584d0,3.6919d0,3.5990d0
     .,3.8881d0,3.4209d0,3.1606d0,3.1938d0,2.9975d0,2.8646d0,3.8138d0
     .,3.7935d0,3.7081d0,3.9155d0,3.5910d0,3.4808d0,3.4886d0,3.3397d0
     .,4.1336d0,4.1122d0,3.9888d0,3.9543d0,3.8917d0,3.5894d0,3.8131d0
     ./)
      r0ab(2171:2240)=(/
     . 3.7635d0,3.7419d0,3.7071d0,3.6880d0,3.6574d0,3.6546d0,3.9375d0
     .,3.6579d0,3.5870d0,3.6361d0,3.5039d0,4.3149d0,4.2978d0,4.1321d0
     .,4.1298d0,3.8164d0,3.7680d0,3.7154d0,3.6858d0,3.6709d0,3.6666d0
     .,3.6517d0,3.8174d0,3.8608d0,4.1805d0,3.9102d0,3.8394d0,3.8968d0
     .,3.7673d0,4.5274d0,4.6682d0,4.3344d0,4.3639d0,4.3384d0,4.3162d0
     .,4.2972d0,4.2779d0,4.2636d0,4.0253d0,4.1168d0,4.1541d0,2.8136d0
     .,3.0951d0,3.4635d0,3.6875d0,3.4987d0,3.5183d0,3.2937d0,3.3580d0
     .,3.1325d0,2.9832d0,3.6078d0,3.8757d0,3.7616d0,3.9222d0,3.6370d0
     .,3.8647d0,3.6256d0,3.4595d0,3.9874d0,4.1938d0,4.0679d0,4.0430d0
     .,3.9781d0,3.3886d0,3.9008d0,3.8463d0,3.8288d0,3.7950d0,3.7790d0
     ./)
      r0ab(2241:2310)=(/
     . 3.7472d0,3.7117d0,3.9371d0,3.6873d0,3.9846d0,3.7709d0,3.6210d0
     .,4.1812d0,4.3750d0,4.2044d0,4.1340d0,3.6459d0,3.5929d0,3.5036d0
     .,3.4577d0,3.4528d0,3.4146d0,3.3904d0,3.9014d0,3.9031d0,4.1443d0
     .,3.8961d0,4.2295d0,4.0227d0,3.8763d0,4.4086d0,4.7097d0,4.4064d0
     .,4.4488d0,4.4243d0,4.4029d0,4.3842d0,4.3655d0,4.3514d0,4.1162d0
     .,4.2205d0,4.1953d0,4.2794d0,2.8032d0,3.0805d0,3.4519d0,3.6700d0
     .,3.4827d0,3.5050d0,3.2799d0,3.3482d0,3.1233d0,2.9747d0,3.5971d0
     .,3.8586d0,3.7461d0,3.9100d0,3.6228d0,3.8535d0,3.6147d0,3.4490d0
     .,3.9764d0,4.1773d0,4.0511d0,4.0270d0,3.9614d0,3.3754d0,3.8836d0
     .,3.8291d0,3.8121d0,3.7780d0,3.7619d0,3.7300d0,3.6965d0,3.9253d0
     ./)
      r0ab(2311:2380)=(/
     . 3.6734d0,3.9733d0,3.7597d0,3.6099d0,4.1683d0,4.3572d0,4.1862d0
     .,4.1153d0,3.6312d0,3.5772d0,3.4881d0,3.4429d0,3.4395d0,3.4009d0
     .,3.3766d0,3.8827d0,3.8868d0,4.1316d0,3.8807d0,4.2164d0,4.0092d0
     .,3.8627d0,4.3936d0,4.6871d0,4.3882d0,4.4316d0,4.4073d0,4.3858d0
     .,4.3672d0,4.3485d0,4.3344d0,4.0984d0,4.2036d0,4.1791d0,4.2622d0
     .,4.2450d0,2.7967d0,3.0689d0,3.4445d0,3.6581d0,3.4717d0,3.4951d0
     .,3.2694d0,3.3397d0,3.1147d0,2.9661d0,3.5898d0,3.8468d0,3.7358d0
     .,3.9014d0,3.6129d0,3.8443d0,3.6054d0,3.4396d0,3.9683d0,4.1656d0
     .,4.0394d0,4.0158d0,3.9498d0,3.3677d0,3.8718d0,3.8164d0,3.8005d0
     .,3.7662d0,3.7500d0,3.7181d0,3.6863d0,3.9170d0,3.6637d0,3.9641d0
     ./)
      r0ab(2381:2450)=(/
     . 3.7503d0,3.6004d0,4.1590d0,4.3448d0,4.1739d0,4.1029d0,3.6224d0
     .,3.5677d0,3.4785d0,3.4314d0,3.4313d0,3.3923d0,3.3680d0,3.8698d0
     .,3.8758d0,4.1229d0,3.8704d0,4.2063d0,3.9987d0,3.8519d0,4.3832d0
     .,4.6728d0,4.3759d0,4.4195d0,4.3952d0,4.3737d0,4.3551d0,4.3364d0
     .,4.3223d0,4.0861d0,4.1911d0,4.1676d0,4.2501d0,4.2329d0,4.2208d0
     .,2.7897d0,3.0636d0,3.4344d0,3.6480d0,3.4626d0,3.4892d0,3.2626d0
     .,3.3344d0,3.1088d0,2.9597d0,3.5804d0,3.8359d0,3.7251d0,3.8940d0
     .,3.6047d0,3.8375d0,3.5990d0,3.4329d0,3.9597d0,4.1542d0,4.0278d0
     .,4.0048d0,3.9390d0,3.3571d0,3.8608d0,3.8056d0,3.7899d0,3.7560d0
     .,3.7400d0,3.7081d0,3.6758d0,3.9095d0,3.6552d0,3.9572d0,3.7436d0
     ./)
      r0ab(2451:2520)=(/
     . 3.5933d0,4.1508d0,4.3337d0,4.1624d0,4.0916d0,3.6126d0,3.5582d0
     .,3.4684d0,3.4212d0,3.4207d0,3.3829d0,3.3586d0,3.8604d0,3.8658d0
     .,4.1156d0,3.8620d0,4.1994d0,3.9917d0,3.8446d0,4.3750d0,4.6617d0
     .,4.3644d0,4.4083d0,4.3840d0,4.3625d0,4.3439d0,4.3253d0,4.3112d0
     .,4.0745d0,4.1807d0,4.1578d0,4.2390d0,4.2218d0,4.2097d0,4.1986d0
     .,2.8395d0,3.0081d0,3.3171d0,3.4878d0,3.5360d0,3.5145d0,3.2809d0
     .,3.3307d0,3.1260d0,2.9940d0,3.4741d0,3.6675d0,3.7832d0,3.6787d0
     .,3.6156d0,3.8041d0,3.5813d0,3.4301d0,3.8480d0,3.9849d0,3.9314d0
     .,3.8405d0,3.8029d0,3.2962d0,3.7104d0,3.6515d0,3.6378d0,3.6020d0
     .,3.5849d0,3.5550d0,3.7494d0,3.6893d0,3.6666d0,3.9170d0,3.7150d0
     ./)
      r0ab(2521:2590)=(/
     . 3.5760d0,4.0268d0,4.1596d0,4.1107d0,3.9995d0,3.5574d0,3.5103d0
     .,3.4163d0,3.3655d0,3.3677d0,3.3243d0,3.2975d0,3.7071d0,3.9047d0
     .,3.8514d0,3.8422d0,3.8022d0,3.9323d0,3.7932d0,4.2343d0,4.4583d0
     .,4.3115d0,4.2457d0,4.2213d0,4.1945d0,4.1756d0,4.1569d0,4.1424d0
     .,4.0620d0,4.0494d0,3.9953d0,4.0694d0,4.0516d0,4.0396d0,4.0280d0
     .,4.0130d0,2.9007d0,2.9674d0,3.8174d0,3.5856d0,3.6486d0,3.5339d0
     .,3.2832d0,3.3154d0,3.1144d0,2.9866d0,3.9618d0,3.8430d0,3.9980d0
     .,3.8134d0,3.6652d0,3.7985d0,3.5756d0,3.4207d0,4.4061d0,4.2817d0
     .,4.1477d0,4.0616d0,3.9979d0,3.6492d0,3.8833d0,3.8027d0,3.7660d0
     .,3.7183d0,3.6954d0,3.6525d0,3.9669d0,3.8371d0,3.7325d0,3.9160d0
     ./)
      r0ab(2591:2660)=(/
     . 3.7156d0,3.5714d0,4.6036d0,4.4620d0,4.3092d0,4.2122d0,3.8478d0
     .,3.7572d0,3.6597d0,3.5969d0,3.5575d0,3.5386d0,3.5153d0,3.7818d0
     .,4.1335d0,4.0153d0,3.9177d0,3.8603d0,3.9365d0,3.7906d0,4.7936d0
     .,4.7410d0,4.5461d0,4.5662d0,4.5340d0,4.5059d0,4.4832d0,4.4604d0
     .,4.4429d0,4.2346d0,4.4204d0,4.3119d0,4.3450d0,4.3193d0,4.3035d0
     .,4.2933d0,4.1582d0,4.2450d0,2.8559d0,2.9050d0,3.8325d0,3.5442d0
     .,3.5077d0,3.4905d0,3.2396d0,3.2720d0,3.0726d0,2.9467d0,3.9644d0
     .,3.8050d0,3.8981d0,3.7762d0,3.6216d0,3.7531d0,3.5297d0,3.3742d0
     .,4.3814d0,4.2818d0,4.1026d0,4.0294d0,3.9640d0,3.6208d0,3.8464d0
     .,3.7648d0,3.7281d0,3.6790d0,3.6542d0,3.6117d0,3.8650d0,3.8010d0
     ./)
      r0ab(2661:2730)=(/
     . 3.6894d0,3.8713d0,3.6699d0,3.5244d0,4.5151d0,4.4517d0,4.2538d0
     .,4.1483d0,3.8641d0,3.7244d0,3.6243d0,3.5589d0,3.5172d0,3.4973d0
     .,3.4715d0,3.7340d0,4.0316d0,3.9958d0,3.8687d0,3.8115d0,3.8862d0
     .,3.7379d0,4.7091d0,4.7156d0,4.5199d0,4.5542d0,4.5230d0,4.4959d0
     .,4.4750d0,4.4529d0,4.4361d0,4.1774d0,4.3774d0,4.2963d0,4.3406d0
     .,4.3159d0,4.3006d0,4.2910d0,4.1008d0,4.1568d0,4.0980d0,2.8110d0
     .,2.8520d0,3.7480d0,3.5105d0,3.4346d0,3.3461d0,3.1971d0,3.2326d0
     .,3.0329d0,2.9070d0,3.8823d0,3.7928d0,3.8264d0,3.7006d0,3.5797d0
     .,3.7141d0,3.4894d0,3.3326d0,4.3048d0,4.2217d0,4.0786d0,3.9900d0
     .,3.9357d0,3.6331d0,3.8333d0,3.7317d0,3.6957d0,3.6460d0,3.6197d0
     ./)
      r0ab(2731:2800)=(/
     . 3.5779d0,3.7909d0,3.7257d0,3.6476d0,3.5729d0,3.6304d0,3.4834d0
     .,4.4368d0,4.3921d0,4.2207d0,4.1133d0,3.8067d0,3.7421d0,3.6140d0
     .,3.5491d0,3.5077d0,3.4887d0,3.4623d0,3.6956d0,3.9568d0,3.8976d0
     .,3.8240d0,3.7684d0,3.8451d0,3.6949d0,4.6318d0,4.6559d0,4.4533d0
     .,4.4956d0,4.4641d0,4.4366d0,4.4155d0,4.3936d0,4.3764d0,4.1302d0
     .,4.3398d0,4.2283d0,4.2796d0,4.2547d0,4.2391d0,4.2296d0,4.0699d0
     .,4.1083d0,4.0319d0,3.9855d0,2.7676d0,2.8078d0,3.6725d0,3.4804d0
     .,3.3775d0,3.2411d0,3.1581d0,3.1983d0,2.9973d0,2.8705d0,3.8070d0
     .,3.7392d0,3.7668d0,3.6263d0,3.5402d0,3.6807d0,3.4545d0,3.2962d0
     .,4.2283d0,4.1698d0,4.0240d0,3.9341d0,3.8711d0,3.5489d0,3.7798d0
     ./)
      r0ab(2801:2870)=(/
     . 3.7000d0,3.6654d0,3.6154d0,3.5882d0,3.5472d0,3.7289d0,3.6510d0
     .,3.6078d0,3.5355d0,3.5963d0,3.4480d0,4.3587d0,4.3390d0,4.1635d0
     .,4.0536d0,3.7193d0,3.6529d0,3.5512d0,3.4837d0,3.4400d0,3.4191d0
     .,3.3891d0,3.6622d0,3.8934d0,3.8235d0,3.7823d0,3.7292d0,3.8106d0
     .,3.6589d0,4.5535d0,4.6013d0,4.3961d0,4.4423d0,4.4109d0,4.3835d0
     .,4.3625d0,4.3407d0,4.3237d0,4.0863d0,4.2835d0,4.1675d0,4.2272d0
     .,4.2025d0,4.1869d0,4.1774d0,4.0126d0,4.0460d0,3.9815d0,3.9340d0
     .,3.8955d0,2.6912d0,2.7604d0,3.6037d0,3.4194d0,3.3094d0,3.1710d0
     .,3.0862d0,3.1789d0,2.9738d0,2.8427d0,3.7378d0,3.6742d0,3.6928d0
     .,3.5512d0,3.4614d0,3.4087d0,3.4201d0,3.2607d0,4.1527d0,4.0977d0
     ./)
      r0ab(2871:2940)=(/
     . 3.9523d0,3.8628d0,3.8002d0,3.4759d0,3.7102d0,3.6466d0,3.6106d0
     .,3.5580d0,3.5282d0,3.4878d0,3.6547d0,3.5763d0,3.5289d0,3.5086d0
     .,3.5593d0,3.4099d0,4.2788d0,4.2624d0,4.0873d0,3.9770d0,3.6407d0
     .,3.5743d0,3.5178d0,3.4753d0,3.3931d0,3.3694d0,3.3339d0,3.6002d0
     .,3.8164d0,3.7478d0,3.7028d0,3.6952d0,3.7669d0,3.6137d0,4.4698d0
     .,4.5488d0,4.3168d0,4.3646d0,4.3338d0,4.3067d0,4.2860d0,4.2645d0
     .,4.2478d0,4.0067d0,4.2349d0,4.0958d0,4.1543d0,4.1302d0,4.1141d0
     .,4.1048d0,3.9410d0,3.9595d0,3.8941d0,3.8465d0,3.8089d0,3.7490d0
     .,2.7895d0,2.5849d0,3.6484d0,3.0162d0,3.1267d0,3.2125d0,3.0043d0
     .,2.9572d0,2.8197d0,2.7261d0,3.7701d0,3.2446d0,3.5239d0,3.4696d0
     ./)
      r0ab(2941:3010)=(/
     . 3.4261d0,3.3508d0,3.1968d0,3.0848d0,4.1496d0,3.6598d0,3.5111d0
     .,3.4199d0,3.3809d0,3.5382d0,3.2572d0,3.2100d0,3.1917d0,3.1519d0
     .,3.1198d0,3.1005d0,3.5071d0,3.5086d0,3.5073d0,3.4509d0,3.3120d0
     .,3.2082d0,4.2611d0,3.8117d0,3.6988d0,3.5646d0,3.6925d0,3.6295d0
     .,3.5383d0,3.4910d0,3.4625d0,3.4233d0,3.4007d0,3.2329d0,3.6723d0
     .,3.6845d0,3.6876d0,3.6197d0,3.4799d0,3.3737d0,4.4341d0,4.0525d0
     .,3.9011d0,3.8945d0,3.8635d0,3.8368d0,3.8153d0,3.7936d0,3.7758d0
     .,3.4944d0,3.4873d0,3.9040d0,3.7110d0,3.6922d0,3.6799d0,3.6724d0
     .,3.5622d0,3.6081d0,3.5426d0,3.4922d0,3.4498d0,3.3984d0,3.4456d0
     .,2.7522d0,2.5524d0,3.5742d0,2.9508d0,3.0751d0,3.0158d0,2.9644d0
     ./)
      r0ab(3011:3080)=(/
     . 2.8338d0,2.7891d0,2.6933d0,3.6926d0,3.1814d0,3.4528d0,3.4186d0
     .,3.3836d0,3.2213d0,3.1626d0,3.0507d0,4.0548d0,3.5312d0,3.4244d0
     .,3.3409d0,3.2810d0,3.4782d0,3.1905d0,3.1494d0,3.1221d0,3.1128d0
     .,3.0853d0,3.0384d0,3.4366d0,3.4562d0,3.4638d0,3.3211d0,3.2762d0
     .,3.1730d0,4.1632d0,3.6825d0,3.5822d0,3.4870d0,3.6325d0,3.5740d0
     .,3.4733d0,3.4247d0,3.3969d0,3.3764d0,3.3525d0,3.1984d0,3.5989d0
     .,3.6299d0,3.6433d0,3.4937d0,3.4417d0,3.3365d0,4.3304d0,3.9242d0
     .,3.7793d0,3.7623d0,3.7327d0,3.7071d0,3.6860d0,3.6650d0,3.6476d0
     .,3.3849d0,3.3534d0,3.8216d0,3.5870d0,3.5695d0,3.5584d0,3.5508d0
     .,3.4856d0,3.5523d0,3.4934d0,3.4464d0,3.4055d0,3.3551d0,3.3888d0
     ./)
      r0ab(3081:3150)=(/
     . 3.3525d0,2.7202d0,2.5183d0,3.4947d0,2.8731d0,3.0198d0,3.1457d0
     .,2.9276d0,2.7826d0,2.7574d0,2.6606d0,3.6090d0,3.0581d0,3.3747d0
     .,3.3677d0,3.3450d0,3.1651d0,3.1259d0,3.0147d0,3.9498d0,3.3857d0
     .,3.2917d0,3.2154d0,3.1604d0,3.4174d0,3.0735d0,3.0342d0,3.0096d0
     .,3.0136d0,2.9855d0,2.9680d0,3.3604d0,3.4037d0,3.4243d0,3.2633d0
     .,3.1810d0,3.1351d0,4.0557d0,3.5368d0,3.4526d0,3.3699d0,3.5707d0
     .,3.5184d0,3.4085d0,3.3595d0,3.3333d0,3.3143d0,3.3041d0,3.1094d0
     .,3.5193d0,3.5745d0,3.6025d0,3.4338d0,3.3448d0,3.2952d0,4.2158d0
     .,3.7802d0,3.6431d0,3.6129d0,3.5853d0,3.5610d0,3.5406d0,3.5204d0
     .,3.5036d0,3.2679d0,3.2162d0,3.7068d0,3.4483d0,3.4323d0,3.4221d0
     ./)
      r0ab(3151:3220)=(/
     . 3.4138d0,3.3652d0,3.4576d0,3.4053d0,3.3618d0,3.3224d0,3.2711d0
     .,3.3326d0,3.2950d0,3.2564d0,2.5315d0,2.6104d0,3.2734d0,3.2299d0
     .,3.1090d0,2.9942d0,2.9159d0,2.8324d0,2.8350d0,2.7216d0,3.3994d0
     .,3.4475d0,3.4354d0,3.3438d0,3.2807d0,3.2169d0,3.2677d0,3.1296d0
     .,3.7493d0,3.8075d0,3.6846d0,3.6104d0,3.5577d0,3.2052d0,3.4803d0
     .,3.4236d0,3.3845d0,3.3640d0,3.3365d0,3.3010d0,3.3938d0,3.3624d0
     .,3.3440d0,3.3132d0,3.4035d0,3.2754d0,3.8701d0,3.9523d0,3.8018d0
     .,3.7149d0,3.3673d0,3.3199d0,3.2483d0,3.2069d0,3.1793d0,3.1558d0
     .,3.1395d0,3.4097d0,3.5410d0,3.5228d0,3.5116d0,3.4921d0,3.4781d0
     .,3.4690d0,4.0420d0,4.1759d0,4.0078d0,4.0450d0,4.0189d0,3.9952d0
     ./)
      r0ab(3221:3290)=(/
     . 3.9770d0,3.9583d0,3.9434d0,3.7217d0,3.8228d0,3.7826d0,3.8640d0
     .,3.8446d0,3.8314d0,3.8225d0,3.6817d0,3.7068d0,3.6555d0,3.6159d0
     .,3.5831d0,3.5257d0,3.2133d0,3.1689d0,3.1196d0,3.3599d0,2.9852d0
     .,2.7881d0,3.5284d0,3.3493d0,3.6958d0,3.3642d0,3.1568d0,3.0055d0
     .,2.9558d0,2.8393d0,3.6287d0,3.5283d0,4.1511d0,3.8259d0,3.6066d0
     .,3.4527d0,3.3480d0,3.2713d0,3.9037d0,3.8361d0,3.8579d0,3.7311d0
     .,3.6575d0,3.5176d0,3.5693d0,3.5157d0,3.4814d0,3.4559d0,3.4445d0
     .,3.4160d0,4.1231d0,3.8543d0,3.6816d0,3.5602d0,3.4798d0,3.4208d0
     .,4.0542d0,4.0139d0,4.0165d0,3.9412d0,3.7698d0,3.6915d0,3.6043d0
     .,3.5639d0,3.5416d0,3.5247d0,3.5153d0,3.5654d0,4.2862d0,4.0437d0
     ./)
      r0ab(3291:3360)=(/
     . 3.8871d0,3.7741d0,3.6985d0,3.6413d0,4.2345d0,4.3663d0,4.3257d0
     .,4.0869d0,4.0612d0,4.0364d0,4.0170d0,3.9978d0,3.9834d0,3.9137d0
     .,3.8825d0,3.8758d0,3.9143d0,3.8976d0,3.8864d0,3.8768d0,3.9190d0
     .,4.1613d0,4.0566d0,3.9784d0,3.9116d0,3.8326d0,3.7122d0,3.6378d0
     .,3.5576d0,3.5457d0,4.3127d0,3.1160d0,2.8482d0,4.0739d0,3.3599d0
     .,3.5698d0,3.5366d0,3.2854d0,3.1039d0,2.9953d0,2.9192d0,4.1432d0
     .,3.5320d0,3.9478d0,4.0231d0,3.7509d0,3.5604d0,3.4340d0,3.3426d0
     .,4.3328d0,3.8288d0,3.7822d0,3.7909d0,3.6907d0,3.6864d0,3.5793d0
     .,3.5221d0,3.4883d0,3.4649d0,3.4514d0,3.4301d0,3.9256d0,4.0596d0
     .,3.8307d0,3.6702d0,3.5651d0,3.4884d0,4.4182d0,4.2516d0,3.9687d0
     ./)
      r0ab(3361:3430)=(/
     . 3.9186d0,3.9485d0,3.8370d0,3.7255d0,3.6744d0,3.6476d0,3.6295d0
     .,3.6193d0,3.5659d0,4.0663d0,4.2309d0,4.0183d0,3.8680d0,3.7672d0
     .,3.6923d0,4.5240d0,4.4834d0,4.1570d0,4.3204d0,4.2993d0,4.2804d0
     .,4.2647d0,4.2481d0,4.2354d0,3.8626d0,3.8448d0,4.2267d0,4.1799d0
     .,4.1670d0,3.8738d0,3.8643d0,3.8796d0,4.0575d0,4.0354d0,3.9365d0
     .,3.8611d0,3.7847d0,3.7388d0,3.6826d0,3.6251d0,3.5492d0,4.0889d0
     .,4.2764d0,3.1416d0,2.8325d0,3.7735d0,3.3787d0,3.4632d0,3.5923d0
     .,3.3214d0,3.1285d0,3.0147d0,2.9366d0,3.8527d0,3.5602d0,3.8131d0
     .,3.8349d0,3.7995d0,3.5919d0,3.4539d0,3.3540d0,4.0654d0,3.8603d0
     .,3.7972d0,3.7358d0,3.7392d0,3.8157d0,3.6055d0,3.5438d0,3.5089d0
     ./)
      r0ab(3431:3500)=(/
     . 3.4853d0,3.4698d0,3.4508d0,3.7882d0,3.8682d0,3.8837d0,3.7055d0
     .,3.5870d0,3.5000d0,4.1573d0,4.0005d0,3.9568d0,3.8936d0,3.9990d0
     .,3.9433d0,3.8172d0,3.7566d0,3.7246d0,3.7033d0,3.6900d0,3.5697d0
     .,3.9183d0,4.0262d0,4.0659d0,3.8969d0,3.7809d0,3.6949d0,4.2765d0
     .,4.2312d0,4.1401d0,4.0815d0,4.0580d0,4.0369d0,4.0194d0,4.0017d0
     .,3.9874d0,3.8312d0,3.8120d0,3.9454d0,3.9210d0,3.9055d0,3.8951d0
     .,3.8866d0,3.8689d0,3.9603d0,3.9109d0,3.9122d0,3.8233d0,3.7438d0
     .,3.7436d0,3.6981d0,3.6555d0,3.5452d0,3.9327d0,4.0658d0,4.1175d0
     .,2.9664d0,2.8209d0,3.5547d0,3.3796d0,3.3985d0,3.3164d0,3.2364d0
     .,3.1956d0,3.0370d0,2.9313d0,3.6425d0,3.5565d0,3.7209d0,3.7108d0
     ./)
      r0ab(3501:3570)=(/
     . 3.6639d0,3.6484d0,3.4745d0,3.3492d0,3.8755d0,4.2457d0,3.7758d0
     .,3.7161d0,3.6693d0,3.6155d0,3.5941d0,3.5643d0,3.5292d0,3.4950d0
     .,3.4720d0,3.4503d0,3.6936d0,3.7392d0,3.7388d0,3.7602d0,3.6078d0
     .,3.4960d0,3.9800d0,4.3518d0,4.2802d0,3.8580d0,3.8056d0,3.7527d0
     .,3.7019d0,3.6615d0,3.5768d0,3.5330d0,3.5038d0,3.5639d0,3.8192d0
     .,3.8883d0,3.9092d0,3.9478d0,3.7995d0,3.6896d0,4.1165d0,4.5232d0
     .,4.4357d0,4.4226d0,4.4031d0,4.3860d0,4.3721d0,4.3580d0,4.3466d0
     .,4.2036d0,4.2037d0,3.8867d0,4.2895d0,4.2766d0,4.2662d0,4.2598d0
     .,3.8408d0,3.9169d0,3.8681d0,3.8250d0,3.7855d0,3.7501d0,3.6753d0
     .,3.5499d0,3.4872d0,3.5401d0,3.8288d0,3.9217d0,3.9538d0,4.0054d0
     ./)
      r0ab(3571:3640)=(/
     . 2.8388d0,2.7890d0,3.4329d0,3.5593d0,3.3488d0,3.2486d0,3.1615d0
     .,3.1000d0,3.0394d0,2.9165d0,3.5267d0,3.7479d0,3.6650d0,3.6263d0
     .,3.5658d0,3.5224d0,3.4762d0,3.3342d0,3.7738d0,4.0333d0,3.9568d0
     .,3.8975d0,3.8521d0,3.4929d0,3.7830d0,3.7409d0,3.7062d0,3.6786d0
     .,3.6471d0,3.6208d0,3.6337d0,3.6519d0,3.6363d0,3.6278d0,3.6110d0
     .,3.4825d0,3.8795d0,4.1448d0,4.0736d0,4.0045d0,3.6843d0,3.6291d0
     .,3.5741d0,3.5312d0,3.4974d0,3.4472d0,3.4034d0,3.7131d0,3.7557d0
     .,3.7966d0,3.8005d0,3.8068d0,3.8015d0,3.6747d0,4.0222d0,4.3207d0
     .,4.2347d0,4.2191d0,4.1990d0,4.1811d0,4.1666d0,4.1521d0,4.1401d0
     .,3.9970d0,3.9943d0,3.9592d0,4.0800d0,4.0664d0,4.0559d0,4.0488d0
     ./)
      r0ab(3641:3710)=(/
     . 3.9882d0,4.0035d0,3.9539d0,3.9138d0,3.8798d0,3.8355d0,3.5359d0
     .,3.4954d0,3.3962d0,3.5339d0,3.7595d0,3.8250d0,3.8408d0,3.8600d0
     .,3.8644d0,2.7412d0,2.7489d0,3.3374d0,3.3950d0,3.3076d0,3.1910d0
     .,3.0961d0,3.0175d0,3.0280d0,2.8929d0,3.4328d0,3.5883d0,3.6227d0
     .,3.5616d0,3.4894d0,3.4241d0,3.3641d0,3.3120d0,3.6815d0,3.8789d0
     .,3.8031d0,3.7413d0,3.6939d0,3.4010d0,3.6225d0,3.5797d0,3.5443d0
     .,3.5139d0,3.4923d0,3.4642d0,3.5860d0,3.5849d0,3.5570d0,3.5257d0
     .,3.4936d0,3.4628d0,3.7874d0,3.9916d0,3.9249d0,3.8530d0,3.5932d0
     .,3.5355d0,3.4757d0,3.4306d0,3.3953d0,3.3646d0,3.3390d0,3.5637d0
     .,3.7053d0,3.7266d0,3.7177d0,3.6996d0,3.6775d0,3.6558d0,3.9331d0
     ./)
      r0ab(3711:3780)=(/
     . 4.1655d0,4.0879d0,4.0681d0,4.0479d0,4.0299d0,4.0152d0,4.0006d0
     .,3.9883d0,3.8500d0,3.8359d0,3.8249d0,3.9269d0,3.9133d0,3.9025d0
     .,3.8948d0,3.8422d0,3.8509d0,3.7990d0,3.7570d0,3.7219d0,3.6762d0
     .,3.4260d0,3.3866d0,3.3425d0,3.5294d0,3.7022d0,3.7497d0,3.7542d0
     .,3.7494d0,3.7370d0,3.7216d0,3.4155d0,3.0522d0,4.2541d0,3.8218d0
     .,4.0438d0,3.5875d0,3.3286d0,3.1682d0,3.0566d0,2.9746d0,4.3627d0
     .,4.0249d0,4.6947d0,4.1718d0,3.8639d0,3.6735d0,3.5435d0,3.4479d0
     .,4.6806d0,4.3485d0,4.2668d0,4.1690d0,4.1061d0,4.1245d0,4.0206d0
     .,3.9765d0,3.9458d0,3.9217d0,3.9075d0,3.8813d0,3.9947d0,4.1989d0
     .,3.9507d0,3.7960d0,3.6925d0,3.6150d0,4.8535d0,4.5642d0,4.4134d0
     ./)
      r0ab(3781:3850)=(/
     . 4.3688d0,4.3396d0,4.2879d0,4.2166d0,4.1888d0,4.1768d0,4.1660d0
     .,4.1608d0,4.0745d0,4.2289d0,4.4863d0,4.2513d0,4.0897d0,3.9876d0
     .,3.9061d0,5.0690d0,5.0446d0,4.6186d0,4.6078d0,4.5780d0,4.5538d0
     .,4.5319d0,4.5101d0,4.4945d0,4.1912d0,4.2315d0,4.5534d0,4.4373d0
     .,4.4224d0,4.4120d0,4.4040d0,4.2634d0,4.7770d0,4.6890d0,4.6107d0
     .,4.5331d0,4.4496d0,4.4082d0,4.3095d0,4.2023d0,4.0501d0,4.2595d0
     .,4.5497d0,4.3056d0,4.1506d0,4.0574d0,3.9725d0,5.0796d0,3.0548d0
     .,3.3206d0,3.8132d0,3.9720d0,3.7675d0,3.7351d0,3.5167d0,3.5274d0
     .,3.3085d0,3.1653d0,3.9500d0,4.1730d0,4.0613d0,4.1493d0,3.8823d0
     .,4.0537d0,3.8200d0,3.6582d0,4.3422d0,4.5111d0,4.3795d0,4.3362d0
     ./)
      r0ab(3851:3920)=(/
     . 4.2751d0,3.7103d0,4.1973d0,4.1385d0,4.1129d0,4.0800d0,4.0647d0
     .,4.0308d0,4.0096d0,4.1619d0,3.9360d0,4.1766d0,3.9705d0,3.8262d0
     .,4.5348d0,4.7025d0,4.5268d0,4.5076d0,3.9562d0,3.9065d0,3.8119d0
     .,3.7605d0,3.7447d0,3.7119d0,3.6916d0,4.1950d0,4.2110d0,4.3843d0
     .,4.1631d0,4.4427d0,4.2463d0,4.1054d0,4.7693d0,5.0649d0,4.7365d0
     .,4.7761d0,4.7498d0,4.7272d0,4.7076d0,4.6877d0,4.6730d0,4.4274d0
     .,4.5473d0,4.5169d0,4.5975d0,4.5793d0,4.5667d0,4.5559d0,4.3804d0
     .,4.6920d0,4.6731d0,4.6142d0,4.5600d0,4.4801d0,4.0149d0,3.8856d0
     .,3.7407d0,4.1545d0,4.2253d0,4.4229d0,4.1923d0,4.5022d0,4.3059d0
     .,4.1591d0,4.7883d0,4.9294d0,3.3850d0,3.4208d0,3.7004d0,3.8800d0
     ./)
      r0ab(3921:3990)=(/
     . 3.9886d0,3.9040d0,3.6719d0,3.6547d0,3.4625d0,3.3370d0,3.8394d0
     .,4.0335d0,4.2373d0,4.3023d0,4.0306d0,4.1408d0,3.9297d0,3.7857d0
     .,4.1907d0,4.3230d0,4.2664d0,4.2173d0,4.1482d0,3.6823d0,4.0711d0
     .,4.0180d0,4.0017d0,3.9747d0,3.9634d0,3.9383d0,4.1993d0,4.3205d0
     .,4.0821d0,4.2547d0,4.0659d0,3.9359d0,4.3952d0,4.5176d0,4.3888d0
     .,4.3607d0,3.9583d0,3.9280d0,3.8390d0,3.7971d0,3.7955d0,3.7674d0
     .,3.7521d0,4.1062d0,4.3633d0,4.2991d0,4.2767d0,4.4857d0,4.3039d0
     .,4.1762d0,4.6197d0,4.8654d0,4.6633d0,4.5878d0,4.5640d0,4.5422d0
     .,4.5231d0,4.5042d0,4.4901d0,4.3282d0,4.3978d0,4.3483d0,4.4202d0
     .,4.4039d0,4.3926d0,4.3807d0,4.2649d0,4.6135d0,4.5605d0,4.5232d0
     ./)
      r0ab(3991:4060)=(/
     . 4.4676d0,4.3948d0,4.0989d0,3.9864d0,3.8596d0,4.0942d0,4.2720d0
     .,4.3270d0,4.3022d0,4.5410d0,4.3576d0,4.2235d0,4.6545d0,4.7447d0
     .,4.7043d0,3.0942d0,3.2075d0,3.5152d0,3.6659d0,3.8289d0,3.7459d0
     .,3.5156d0,3.5197d0,3.3290d0,3.2069d0,3.6702d0,3.8448d0,4.0340d0
     .,3.9509d0,3.8585d0,3.9894d0,3.7787d0,3.6365d0,4.1425d0,4.1618d0
     .,4.0940d0,4.0466d0,3.9941d0,3.5426d0,3.8952d0,3.8327d0,3.8126d0
     .,3.7796d0,3.7635d0,3.7356d0,4.0047d0,3.9655d0,3.9116d0,4.1010d0
     .,3.9102d0,3.7800d0,4.2964d0,4.3330d0,4.2622d0,4.2254d0,3.8195d0
     .,3.7560d0,3.6513d0,3.5941d0,3.5810d0,3.5420d0,3.5178d0,3.8861d0
     .,4.1459d0,4.1147d0,4.0772d0,4.3120d0,4.1207d0,3.9900d0,4.4733d0
     ./)
      r0ab(4061:4130)=(/
     . 4.6157d0,4.4580d0,4.4194d0,4.3954d0,4.3739d0,4.3531d0,4.3343d0
     .,4.3196d0,4.2140d0,4.2339d0,4.1738d0,4.2458d0,4.2278d0,4.2158d0
     .,4.2039d0,4.1658d0,4.3595d0,4.2857d0,4.2444d0,4.1855d0,4.1122d0
     .,3.7839d0,3.6879d0,3.5816d0,3.8633d0,4.1585d0,4.1402d0,4.1036d0
     .,4.3694d0,4.1735d0,4.0368d0,4.5095d0,4.5538d0,4.5240d0,4.4252d0
     .,3.0187d0,3.1918d0,3.5127d0,3.6875d0,3.7404d0,3.6943d0,3.4702d0
     .,3.4888d0,3.2914d0,3.1643d0,3.6669d0,3.8724d0,3.9940d0,4.0816d0
     .,3.8054d0,3.9661d0,3.7492d0,3.6024d0,4.0428d0,4.1951d0,4.1466d0
     .,4.0515d0,4.0075d0,3.5020d0,3.9158d0,3.8546d0,3.8342d0,3.8008d0
     .,3.7845d0,3.7549d0,3.9602d0,3.8872d0,3.8564d0,4.0793d0,3.8835d0
     ./)
      r0ab(4131:4200)=(/
     . 3.7495d0,4.2213d0,4.3704d0,4.3300d0,4.2121d0,3.7643d0,3.7130d0
     .,3.6144d0,3.5599d0,3.5474d0,3.5093d0,3.4853d0,3.9075d0,4.1115d0
     .,4.0473d0,4.0318d0,4.2999d0,4.1050d0,3.9710d0,4.4320d0,4.6706d0
     .,4.5273d0,4.4581d0,4.4332d0,4.4064d0,4.3873d0,4.3684d0,4.3537d0
     .,4.2728d0,4.2549d0,4.2032d0,4.2794d0,4.2613d0,4.2491d0,4.2375d0
     .,4.2322d0,4.3665d0,4.3061d0,4.2714d0,4.2155d0,4.1416d0,3.7660d0
     .,3.6628d0,3.5476d0,3.8790d0,4.1233d0,4.0738d0,4.0575d0,4.3575d0
     .,4.1586d0,4.0183d0,4.4593d0,4.5927d0,4.4865d0,4.3813d0,4.4594d0
     .,2.9875d0,3.1674d0,3.4971d0,3.6715d0,3.7114d0,3.6692d0,3.4446d0
     .,3.4676d0,3.2685d0,3.1405d0,3.6546d0,3.8579d0,3.9637d0,4.0581d0
     ./)
      r0ab(4201:4270)=(/
     . 3.7796d0,3.9463d0,3.7275d0,3.5792d0,4.0295d0,4.1824d0,4.1247d0
     .,4.0357d0,3.9926d0,3.4827d0,3.9007d0,3.8392d0,3.8191d0,3.7851d0
     .,3.7687d0,3.7387d0,3.9290d0,3.8606d0,3.8306d0,4.0601d0,3.8625d0
     .,3.7269d0,4.2062d0,4.3566d0,4.3022d0,4.1929d0,3.7401d0,3.6888d0
     .,3.5900d0,3.5350d0,3.5226d0,3.4838d0,3.4594d0,3.8888d0,4.0813d0
     .,4.0209d0,4.0059d0,4.2810d0,4.0843d0,3.9486d0,4.4162d0,4.6542d0
     .,4.5005d0,4.4444d0,4.4196d0,4.3933d0,4.3741d0,4.3552d0,4.3406d0
     .,4.2484d0,4.2413d0,4.1907d0,4.2656d0,4.2474d0,4.2352d0,4.2236d0
     .,4.2068d0,4.3410d0,4.2817d0,4.2479d0,4.1921d0,4.1182d0,3.7346d0
     .,3.6314d0,3.5168d0,3.8582d0,4.0927d0,4.0469d0,4.0313d0,4.3391d0
     ./)
      r0ab(4271:4340)=(/
     . 4.1381d0,3.9962d0,4.4429d0,4.5787d0,4.4731d0,4.3588d0,4.4270d0
     .,4.3957d0,2.9659d0,3.1442d0,3.4795d0,3.6503d0,3.6814d0,3.6476d0
     .,3.4222d0,3.4491d0,3.2494d0,3.1209d0,3.6324d0,3.8375d0,3.9397d0
     .,3.8311d0,3.7581d0,3.9274d0,3.7085d0,3.5598d0,4.0080d0,4.1641d0
     .,4.1057d0,4.0158d0,3.9726d0,3.4667d0,3.8802d0,3.8188d0,3.7989d0
     .,3.7644d0,3.7474d0,3.7173d0,3.9049d0,3.8424d0,3.8095d0,4.0412d0
     .,3.8436d0,3.7077d0,4.1837d0,4.3366d0,4.2816d0,4.1686d0,3.7293d0
     .,3.6709d0,3.5700d0,3.5153d0,3.5039d0,3.4684d0,3.4437d0,3.8663d0
     .,4.0575d0,4.0020d0,3.9842d0,4.2612d0,4.0643d0,3.9285d0,4.3928d0
     .,4.6308d0,4.4799d0,4.4244d0,4.3996d0,4.3737d0,4.3547d0,4.3358d0
     ./)
      r0ab(4341:4410)=(/
     . 4.3212d0,4.2275d0,4.2216d0,4.1676d0,4.2465d0,4.2283d0,4.2161d0
     .,4.2045d0,4.1841d0,4.3135d0,4.2562d0,4.2226d0,4.1667d0,4.0932d0
     .,3.7134d0,3.6109d0,3.4962d0,3.8352d0,4.0688d0,4.0281d0,4.0099d0
     .,4.3199d0,4.1188d0,3.9768d0,4.4192d0,4.5577d0,4.4516d0,4.3365d0
     .,4.4058d0,4.3745d0,4.3539d0,2.8763d0,3.1294d0,3.5598d0,3.7465d0
     .,3.5659d0,3.5816d0,3.3599d0,3.4024d0,3.1877d0,3.0484d0,3.7009d0
     .,3.9451d0,3.8465d0,3.9873d0,3.7079d0,3.9083d0,3.6756d0,3.5150d0
     .,4.0829d0,4.2780d0,4.1511d0,4.1260d0,4.0571d0,3.4865d0,3.9744d0
     .,3.9150d0,3.8930d0,3.8578d0,3.8402d0,3.8073d0,3.7977d0,4.0036d0
     .,3.7604d0,4.0288d0,3.8210d0,3.6757d0,4.2646d0,4.4558d0,4.2862d0
     ./)
      r0ab(4411:4465)=(/
     . 4.2122d0,3.7088d0,3.6729d0,3.5800d0,3.5276d0,3.5165d0,3.4783d0
     .,3.4539d0,3.9553d0,3.9818d0,4.2040d0,3.9604d0,4.2718d0,4.0689d0
     .,3.9253d0,4.4869d0,4.7792d0,4.4918d0,4.5342d0,4.5090d0,4.4868d0
     .,4.4680d0,4.4486d0,4.4341d0,4.2023d0,4.3122d0,4.2710d0,4.3587d0
     .,4.3407d0,4.3281d0,4.3174d0,4.1499d0,4.3940d0,4.3895d0,4.3260d0
     .,4.2725d0,4.1961d0,3.7361d0,3.6193d0,3.4916d0,3.9115d0,3.9914d0
     .,3.9809d0,3.9866d0,4.3329d0,4.1276d0,3.9782d0,4.5097d0,4.6769d0
     .,4.5158d0,4.3291d0,4.3609d0,4.3462d0,4.3265d0,4.4341d0
     ./)

      k=0
      DO i=1,max_elem
         DO j=1,i
            k=k+1
            r(i,j)=r0ab(k)/autoang
            r(j,i)=r0ab(k)/autoang
         ENDDO
      ENDDO

      END SUBROUTINE setr0ab


