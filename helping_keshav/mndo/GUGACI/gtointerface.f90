!
! This file contains interface routines for calling module subroutines
! from outside of this directory.
!

Subroutine CalcPropInt(S, R, Del, RxDel, Z, First, Last, Coord, zs, zp, zd, ns, np, nd, LMZ, natoms, naos, stdout, ecp)
  Use gaussoei, Only: CalcPropOEI
  Implicit None
  Integer                                  :: LMZ, natoms, naos
  Double Precision, Dimension(naos,naos)   :: S
  Double Precision, Dimension(naos,naos,3) :: R, Del, RxDel
  Integer,          Dimension(natoms)      :: Z, First, Last
  Double Precision, Dimension(3,natoms)    :: Coord
  Double Precision, Dimension(LMZ)         :: zs, zp, zd
  Integer,          Dimension(LMZ)         :: ns, np, nd
  Integer                                  :: stdout
  Logical                                  :: ecp
  ! end of dummy parameters
  Integer,          Dimension(LMZ)         :: ns0, np0, nd0

  ns0 = 0
  np0 = 0
  nd0 = 0
                                               ns0(Z(1:natoms)) = ns(Z(1:natoms))
  Where (Last(1:natoms)-First(1:natoms) >= 1)  np0(Z(1:natoms)) = np(Z(1:natoms))
  Where (Last(1:natoms)-First(1:natoms) >= 4)  nd0(Z(1:natoms)) = nd(Z(1:natoms))

  Call CalcPropOEI(S, R, Del, RxDel, Z, First, Last, Coord, zs, zp, zd, ns0, np0, nd0, natoms, stdout, ecp)

  Return
End Subroutine CalcPropInt
