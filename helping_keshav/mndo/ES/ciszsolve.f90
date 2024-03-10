!
! created by jie liu on 09/16
!
module cisz_solve

use es_global
use cis_solve

private

integer ndim

public z_solve

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve Ax=b equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_solve(H,nq)
  implicit none
  integer ierror,nq,iprt
  real*8 H(*)
  real*8, dimension(:), allocatable :: sAs,Zs

  ndim = nq
  allocate(sAs(ns*ns),Zs(ns*nz))
  call matmult(Tv,At,sAs,ns,ns,ndim,ndim,ndim,ns,2)
  call matmult(Tv,Zr,Zs,ns,nz,ndim,ndim,ndim,ns,2)
  if(iprtz.eq.6) then
    write(6,*) "A in subspace"
    call matprnt(sAs,ns,ns,6)
    write(6,*) "B in subspace"
    call matprnt(Zs,ns,nz,6)
  endif

  ! Diagonalization in subspace
  ierror = 0
  call axeqb(sAs,Zs,ns,nz,ns,ns,ierror)
  if(ierror.ne.0) then
    write(NB6,*) "CIS DIAGON fails"
    stop
  endif
  if(iprtz.eq.6) then
    write(6,*) "solution in subspace"
    call matprnt(Zs,ns,nz,6)
  endif

  ! Check for convergence
  call z_chkcnv(Zs,Zr,H)

  deallocate(sAs)
  deallocate(Zs)
end subroutine z_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check convergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_chkcnv(Zs,B,H)
  implicit none
  integer i,is,m
  real*8 Zs(ns,*),B(ndim,*),H(*)
  real*8, dimension(:), allocatable :: R
  real*8 dot, tol, maxd

  allocate(R(ndim))

  maxd = 0d0
  tol = 1d1**(1-nrmz)
  tol = ndim * tol *tol
  is = ns
  do m=1,nz
    dot = 0d0
    call veccopy(R,B(1,m),ndim)
    do i=1,ns
      call vecadd2(R,-Zs(i,m),At(1,i),R,ndim)
    enddo
    do i=1,ndim
      if(H(i).eq.0.d0) then
        R(i) = 0.d0
      else 
        R(i) = R(i) / H(i)
      endif
      dot = dot + R(i)**2
    enddo
    if(dot.gt.maxd) maxd = dot
    if(dot.gt.tol) call ssexpand(R,is,ndim)
  enddo

  maxd = sqrt(maxd/ndim)
  nl = is - ns
  ns = is
  if(nl.eq.0) then
    if(iprtz.ge.2) write(NB6,*) niter, nl, maxd," iterz converged"
    if(.not.sasfcis) then
      call matmult(Tv,Zs,Zv,ndim,nz,ns,ndim,ns,ndim,1)
    else
      call matmult(Tv,Zs,Zr,ndim,nz,ns,ndim,ns,ndim,1)
    endif
  else if(ns.gt.maxdav) then
    write(NB6,*) "Max Subspace is reached in iterz"
    stop
  else
    if(iprtz.ge.5) write(NB6,*) niter, nl, maxd
  endif

  return
end subroutine z_chkcnv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ax=b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine axeqb(AA,BB,N,NRHS,LDA,LDB,INFO)
  implicit none
  integer N,NRHS,LDA,LDB,INFO
  real*8 AA(*),BB(*)
  integer, dimension(:), allocatable :: IPIV

  allocate(IPIV(N))

  call dgesv( N, NRHS, AA, LDA, IPIV, BB, LDB, INFO )

  deallocate(IPIV)
  return
end subroutine axeqb

end module cisz_solve
