
module cissym

use es_global

public 

integer numlab
integer, dimension(:), allocatable :: symlabel

contains

subroutine orbsym(V,nsym)
  implicit none
  real*8  V(*)
  logical uhf
  integer isub, ilab, jlab, nsym(*)
  integer, dimension(:), allocatable :: isym
  integer IRREPN,IRREPA,NUMAT
  character*3 NGROUP, grplab
  character*4 IRREP
  common /SYMLAB/ NGROUP(7),IRREP(26)
  common /SYMLIM/ IRREPA(7),IRREPN(7)
  common /UHF   / UHF
  common /ATOMS / NUMAT

  allocate(isym(3*NUMAT))  

  call mosym(V(jCa),V(jEa),LM2,LM3,NOrb,nsym,isym,isub,0)
  if(isub.gt.0 .and. ncisym.ge.0) then
    numlab = irrepn(isub)
    ilab   = irrepa(isub)+1
    jlab   = irrepa(isub)+numlab
    grplab = ngroup(isub)
    labels(1:numlab) = irrep(ilab:jlab)
  else
    numlab = 1
    labels(1) = 'A   '
    grplab = 'C1 '
    nsym(1:LM3)   = 1
  endif
  if(numlab.lt.8) labels(numlab+1:8) = 'XXX '

  if(uhf) then
    call mosym(V(jCb),V(jEb),LM2,LM3,NOrb,nsym(LM3+1),isym,isub,iprtci)
  endif
 
  deallocate(isym)
  return
end subroutine orbsym

subroutine ovsym(V)
  implicit none
  real*8  V(*),dmax
  integer i,j,k,m,n,s1,s2,s1s2
  integer, dimension(:), allocatable :: nsym
  integer, dimension(8,8) :: MultiplicationTable
  data MultiplicationTable /  1, 2, 3, 4, 5, 6, 7, 8,  &
                           &  2, 1, 4, 3, 6, 5, 8, 7,  &
                           &  3, 4, 1, 2, 7, 8, 5, 6,  &
                           &  4, 3, 2, 1, 8, 7, 6, 5,  &
                           &  5, 6, 7, 8, 1, 2, 3, 4,  &
                           &  6, 5, 8, 7, 2, 1, 4, 3,  &
                           &  7, 8, 5, 6, 3, 4, 1, 2,  &
                           &  8, 7, 6, 5, 4, 3, 2, 1   /

  allocate(jsym(nroots))
  if(.not.is_sym) then 
    jsym(:) = 1
    labels(1) = 'A   '
    return
  endif

  allocate(symlabel(NOV))
  allocate(nsym(LM3*nden))

  call orbsym(V,nsym)
  k = 1
  if(sasfcis) then
    s1 = nsym(NOb+1)
    s2 = nsym(NOb+2)
    s1s2 = MultiplicationTable(s1,s2)
    symlabel(1) = MultiplicationTable(s1,s2)
    symlabel(2) = MultiplicationTable(s1,s2)
    symlabel(3) = MultiplicationTable(s1,s1)
    do i=1,NOb
      m = nsym(i)
      symlabel(3+i) = MultiplicationTable(m,s1)
      symlabel(3+NOb+i) = MultiplicationTable(m,s2)
    enddo
    do i=1,NVa
      m = nsym(i+NOa)
      symlabel(3+NOb*2+i) = MultiplicationTable(m,s1)
      symlabel(3+NOb*2+NVa+i) = MultiplicationTable(m,s2)
    enddo
    do j=1,NOb
      do i=1,NVa
        m = nsym(i+NOa)
        n = nsym(j)
        symlabel(3+NOb*2+NVa*2+k) = MultiplicationTable(m,n)
        symlabel(3+NOb*2+NVa*2+NVa*NOb+k) = MultiplicationTable(m,n)
        if(sf_xcis) then
          symlabel(3+NOb*2+NVa*2+NVa*NOb*2+k) = MultiplicationTable(m,n)
          symlabel(3+NOb*2+NVa*2+NVa*NOb*3+k) = MultiplicationTable(m,n)
        endif
        k = k + 1
      enddo
    enddo
    do i=1,3+NOb*2+NVa*2+NVa*NOb*2
      symlabel(i) = MultiplicationTable(s1s2,symlabel(i))
    enddo
  else
    do i=1,NOa  
      do j=1,NVa
        m = nsym(i)
        n = nsym(j+NOa)
        symlabel(k) = MultiplicationTable(m,n)
        k = k + 1
      enddo
    enddo
    if(nden.eq.2) then
      do i=1,NOb  
        do j=1,NVb
          m = nsym(i+NOrb)
          n = nsym(j+NOa+NOrb)
          symlabel(k) = MultiplicationTable(m,n)
          k = k + 1
        enddo
      enddo
    else 
      symlabel(NOVa+1:NOV) = symlabel(1:NOVa)
    endif
  endif
  
  do m=1,nroots
    dmax = 0d0
    k = 1
    do j = 1,NOV
      if(abs(Xv(j,m)).gt.dmax) then
        k = j
        dmax = abs(Xv(j,m))
      endif
    enddo
    jsym(m) = symlabel(k)
  enddo
  if(.not.is_sf) jsym(1) = 1
 
  if(iprtci.ge.5) then
    do m=1,nroots
      write(nb6,*) "State:",m,"Symmetry:",labels(jsym(m))
    enddo
  endif

  deallocate(nsym)
  deallocate(symlabel)
  return
end subroutine ovsym

end module cissym
