module constraint
  use global
  use matrixlib
  implicit none

!------------------------------------------------------------------------------
! General:
!
! - Reference:
!   Implements J. Baker, J. Chem. Phys. 105 (1996) 192.
!
! - Procedure:
!   1. (done in <interface>.f90:hdlc_get_params)
!      Read in constraints (search.f90:search)
!   2. assign_cons
!      Assign constraints to particular residues and renumber atoms
!      (search.f90:search, after dlc_connect -> assign_cons)
!   3. ci_cons
!      Check in required connections (hdlclib.f90:hdlc_create after
!      hdlc_con_checkin -> ci_cons)
!   4. ck_angle_cons
!      Check in required bends and rots (hdlclib.f90:hdlc_create before
!      end if ! (associated(con)) -> ck_angle_cons)
!   5. ortho_cons
!      Project and orthogonalise active space to constraints
!      (dlc_manager.f:dlc_manager after dlc_make_ut -> ortho_cons)
!   6. split_cons
!      Separate between active and constrained coordinate values
!      (prfossa.f: optef_main after dlc_get_dlc_values and
!      dlc_get_dlc_gradients -> split_cons)
!   7. rest_cons
!      Restore Ut from constrained Ut and unprojected constraints
!      (prfossa.f: optef_main -> rest_cons)
!
! - Data to be passed between steps:
!        context scope  form symbol name description
!   1-2: global  local  a)               constraints specification
!   2-3: residue arg    b)               constraints specification
!   3-5: residue object c)               int. coord. seq. nbrs. of constr.
!   4-5: residue object c)               do.
!   5-6: residue object --   Ut,ncons    transposed V matrix, nbrs. of cns.
!   6-7: residue object d)               constraints values, number
!
! - Form specification:
!   a) iconstr(4,nconstr): i, j, k, l
!      vconstr(nconstr): value
!      atom sequence numbers global
!      k=l=0 -> bond, l=0 -> bond angle, all >0 -> torsion
!      type: array, double(5,ncnstr)
!   b) cns(7,ncons): i, j, k, l, itype, iseq, value; matrix
!      itype = 1: bond; i, j specified
!      itype = 2: bond angle; i, j, k specified
!      itype = 3: torsion; i, j, k, l specified
!      iseq is undefined: see c)
!      value is defined by (see d):
!        sequence number (mcnstr.eq.0)
!        target value (mcnstr.eq.1)
!      i, j, k, l are relative to the residue, ncons is stored per residue
!      object, as is icons(7,ncons) as well
!      type: matrix, (7,N), double
!   c) see b) but iseq points to the internal coordinate to be constrained
!      types: icons: array, (6,N), integer, vcons: array, (N), double
!   d) see b) but ivalue points to the target value in constraints(ivalue)
!      if extended Lagrangean is chosen, ivalue is kept from before
!   constraints(nconsvals) is stored globally
!
! - Notes:
!   * U/natoms <-> V/natoms/ncons/ (C/C_proj)?
!   * ncons has to be subtracted from n6 at each occurrence of if(internal)
!   * Correspondence between 'rows' and 'columns' by Alex Turner, Fortran
!     and mathematics is straightforward: the place of the index is always
!     the same, no matter whether the faster index or not.
!     First index: row, second index: column; a(row,column)
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
! subroutine assign_cons
!
! Arguments:
! cns(7,ncns):       matrix of constraints of one residue, format b (out)
! iconstr(4,ncnstr): constraints as read in, format a (in)
! nconstr:           number of constraints provided (in)
! mconstr:           type of constraints (in)
! vconstr(nconstr):  target or actual values (in)
! start:             first sequence number in the residue (in)
! finish:            last sequence number in the residue (in)
!
! Description: step 2 of the procedure described above
!------------------------------------------------------------------------------

  subroutine assign_cons (cns, iconstr, nconstr, vconstr, mconstr, start, &
       finish, internal)

! args
    logical internal
    integer nconstr, mconstr, start, finish
    integer, dimension(4,nconstr) :: iconstr
    real(kind=8), dimension(nconstr) :: vconstr
    type(matrix), pointer :: cns

! local vars
    logical li, lj, lk, ll
    integer i, iatom, j, k, kc, l, m, mcomp, ncns, type
    real(kind=8), dimension(:,:), allocatable :: icns_tmp

! begin, allocate three times more space for the Cartesian components
    if (ctrl%printl.ge.3) write (stdout, '(3X,A,/)') 'Looking for constraints'
    ncns = 0
    allocate (icns_tmp(7,nconstr*3))

! loop over constraints and check if they affect the current residue
    do kc = 1,nconstr
       i = iconstr(1,kc)
       j = iconstr(2,kc)
       k = iconstr(3,kc)
       l = iconstr(4,kc)

       li = (i.ge.start .and. i.le.finish)
       lj = (j.ge.start .and. j.le.finish)
       lk = (k.ge.start .and. k.le.finish)
       ll = (l.ge.start .and. l.le.finish)
       if (li .or. lj .or. lk .or. ll) then
          if (li .and. lj .and. lk .and. ll) then
             type = 3
          elseif ((l.eq.0) .and. li .and. lj .and. lk) then
             type = 2
          elseif ((l.eq.0) .and. (k.eq.0) .and. li .and. lj) then
             type = 1
          elseif ((l.le.0) .and. (k.le.0) .and. (j.le.0) .and. li) then
             type = 4
          else
             write (stdout, '(A,I5,A,I5,A,I5,A,I5,A)') 'Constraint ', i, &
                  ' - ', j, ' - ', k, ' - ', l, ' crosses residue boundary!'
             call hdlc_errflag ('Constraints error', 'stop')
          end if

! found an internal constraint
          if (type .le. 3) then
             ncns = ncns + 1
             if (ctrl%printl.ge.1) then
                write (stdout,'(3X,A,4I5,/)') &
                     'Found constraint between atoms ', &
                     (iconstr(m,kc), m = 1,type+1)
             end if
             do m = 1,4
                if (iconstr(m,kc).ne.0.0) then
                   icns_tmp(m,ncns) = real ((iconstr(m,kc)-(start-1)), 8)
                else
                   icns_tmp(m,ncns) = 0.0_8
                end if
             end do
             icns_tmp(5,ncns) = real (type, 8)
             icns_tmp(6,ncns) = 0.0_8
             if (mconstr.ge.1) then
                icns_tmp(7,ncns) = vconstr(kc)
             else
                icns_tmp(7,ncns) = 0.0_8
             end if

! found Cartesian constraints
          else ! (type .le. 3)
             if (internal) then
                write (stdout, '(A,I5,A)') 'Constraint ', i, &
                     ': no Cartesian components can be constrained with DLC'
                call hdlc_errflag ('Constraints error', 'stop')
             end if
             ncns = ncns + 1
             iatom = 1 + int((iconstr(1,kc)-1)/3)
             mcomp = iconstr(1,kc)-3*(iatom-1)
             if (ctrl%printl.ge.1) then
                write (stdout,'(3X,A,4I5,/)') &
                     'Found constrained Cartesian, atom ', iatom, &
                     ', component ', mcomp
             end if
             icns_tmp(1,ncns) = real ((mcomp+3*(iatom-start)), 8)
             do m = 2,4
                icns_tmp(m,ncns) = 0.0_8
             end do
             icns_tmp(5,ncns) = real (type, 8)
             icns_tmp(6,ncns) = 0.0_8
             if (mconstr.ge.1) then
                icns_tmp(7,ncns) = vconstr(kc)
             else
                icns_tmp(7,ncns) = 0.0_8
             end if
          end if ! (type .le. 3) ... else
       end if ! (li .or. lj .or. lk .or. ll)
    end do ! kc = 1,nconstr
    cns => matrix_create (7, ncns, 'constraints')
    do kc = 1,ncns
       i = matrix_set_column (cns, icns_tmp(1,kc), kc)
    end do
    deallocate (icns_tmp)
  end subroutine assign_cons

!------------------------------------------------------------------------------
! subroutine ci_cons
!
! Description:
! Step 3 of the procedure described above. The constraints provided in cns are
! copied to icons, and all connections not yet occurring in iconn are added.
! For bonds, the bond sequence number is stored in icons(iseq,*).
! The target value of the constraints is copied to vcons. This only matters if
! the extended Lagrangean method is chosen.
!
! Arguments:
! cns   matrix(7,ncons) (in)  constraint specification for the residue
! ncons                 (in)  number of constraints provided in cns
! icons (6,ncons)       (out) constraint specification copied from cns
! vcons (ncons)         (out) constraint values copied from cns
! nconn                 (i/o) number of connections between atoms in the res.
! iconn (2,nconn)       (i/o) connections list, new pointer returned
!------------------------------------------------------------------------------

  subroutine ci_cons (cns, ncons, icons, vcons, nconn, iconn)

! args
    type(matrix), pointer :: cns
    integer ncons, nconn
    integer, pointer, dimension(:,:) :: icons, iconn
    real (kind=8), pointer, dimension(:) :: vcons

! local vars
    integer, pointer, dimension(:,:) :: iconn_tmp, iconn_new
    integer el, k, kc
    integer ibnd, jbnd, kconn, kibnd, kjbnd, ktype, nconn_new

! begin, allocate space to hold potential new connections
    nconn_new = nconn
    allocate (iconn_tmp(2,ncons))

! copy the constraints
    do kc = 1,ncons
       do el = 1,6
          icons(el,kc) = int(cns%data(el,kc))
       end do
       vcons(kc) = cns%data(7,kc)

! check all constraints including angles and torsions
       ktype = icons(5,kc)
       jbnd = icons(1,kc)
       if (ctrl%printl.ge.2) then
          write (stdout,'(5X,A,I1,A,4I5)') 'Considering constraint (type ', &
               ktype, '): ', (icons(k,kc), k = 1,ktype+1)
       end if
       do k = 1,ktype
          ibnd = jbnd
          jbnd = icons(1+k,kc)

! check if the connection ibnd - jbnd already occurs
          do kconn = 1,nconn
             kibnd = iconn(1,kconn)
             if (ibnd.eq.kibnd .or. jbnd.eq.kibnd) then
                kjbnd = iconn(2,kconn)
                if (ibnd.eq.kjbnd .or. jbnd.eq.kjbnd) then

! connection is already defined - store its sequence number and break the loop
                   if (ktype.eq.1) then
                      icons(6,kc) = kconn
                      if (ctrl%printl.ge.2) then
                         write (stdout,'(7X,A,I4)') '... constraining bond ', &
                              kconn
                      endif
                   endif
                   goto 10
                endif
             endif
          end do ! kconn = 1,nconn

! connection is not yet defined - do it now
          nconn_new = nconn_new + 1
          iconn_tmp(1,nconn_new-nconn) = ibnd
          iconn_tmp(2,nconn_new-nconn) = jbnd
          if (ktype.eq.1) then
             icons(6,kc) = nconn_new
             if (ctrl%printl.ge.2) then
                write (stdout,'(9X,A,I4,A,I4,A,I4,/)') '... adding stretch ', &
                     ibnd, ' - ', jbnd, ' to constrain it as bond ', nconn_new
             endif
          else
             if (ctrl%printl.ge.2) then
                write (stdout,'(9X,A,I4,A,I4,/)') '... adding stretch ', &
                     ibnd, ' - ', jbnd
             endif
          endif

! jump here if the connection is already defined
10        continue
       end do ! do k = 1,ktype
    end do ! kc = 1,ncons

! add new connections to the old ones and reset pointers
    if (nconn_new .gt. nconn) then
       allocate (iconn_new(2,nconn_new))
       do k = 1,nconn
          iconn_new(1,k) = iconn(1,k)
          iconn_new(2,k) = iconn(2,k)
       end do
       do k = 1,nconn_new-nconn
          iconn_new(1,k+nconn) = iconn_tmp(1,k)
          iconn_new(2,k+nconn) = iconn_tmp(2,k)
       end do
       deallocate (iconn)
       iconn => iconn_new
    endif
    deallocate (iconn_tmp)
    nconn = nconn_new
  end subroutine ci_cons

!------------------------------------------------------------------------------
! subroutine ck_angle_cons
!
! Description:
! Step 4 of the procedure described above. The angles to be constrained must
! have been found already. Check their sequence numbers and store them.
! Only bond angles are scanned (ibend(4,*).eq.0), no impropers or linears.
!
! Arguments:
! ncons           (in)  number of constraints
! icons (6,ncons) (i/o) constraint specification
! nconn           (in)  number of connections between atoms in the res.
! nbend           (i/o) number of bends
! ibend (4,nbend) (i/o) bends list
! nrots           (i/o) number of torsions
! irots (4,nrots) (i/o) torsions list
!------------------------------------------------------------------------------

  subroutine ck_angle_cons (ncons, icons, nconn, ibend, nbend, irots, nrots)

! args
    integer ncons, nconn, nbend, nrots
    integer ibend(4,nbend), irots(4,nrots)
    integer, pointer, dimension(:,:) :: icons

! local vars
    integer k, kangle, kc, ktype, ni

! begin
    do kc = 1,ncons
       ktype = icons(5,kc)

! bends
       if (ktype.eq.2) then

          do kangle = 1,nbend
             if (ibend(4,kangle) .eq. 0) then
                if (ibend(2,kangle) .eq. icons(2,kc)) then
                   if ((ibend(1,kangle) .eq. icons(1,kc) .and. &
                        ibend(3,kangle) .eq. icons(3,kc)) .or. &
                        (ibend(1,kangle) .eq. icons(3,kc) .and. &
                        ibend(3,kangle) .eq. icons(1,kc))) then
                      goto 10
                   end if
                end if
             end if
          end do

! angle not found
          write (stdout,'(A,I5,A,I5,A,I5,A)') 'Constrained angle ', &
               icons(1,kc), ' - ', icons(2,kc), ' - ', icons(3,kc), &
               ' not found!'
          call hdlc_errflag ('Constraints error', 'abort')

! angle found
10        continue
          ni = kangle
          do k = 1,kangle
             if (ibend(4,k).lt.0 .and. ibend(1,k).ne.0) ni = ni + 1
          end do
          if (ctrl%printl.ge.2) write (stdout,'(5X,A,I3,A,4I5)') &
               'Constraining angle ', ni, ': ', (icons(k,kc), k = 1,3)
          icons(6,kc) = ni + nconn

! torsions
       elseif (ktype.eq.3) then

          do kangle = 1,nrots

             if (irots(2,kangle) .eq. icons(2,kc)) then
                if (irots(1,kangle) .eq. icons(1,kc) .and. &
                     irots(3,kangle) .eq. icons(3,kc) .and. &
                     irots(4,kangle) .eq. icons(4,kc)) goto 20
             elseif (irots(2,kangle) .eq. icons(3,kc)) then
                if (irots(1,kangle) .eq. icons(4,kc) .and. &
                     irots(3,kangle) .eq. icons(2,kc) .and. &
                     irots(4,kangle) .eq. icons(1,kc)) goto 20
             endif
          end do

! torsion not found
          write (stdout,'(A,I5,A,I5,A,I5,A,I5,A)') 'Constrained torsion ', &
               icons(1,kc), ' - ', icons(2,kc), ' - ', icons(3,kc), ' - ', &
               icons(4,kc), ' not found!'
          call hdlc_errflag ('Constraints error','abort')

! torsion found - add stretches, bends and linears to the sequence number
20        continue
          ni = kangle
          do k = 1,nbend
             if(ibend(4,k).lt.0.and.ibend(1,k).ne.0) ni = ni + 1
          end do
          if (ctrl%printl.ge.2) write (stdout,'(5X,A,I3,A,4I5)') &
               'Constraining torsion ', kangle, ': ', (icons(k,kc), k = 1,4)
          icons(6,kc) = ni + nconn + nbend

! Cartesians - add stretches, bends and linears to the sequence number
       elseif (ktype.eq.4) then ! (ktype.eq.2) .. elseif (ktype.eq.3)
          ni = icons(1,kc)
          kangle = ni
          do k = 1,nbend
             if(ibend(4,k).lt.0.and.ibend(1,k).ne.0) ni = ni + 1
          end do
          if (ctrl%printl.ge.2) write (stdout,'(5X,A,I3)') &
               'Constraining Cartesian ', kangle
          icons(6,kc) = ni + nconn + nbend
       end if ! (ktype.eq.2) .. elseif (ktype.eq.4)
    end do ! kc = 1,ncons
  end subroutine ck_angle_cons

!------------------------------------------------------------------------------
! subroutine ortho_cons
!
! Implements the method for adding constraints to natural DLC
! described by Baker et al.
! based on Schmidt orthogonalisation of the column vectors of U to the
! projected constraint vectors C
!
! The np-dimensional space of the internal, primitive, local coordinates is
! represented as in 
!
! Input:
!   Ut(m,np): Transpose of the nonredundant eigenspace of G
!     m rows: nonredundant eigenvectors of G represented using
!     np columns: primitive, local, internal coordinates
!   nc: number of internal coordinates to be constrained
!
! Output:
!   Ut(m-nc,np): Transpose of the active, nonredundant eigenspace of G
!   Ut(m-nc..m,np): Transpose of the constrained eigenspace of G
!     m rows: see above, orthogonal to the space {C(k), k=1..nc} projected
!
! Temporarily used:
!   C(np,nc): Constraints
!     np rows: primitive, local, internal coordinates
!     nc columns: constraint vectors represented using
!   np: number of primitive internal coordinates
!   m: dimension of the nonredundant eigenspace of G
!
! Description: step 5 of the procedure described above
!------------------------------------------------------------------------------


  subroutine ortho_cons (ut_mat, nc, icons)

! args
    integer nc, icons(7,*)
    type(matrix), pointer :: ut_mat

! local vars
    integer j, k, l, m, np
    type(matrix), pointer :: c_mat, v_mat
    real(kind=8), dimension(:), allocatable :: work

! begin, find dimensions, allocate space, np is always .gt. m
    m = matrix_dimension (ut_mat, 1)
    np = matrix_dimension (ut_mat, 2)
    c_mat => matrix_create (np, nc, 'C matrix')
    v_mat => matrix_create (np, m+nc, 'V matrix')
    allocate (work(np+m))

! set constraints matrix and project constrained coordinates
    call gen_cons (c_mat%data, icons, nc, np)
    call proj_cons (c_mat%data, ut_mat%data, work, m, nc, np)

! now we have Ut and Cp, start orthogonalisation
    call merge_cons (v_mat, c_mat, ut_mat, work, m, nc, np)
    call ortho_mat (v_mat%data, work, m, nc, np)

! move the constraints to the end of the non-zero vectors and transpose
    call move_cons (v_mat, ut_mat, work, m, nc, np)

! replace projected by unprojected constraints
    call gen_cons (c_mat%data, icons, nc, np)
    call unproj_cons (ut_mat, c_mat, work, m, nc, np)
    j = matrix_destroy (c_mat)
    j = matrix_destroy (v_mat)
    deallocate (work)

! end ortho_cons
  end subroutine ortho_cons

!------------------------------------------------------------------------------
! subroutine split_cons
!
! Is called by the 'method routine' hdlc_split_cons of hdlc_manager, but is
! implemented separately.
! matrix_set copies the data rather than setting the handle to the data
!
! Arguments:
! hdlc:    matrix of full HDLC coordinates (in)
!          matrix of active HDLC coordinates (out)
! lstore:  constrained HDLC coordinates are stored to vcons if (lstore) (in)
! natom:   number of atoms in the residue (in)
! ncons:   number of constraints in the residue (in)
! vcons:   (ncons) values of the constrained coordinates (out)
!
! Description: step 6 of the procedure described above
!------------------------------------------------------------------------------

  subroutine split_cons (hdlc, lstore, natom, ncons, vcons)

! args
    logical lstore
    integer natom, ncons
    real(kind=8), dimension(ncons) :: vcons
    type(matrix), pointer :: hdlc

! local vars
    integer idum, ncoords, nactive
    real(kind=8), dimension(:), allocatable :: dhdlc

! begin
    ncoords = 3 * natom
    nactive = ncoords - ncons

! move all data to temp array
    allocate (dhdlc(ncoords))
    idum = matrix_get (hdlc, dhdlc)
    idum = matrix_destroy (hdlc)

! store data to new HDLC matrix and residue%vcons if requested
    hdlc => matrix_create (nactive, 1, 'HDLC active')
    idum = matrix_set (hdlc, dhdlc)
    if (lstore) then
       call copy_coords (vcons, dhdlc(nactive+1), ncons)
    end if
    deallocate (dhdlc)

! end split_cons
  end subroutine split_cons

!------------------------------------------------------------------------------
! subroutine rest_cons
!
! Is called by the 'method routine' hdlc_rest_cons of hdlc_manager, but is
! implemented separately.
! matrix_set copies the data rather than setting the handle to the data
!
! Arguments:
! hdlc:    matrix of active HDLC coordinates (in)
!          matrix of full HDLC coordinates (out)
! natom:   number of atoms in the residue (in)
! ncons:   number of constraints in the residue (in)
! vcons:   (ncons) values of the constrained coordinates (in)
!
! Description: step 7 of the procedure described above
!------------------------------------------------------------------------------

  subroutine rest_cons (hdlc, natom, ncons, vcons)

! args
    integer natom, ncons
    real(kind=8), dimension(ncons) :: vcons
    type(matrix), pointer :: hdlc

! local vars
    integer idum, ncoords, nactive
    real(kind=8), dimension(:), allocatable :: dhdlc

! begin
    ncoords = 3 * natom
    nactive = ncoords - ncons

! move all data from old HDLC matrix and vcons to temp array
    allocate (dhdlc(ncoords))
    idum = matrix_get (hdlc, dhdlc)
    idum = matrix_destroy (hdlc)
    call copy_coords (dhdlc(nactive+1), vcons, ncons)

! store data to new HDLC
    hdlc => matrix_create (ncoords, 1, 'HDLC all')
    idum = matrix_set (hdlc, dhdlc)
    deallocate (dhdlc)

! end split_cons
  end subroutine rest_cons

!------------------------------------------------------------------------------
! Helpers
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This tiny helper is self explaining
!------------------------------------------------------------------------------

  subroutine copy_coords (target, source, length)

! args
    integer length
    real(kind=8), dimension(length) :: target, source

! local vars
    integer i

! begin
    do i = 1,length
       target(i) = source(i)
    end do

! end copy_coords
  end subroutine copy_coords

!------------------------------------------------------------------------------
! subroutine gen_cons
!
! Generates a constraints matrix in the primitive internal basis from the
! specification in icons(iseq,*)
!
! Arguments:
! cdat(np,nc): constraints matrix (out)
! icons(6,nc): constraints specification (in)
! nc:          number of constraints (in)
! np:          dimension of the primitive space
!------------------------------------------------------------------------------

  subroutine gen_cons (cdat, icons, nc, np)

! args
    integer nc, np, icons(6,nc)
    real(kind=8), dimension(np,nc) :: cdat

! local vars
    integer ip, ic

! begin
    do ic = 1,nc
       do ip = 1,np
          cdat(ip,ic) = 0.0_8
       end do
       cdat(icons(6,ic),ic) = 1.0_8
    end do

! end gen_cons
  end subroutine gen_cons

!------------------------------------------------------------------------------
! subroutine proj_cons
!
! Projects the vectors in the matrix c into the space spanned by utmat
!
! Arguments:
! cdat(np,nc): unprojected constraints matrix (in) / projected matrix (out)
! utdat(m,np): m-dimensional space of vectors of dimension np (in)
! work(m):     scratch array, used for dp_{ic,j} (in)
! m:           dimension of the space spanned by utmat (nonredundant) (in)
! nc:          number of constraints (in)
! np:          dimension of the space in which utmat is represented (in)
!------------------------------------------------------------------------------

  subroutine proj_cons (cdat, utdat, work, m, nc, np)

! args
    integer nc, np, m
    real(kind=8), dimension(m) :: work
    real(kind=8), dimension(np,nc) :: cdat
    real(kind=8), dimension(m,np) :: utdat

! local vars
    integer ic, ip, j

! begin, dp_{ic,j} = <C_ic|U_j>
    do ic = 1,nc
       do j = 1,m
          work(j) = 0.0_8
          do ip = 1,np
             work(j) = work(j) + utdat(j,ip)*cdat(ip,ic)
          end do
       end do

! C_ic = 0
       do ip = 1,np
          cdat(ip,ic) = 0.0_8
       end do

! C_ic = sum_j dp_{ic,j}*U_j
       do j = 1,m
          do ip = 1,np
             cdat(ip,ic) = cdat(ip,ic) + work(j)*utdat(j,ip)
          end do
       end do
    end do ! ic = 1,nc

! end proj_cons
  end subroutine proj_cons

!------------------------------------------------------------------------------
! subroutine merge_cons
!
! Composes a new matrix V out of [C,Ut]
!
! Arguments:
! v_mat:  matrix (np,m+nc) of [cmat,utmat transposed] (out)
! c_mat:  matrix (np,nc) of projected constraints (in)
! ut_mat: matrix (m,np) of U transposed (in)
! work:   work array(np) (in)
! m:      dimension of the space spanned by utmat (nonredundant) (in)
! nc:     number of constraints (in)
! np:     dimension of the space in which utmat is represented (in)
!------------------------------------------------------------------------------

  subroutine merge_cons (v_mat, c_mat, ut_mat, work, m, nc, np)

! args
    integer m, nc, np
    type(matrix), pointer :: v_mat, c_mat, ut_mat
    real(kind=8), dimension(np) :: work

! local vars
    integer i, idum, j

! begin
    do i = 1,nc
       idum = matrix_get_column (c_mat, work, i)
       idum = matrix_set_column (v_mat, work, i)
    end do
    j = nc
    do i = 1,m
       j = j + 1
       idum = matrix_get_row (ut_mat, work, i)
       idum = matrix_set_column (v_mat, work, j)
    end do

! end merge_cons
  end subroutine merge_cons

!------------------------------------------------------------------------------
! subroutine ortho_mat
!
! Applies Schmidt orthogonalisation to the columns of vmat
! Taken from mankyopt: orthog.F
!
! Arguments:
! v_dat(np,nc+m): [cmat,utmat transposed] (in), orthogonalised (out)
! work(np):       work array (in)
! m, nc, np:      see above (in)
!------------------------------------------------------------------------------

  subroutine ortho_mat (v_dat, work, m, nc, np)

! args
    integer m, nc, np
    real(kind=8), dimension(np,nc+m) :: v_dat
    real(kind=8), dimension(np) :: work

! local vars
    integer i, idum, j, k, nelem, nvec
    real(kind=8) dnorm, scapro, tol

! data
    data tol/1.0d-10/

! begin, orthogonalise vectors i = 2,nvec
    nvec = m + nc
    nelem = np
    do i = 1,nvec
       do j = 1,nelem
          work(j) = 0.0_8
       end do

! make vector i orthogonal to vectors  k = 1,i-1
       do k = 1,i-1
          scapro = 0.0_8
          do j = 1,nelem
             scapro = scapro + v_dat(j,i)*v_dat(j,k)
          end do
          do j = 1,nelem
             work(j) = work(j) + v_dat(j,k)*scapro
          end do
       end do

! subtract the collinear vector to make vector i orthogonal to k = 1,i-1
       do j = 1,nelem
          v_dat(j,i) = v_dat(j,i) - work(j)
       end do

! normalise vector i
       dnorm = 0.0_8
       do j = 1,nelem
          dnorm = dnorm + v_dat(j,i)*v_dat(j,i)
       end do
       if (abs(dnorm).lt.tol)then
          dnorm = 0.0_8
       else
          dnorm = 1.0_8/dsqrt(dnorm)
       end if
       do j = 1,nelem
          v_dat(j,i) = v_dat(j,i)*dnorm
       end do
    end do

! report resulting matrix
!   do i = 1,nvec
!      work(i) = 0.0_8
!      do j = 1,nelem
!         work(i) = work(i) + v_dat(j,i)*v_dat(j,i)
!      end do
!   end do

! end ortho_mat
  end subroutine ortho_mat

!------------------------------------------------------------------------------
! subroutine move_cons
!
! Moves the constraints behind the active space in the V matrix, transposes
! V and stores it in Ut
!
! Arguments: see merge_cons
!------------------------------------------------------------------------------

  subroutine move_cons (v_mat, ut_mat, work, m, nc, np)

! args
    integer m, nc, np
    type(matrix), pointer :: v_mat, ut_mat
    real(kind=8), dimension(np) :: work

! local vars
    integer is, it, i
    real(kind=8) dnorm, tol

! data
    data tol /1.0d-10/

! begin
    do is = 1,nc
       i = matrix_get_column (v_mat, work, is)
       i = matrix_set_row (ut_mat, work, m-nc+is)
    end do
    it = 0
    do is = nc+1,nc+m
       i = matrix_get_column (v_mat, work, is)
       dnorm = 0.0_8
       do i = 1,np
          dnorm = dnorm + work(i)*work(i)
       end do
       if (dnorm .gt. tol) then
          it = it + 1
          if (it .gt. m-nc) then
             write(stdout,'(A,I5)') 'Too many active vectors, required: ', m-nc
             call hdlc_errflag ('Constraints error', 'abort')
          else
             i = matrix_set_row (ut_mat, work, it)
          end if
       end if
    end do

! end move_cons
  end subroutine move_cons

!------------------------------------------------------------------------------
! subroutine unproj_cons
!
! Replaces projected constraints by the unprojected ones in the matrix
! U transposed
!
! Arguments: see merge_cons
!------------------------------------------------------------------------------

  subroutine unproj_cons (ut_mat, c_mat, work, m, nc, np)

! args
    integer m, nc, np
    type(matrix), pointer :: ut_mat, c_mat
    real(kind=8), dimension(np) :: work

! local vars
    integer is, i

! begin
    do is = 1,nc
       i = matrix_get_column (c_mat, work, is)
       i = matrix_set_row (ut_mat, work, m-nc+is)
    end do

! end unproj_cons
  end subroutine unproj_cons

end module constraint
