module primitive
  use global
  use matrixlib
  implicit none

contains

!------------------------------------------------------------------------------
! subroutine connect_prim
!
! Automatically constructs a connection pattern for the system
!
! Method
! ======
!
! This method goes in two stages.  First a simple covalent radius method is
! employed, second the shortest branched path that internconnects every atom
! is found. The result is the fusion of these two connection patterns.
!
! The covalent interconnection method uses a fixed radius for each row of the
! periodic table, this is a bit crude. A connection is made when the square of
! the sum of the atomic radii is greater than 0.8 times the square of the
! inter atomic distance.
!
! The minimum distance method can be described as a two bucket algorithm.
! The is a 'known' bucket and an 'unknown' bucket.
!
! Initially, one atom is placed in the known and the rest in unknown.
! The algorithm cycles by finding the mimium distance between any atom in
! the known and unknown buckets.  For each cycle, a connection is made
! between the atoms found to have this shortest distance, and the atom
! of the connection that is in the unknows, and the unknown atom is placed
! into the known bucket.
!
! Once all the atoms are in the known bucket, the connection pattern is 
! finished.
!
! Parameters
! =========
!
! natom               (in)  number of atoms, integer
! types   (natom)     (in)  integer array of nuclear charges
! xyz matrix(3,natom) (in)  cartesian coordiantes in AU
! nconn               (out) number of connections found
! iconn   (2,nconn)   (out) pointer to connections list
!------------------------------------------------------------------------------

  subroutine connect_prim (natom, types, nconn, iconn, xyz)

! args
    integer natom, nconn
    integer, dimension(natom) :: types
    integer, pointer, dimension(:,:) :: iconn
    type(matrix), pointer :: xyz

! local vars
    logical, allocatable, dimension(:) :: lconn
    integer i, ii, idum, j, jj, k, maxcon, nico
    real (kind=8) dx, dy, dz, rad, radk, rr, s1, x, y, z
    real (kind=8), allocatable, dimension(:,:) :: xcoords
    logical tmp

! begin, get coordinates and allocate temporary space
    allocate (xcoords(3,natom))
    idum = matrix_get (xyz, xcoords)
    nconn = 0
    maxcon = natom
    allocate (iconn(2,maxcon))

! loop over primary atom i, get covalent radius and coordinates
    do i = 1,natom
       rad = covrad(types(i))
       if (ctrl%printl.ge.4) write (stdout,'(5x,a,i4,a,i2,a,f10.4)') &
            'Atom ', i, ', type: ', types(i), ', radius: ', rad
       x = xcoords(1,i)
       y = xcoords(2,i)
       z = xcoords(3,i)

! loop over secondary atom k, get threshold distance
       do k = 1,i-1
          radk = rad + covrad(types(k))
          dx = x-xcoords(1,k)
          dy = y-xcoords(2,k)
          dz = z-xcoords(3,k)
          rr = 0.8_8 * (dx*dx + dy*dy + dz*dz)

! atoms i and k are covalently bound, connect them
          if (rr .lt. radk*radk) then
             nconn = nconn + 1

! allocate more space if required
             if (nconn.gt.maxcon) then
                call conn_grow (iconn, maxcon, maxcon+natom, maxcon)
                maxcon = maxcon + natom
             endif
             iconn(1,nconn) = i
             iconn(2,nconn) = k
             if (ctrl%printl.ge.4) write &
                  (stdout,'(5x,a,i5,a,i5,a,f10.4,a,f10.4)') &
                  'Covalent bond ', i, ' - ', k, ': ', sqrt (rr/0.8_8), &
                  ', cov. distance: ', radk
          endif
       end do ! k = 1,i-1
    end do ! i = 1,natom

! done covalent connections
! now insert connections from the branched minimum distance path

! allocate space to hold information which atom is already connected
    nico = 0
    allocate (lconn(natom))
    do i = 2,natom
       lconn(i) = .false.
    end do
    lconn(1) = .true.

! jump here if more connections need to be made
10  continue
    s1 = -1.0_8
    do i = 1,natom
       if (lconn(i)) then
          do j = 1,natom
             if (.not. lconn(j)) then
                dx = xcoords(1,i) - xcoords(1,j)
                dy = xcoords(2,i) - xcoords(2,j)
                dz = xcoords(3,i) - xcoords(3,j)
                rr = dx*dx + dy*dy + dz*dz
                if (rr.lt.s1 .or. s1.eq.-1.0_8) then
                   s1 = rr
                   jj = j
                   ii = i
                end if
             end if ! if (.not. lconn(j))
          end do ! j = 1,natom
       end if ! if (lconn(i))
    end do ! i = 1,natom

! check if shortest distance is not yet included in the connections list
    lconn(jj) = .true.
    nico = nico + 1
    tmp = not_connected (ii, jj, iconn, nconn)
    if (tmp) then
       if (ctrl%printl.ge.2) then
          write(stdout,'(5x,a,i5,i5)') 'Adding inter fragment connction ', &
               jj, ii
       end if
       nconn = nconn + 1

! allocate more space if required
       if (nconn.gt.maxcon) then
          call conn_grow (iconn, maxcon, maxcon+natom, maxcon)
          maxcon = maxcon + natom
       endif
       iconn(1,nconn) = jj
       iconn(2,nconn) = ii
    end if ! (not_connected (ii, jj, iconn, nconn))

! check if all done
    if (nico.ne.natom-1) goto 10

! free flags array and tidy up size of returned memory
    deallocate (lconn)
    call conn_grow (iconn, maxcon, nconn, nconn)
    deallocate (xcoords)

  end subroutine connect_prim

!//////////////////////////////////////////////////////////////////////////////
! Used by connect_prim only: Returns covalent radius from nuclear charge
! Currently: fixed radius for all atoms of a row
!//////////////////////////////////////////////////////////////////////////////

  function covrad (nuccg)
    real (kind=8) covrad
! args
    integer nuccg
! local vars
    integer i, irelem
    real (kind=8) radius
! begin
    irelem = 0
    do i = 1,ctrl%nrad
       if (ctrl%irad(i).eq.nuccg) irelem = i
    end do
    if (irelem.eq.0) then
       if (nuccg.le.2) then
          radius = 1.0_8 ! H..He
       elseif (nuccg.le. 10) then
          radius = 1.6_8 ! Li..Ne
       elseif (nuccg.le.18) then
          radius = 2.5_8 ! Na..Ar
       elseif (nuccg.le.36) then
          radius = 3.0_8 ! K..Kr
       elseif (nuccg.le.54) then
          radius = 4.0_8 ! Rb..Xe
       else
          radius = 5.0_8 ! Cs....
       endif
    else
       radius = ctrl%vrad(irelem)
    end if
    covrad = radius
  end function covrad

!//////////////////////////////////////////////////////////////////////////////
! Used by connect_prim only: Check if ii - jj is not contained in iconn
!//////////////////////////////////////////////////////////////////////////////

  function not_connected (ii, jj, iconn, nconn)
! retval
    logical not_connected
! args
    integer ii, jj, nconn
    integer iconn(2,nconn)
! local vars
    integer i
! begin
    not_connected = .true.
    do i = 1,nconn
       if ((iconn(1,i).eq.ii .and. iconn(2,i).eq.jj) .or. &
            (iconn(1,i).eq.jj .and. iconn(2,i).eq.ii)) then
          not_connected = .false.
          return
       endif
    end do
  end function not_connected

!//////////////////////////////////////////////////////////////////////////////
! Used by connect_prim only: Copy the connectivity into a larger list
!//////////////////////////////////////////////////////////////////////////////

  subroutine conn_grow (iconn, old_size, new_size, copy_size)
! args
    integer old_size, new_size, copy_size
    integer, pointer, dimension(:,:) :: iconn
! local vars
    integer, pointer, dimension(:,:) :: iconn_new
    integer i
! begin
    allocate (iconn_new(2,new_size))
    do i = 1,copy_size
       iconn_new(1,i) = iconn(1,i)
       iconn_new(2,i) = iconn(2,i)
    end do
    deallocate (iconn)
    iconn => iconn_new
  end subroutine conn_grow

!//////////////////////////////////////////////////////////////////////////////
! end routines used by connect_prim
!//////////////////////////////////////////////////////////////////////////////

!------------------------------------------------------------------------------
! subroutine connect_all
!
! Alternative to connect_prim and valcoor. Every atom is connected to every
! atom. No angular coordinates required.
!
! Arguments
! natom             (in)  number of atoms, integer
! nconn             (out) number of connections found
! iconn   (2,nconn) (out) pointer to connections list
!------------------------------------------------------------------------------

  subroutine connect_all (natom, nconn, iconn)

! args
    integer natom, nconn
    integer, pointer, dimension(:,:) :: iconn

! local vars
    integer i, j

! begin
    allocate (iconn(2,natom*(natom-1)/2))
    nconn = 0
    do i = 1,natom
       do j = 1,i-1
          nconn = nconn + 1
          iconn(1,nconn) = i
          iconn(2,nconn) = j
       end do
    end do
    if (ctrl%printl.ge.1) then
       write(stdout,'(/,a,/)') 'Generating total connection'
       write(stdout,'(5x,a,i5,a,/)') 'System has ', nconn, ' connections'
    end if
  end subroutine connect_all

!------------------------------------------------------------------------------
! subroutine ci_conn
!
! Check in user connections and create a connections (integer) matrix. Note
! that all sequence numbers are relative to the start of the residue at exit.
!
! Arguments:
! con    int_matrix(2,*) (out) connections matrix, relative
! nconn                  (in)  number of calculated connections
! iconn  (2,nconn)       (in)  calculated connections, relative
! nincon                 (in)  number of user connections
! incon  (2,nincon)      (in)  user connections, absolute
! start                  (in)  first atom of the residue
! finish                 (in)  last atom of the residue
!------------------------------------------------------------------------------

  subroutine ci_conn (con, nconn, iconn, nincon, incon, start, finish)

! args
    integer nconn, nincon, start, finish
    integer, dimension(2,nconn) :: iconn
    integer, dimension(2,nincon) :: incon
    type(int_matrix), pointer :: con

! local vars
    integer i, j, ninres
    integer, dimension(:,:), allocatable :: incon_tmp

! begin, allocate scratch
    allocate (incon_tmp(2,nincon))

! check if input connections are in this residue
    ninres = 0
    do j = 1,nincon
       if (incon(1,j).ge.start .and. incon(1,j).le.finish) then
          if (incon(2,j).ge.start .and. incon(2,j).le.finish) then
             if (ctrl%printl.ge.1) then
                write (stdout,'(5x,a,i4,a,i4)') 'Adding user stretch', &
                     ctrl%incon(1,j), ' - ', ctrl%incon(2,j)
             end if
             ninres = ninres + 1
             incon_tmp(1,ninres) = incon(1,j) - (start-1)
             incon_tmp(2,ninres) = incon(2,j) - (start-1)
          end if
       else
          if (ctrl%printl.ge.-1) write (stdout,'(a,i4,a,i4,a,/)') &
               'Warning: user stretch ', incon(1,j), ' - ', incon(2,j), &
               'crosses the residue boundary'
       end if
    end do

! create the matrix and check in the calculated connections
    con => int_matrix_create (2, nconn+ninres, 'connections')
    do i = 1,nconn
       j = int_matrix_set_column (con, iconn(1,i), i)
    end do

! check in the user connections
    do i = 1,ninres
       j = int_matrix_set_column (con, incon_tmp(1,i), i+nconn)
    end do
    deallocate (incon_tmp)
  end subroutine ci_conn

!------------------------------------------------------------------------------
! subroutine valcoor
!
! Alexander J Turner Oct 1997
!
! This algorithm will make a fully redundant set of interanal coordinates.
! It defines rotation about bonds but not dihedral angles.
!
! The rotational orientation of atoms should be handled using
! rotation coordinates along the connection axes.
!
! Oct. 1997
! Changed to handle Quasi style matrix objects and to generate
! redundant co-ordinates. This is a derivative of valcoor.f in
! GRACE - I.H.Williams and A.J.Turner - University of Bath
! U.K. 1997
!
! June 1999
! Rewritten object and memory management - Zurich, SB
!
! On input
! ========
! natom : number of atoms
! nconn : Number of connections
! x,y,z : simple cartesian arrays
! iconn : 2 by n array - integer: column one is connected to column two
!
! On output
! ========
! ibend : 4 by nbend  - pointer to integer array
! irots : 3 by nnrots - pointer to integer array
! x,y,z,iconn - unchanged
!
! Definition of ibend and irots
! =============================
! ibend : bends are between ibend(1,i) and ibend(2,i) about ibend(3,1)
!     ibend(1,i) is zero - no bend defined
!       : linear bend ibend(1,i)-ibend(2,i)-ibend(3,i) 
!     ibend(4,i)=-2 => relative to xy plane
!     ibend(4,i)=-1 => relative to yz plane
!
! irots : Diherdal angle 1-2-3 w.r.t 2-3-4 , auto generated dihedrals
!     will either be proper or inproper 
!     *** FOR HISTORICAL REASONS - THESE ARE OFTEN REFERRED ***
!               TO AS ROTATIONS IN THIS CODE
!  
! nrots : Number of auto-dihedrals
! nbend : Number of bends
!
! Comments:
! =========
! The algorithm uses dynamic storage for valence information
! The scratch storage is in four dynamically allocated parts
!   ivale - an array giving the valence of each atom
!   ipvle - an array of pointers into ipple
!      ipple(ipvle+i-1) points to the start of the record of
!      atom numbers in the valence of atom i
!   ipple  - the array used to store
!      the atoms in an atoms valence shell, the pointers
!      stored from ipvle point into this array 
!   ipcvl  - counters for atomic valence collection
!
! Arguments:
! ==========
! See the comments to the fields of hdlc_obj (hdlclib.f90)
! in:  nconn, iconn, natom, x, y, z
! out: nbend, ibend, nrots, irots
!------------------------------------------------------------------------------

  subroutine valcoor (nconn, iconn, natom, x, y, z, nbend, ibend, nrots, irots)

! args
    integer nconn, natom, nbend, nrots
    integer, pointer, dimension(:,:) :: iconn
    integer, pointer, dimension(:,:) :: ibend
    integer, pointer, dimension(:,:) :: irots
    real(kind=8), pointer, dimension(:) :: x, y, z

! local vars
    logical notline
    integer i, ib, ib1, ib2, ib3, icent, ileft, ip, ibp, j, k, l
    integer maxbend, maxrots
    integer, dimension(:), allocatable :: ivale, ipvle, ipcvl, ipple, index
    real(kind=8) :: angle, tx, ty, tz, ax, ay, az, ex, ey, ez, dx, dy, dz, r

! begin
    if (ctrl%printl.ge.3) then
       write (stdout,'(5X,A,/)') 'Generating primitive internal coordinates'
    end if

! get some scratch space
    allocate (ivale(natom))
    allocate (ipvle(natom))
    allocate (ipcvl(natom))

!//////////////////////////////////////////////////////////////////////////////
! compute the valence of each atom
!
! ivale stores the total valence of each atom
! ipvle(i) points to the start of the valence of atom i in ipcvl(*)
!//////////////////////////////////////////////////////////////////////////////

    do i = 1,natom
       ivale(i)=0
    end do
    do i = 1,nconn
       j = iconn(1,i)
       ivale(j) = ivale(j) + 1
       j = iconn(2,i)
       ivale(j) = ivale(j) + 1
    end do

! get the number of bends, ib counts the total number of valence records
    nbend = 0
    ib = 0
    do i = 1,natom
       ib = ib + ivale(i)
       do j = 1,ivale(i)-1
          nbend = nbend + j
       end do
    end do
    maxbend = nbend
    allocate (ipple(ib))
    allocate (ibend(4,maxbend))
    ipvle(1) = 1

! set pointers for valence records (formerly valcoor2)
    nbend = 0
    do i = 2,natom
       ipvle(i) = ipvle(i-1)+ivale(i-1)
    end do

! set counters for atomic valence collection to -1
    do i = 1,natom
       ipcvl(i) = -1
    end do

! record the valence atoms around each atom from connection list
    do i = 1,nconn      
       j = iconn(1,i)
       k = iconn(2,i)

! increment record counter for atom j
       ipcvl(j) = ipcvl(j) + 1

! acquire pointer to records for j
       l = ipvle(j) + ipcvl(j)

! set that element to the atom to which it is connected
       ipple(l) = k

! repeat for other atom
       j = iconn(2,i)
       k = iconn(1,i)

! increment record counter for j
       ipcvl(j) = ipcvl(j) + 1

! acquire pointer to records for j
       l = ipvle(j) + ipcvl(j)

! set that element to the atom to which it is connected
       ipple(l) = k
    end do

! print the valence
    if (ctrl%printl.ge.3) then
       write(stdout,'(5X,A,/)') 'Valence of atomic centres' 
       do  i=1,natom
          j=ipvle(i) 
          k=ivale(i) 
          write(stdout,'(5x,A,I5,A,8(1X,I5))') 'Atom no ', i, &
               ' has valence =', (ipple(j+l),l=0,k-1)
       end do
       write(stdout,'(/)')
    end if

!//////////////////////////////////////////////////////////////////////////////
! Generate bends
!
! Loop scanning atomic centres. It looks for bends and also puts in linear
! bends and wags etc.
!
! I scans each atom
! J scans the valence atoms around the I'th atom
! K scans the 'other' valence atoms around the I'th atom so bends are made
!
! Other variables used
!
! ib     : The number of bends allocated about a given centre
! ib1    : The valence atom from which bends are coming
! ib2    : The valence atom to which bends are going
! ib3    : The waging atom ib3 wags in plane of ib1-i-ib2
! ip     : Pointer to valence records for I'th atom
! ibp    : Pointer to the bend record that is being written 
! nbend  : Number of allocated bends in total
! maxbend: Ammount of space allocated for bends
!//////////////////////////////////////////////////////////////////////////////

    do i = 1,natom
       ib = 0
       if (ctrl%printl.ge.4) then
          write (stdout,'(7X,A,I5)') 'Scanning bends around ', i
       end if

! pointer to start of valence info of atom i and loop through atoms j
       ip = ipvle(i)
       do j = 1,ivale(i)-1

! set ib1 to the atom connected to atom i from which bends are to be made
          ib1 = ipple(ip+j-1)
          if (ctrl%printl.ge.4) then
             write (stdout,'(7X,A,I5)') 'Making bends from atom ', ib1
          end if

! loop through all atoms higher in order than atom j and make bends
          do k=j+1,ivale(i)
             nbend = nbend + 1

! detect a potential problem in the code
             if(nbend.gt.maxbend) then
                call hdlc_errflag &
                     ('Insufficient bends allocated in valcoor - code error', &
                     'abort')
             end if

! increment bond counter
             ib = ib + 1
 
! set ib2 to the atom connected to atom i to which bends are to be made
             ib2 = ipple(ip+k-1)

! ibp is a pointer to the present bend record: make bend
             ibp = nbend 
             ibend(1,ibp) = ib1
             ibend(2,ibp) = i
             ibend(3,ibp) = ib2
             ibend(4,ibp) = 0

!//////////////////////////////////////////////////////////////////////////////
! check for linearity
!//////////////////////////////////////////////////////////////////////////////

             angle=vangled(x(ib1),x(i),x(ib2),y(ib1),y(i),y(ib2), &
                  z(ib1),z(i),z(ib2))
             if (ctrl%printl.ge.4) then
                write (stdout,'(7X,A5,I5,A2,3I5,A8,F12.7)') &
                     'Bend ', nbend, ': ', ib1, i, ib2, ' angle: ', angle
             end if

! make an l_function if the centre is bi-valent & > 179 deg
             if (angle.ge.175.0_8) then
                if (ivale(i).eq.2) then
                   if (ctrl%printl.ge.2) then
                      write(stdout,'(7X,A,I5)') 'Making linear about ', i
                   end if

! test linear function for being perp to xy plane:
! translate one point in z
! take angle between point and other two points
! must be 5<angle<175 or yz plane signaled (istra(nline,4)=-1)

                   angle=vangled(x(ib1),x(ib2),x(ib2),y(ib1),y(ib2),y(ib2), &
                        z(ib1),z(ib2),z(ib2)+1.0_8)
                   if (ctrl%printl.ge.3) then
                      write (stdout,'(7X,A,F8.3)') 'Tester angle = ', angle
                   end if
                   if(angle.le.5.0_8.or.angle.ge.175.0_8) then
                      ibend(4,ibp)=-2
                   else
                      ibend(4,ibp)=-1
                   end if

! note: removed if (.not.internal) then (remove the bend) endif
                else ! if (ivale(i).eq.2) then ...
                   if (ctrl%printl.ge.3) write (stdout,'(7X,A)') &
                        'Bend ignored, linear non bivalent'
                   ibend(1,ibp) = 0
                   nbend = nbend - 1
                end if ! if (ivale(i).eq.2) then ... else
             end if ! if (angle.ge.175.0_8)
          end do ! do k=j+1,ivale(i)
       end do ! do j = 1,ivale(i)-1

!//////////////////////////////////////////////////////////////////////////////
! check for planar systems
!//////////////////////////////////////////////////////////////////////////////

       if (ivale(i).gt.2) then

! acquire pointer to start of valence info for this atom
          ip=ipvle(i)
          tx=1.0_8
          ty=1.0_8
          tz=1.0_8

! loop through all valence atoms dotting the perps to the normalised vectors
          do j = 1,ivale(i)-1

! set ib1 to the atom connected to atom i from which bends are to be made
             ib1 = ipple(ip+j-1)
             ax = x(ib1) - x(i)
             ay = y(ib1) - y(i)
             az = z(ib1) - z(i)
             r = sqrt(ax*ax + ay*ay + az*az)
             ax = ax/r
             ay = ay/r
             az = az/r
             ib1 = ipple(ip+j)
             ex = x(ib1) - x(i)
             ey = y(ib1) - y(i)
             ez = z(ib1) - z(i)
             r = sqrt (ex*ex + ey*ey + ez*ez)
             ex = ex/r
             ey = ey/r
             ez = ez/r
             dx = (Ay*Ez - Az*Ey)
             dy = (Az*Ex - Ax*Ez)
             dz = (Ax*Ey - Ay*Ex)
             r = sqrt (dx*dx + dy*dy + dz*dz)
             tx = tx*dx/r
             ty = ty*dy/r
             tz = tz*dz/r
          end do
          if (ctrl%printl.ge.3) then
             write (stdout,'(7X,A,I5,A,F12.8)') 'Dot-product about ', i, &
                  ' is ', tx+ty+tz
          end if

! if dot greater than 5% off 1 put in dihedral instead of bend
          if (abs(tx+ty+tz).gt.0.95_8) then
             if (ctrl%printl.ge.2) then
                write (stdout,'(7X,A)')'System near planar - making improper'
             end if

! construct improper from last bend - i.e. nbend
             k = -1
             do j = 1,ivale(i)
                ib1 = ipple(ip+j-1)
                if (ib1.ne.ibend(1,nbend) .and. ib1.ne.ibend(3,nbend)) &
                     k = ib1
             end do
             if (k.eq.-1) then
                call hdlc_errflag ('Failed to contruct internals','abort')
             end if
             ibend(4,nbend) = k

! swap 1 and 3 if 2-3-4 is co-linear
             dx = x(ibend(2,nbend)) - x(ibend(3,nbend))
             dy = y(ibend(2,nbend)) - y(ibend(3,nbend))
             dz = z(ibend(2,nbend)) - z(ibend(3,nbend))
             r = sqrt(dx*dx + dy*dy + dz*dz)
             tx = dx/r
             ty = dy/r
             tz = dz/r
             dx = x(ibend(3,nbend)) - x(ibend(4,nbend))
             dy = y(ibend(3,nbend)) - y(ibend(4,nbend))
             dz = z(ibend(3,nbend)) - z(ibend(4,nbend))
             r = sqrt(dx*dx + dy*dy + dz*dz)
             tx = tx*dx/r
             ty = ty*dy/r
             tz = tz*dz/r
             if (tx+ty+tz.gt.0.095_8) then
                k = ibend(3,nbend)
                ibend(3,nbend) = ibend(1,nbend)
                ibend(1,nbend) = k
             end if
          end if ! (abs(tx+ty+tz).gt.0.95_8)
       end if ! (ivale(i).gt.2)
    end do ! i = 1,natom (end of valcoor2)

!//////////////////////////////////////////////////////////////////////////////
! generate dihedrals and improper dihedrals
!//////////////////////////////////////////////////////////////////////////////

! initialise (formerly valcoor3)
    deallocate (ipple)
    allocate (index(natom))
    maxrots = nbend*3
    allocate (irots(4,maxrots))
    if (ctrl%printl.ge.5) write (stdout,'(7X,A,I5)') 'Guessing maxrots at ', maxrots
    do i = 1,natom
       index(i) = -1
    end do

! set index(i) to point to records of bends about i
    k = 0
    do i = 1,nbend
       if (ibend(2,i).ne.k) then
          k = ibend(2,i)
          if(ibend(1,i).ne.0) index(k) = i
       end if
    end do

! sort all bends so if one end is univalent and the other bivalent, the
! bivalent is the far end
    do i = 1,nbend
       if (index(ibend(3,i)).eq.-1 .and. ibend(4,i).ge.0) then
          k = ibend(3,i) 
          ibend(3,i) = ibend(1,i)
          ibend(1,i) = k
       end if
    end do
    do i = 1,natom
       index(i) = -1
    end do
    k = -1
    do i = 1,nbend
       if (ibend(2,i).ne.k) then
          k = ibend(2,i)
          if (ibend(1,i).ne.0) index(k) = i
       end if
    end do
!WT write (stdout,'(1X,A)') ' break 8 '
!WT do i = 1,nbend
!WT    write (stdout,'(1X,A,5I5)') ' i,ibend(1-4,i) ', i, (ibend(j,i),j=1,4)
!WT end do
!WT do i = 1,natom
!WT    write (stdout,'(1X,A,2I5)') ' i,index(i) ', i, index(i)
!WT end do

! main loop through bends constructing dihedrals
    nrots = 0
    do i = 1,nbend
       if (ibend(4,i).lt.0) goto 300 ! skip if linear
!WT write (stdout,'(1X,A,2I5)') ' nbend,i ', nbend,i

! flag this is a new bend and not due to expansion of a dihedral
       notline = .true.

! k points to the bend record of the far end of bend i
       ileft = ibend(2,i)
       icent = ibend(3,i)
       k = index (icent)

! follow linears until a real bend (loop point for a 'do until' type loop)
100    continue

! if the atom ibend(3,i) is univalent, you cannot make a dihedral -> break
          if (k.eq.-1) goto 300

! expand the dihedral (follow the bond) if atom is linear bivalent
          if (ibend(4,k).lt.0) then

! test for three membered ring, break if found
             if ((ibend(3,i).eq.ibend(3,k) .or. ibend(1,i).eq.ibend(1,k) .or. &
                  ibend(1,i).eq.ibend(3,k) .or. ibend(3,i).eq.ibend(1,k)) &
                  .and. notline) goto 300 
             notline=.false.

! test which way around the bend is
          ileft = ibend(2,k)
          if (ileft.ne.max(ibend(1,k),ibend(3,k))) then
             icent = max(ibend(1,k),ibend(3,k))
             k = index(icent)
             goto 100
          else
             icent = min(ibend(1,k),ibend(3,k))
             k = index(icent)
             goto 100
          end if
          end if ! (ibend(4,k).lt.0)

! cycle over all bends around icent
200    do while (ibend(2,k) .eq. icent)

! test for three membered ring, break if found
          if ((ibend(3,i).eq.ibend(3,k) .or. ibend(1,i).eq.ibend(1,k) .or. &
               ibend(1,i).eq.ibend(3,k) .or. ibend(3,i).eq.ibend(1,k)) .and. &
               notline) goto 300 

! test added by Martin Graf
! test for identical central atom on creating dihedral MG
          if (ibend(2,k).eq.ibend(2,i)) goto 300

! test if bend k is valid and matching bend i
          if (ibend(1,k).ne.0 .and. ibend(4,k).eq.0 .and. &
              (ibend(1,k).eq.ileft .or. ibend(3,k).eq.ileft)) then
             nrots=nrots+1

! allocate more space if required
             if (nrots.gt.maxrots) then
                call rots_grow (irots, maxrots, maxrots+nbend*9, maxrots)
                maxrots = maxrots + nbend*9
             end if

! store the dihedral (formerly valcoor4)
             if (ibend(1,k).eq.ileft) then
                irots(4,nrots) = ibend(3,k)
             else
                irots(4,nrots) = ibend(1,k)
             end if
             do j = 1,3
                irots(j,nrots) = ibend(j,i)
             end do

! to ensure identical dihedrals are not included
             do j = 1,nrots-1
                if (irots(1,j).eq.irots(4,nrots)) then
                   if(irots(2,j).eq.irots(3,nrots) .and. &
                        irots(3,j).eq.irots(2,nrots)) then
                      nrots = nrots-1
                   end if
                end if
             end do
          end if ! valid and matching bend

! end cycle over all bends around icent
          k = k + 1
          if (k > maxbend) exit
       end do

! loss of indentation due to end of loop to 200
300    continue
    end do ! do i = 1,nbend (end of valcoor3)

! deallocate scratch space
    deallocate (ivale)
    deallocate (ipvle)
    deallocate (ipcvl)

  end subroutine valcoor

!//////////////////////////////////////////////////////////////////////////////
! Used by valcoor only: Copy the torsions list into a larger list
!//////////////////////////////////////////////////////////////////////////////

  subroutine rots_grow (irots, old_size, new_size, copy_size)
! args
    integer old_size, new_size, copy_size
    integer, pointer, dimension(:,:) :: irots
! local vars
    integer, pointer, dimension(:,:) :: irots_new
    integer i
! begin
    allocate (irots_new(4,new_size))
    do i = 1,copy_size
       irots_new(1,i) = irots(1,i)
       irots_new(2,i) = irots(2,i)
       irots_new(3,i) = irots(3,i)
       irots_new(4,i) = irots(4,i)
    end do
    deallocate (irots)
    irots => irots_new
  end subroutine rots_grow

!//////////////////////////////////////////////////////////////////////////////
! end routines used by valcoor
!//////////////////////////////////////////////////////////////////////////////

!------------------------------------------------------------------------------
! subroutine valcoor_print
!
! Arguments:
! See the comments to the fields of hdlc_obj (hdlclib.f90)
! in:  nconn, iconn, natom, x, y, z, nbend, ibend, nrots, irots
!------------------------------------------------------------------------------

  subroutine valcoor_print (nconn, iconn, nbend, ibend, nrots, irots, &
       x, y, z, natom)

! args
    integer nconn, nbend, nrots, natom
    integer iconn(2,nconn), ibend(4,nbend), irots(4,nrots)
    real(kind=8) :: x(natom), y(natom), z(natom)

! local vars
    integer i, j, k
    real(kind=8) :: dx, dy, dz, r

! begin
    if (ctrl%printl.ge.2) then
       write (stdout,'(/,3X,A)') 'The Cartesian coordinates'
       do i = 1,natom
          write (stdout,'(5X,I5,1X,F20.14,1X,F20.14,1X,F20.14)') &
               i, x(i), y(i), z(i)
       end do
    end if

    if (ctrl%printl.ge.1) then
       write (stdout,'(/,3X,A)') 'The primitive internal coordinates'
       k = 0

       if (ctrl%printl.ge.2) write (stdout,'(/)')
       write (stdout,'(5X,A,I5,A)') 'The system has ', nconn, ' stretches'
       if (ctrl%printl.ge.2) then
          do i = 1,nconn
             k = k + 1
             dx = x(iconn(1,i)) - x(iconn(2,i))
             dy = y(iconn(1,i)) - y(iconn(2,i))
             dz = z(iconn(1,i)) - z(iconn(2,i))
             r = sqrt (dx*dx + dy*dy + dz*dz)
             write (stdout,'(5X,i5,A,2(I5,1X),12X,A,F13.8)') k, ' Stre ', &
                  (iconn(j,i), j=1,2), ' = ', r
          end do
          write (stdout,'(/)')
       end if ! (ctrl%printl.ge.2)

       write (stdout,'(5X,A,I5,A)') 'The system has ', nbend, &
            ' bends and impropers'
       if (ctrl%printl.ge.2) then
          do i = 1,nbend
             k = k + 1
             if (ibend(4,i).eq.0) then
                r = vangled(x(ibend(1,i)),x(ibend(2,i)),x(ibend(3,i)), &
                     y(ibend(1,i)),y(ibend(2,i)),y(ibend(3,i)), &
                     z(ibend(1,i)),z(ibend(2,i)),z(ibend(3,i)))
                write (stdout,'(5X,I5,A,3(I5,1X),6X,A,F13.8)') k, ' Bend ', &
                     (ibend(j,i), j=1,3), ' = ', r
             elseif (ibend(4,i).eq.-1) then
                write (stdout,'(5X,I5,A,3(I5,1X),A)') k, ' Line ', &
                     (ibend(j,i), j=1,3), ' wrt xy'
             elseif(ibend(4,i).eq.-2) then
                write (stdout,'(5X,I5,A,3(I5,1X),A)') k, ' Line ', &
                     (ibend(j,i), j=1,3), ' wrt yz'
             else
                r = vdihedrald(x(ibend(1,i)),y(ibend(1,i)),z(ibend(1,i)), &
                     x(ibend(2,i)),y(ibend(2,i)),z(ibend(2,i)), &
                     x(ibend(3,i)),y(ibend(3,i)),z(ibend(3,i)), &
                     x(ibend(4,i)),y(ibend(4,i)),z(ibend(4,i)))
                write (stdout,'(5X,I5,A,4(I5,1X),A,F13.8)') k, ' Impr ', &
                     (ibend(j,i), j=1,4), ' = ', r
             end if ! (ibend(4,i).eq.0)
          end do ! i = 1,nbend
          write (stdout,'(/)')
       end if ! (ctrl%printl.ge.2)

       write (stdout,'(5X,A,I5,A)') 'The system has ', nrots, ' dihedrals'
       if (ctrl%printl.ge.2) then
          do i = 1,nrots
             k = k + 1
             r = vdihedrald(x(irots(1,i)),y(irots(1,i)),z(irots(1,i)), &
                  x(irots(2,i)),y(irots(2,i)),z(irots(2,i)), &
                  x(irots(3,i)),y(irots(3,i)),z(irots(3,i)), &
                  x(irots(4,i)),y(irots(4,i)),z(irots(4,i)))
             write (stdout,'(5X,I5,A,4(I5,1X),A,F13.8)') k, ' Dihe ', &
                  (irots(j,i), j=1,4), ' = ', r
          end do
          write (stdout,'(/)')
       end if ! (ctrl%printl.ge.2)
    end if ! (ctrl%printl.ge.1)
    call flushout()
  end subroutine valcoor_print

!------------------------------------------------------------------------------
! Utilities by Alex Turner (originally dlc_angle_utils)
!------------------------------------------------------------------------------

  subroutine rots_dlc(noint,ii,b,ib,c)
!
! Alexander J Turner
!
! Diagram:
!
!   p
!   |\ A
! C | \
!   s--q-----r
!    D   E
!
!  What it does:
!
! Computes the cartesian force on  A
! given a unit moment about axis   E
!
! It is intended for the generation of B matrix elements
!
! Method:
!
! Force on p is perpendicular to plane r,q,p
! Force on p is proportional to |C|
!
! C = A - (A.E x E/|E|)
!
! Direction perpendicular to p = ((AyEz-AzEy),(AzEx-AxEz),(AxEy-AyEx))
! Call this vector Z
!
! P.S. The use of case in this subroutine is just for human reading
!      Uppercase => vector , lowercase => scalar or point
!
    integer i,ib(4,*),iaind,jaind,kaind,m
    integer noint,ii
!
    real(kind=8) :: b(3,4,*),c(*),dzero,done,eps
    real(kind=8) :: dx,dy,dz
    real(kind=8) :: Ex,Ey,Ez
    real(kind=8) :: Tx,Ty,Tz
    real(kind=8) :: Ax,Ay,Az
    real(kind=8) :: Cx,Cy,Cz
    real(kind=8) :: Zx,Zy,Zz
    real(kind=8) :: magAE,magE,fact,magC
    parameter(eps=1.0e-6_8)
!
    data dzero/0.0_8/,done/1.0_8/
!
! Set up position in B matrix pointers
! Zero bmat for premeture return
!
    ib(1,noint)=ii
    ib(2,noint)=0
    ib(3,noint)=0
    ib(4,noint)=0
    b(1,1,noint)=dzero
    b(2,1,noint)=dzero
    b(3,1,noint)=dzero
!
! Get Vectors
!
    Ex=c(4)-c(7)
    Ey=c(5)-c(8)
    Ez=c(6)-c(9)

    Ax=c(1)-c(4)
    Ay=c(2)-c(5)
    Az=c(3)-c(6)
!
! Get A.E / |E|
!
    magE =sqrt(Ex*Ex+Ey*Ey+Ez*Ez)
    magAE=Ax*Ex+Ay*Ey+Az*Ez
    fact =-1.0_8 * magAE/magE
!
! If E is small - return nothing
!
    if(abs(magE).lt.eps) return
!
! Scale E by factor
!
    Tx=Ex*fact
    Ty=Ey*fact
    Tz=Ez*fact
!
! Get C
!
    Cx=Ax-Tx
    Cy=Ay-Ty
    Cz=Az-Tz
!
! Get magnitude of C - and thus force
!
    magC=sqrt(Cx*Cx+Cy*Cy+Cz*Cz)
!
! If C is small - angle linear - zero force, or numerical problems can occur
!
    if(magC.lt.eps) return
!
! Get Z
!
    Zx=(Ay*Ez-Az*Ey)
    Zy=(Az*Ex-Ax*Ez)
    Zz=(Ax*Ey-Ay*Ex)
!
! Normalise and times by magC
!
    fact=magC/sqrt(Zx*Zx+Zy*Zy+Zz*Zz)
    Zx=Zx*fact
    Zy=Zy*fact
    Zz=Zz*fact
!
! Put it into B
!
    b(1,1,noint)=Zx
    b(2,1,noint)=Zy
    b(3,1,noint)=Zz
!
! All done
!
    return
  end subroutine rots_dlc
!
! ------------------------------------------------------------
!
  function dlc_free_rot ( &
       pxa, pxb, pxc,  pya, pyb, pyc,  pza, pzb, pzc, &
       poxa,poxb,poxc, poya,poyb,poyc, poza,pozb,pozc)
    real(kind=8) :: dlc_free_rot
!
! Alexander J Turner
! Feb 1998
!
! Compute the magnitude of the rotation of an atom about the axis of a bond,
! given that the vector of the bond is constant
!
    real(kind=8) :: xa,xb,xc,ya,yb,yc,za,zb,zc, &
         oxa,oxb,oxc,oya,oyb,oyc,oza,ozb,ozc,dx,dy,dz, &
         s,pxa,pya,pza,pxb,pyb,pzb,pxc,pyc,pzc, &
         poxa,poya,poza,poxb,poyb,pozb,poxc,poyc,pozc
    real(kind=8) :: dzero,theta1,theta2,theta3,theta4
    parameter(dzero=0.0_8)

!
! Move by refernece pxa ... to by value oxa ....
!
    xa=pxa
    xb=pxb
    xc=pxc
    ya=pya
    yb=pyb
    yc=pyc
    za=pza
    zb=pzb
    zc=pzc
    oxa=poxa
    oxb=poxb
    oxc=poxc
    oya=poya
    oyb=poyb
    oyc=poyc
    oza=poza
    ozb=pozb
    ozc=pozc
!
! move A,B,C by B->origin
!
    dx=-oxb
    dy=-oyb
    dz=-ozb

    oxa=oxa-dx
    oya=oya-dy
    oza=oza-dz

    oxc=oxc-dx
    oyc=oyc-dy
    ozc=ozc-dz

    oxb=dzero
    oyb=dzero
    oyb=dzero

    xa=xa-dx
    ya=ya-dy
    za=za-dz

    xc=xc-dx
    yc=yc-dy
    zc=zc-dz

    xb=dzero
    yb=dzero
    yb=dzero
!
! Translate A' B'->C'
!
! oxa=oxa-xc
! oya=oya-yc
! oza=oza-zc
!
    dlc_free_rot=-vdihedral(xa,ya,za, oxc,oyc,ozc, oxb,oyb,ozb, oxa,oya,oza)
!      if(abs(dlc_free_rot).gt.1.0d-1) then
!         write(stdout,*)
!         write(stdout,*)'Dihed ',dlc_free_rot
!         write(stdout,*)'XY Theta: ',theta1
!         write(stdout,*)'XZ Theta: ',theta2
!         write(stdout,*)'Old c ',oxc,oyc,ozc
!         write(stdout,*)'New c',xc,yc,zc
!         write(stdout,*)'Old b ',oxb,oyb,ozb
!         write(stdout,*)'New b',xb,yb,zb
!         write(stdout,*)'Old a ',oxa,oya,oza
!         write(stdout,*)'New a',xa,ya,za
!      endif
!
    return
  end function dlc_free_rot
!
! ------------------------------------------------------------
!
  function vdihedral(xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl)
    real(kind=8) :: vdihedral
    real(kind=8) :: XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,XL,YL,ZL
    real(kind=8) :: FX,FY,FZ,FR2,FR,HX,HY,HZ,HR2,HR,EX,EY,EZ,ER2,ER
    real(kind=8) :: RIK,TH1,TH2,TH3,PHI,GX,GY,GZ,GR2,GR,GRR,GXR,GYR,GZR
    real(kind=8) :: FRR,FXR,FYR,FZR,CST,HRR,HXR,HYR,HZR,ERR,EXR,EYR,EZR
    real(kind=8) :: AX,AY,AZ,BX,BY,BZ,RA2,RB2,RA,RB,RAR,RBR,AXR,AYR,AZR
    real(kind=8) :: BXR,BYR,BZR,CX,CY,CZ
    real(kind=8) :: cut,one,cnv
    parameter(cnv=0.017453292519943295_8)
    parameter(cut=1.0e-8_8)
    parameter(one=1.0_8)
    vdihedral=0.0_8
!
! BOND DISTANCES
!
    FX=XI-XJ
    FY=YI-YJ
    FZ=ZI-ZJ
    FR2=FX*FX + FY*FY + FZ*FZ
    if(fr2.lt.cut) return
!
    HX=XL-XK
    HY=YL-YK
    HZ=ZL-ZK
    HR2=HX*HX + HY*HY + HZ*HZ
    if(hr2.lt.cut) return
!
    EX=XI-XK
    EY=YI-YK
    EZ=ZI-ZK
    ER2=EX*EX + EY*EY + EZ*EZ
    if(er2.lt.cut) return
!
    GX=XJ-XK
    GY=YJ-YK
    GZ=ZJ-ZK
    GR2=GX*GX+GY*GY+GZ*GZ
    if(gr2.lt.cut) return
!
! DIHEDRAL VALUES
!
! AX perp to F-G
! BX perp to H-G
!
    AX=FY*GZ-FZ*GY
    AY=FZ*GX-FX*GZ
    AZ=FX*GY-FY*GX
    BX=HY*GZ-HZ*GY
    BY=HZ*GX-HX*GZ
    BZ=HX*GY-HY*GX
!
    RA2=AX*AX+AY*AY+AZ*AZ
    RB2=BX*BX+BY*BY+BZ*BZ
    RA=SQRT(RA2)
    RB=SQRT(RB2)
    IF(RA.LE.0.01_8) THEN
!      write(stdout,'(a)')'Warning: Dihedral near linear'
!      write(stdout,'(4(x,g13.5))')xi,yi,zi
!      write(stdout,'(4(x,g13.5))')xj,yj,zj
!      write(stdout,'(4(x,g13.5))')xk,yk,zk
!      write(stdout,'(4(x,g13.5))')xl,yl,zl
!      RA=0.01_8
       return
    ENDIF
    RAR=1.0_8/RA
    IF(RB.LE.0.01_8) THEN
!       write(stdout,'(a)')'Warning: Dihedral near linear'
!       RB=0.01_8
       return
    ENDIF
    RBR=1.0_8/RB
!
! Normalise
!
    AXR=AX*RAR
    AYR=AY*RAR
    AZR=AZ*RAR
    BXR=BX*RBR
    BYR=BY*RBR
    BZR=BZ*RBR
!
    CST=AXR*BXR+AYR*BYR+AZR*BZR
    IF(ABS(CST).GE.1.0_8) CST=SIGN(ONE,CST)
    PHI=ACOS(CST)
    CX=AYR*BZR-AZR*BYR
    CY=AZR*BXR-AXR*BZR
    CZ=AXR*BYR-AYR*BXR
    IF(GX*CX+GY*CY+GZ*CZ.GT.0.0_8) PHI=-PHI
    vdihedral=phi
!   vdihedral=0.0_8
!   write(stdout,*)phi
  end function vdihedral
!
! ------------------------------------------------------------
!
  function vangled (x1,x2,x3,y1,y2,y3,z1,z2,z3)
    real(kind=8) :: vangled
    real(kind=8) :: d21,d23,d13,dacosd
    real(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(kind=8) :: two,eps,pi,sd21,sd23,sd13,cnv
    parameter(two=2.0_8)
    parameter(eps=1.0e-6_8)
!
! Alexander J Turner - June 1997
! Used by valence, it returns the angle 1 and 3 about 2
! by cosine law in degrees
!
    d21=(x1-x2)**two + (y1-y2)**two + (z1-z2)**two
    d23=(x3-x2)**two + (y3-y2)**two + (z3-z2)**two
    d13=(x1-x3)**two + (y1-y3)**two + (z1-z3)**two
!
! Allow for angle near 0 deg not causing numerical problems
!
    if(d21.lt.eps.or.d23.lt.eps.or.d21.lt.eps) then
       vangled=0.0_8
       return
    endif
!
! Allow for angle near 180 deg not causing numberical problems
!
    sd21=sqrt(d21)
    sd23=sqrt(d23)
    sd13=sqrt(d13)
!
    if(abs(sd13-sd21-sd23).lt.eps) then
       vangled=180.0_8
    else
       vangled=57.295779513082316_8*acos((d21 + d23 - d13)/(2.0_8*sd21*sd23))
    endif
!
    return
  end function vangled
!
! ------------------------------------------------------------
!
  function vangle (x1,x2,x3,y1,y2,y3,z1,z2,z3)
    real(kind=8) :: vangle
    real(kind=8) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(kind=8) :: cnv
    parameter(cnv=0.17453292519943295e-1_8)
!
    vangle=vangled(x1,x2,x3,y1,y2,y3,z1,z2,z3)
    vangle=cnv * vangle
    return
  end function vangle
!
! ------------------------------------------------------------
!
  function vdihedrald(xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl)
    real(kind=8) vdihedrald
    real(kind=8) :: XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,XL,YL,ZL
    real(kind=8) :: cnv
    parameter(cnv=57.295779513082316_8)
!
    vdihedrald=vdihedral(xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl)
    vdihedrald=vdihedrald*cnv
    return
  end function vdihedrald
end module primitive
