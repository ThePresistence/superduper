module random
! for langevin dynamics calling 
  implicit none
  private ran_state

  integer, save :: n_seed = 0
  integer, save :: n_seed0 = 0
  integer, save, dimension(100) :: ran_state, ran_state0
  double precision, save, dimension(10000) :: ran_all
  integer, save :: n_ran = 0

contains
  
  subroutine random_gauss2d(X, Y)
    implicit none
    double precision, intent(out) :: X, Y
    double precision :: PI, R1, R2, R
    ! define pi value
    PI=4.D0*DATAN(1.D0)
    call random_number(R1)
    call random_number(R2)
    R1 = DSQRT(-2.0*DLOG(R1))
    R2 = 2.0 * PI * R2
    X  = R1*COS(R2)
    Y  = R1*SIN(R2)
  end subroutine random_gauss2d

  subroutine get_random_state(iflag)
    implicit none
    integer :: iflag
    if (iflag == 0) then
       !call random_seed(size=n_seed)
       call random_seed(get=ran_state)    
    else if (iflag == 1) then
       call random_seed(size=n_seed)
       call random_seed(get=ran_state0)   
    else
       write(*,*) "no such option"
       stop
    endif
  end subroutine get_random_state

  subroutine set_random_state(iflag)
    implicit none
    integer :: iflag
    if (iflag == 0) then
       call random_seed(put=ran_state(1:n_seed))
    else if (iflag == 1) then
       call random_seed(put=ran_state0(1:n_seed))
    else
       write(*,*) "no such option"
       stop
    endif
  end subroutine set_random_state


  subroutine set_random_state2(seed, n_seed)
    implicit none
    integer :: n_seed
    integer, dimension(:) :: seed
    call random_seed(put=seed(1:n_seed))
  end subroutine set_random_state2


  ! in which \mu is zero, and \deta is 1.
  subroutine ran_num(ran, nr)
    implicit none
    integer :: nr
    double precision, dimension(:) :: ran
    integer :: i
    double precision :: num

    call get_random_state(1)
    ! set seed state
    call set_random_state(0)
    
    ! get one gaussian random number
    do i = 1, nr
       call get_std_normal_one(num)
       ran(i) = num
    enddo

    ! get current state
    call get_random_state(0)
    ! eliminate the effect of this subroutine
    call set_random_state(1)

  end subroutine ran_num



  subroutine get_std_normal_one(num)
    implicit none
    integer, save :: iset = 0
    double precision, save :: saveY
    double precision, intent(out) :: num
    double precision :: X, Y
    if (iset .eq. 0) then
       call random_gauss2d(X, Y)
       num = X
       saveY = Y
       iset = 1
    else
       num = saveY
       iset = 0
    endif

  end subroutine get_std_normal_one



  ! generate all value at once & store in memory
  subroutine get_std_normal_all(n_ran)
    implicit none
    integer :: i
    integer :: n_ran
    double precision :: X, Y
 
    call get_random_state(0)
    do i = 1, n_ran/2+1
       call random_gauss2d(X, Y)
       ran_all(i) = X
       ran_all(i+1) = Y
    enddo

    call set_random_state(0)

  end subroutine get_std_normal_all


  ! generate all value at once & store in memory
  subroutine dump_std_normal_all(n_ran)
    implicit none
    integer :: i
    integer :: n_ran
    double precision :: X, Y
    !allocate(ran_all(1:n_ran))

    call get_random_state(0)
    write(*,*) n_ran
    do i = 1, n_ran/2+1
       call random_gauss2d(X, Y)
       write(*,*) X
       write(*,*) Y
    enddo

    call set_random_state(0)

  end subroutine dump_std_normal_all


end module random

