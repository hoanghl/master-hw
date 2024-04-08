!--------------------------------------------------------!
!                                                        !
! 2D molecular dynamics simulation code for course       !
! MATR326 Tools for high performance computing           !
! Antti Kuronen, University of Helsinki                  !
!                                                        !
! Atoms set in a square lattice (lattice constant 1)     !
! and interact with harmonic (E~k2×Δr) and anharmonic    !
! (E~k3×Δr) potentials. All masses are one.              !
!                                                        !
! Try e.g.                                               !
! ./a.out 100 0.001 100000 0.5 100 0                     !
! This should result in the following average energies   !
! (skipping the initial transient)                       !
!   <Etot>=209.77                                        !
!   <Epot>=103.53                                        !
!   <Ekin>=106.23                                        !
!                                                        !
!--------------------------------------------------------!

module sizes
  integer,parameter :: rk0=selected_real_kind(5,10)    ! single precision
  integer,parameter :: rk=selected_real_kind(10,40)    ! double precision
  integer,parameter :: rk2=selected_real_kind(30,200)  ! quadruple precision
  integer,parameter :: MAXBUF=200
end module sizes

! --------------------------------------------------------------------------------------

program md2d
  use sizes
  implicit none  
  
  type :: vector
     real(rk) :: x
     real(rk) :: y
  end type vector

  integer,parameter :: NNEIGH=4
  
  type(vector),allocatable :: x(:)   ! atom positions
  type(vector),allocatable :: v(:)   !      velocities
  type(vector),allocatable :: v0(:)  !      previous veloocities (leap frog needs them)
  type(vector),allocatable :: a(:)   !      accelerations
  real(rk),allocatable :: ep(:)      !      potential energies
  real(rk),allocatable :: ek(:)      !      kinetic energies
  integer,allocatable :: neigh(:,:)  !      neghbor list

  real(rk) :: epsum,eksum            ! system energies

  real(rk) :: dt                     ! time step 
  real(rk) :: vsc                    ! mean initial velocity
  real(rk) :: box                    ! system size
  real(rk) :: k2                     ! harmonic 'spring constant'
  real(rk) :: k3                     ! anharmonic 'spring constant'
  integer :: nuc                     ! number of unitcells in one direction
  integer :: nat                     ! number of atoms
  integer :: maxt                    ! number of time steps simulated
  integer :: eout                    ! energy output interval
  integer :: cout                    ! coordinate output interval (lot of data, beware!)

  ! Misc variables
  integer :: t0,t1,wct,rate
  integer :: i,n,ia
  character(len=MAXBUF) :: arg,cmdline,exe
  real(rk) :: vx,vy
  logical :: debug=.true.

  ! Get number of atoms, time step and simulation length from command line

  ia=command_argument_count()
  if (ia<4.or.ia>8) then
     call get_command_argument(0,arg)
     write(6,'(a,a,a)')   'usage: ',trim(arg),' nat dt maxt vsc [eout [cout [k2 [k3]]]]'
     write(6,'(a)')       '    nuc  = number of ''unit cells'' in one direction'
     write(6,'(a)')       '    dt   = time step'
     write(6,'(a)')       '    maxt = number of time steps in simulation'
     write(6,'(a)')       '    vsc  = mean velocity of atoms in the beginning (''temperature'')'
     write(6,'(a)')       '    eout = interval for printing energies to stdout'
     write(6,'(a)')       '    cout = interval for printing coordinates to ''fort.10'''
     write(6,'(a)')       '    k2   = harmonic potential ''spring constant'''
     write(6,'(a)')       '    k3   = anrmonic potential ''spring constant'''
     stop
  end if

  ! Some default values in case they are not given on command line

  cout=0
  eout=1
  k2=1.0
  k3=0.1
  
  cmdline=''
  call get_command_argument(0,arg)                      ; cmdline=trim(cmdline)//' '//arg
  call get_command_argument(1,arg); read(arg,*) nuc     ; cmdline=trim(cmdline)//' '//arg
  call get_command_argument(2,arg); read(arg,*) dt      ; cmdline=trim(cmdline)//' '//arg
  call get_command_argument(3,arg); read(arg,*) maxt    ; cmdline=trim(cmdline)//' '//arg
  call get_command_argument(4,arg); read(arg,*) vsc     ; cmdline=trim(cmdline)//' '//arg
  if (ia>4) then
     call get_command_argument(5,arg); read(arg,*) eout ; cmdline=trim(cmdline)//' '//arg
  end if
  if (ia>5) then
     call get_command_argument(6,arg); read(arg,*) cout ; cmdline=trim(cmdline)//' '//arg
  end if
  if (ia>6) then
     call get_command_argument(7,arg); read(arg,*) k2   ; cmdline=trim(cmdline)//' '//arg
  end if
  if (ia>7) then
     call get_command_argument(8,arg); read(arg,*) k3   ; cmdline=trim(cmdline)//' '//arg
  end if

  print '(a,a)','# Command line: ',trim(cmdline)
  print '(a)','# Columns: time, total energy, potential energy, kinetic energy'

  nat=nuc*nuc  ! Input parameter nuc is the number of atoms in one row
  allocate(x(nat),v(nat),v0(nat),a(nat),ep(nat),ek(nat),neigh(nat,NNEIGH))

  ! Initialize atoms positions and give them random velocities

  box=nuc
  x%x=[(i/nuc,i=0,nat-1)]
  x%y=[(mod(i,nuc),i=0,nat-1)]
  call random_number(v%x)
  call random_number(v%y)
  ! Scale the velocities to vsc*[-½,½]
  v%x=vsc*(v%x-0.5) 
  v%y=vsc*(v%y-0.5)
  ! Remove center of mass velocity
  v%x=(v%x-sum(v%x)/nat) 
  v%y=(v%y-sum(v%y)/nat)

  ! Create the (fixed) neighbor list

  call get_neigh_list(neigh)
  
  ! If the user wants to calculate initial energy and print initial coords

  if (cout>0) then
     do i=1,nat
        call accel(i,ep(i),a(i))
     end do
     n=0
     call printcoords()
  end if


  ! Simulation proper

  call system_clock(count_rate=rate)
  call system_clock(count=t0)

  time_loop: do n=1,maxt

     v0=v

     atom_loop1: do i=1,nat
        ! New potential energy and acceleration
        call accel(i,ep(i),a(i))
     end do atom_loop1

     atom_loop2: do i=1,nat
        ! Leap frog integration algorithm: update position and velocity
        v(i)%x=v(i)%x+dt*a(i)%x
        v(i)%y=v(i)%y+dt*a(i)%y
        x(i)%x=x(i)%x+dt*v(i)%x
        x(i)%y=x(i)%y+dt*v(i)%y
        ! Check periodic boundary conditions
        if (x(i)%x<0.0 ) x(i)%x=x(i)%x+box
        if (x(i)%y<0.0 ) x(i)%y=x(i)%y+box
        if (x(i)%x>=box) x(i)%x=x(i)%x-box
        if (x(i)%y>=box) x(i)%y=x(i)%y-box
        ! Calculate kinetic energy (note: mass=1)
        vx=(v0(i)%x+v(i)%x)/2.0
        vy=(v0(i)%y+v(i)%y)/2.0
        ek(i)=1.0/2.0*(vx**2+vy**2)
     end do atom_loop2

     ! Calculate and print total potential end kinetic energies
     ! and their sum that should be conserved.
     epsum=sum(ep)
     eksum=sum(ek)
     if (eout>0) then
        if (mod(n,eout)==0) print '(4g20.10)',dt*n,epsum+eksum,epsum,eksum
     end if
     if (cout>0) then
        if (mod(n,cout)==0) call printcoords()
     end if

  end do time_loop

  call system_clock(count=t1)
  write(0,'(a,g16.8,a)') 'Wall clock time: ',real(t1-t0,rk)/rate,' seconds'

  stop

contains


  ! --------------------------------------------------------------------------------------
  
  subroutine get_neigh_list(neigh)
    !
    ! Setup the neighbor list array:
    ! neigh(i,1:4) = atom indices of the four neighbors of atom i
    ! Numbering of the neighbors:
    !
    !        (4)
    !         |
    !        ___
    !  (1)---|i|---(2)
    !        ---
    !         |
    !        (3)
    !

    implicit none
    integer,intent(out) :: neigh(nat,NNEIGH)
    integer :: i,j,k,jle,kle,jri,kri,jlo,klo,jup,kup

    do i=1,nat
       ! 2D indices ('lattice indices') of atom i
       j=mod(i-1,nuc)+1
       k=(i-1)/nuc+1
       ! 2D indices of neighbors with pbc's included
       jle=j-1; if (jle<1) jle=nuc
       kle=k
       jri=j+1; if (jri>nuc) jri=1
       kri=k
       jlo=j
       klo=k-1; if (klo<1) klo=nuc
       jup=j
       kup=k+1; if (kup>nuc) kup=1
       ! Convert these to 1D indices
       neigh(i,1)=jle+(kle-1)*nuc
       neigh(i,2)=jri+(kri-1)*nuc
       neigh(i,3)=jlo+(klo-1)*nuc
       neigh(i,4)=jup+(kup-1)*nuc
    end do

    return
  end subroutine get_neigh_list
  
  ! --------------------------------------------------------------------------------------

  subroutine accel(i,u,a)
    ! Calculate the potential energy u 
    ! and acceleration a of atom i.
    integer,intent(in) :: i
    real(rk),intent(out) :: u
    type(vector),intent(out) :: a
    real(rk),parameter :: d=1.0  !,k2=1.0,k3=0.4
    integer :: in,j
    real(rk) :: dx,dy,r,u2,u3,fx2,fy2,fx3,fy3

    u=0
    u2=0
    u3=0
    a%x=0
    a%y=0
    
    neighloop: do j=1,NNEIGH
       in=neigh(i,j)
       dx=x(in)%x-x(i)%x
       if (dx<-box/2.0) dx=dx+box
       if (dx>=box/2.0) dx=dx-box
       dy=x(in)%y-x(i)%y
       if (dy<-box/2.0) dy=dy+box
       if (dy>=box/2.0) dy=dy-box
       r=sqrt(dx**2+dy**2)
       u2=u2+1.0/2.0*k2*(r-d)**2
       u3=u3+1.0/3.0*k3*(r-d)**3
       fx2=k2*(r-d)*dx/r
       fx3=k3*(r-d)**2*dx/r
       fy2=k2*(r-d)*dy/r
       fy3=k3*(r-d)**2*dy/r
       a%x=a%x+fx2+fx3
       a%y=a%y+fy2+fy3
    end do neighloop

    u=(u2+u3)/2.0 ! Remember the factor 1/2 due to 'double counting'

    return
  end subroutine accel

  ! --------------------------------------------------------------------------------------

  subroutine printcoords()
    integer :: ia
    real(rk),parameter :: xsc=1.0
    write(10,*) nat
    write(10,'(a,x,i0,x,i0,x,a,3f14.4)') 'Frame number ',n,n,' fs boxsize',xsc*box,xsc*box,10.0
    do ia=1,nat
       write(10,'(a,x,i0,x,8g20.10)') 'Fe',ia,xsc*(x(ia)%x-box/2),xsc*(x(ia)%y-box/2),0.0,v(ia)%x,v(ia)%y,0.0,ep(ia),ek(ia)
    end do
    return
  end subroutine printcoords

end program md2d
