! $Id: solve.f90 17 2008-09-16 07:59:06Z aakurone $


subroutine Solve1(x,v,a,N,size,periodic,deltat)
  USE GLOBAL
  implicit none

  type(vector)  :: x(*),v(*),a(*)
  integer       :: N
  real(double)  :: deltat
  type(vector)  :: size
  type(ivector) :: periodic

  real(double)  :: halfdeltatsq
  type(vector)  :: half

  integer :: i

  real(double) :: dxmax,dx
  integer :: idxmax
  
  halfdeltatsq=0.5*deltat*deltat

  half%x=size%x/2.0
  half%y=size%y/2.0
  half%z=size%z/2.0

  dxmax=-100.0

  do i=1,N

     ! First part of velocity Verlet solution
     
     x(i)%x=x(i)%x+deltat*v(i)%x+halfdeltatsq*a(i)%x
     x(i)%y=x(i)%y+deltat*v(i)%y+halfdeltatsq*a(i)%y
     x(i)%z=x(i)%z+deltat*v(i)%z+halfdeltatsq*a(i)%z

     dx=sqrt( (deltat*v(i)%x+halfdeltatsq*a(i)%x)**2 + & 
          &   (deltat*v(i)%y+halfdeltatsq*a(i)%y)**2 + &
          &   (deltat*v(i)%z+halfdeltatsq*a(i)%z)**2)
     if (dx>dxmax) then
        dxmax=dx
        idxmax=i
     end if

     v(i)%x=v(i)%x+0.5*deltat*a(i)%x
     v(i)%y=v(i)%y+0.5*deltat*a(i)%y
     v(i)%z=v(i)%z+0.5*deltat*a(i)%z
     
     ! Enforce periodics if needed
     if (periodic%x==1 .and. x(i)%x < -half%x) x(i)%x=x(i)%x+size%x
     if (periodic%x==1 .and. x(i)%x >= half%x) x(i)%x=x(i)%x-size%x
     
     if (periodic%y==1 .and. x(i)%y < -half%y) x(i)%y=x(i)%y+size%y
     if (periodic%y==1 .and. x(i)%y >= half%y) x(i)%y=x(i)%y-size%y
     
     if (periodic%z==1 .and. x(i)%z < -half%z) x(i)%z=x(i)%z+size%z
     if (periodic%z==1 .and. x(i)%z >= half%z) x(i)%z=x(i)%z-size%z

  enddo

end subroutine Solve1


subroutine Solve2(v,a,N,deltat,btctau,T,desiredT)
  USE GLOBAL
  implicit none

  type(vector) :: v(*),a(*)
  integer      :: N
  real(double) :: deltat
  real(double) :: T,desiredT,btctau

  real(double) :: btclambda
  
  integer :: i

  ! Compute btc lambda. See J. Chem. Phys. 81 (1984) 3684.
  ! Note the factor 2 inside the sqrt. For details see the lecture notes.
  if (btctau > 0.0 .and. T > 0.0) then
     btclambda = sqrt(1+2.0*deltat/btctau*(desiredT/T-1.0))
  else
     btclambda=1.0;
  endif

  do i=1,N

     ! Second part of velocity Verlet solution
     
     v(i)%x=v(i)%x+0.5*deltat*a(i)%x
     v(i)%y=v(i)%y+0.5*deltat*a(i)%y
     v(i)%z=v(i)%z+0.5*deltat*a(i)%z

     v(i)%x=v(i)%x*btclambda
     v(i)%y=v(i)%y*btclambda
     v(i)%z=v(i)%z*btclambda

  enddo

end subroutine Solve2

