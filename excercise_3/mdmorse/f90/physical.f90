! $Id: physical.f90 19 2008-09-16 13:03:49Z aakurone $

subroutine SetTemperature(v,N,m,T,seed)
  USE GLOBAL
  implicit none

  !
  ! Subroutine sets the temperature of the atoms in the system
  ! to T assuming a Maxwell-Boltzmann velocity distribution
  !

  type(vector) :: v(*)
  integer      :: N
  real(double) :: m,T
  integer      :: seed

  type(vector) :: vtot
  real(double) :: vxrms,gaussianrand
  integer :: i
  
  ! RMS v in one dimension. See e.g. Mandl p. 191
  
  vxrms=sqrt(kB*T/(m*u))/vunit

  vtot%x=0.0; vtot%y=0.0; vtot%z=0.0; 
  do i=1,N

     v(i)%x=vxrms*gaussianrand(seed)
     v(i)%y=vxrms*gaussianrand(seed)
     v(i)%z=vxrms*gaussianrand(seed)

     vtot%x=vtot%x+v(i)%x
     vtot%y=vtot%y+v(i)%y
     vtot%z=vtot%z+v(i)%z
     
  enddo

  vtot%x=vtot%x/N
  vtot%y=vtot%y/N
  vtot%z=vtot%z/N

  ! Subtract total velocity of atoms 
  do i=1,N
     v(i)%x=v(i)%x-vtot%x
     v(i)%y=v(i)%y-vtot%y
     v(i)%z=v(i)%z-vtot%z
  enddo
  
  return

end subroutine SetTemperature


subroutine GetTemperature(v,N,m,T)
  USE global
  implicit none

  !
  ! Subroutine gets the temperature of the atoms in the system
  !

  type(vector) :: v(*)
  integer :: N
  real(double) :: T,m

  real(double) :: Ekinsum,gaussianrand
  integer :: i

  ! Get sum of kinetic energies
  Ekinsum=0
  do i=1,N

     Ekinsum=Ekinsum+v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z

  enddo
  ! Do necessary unit transformations to get E in J
  Ekinsum=0.5*m*u*Ekinsum*vunit*vunit
  
  ! and get temperature using E=3/2 N k T
  T = Ekinsum/(1.5*N*kB)
  
  return

end subroutine GetTemperature

subroutine GetEnergies(v,N,m,Ekinsum,Epotsum,Etot,Ekin,Epot)
  USE GLOBAL
  implicit none

  type(vector) :: v(*)
  real(double) :: Ekinsum,Epotsum,Etot,m,Ekin(*),Epot(*)
  integer :: N
  
  real(double) :: help1
  integer :: i
  
  Ekinsum=0.0; Epotsum=0.0;
  help1=0.5*m*u*vunit*vunit/e;
  do i=1,N
     Ekin(i)=help1*(v(i)%x*v(i)%x+v(i)%y*v(i)%y+v(i)%z*v(i)%z)
     Ekinsum=Ekinsum+Ekin(i)
     Epotsum=Epotsum+Epot(i)
  enddo
  Etot=Ekinsum+Epotsum

  return

end subroutine GetEnergies
  

!---------------------------------------------------------------------------

function uniformrand(seed)
  use global
  implicit none

  ! -----------------------------------------------------------
  ! Park-Miller "minimal" Random number generator
  ! uniform distribution ]0,1[ . See Numerical Recipes ch. 7.0
  ! -----------------------------------------------------------

  real(double) :: uniformrand
  integer :: seed
  integer :: IA,IM,IQ,IR,MASK
  real(double) :: AM
  integer :: k
  
  parameter (IA=16807,IM=2147483647,AM=1.0/IM,IQ=127773,IR=2836,MASK=123459876)

  seed=ieor(seed,MASK)
  k=seed/IQ
  seed=IA*(seed-k*IQ)-IR*k
  if (seed < 0) seed=seed+IM
  uniformrand=AM*seed
  seed=ieor(seed,MASK)
  
  return  
end function uniformrand



function gaussianrand(seed)
  use global
  implicit none
  
  ! ---------------------------------------------------------
  ! Random numbers with normal (Gaussian) distribution.
  ! Mean is 0 and standard deviation is 1
  ! See W.H.Press et al., Numerical Recipes 1st ed., page 203
  ! ---------------------------------------------------------
  
  real(double) :: gaussianrand
  integer :: seed
  real(double) :: fac,v1,v2,r
  real(double), save :: gset
  integer, save :: iset=0
  real(double) :: uniformrand
  
  
  
  if (iset==0) then
1    v1 = 2.*uniformrand(seed)-1.
     v2 = 2.*uniformrand(seed)-1.
     r = v1*v1+v2*v2
     if (r>=1.0) goto 1
     fac = sqrt(-2.0*log(r)/r)
     gset = v1*fac
     gaussianrand = v2*fac
     iset = 1
  else
     gaussianrand = gset
     iset = 0
  endif
  
  return
end function gaussianrand


