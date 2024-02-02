! $Id: modules.f90 21 2008-10-13 10:10:15Z aakurone $

module global
  implicit none

  ! This gives us double precision.
  integer, parameter :: double=SELECTED_REAL_KIND(P=14,R=100)
  
  ! Define a vector type
  type vector
     real(double) :: x,y,z
  end type vector

  ! Define an integer vector type
  type ivector
     integer :: x,y,z
  end type ivector

  ! Define some needed unit conversion factors
  !
  ! To get internal units from SI, divide by one of these
  ! To get SI from internal units, multiply by one of these
  !
  real(double), parameter :: pi=3.14159265358979
  real(double), parameter :: e=1.6021892d-19
  real(double), parameter :: u=1.660565d-27
  real(double), parameter :: kB=1.380662d-23
  real(double), parameter :: lunit=1e-10
  real(double), parameter :: tunit=1e-15
  real(double), parameter :: vunit=lunit/tunit
  real(double), parameter :: aunit=vunit/tunit
  
  ! Set array sizes

  integer, parameter :: MAXAT=10000
  integer, parameter :: MAXNEIGHBOURS=100

  ! Do we read atom velocities from the xyz file?
  logical :: vfromxyz
  ! Is this a scattering study
  logical :: scattering
  
end module global
