! $Id: neighbourlist.f90 17 2008-09-16 07:59:06Z aakurone $

subroutine UpdateNeighbourlist(x,N,size,periodic,neighbourlist,rskincut)
  USE GLOBAL
  implicit none

  ! Subroutine constructs a Verlet neighbour list in the
  ! following format:
  !
  ! nngbr1          Number of neighbours for atom 1
  ! ingbr1,1        Index of neighbour 1 of atom 1
  ! ingbr2,1        Index of neighbour 2 of atom 1
  ! ....
  ! ingbrnngbr1,1   Index of neighbour nngbr1 of atom 1
  ! nngbr2          Number of neighbours for atom 2
  ! ...

  
  type(vector)  :: x(*)
  integer       :: N,neighbourlist(*)
  type(vector)  :: size
  type(ivector) :: periodic
  real(double)  :: rskincut

  real(double)  :: dx,dy,dz,rsq,rskincutsq
  type(vector)  :: half
  integer       :: i,j,nneighboursi,startofineighbourlist,nneighbourstot
  
  half%x=size%x/2.0
  half%y=size%y/2.0
  half%z=size%z/2.0
  rskincutsq=rskincut*rskincut;

  nneighbourstot=0
  startofineighbourlist=1
  do i=1,N

     nneighboursi=0
     do j=1,N
	
	if (i==j) cycle
	
	dx=x(j)%x-x(i)%x;
	if (periodic%x==1 .and. dx>half%x ) dx=dx-size%x;
	if (periodic%x==1 .and. dx<-half%x) dx=dx+size%x;

	dy=x(j)%y-x(i)%y;
	if (periodic%y==1 .and. dy>half%y ) dy=dy-size%y;
	if (periodic%y==1 .and. dy<-half%y) dy=dy+size%y;

	dz=x(j)%z-x(i)%z;
	if (periodic%z==1 .and. dz>half%z ) dz=dz-size%z;
	if (periodic%z==1 .and. dz<-half%z) dz=dz+size%z;

	rsq=dx*dx+dy*dy+dz*dz

	if (rsq <= rskincutsq) then
	   !print *,'  ',j,x(j)%x,x(j)%y,x(j)%z
	   !print *,'accepted'
	   nneighboursi=nneighboursi+1
	   nneighbourstot=nneighbourstot+1
	   if (startofineighbourlist+nneighboursi > MAXAT*MAXNEIGHBOURS) then
	      STOP 'Neighbour list overflow, increase MAXNEIGHBOURS'
	   endif
	   neighbourlist(startofineighbourlist+nneighboursi)=j
	endif
		
     enddo

     ! Write in number of i's neighbours into list
     neighbourlist(startofineighbourlist)=nneighboursi
     ! Set starting position for next atom
     startofineighbourlist=startofineighbourlist+nneighboursi+1
     
  enddo

  !print '(A,G13.5,A)','Neighbour list update found ',    &
  !     1.0*nneighbourstot/N,' neighbours per atom'
  
end subroutine UpdateNeighbourlist


