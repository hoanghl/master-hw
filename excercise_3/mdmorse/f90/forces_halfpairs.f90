! $Id: forces_halfpairs.f90 11 2006-11-02 14:14:45Z akuronen $

subroutine GetForces(x,a,N,size,periodic,m,neighbourlist,rpotcut, &
     morseD,morsealpha,morser0,Epot,virial)
  
  USE GLOBAL
  implicit none

  type(vector)  :: x(*),a(*)
  integer       :: N,neighbourlist(*)
  real(double)  :: m,rpotcut
  type(vector)  :: size
  type(ivector) :: periodic
  real(double)  :: morseD,morsealpha,morser0,Epot(*)
  real(double)  :: virial

  type(vector)  :: half
  integer       :: i,j,jj,startofineighbourlist,nneighboursi
  real(double)  :: r,rsq,rpotcutsq,dx,dy,dz,help1,help2,V,dVdr
  
  half%x=size%x/2.0
  half%y=size%y/2.0
  half%z=size%z/2.0
  rpotcutsq=rpotcut*rpotcut

  do i=1,N
     Epot(i)=0.0;
     a(i)%x=0.0;
     a(i)%y=0.0;
     a(i)%z=0.0;
  enddo

  virial=0.0;
  startofineighbourlist=1
  do i=1,N-1
     nneighboursi=neighbourlist(startofineighbourlist)
     !print *,i,nneighboursi     
     
     ! Loop over neighbours of i in neighbourlist
     do jj=1,nneighboursi
	! Retrieve real atom index
	j=neighbourlist(startofineighbourlist+jj)

	! Get distance
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

	!print *,i,j,rsq
	
	if (rsq <= rpotcutsq) then

	   ! i-j interaction within cutoff, evaluate it
  
	   r=sqrt(rsq)

	   ! Utilize the fact that one of the Morse terms is the
	   ! square of the other !
	   help1=exp(-morsealpha*(r-morser0));
	   help2=help1*help1

	   ! Calculate potential energy
           V=0.5*morseD*(help2-2.0*help1)
	   Epot(i)=Epot(i)+V
	   Epot(j)=Epot(j)+V

	   ! Calculate forces. Be careful with signs and the factor 1/2 !

	   ! First get dV/dr in units of eV/Å
           dVdr=-morseD*(help2*(-2.0*morsealpha)-2.0*help1*(-morsealpha))

	   ! Get virial = Sum (r·f) for pressure calculation in units of eV 
	   virial=virial+dVdr*(dx/r*dx+dy/r*dy+dz/r*dz)

	   ! then do transformations to get dV/dr in SI units
	   dVdr=dVdr*e/lunit
	   ! Then get a=F/m in units of Å/fs^2
	   dVdr=(dVdr/(m*u))/aunit/r

	   ! Now project on each vector

	   a(i)%x=a(i)%x-dVdr*dx
	   a(j)%x=a(j)%x+dVdr*dx
	   a(i)%y=a(i)%y-dVdr*dy
	   a(j)%y=a(j)%y+dVdr*dy
	   a(i)%z=a(i)%z-dVdr*dz
	   a(j)%z=a(j)%z+dVdr*dz
	   
	endif

     enddo ! End of loop over neighbours jj

     ! Set starting position for next atom
     startofineighbourlist=startofineighbourlist+nneighboursi+1

  enddo ! End of loop over atoms i
  


end subroutine GetForces


