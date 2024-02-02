! $Id: forces.f90 29 2009-07-23 12:38:11Z aakurone $

subroutine GetForces(x,a,N,size,periodic,m,neighbourlist,rpotcut, &
     morseD,morsealpha,morser0,Epot,virial,Wat)
  
  USE GLOBAL
  implicit none

  type(vector)  :: x(*),a(*)
  integer       :: N,neighbourlist(*)
  real(double)  :: m,rpotcut
  type(vector)  :: size
  type(ivector) :: periodic
  real(double)  :: morseD,morsealpha,morser0,Epot(*),Wat(*)
  real(double)  :: virial

  type(vector)  :: half
  integer       :: i,j,jj,startofineighbourlist,nneighboursi
  real(double)  :: r,rsq,rpotcutsq,dx,dy,dz,help1,help2,V,dVdr

  logical,save :: firsttime=.true.
  real(double),save :: vrcut,dvrcut,etmp
  integer,save :: counter=0
  logical,save :: DBENEX=.false.
  real(double) :: w1
  
  if (firsttime) then
     firsttime=.false.
     etmp=exp(-morsealpha*(rpotcut-morser0))
     vrcut=morseD*(etmp**2-2.0*etmp)
     dvrcut=2.0*morseD*morsealpha*(-etmp**2+etmp)
     write(6,'(/,a,2g20.10,/)') 'Energy and force shift-tilt terms:',vrcut,dvrcut
  end if
  

  half%x=size%x/2.0
  half%y=size%y/2.0
  half%z=size%z/2.0
  rpotcutsq=rpotcut*rpotcut

  do i=1,N     
     Epot(i)=0.0;
     Wat(i)=0.0
     a(i)%x=0.0;
     a(i)%y=0.0;
     a(i)%z=0.0;
  enddo

  virial=0.0;
  startofineighbourlist=1
  do i=1,N
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

	if (rsq <= rpotcutsq) then

	   ! i-j interaction within cutoff, evaluate it
  
	   r=sqrt(rsq)

	   ! Utilize the fact that one of the Morse terms is the
	   ! square of the other !
	   help1=exp(-morsealpha*(r-morser0));
	   help2=help1*help1

	   ! Calculate potential energy
           ! The two factors of 1/2 are due to
           !    1) division of the bond energy to two atoms
           !    2) double counting in neighbor list
	   V=0.5*0.5*(morseD*(help2-2.0*help1)-dvrcut*(r-rpotcut)-vrcut)
	   Epot(i)=Epot(i)+V
	   Epot(j)=Epot(j)+V

	   ! Calculate forces. Be careful with signs and the factor 1/2 !

	   ! First get dV/dr in units of eV/Å
           ! The factor 1/2 is due to double counting in neighbor list
           dVdr=(morseD*(help2*(-2.0*morsealpha)-2.0*help1*(-morsealpha))-dvrcut)/2.0

          if (DBENEX) then
             ! Debug output to compare with the mmc results (thus the factors 4 and 2).
             counter=counter+1
             if (mod(counter,10000)==0) then
                write(100,*) r,4.0*V
                write(101,*) r,2.0*dVdr
             end if
          end if

	   ! Get virial = Sum (r·f) for pressure calculation in units of eV 
           w1=-dVdr*r
	   virial=virial+w1
           Wat(i)=Wat(i)+w1

	   ! then do transformations to get dV/dr in SI units
	   dVdr=dVdr*e/lunit
	   ! Then get a=F/m in units of Å/fs^2
	   dVdr=(dVdr/(m*u))/aunit/r

	   ! Now project on each vector

	   a(i)%x=a(i)%x+dVdr*dx
	   a(j)%x=a(j)%x-dVdr*dx
	   a(i)%y=a(i)%y+dVdr*dy
	   a(j)%y=a(j)%y-dVdr*dy
	   a(i)%z=a(i)%z+dVdr*dz
	   a(j)%z=a(j)%z-dVdr*dz
	   
	endif

     enddo ! End of loop over neighbours jj


     ! Set starting position for next atom
     startofineighbourlist=startofineighbourlist+nneighboursi+1

  enddo ! End of loop over atoms i
  

end subroutine GetForces


