! $Id: main.f90 26 2009-07-23 12:11:52Z aakurone $

!
! Simple MD program made for demonstration purposes
! 
!
! Internal units:
!
!   Length:       Å  = 1e-10 m
!   Times:        fs = 1e-15 s
!   Velocity:     Å/fs
!   Accelration:  Å/fs^2
!   Forces:       eV/Å
!   Energy:       eV
!
! Simulation cell centered at the origin.
!

program mdmorse
  USE GLOBAL
  implicit none
  
  ! Arrays for atom data:
  !
  ! x,v,a    : atom positions, velocities and accelerations
  ! Epot     : potential energy for each atom
  ! Ekin     : kinetic energy     - " -
  ! Wat      : virial             - " -
  ! atomname : simply the atom name

  type(vector)  :: x(MAXAT),v(MAXAT),a(MAXAT)
  real(double)  :: Epot(MAXAT),Ekin(MAXAT),Wat(MAXAT)
  character(5)  :: atomname(MAXAT)
  integer       :: neighbourlist(MAXAT*MAXNEIGHBOURS)
  
  ! Other parameters and variables:
  !
  ! time                 : current time
  ! tmax                 : maximum time
  ! deltat               : time step
  ! itime                : integer number of the current time step
  ! T                    : temperature (in K) 
  ! m                    : atom mass
  ! Ekinsum, Epotsum     : system kinetic and potential energies
  ! Etot                 : total energy for all atoms
  ! N                    : number of atoms
  ! nneighbourlistupdate : interval between neighbour list updates
  ! noutput              : interval between 'ec' line output
  ! nmovieoutput         : interval between atom position output
  ! rpotcut              : potential cutoff radius
  ! rskincut             : neighbour list cutoff radius
  ! size                 : simulation cell size
  ! periodic             : periodic boundaries will be used (1=yes, 0=no)
  ! seed                 : random number generator seed number
  !
  real(double)  :: time,tmax,deltat,T,Ekinsum,Epotsum,Etot,m
  !
  integer       :: N,nneighbourlistupdate,itime,noutput,nmovieoutput
  real(double)  :: rpotcut,rskincut
  type(vector)  :: size
  type(ivector) :: periodic
  integer       :: seed

  ! Morse potential parameters
  real(double) :: morseD,morsealpha,morser0

  ! Parameters for Berendsen temperature and pressure control 
  ! [J. Chem. Phys.81 (1984) 3684]
  !
  ! btctau   : Berendsen temperature control tau, 
  ! initialT : initial temperature
  ! desiredT : desired temperature
  ! bpctau   : Berendsen pressure control tau
  ! bpcbeta  : 1/B
  ! P        : instantaneous pressure
  ! desiredP : desired pressure
  real(double) :: btctau,initialT,desiredT
  real(double) :: bpctau,bpcbeta,P,desiredP,virial,mu
  real(double) :: aoutt0
  
  integer :: i

  print *,''
  print *,'--------------- mdmorse08 1.0f --------------------'
  print *,''
  
  !
  ! Initialize the simulation
  !
  call ReadParams(initialT,m,tmax,deltat,size,periodic,  &
       nneighbourlistupdate,rpotcut,rskincut,  &
       morseD,morsealpha,morser0,nmovieoutput,    &
       btctau,desiredT,bpctau,bpcbeta,desiredP,seed,aoutt0,noutput)

  call ReadAtoms(x,v,atomname,N)
  
  ! Factor 2 because of energy equipartition theorem
  if (.not.vfromxyz) call SetTemperature(v,N,m,2.0*initialT,seed)
  ! Check that we got it right
  call GetTemperature(v,N,m,T)
  print *
  print *,'Initial atom temperature is ',T

  ! Initialize accelerations to 0
  a = vector(0.0, 0.0, 0.0)
  time=0.0
  itime=0
  P=0.0

  ! Initial outputs
  call UpdateNeighbourlist(x,N,size,periodic,neighbourlist,rskincut)
  call GetForces(x,a,N,size,periodic,m,neighbourlist,rpotcut, &
       morseD,morsealpha,morser0,Epot,virial,Wat)
  call GetTemperature(v,N,m,T)
  call GetEnergies(v,N,m,Ekinsum,Epotsum,Etot,Ekin,Epot)
  P=(N*kB*T+1.0/3.0*virial*e)/(size%x*size%y*size%z*1e-30)/1e8;
  if (noutput>0) print '(A,g18.8,1X,g18.8,5(1X,g18.8))','ec ',time,T,Ekinsum/N,Epotsum/N,Etot/N,P,virial
  if (nmovieoutput>0 .and. aoutt0<=0.0) then
     print '(A,F10.3)','Outputting atom movie at t = ',time
     call WriteAtoms(x,atomname,Ekin,Epot,N,time,size,Wat,itime)
  endif


  ! Start main loop over times

  do 

     if (mod(itime,nneighbourlistupdate)==0) then
	call UpdateNeighbourlist(x,N,size,periodic,neighbourlist,rskincut)
     endif
     
     call Solve1(x,v,a,N,size,periodic,deltat)
     
     call GetForces(x,a,N,size,periodic,m,neighbourlist,rpotcut, &
	  morseD,morsealpha,morser0,Epot,virial,Wat)

     call GetTemperature(v,N,m,T)
     
     call Solve2(v,a,N,deltat,btctau,T,desiredT)
     
     ! Update time
     time=time+deltat
     
     ! Actual solution done, now compute result quantities

     call GetEnergies(v,N,m,Ekinsum,Epotsum,Etot,Ekin,Epot)

     ! Compute pressure in units of kbar
     P=(N*kB*T+1.0/3.0*virial*e)/(size%x*size%y*size%z*1e-30)/1e8;

     ! Output pressure in units of kbar
     if (noutput>0) then
	if (mod(itime,noutput)==0) then
	   print '(A,g18.8,1X,g18.8,5(1X,g18.8))','ec ',time,T,Ekinsum/N,Epotsum/N,Etot/N,P,virial
	end if
     end if

     if (nmovieoutput>0 .and. time>aoutt0) then
        if (mod(itime,nmovieoutput)==0) then
           print '(A,F10.3)','Outputting atom movie at t = ',time
           call WriteAtoms(x,atomname,Ekin,Epot,N,time,size,Wat,itime)
        endif
     end if
     if (scattering) then
        call scatteringoutput(time,x,v,Epot,Ekin)
     end if

     ! Do Berendsen pressure control
     if (bpctau > 0.0) then
	mu=(1.0-bpcbeta*deltat/bpctau*(desiredP-P))**(1.0/3.0)

	size%x=size%x*mu
	size%y=size%y*mu
	size%z=size%z*mu
	do i=1,N
	   x(i)%x=x(i)%x*mu
	   x(i)%y=x(i)%y*mu
	   x(i)%z=x(i)%z*mu
	enddo
	if (mod(itime,10*noutput)==1) then
	   print '(A,F16.3,3(1X,F16.6),3(1X,F16.5))','bpc ',time,size%x,size%y,size%z, &
		size%x*size%y*size%z,P,mu
	endif
     endif
	   
     ! Check whether we are done
     if (time >= tmax) exit     

     itime=itime+1

  end do
  

  ! Always output final energy and atom coordinates using high precision.
  print '(A,g18.8,1X,g18.8,5(1X,g18.8))','ec ',time,T,Ekinsum/N,Epotsum/N,Etot/N,P,virial
  !print '(A,F12.5,1X,F12.5,4(1X,F14.7))','ec ',time,T,Ekinsum/N,Epotsum/N,Etot/N,P

  call WriteAtoms(x,atomname,Ekin,Epot,N,time,size,Wat,itime)

  ! Close movie file 
  close(31)

end program mdmorse
  
