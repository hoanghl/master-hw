! $Id: inout.f90 26 2009-07-23 12:11:52Z aakurone $

! Input-output subroutines for the mdmorse code

subroutine ReadParams(initialT,m,tmax,deltat,size,periodic,  &
     nneighbourlistupdate,rpotcut,rskincut,  &
     morseD,morsealpha,morser0,nmovieoutput, &
     btctau,desiredT,bpctau,bpcbeta,desiredP,seed,aoutt0,noutput)
     
  USE GLOBAL
  implicit none

  real(double)  :: initialT,m,tmax
  type(vector)  :: size
  type(ivector) :: periodic
  real(double)  :: deltat,rpotcut,rskincut,morseD,morsealpha,morser0
  integer       :: nneighbourlistupdate,nmovieoutput,seed,noutput
  real(double)  :: btctau,desiredT,bpctau,bpcbeta,desiredP,aoutt0
  
  integer :: i,ios
  character(80) :: line
  character(40) :: string
  real(double) :: x
  
  open(20,file='mdmorse.in',status='old',iostat=ios)
   
  if (ios /= 0) STOP 'Open error for file mdmorse.in'

  ! Initialize variables to sensible (?) default values
  initialT=0.0; m=63.546;
  size=vector(0.0,0.0,0.0)

  deltat=1.0;
  nneighbourlistupdate=5;
  rpotcut=4.02; rskincut=5.0;
  ! Initialize Morse parameters for Cu
  ! Values from Girifalco and Weizer, Phys. Rev. 114 (1959) 687.
  morseD=0.3429; morsealpha=1.3588; morser0=2.866;

  btctau=0.0; desiredT=0.0;
  bpctau=0.0; bpcbeta=1e-3; desiredP=0.0;
  
  seed=23267761
  vfromxyz=.false.
  scattering=.false.

  aoutt0=-1.0
  noutput=1

  
  ! Loop through file and find all variables
  i=0
  do
     i=i+1

     ! First read in line in whole.
     read(20,fmt='(A80)',iostat=ios) line
     ! Check for probable end-of-file
     if (ios < 0) exit

     ! If not parameter line, cycle
     if (line(1:1) /= '$') cycle
     
     ! Attempt to parse line into command and value using internal formatted read.
     read(unit=line,fmt=*,iostat=ios) string,x
     if (ios > 0) then
        print *,'ERROR: Data read in error on line',i
        stop 'Invalid parameter file'
     endif

     ! From here on line should be in variable format
     
     ! Look for variable-defining strings
     if (string=='$initialT') then
	initialT=x
     else if (string=='$desiredT') then
	desiredT=x
     else if (string=='$btctau') then
	btctau=x
     else if (string=='$mass') then
	m=x
     else if (string=='$xsize') then
	size%x=x
     else if (string=='$ysize') then
	size%y=x 
     else if (string=='$zsize') then
	size%z=x
     else if (string=='$periodicx') then
	periodic%x=int(x+0.5)
     else if (string=='$periodicy') then
	periodic%y=int(x+0.5)
     else if (string=='$periodicz') then
	 periodic%z=int(x+0.5)
     else if (string=='$seed') then
	seed=int(x+0.5)
     else if (string=='$tmax') then
	tmax=x
     else if (string=='$deltat') then
	deltat= x
     else if (string=='$nupdate') then
	nneighbourlistupdate=int(x+0.5)
     else if (string=='$rpotcut') then
	rpotcut=x 
     else if (string=='$rskincut') then
	rskincut=x     
     else if (string=='$morseD') then
	morseD=x
     else if (string=='$morsealpha') then
	morsealpha=x
     else if (string=='$morser0') then
	 morser0=x
     else if (string=='$bpctau') then
	 bpctau=x
     else if (string=='$bpcbeta') then
	 bpcbeta=x
     else if (string=='$desiredP') then
	 desiredP=x
     else if (string=='$aoutt0') then
	 aoutt0=x
     else if (string=='$noutput') then
	noutput=int(x+0.5)
     else if (string=='$nmovieoutput') then
	nmovieoutput=int(x+0.5)
     else if (string=='$scattering') then
	scattering=int(x)>0
     else
	print '(A,A)','Unknown parameter',string
	stop 'Parameter read in error'
     endif

     print '(A,A16,A,G13.6)','Read in parameter ',string,' value',x
    
  enddo

  if (initialT<0.0) then
     vfromxyz=.true.
     print '(a)','Atom velocities read from the xyz file.'
  end if
  if (scattering) then
     print '(a)','Scattering study.'
  end if

  close(20)

  ! Print out some general information about read results:

  print '(A,3I2)','Using periodics (1=on, 0=off)',     &
       periodic%x,periodic%y,periodic%z

  print '(A,3F12.6)','Morse potential parameters: D alpha r0', &
       morseD,morsealpha,morser0

  print '(A,I0,A)','Output selected every ',noutput,' steps'
  print '(A,I0,A)','Movie output selected every ',nmovieoutput,' steps'

  if (btctau > 0.0) then 
     print '(/,A,2F10.3,/)','Doing Berendsen temperature control with tau T', &
	  btctau,desiredT
  endif
  
  if (bpctau > 0.0) then 
     print '(/,A,2F10.3,/)','Doing Berendsen pressure control with tau beta', &
	  bpctau,bpcbeta
  endif
  
  return
  
end subroutine ReadParams


subroutine ReadAtoms(x,v,atomname,N)
  USE GLOBAL
  implicit none

  integer :: N
  type(vector) :: x(*),v(*)
  character(5) atomname(*)
  
  character(80) :: comment
  integer :: i,ios

  open(30,file='atoms.in',status='old',iostat=ios)

  if (ios/=0) STOP 'File atoms.in open failed'

  read(30,fmt=*) N

  if (N > MAXAT) STOP 'Error: MAXAT overflow, increase it.'

  read(30,fmt='(A80)') comment
  print '(A,I7,A,A40)','Reading in ',N,' atoms described as ',comment
  
  do i=1,N

     if (vfromxyz) then
        read(30,fmt=*,iostat=ios) atomname(i),x(i)%x,x(i)%y,x(i)%z,v(i)%x,v(i)%y,v(i)%z
     else
        read(30,fmt=*,iostat=ios) atomname(i),x(i)%x,x(i)%y,x(i)%z
     end if

     if (ios/=0) then
	print *,'Atom read in error for atom',i
	print *,atomname(i),x(i)%x,x(i)%y,x(i)%z
	STOP ' Atom read in error'
     endif

  enddo
  
  
end subroutine ReadAtoms



subroutine WriteAtoms(x,atomname,Ekin,Epot,N,time,size,Wat,istep)
  USE GLOBAL
  implicit none

  type(vector) :: x(*)
  character(5) :: atomname(*)
  integer :: N,istep
  real(double) :: Ekin(*),Epot(*),time,Wat(*)
  type(vector) :: size
  
  character(80) :: comment
  integer :: i,ios

  logical, save :: firsttime=.true.

  if (firsttime) then
     open(31,file='atoms.out',status='replace',iostat=ios)
     if (ios/=0) STOP 'File atoms.out open failed'
     firsttime=.false.
  endif
     

  write(31,fmt='(I10)') N

  !write(31,fmt='(A,F11.3,A,3F12.4)') 'mdmorse atom output at time ',time,' fs boxsize ',size%x,size%y,size%z
  write(31,fmt='('' Frame number '',i7,g14.5,'' fs boxsize'',3(F11.6,1X))') istep,time,size%x,size%y,size%z
  
  do i=1,N
     !write(31,fmt='(A2,1X,3g18.8,2g16.6)') atomname(i),x(i)%x,x(i)%y,x(i)%z,Epot(i),Ekin(i)
     write(31,fmt='(A2,1X,3g18.8,2g16.6)') atomname(i),x(i)%x,x(i)%y,x(i)%z,Epot(i),Wat(i)
  enddo
  
end subroutine WriteAtoms


subroutine scatteringoutput(time,x,v,Epot,Ekin)
  USE GLOBAL
  implicit none
  type(vector) :: x(*),v(*)
  real(double) :: time,Epot(*),Ekin(*)
  character(len=100),save :: form='(4g18.8)'

  write(100,fmt=form) time,x(1)%x,x(1)%y,x(1)%z
  write(200,fmt=form) time,x(2)%x,x(2)%y,x(2)%z

  write(101,fmt=form) time,v(1)%x,v(1)%y,v(1)%z
  write(201,fmt=form) time,v(2)%x,v(2)%y,v(2)%z

  write(102,fmt=form) time,Epot(1),Ekin(1)
  write(202,fmt=form) time,Epot(2),Ekin(2)

  return
end subroutine scatteringoutput
