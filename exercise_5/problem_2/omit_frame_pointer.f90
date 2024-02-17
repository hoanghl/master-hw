

module funcmod
  integer,parameter :: rk=8
contains

  real(rk) function f1(x)
    implicit none
    real(rk),intent(in) :: x
    f1=x/1e10
    return
  end function f1
  
  real(rk) function f2(x)
    implicit none
    real(rk),intent(in) :: x
    f2=f1(x)+f1(x)
    return
  end function f2
  
  real(rk) function f3(x)
    implicit none
    real(rk),intent(in) :: x
    f3=f2(x)+f2(x)
    return
  end function f3
  
  real(rk) function f4(x)
    implicit none
    real(rk),intent(in) :: x
    f4=f3(x-1.0)*f3(x+1.0)
    return
  end function f4
  
  real(rk) function f5(x)
    implicit none
    real(rk),intent(in) :: x
    f5=f4(x*2.0)*f3(x/2.0)
    return
  end function f5
  
  real(rk) function f6(x)
    implicit none
    real(rk),intent(in) :: x
    f6=2.0*f5(x)
    return
  end function f6
end module funcmod
  
program omit_frame_pointer
  use funcmod
  implicit none
  integer :: n,i
  real(rk),allocatable :: a(:)
  real(rk) :: t1,t2

  n=10000000
  allocate(a(n))

  call cpu_time(t1)
  do i=1,n
     a(i)=f6(0.1_rk*n)
  end do
  call cpu_time(t2)

  print '(2f20.10)',sum(a),t2-t1

end program omit_frame_pointer

