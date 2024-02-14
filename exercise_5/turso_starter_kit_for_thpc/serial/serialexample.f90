program serialexample
  implicit none
  character(len=200) :: host

  call get_environment_variable('HOSTNAME',host)
  print '(a,a)','host: ',trim(host)

end program serialexample
