      program test
      use iso_c_binding
      implicit none

      integer(c_int) L,V
      real(c_double),allocatable,dimension(:) :: out,in
      real(c_double):: sigma
      
      L=16
      V=16
      allocate(out(L))
      allocate(in(L))
      sigma=0.001

      in=1
      
      call inverse_laplace_operator(out,in,sigma,L,V)
      
      write(*,*) out
      
      stop
      end

      
