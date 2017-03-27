      program test
      use iso_c_binding
      implicit none
      
      integer(c_int) n
      real(c_double),allocatable,dimension(:) :: mat
      real(c_double):: cond_numb
      
      n=2
      allocate(mat(n*n))
      cond_numb=0.001
      
      call fill_defpos_symm_matrix(mat,cond_numb,n)
      
      write(*,*) mat
      
      stop
      end
