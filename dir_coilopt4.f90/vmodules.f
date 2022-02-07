      module times_stuff
      use cparms    
      real(rprec) :: strt_time, fin_time
      end module times_stuff

      module Vname
      use cparms
#if defined(RISC)
      integer, parameter :: num_processors    = 1
      integer, parameter :: num_levmar_params = 1
#else
      integer, parameter :: num_processors    = 8
      integer, parameter :: num_levmar_params = 10
#endif      
      character :: input_data_file*120
      character :: for05_file*120, for05_file_new*126
      character :: extension*120
      real(rprec), parameter :: min_chisq = 1000000
      real(rprec) :: chisq_min = min_chisq
      logical :: lbatch
      end module Vname
 
 
      module Vwire
      integer :: nwire, nwire1
      end module Vwire

      module coiltypes
      use cparms
 
         type modularcoil
            real(rprec), dimension(:), pointer :: rhoc, rhos,
     1         phic, phis, rho, phi
            real(rprec) :: current
         end type modularcoil      
 
         type saddlecoil
            real(rprec), dimension(:), pointer :: 
     1         v_c, v_s, u_c, u_s, phi, theta
            real(rprec) :: current
         end type saddlecoil      
 
         type vfcoil
            real(rprec) :: radius, height, current
         end type vfcoil      
 
      end module coiltypes

      module coils
      use coiltypes
 
! This module declares an allocatable array of coils for each type of coil.
! The actual size allocation is done in the initialization routine. 
 
         type(modularcoil), dimension(:), allocatable :: modular
         type(saddlecoil),  dimension(:), allocatable :: saddle
         type(vfcoil),      dimension(:), allocatable :: vertical
          
      end module coils
 
      module mpi_params                                      !mpi stuff
      integer :: myid=0, numprocs                            !mpi stuff
      integer, parameter :: master = 0                       !mpi stuff
      end module mpi_params                                  !mpi stuff
