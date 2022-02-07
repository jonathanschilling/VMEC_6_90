      program coilopt4
      use boundary
      use bnorm_mod
      use modular_coils
      use saddle_coils
      use saddle_surface
      use bcoils_mod
      use vf_coils
      use Vcoilpts
      use Vname
      use control_mod
      use system_mod
      use safe_open_mod
      use mpi_params                                         !mpi stuff
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: nxc=1000
      integer :: istat, ifun, nvar, index_blank, numargs, numchars
      integer :: iargc, getarg
      integer :: i, j, n, iunit
      real(rprec) :: xc(nxc)
      character*120  :: arg1, arg2
!     For B-normal from BNORM code
      character(len=3) :: decision = 'no'
      integer :: nvariables
      integer :: ncl
      integer :: ierr                                        !mpi stuff
      real(rprec) :: start_time, finish_time, del_time, del_time_max
      real(rprec) :: crv
      real(rprec) :: u0, v0
      real(rprec) :: x, y, z
      real(rprec) :: tx, ty, tz
      real(rprec) :: nx, ny, nz
      logical :: lplasma, isthere, isfor05there
      character*200 :: scratch_dir, wout_file, temp
      integer :: k, ierr_mpi
!-----------------------------------------------
 
!     Call MPI Initialization routines:
 
#ifdef MPI_OPT
      call MPI_INIT( ierr )                                  !mpi stuff
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )       !mpi stuff
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )   !mpi stuff
#endif
 
      call second0(start_time)
!-----------------------------------------------
 
!     CALLING SYNTAX:
 
!     XCOIL EXTENSION_NAME
 
!     THE FILES <FOR05.EXTENSION_NAME>, <WOUT.EXTENSION_NAME>
!     MUST EXIST IN PATH. IF THEY DO NOT, THEN THIS IS BEING
!     RUN IN STELLOPT OPTIMIZER, AND THE NAMELIST IS LOOKED FOR
!     IN THE <INPUT_DATA_FILE>
 
!     READ COMMAND LINE EXTENSION
 
!     numargs = iargc()
!     numchars = getarg(1,arg1)
      call getcarg (1, arg1, numargs)
      if( numargs.le.0 )then
         if( myid .eq. master ) then                         !mpi stuff
            print *,' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         end if     !if( myid .eq. master )                  !mpi stuff
        stop 01
      endif
  
      index_blank = index(arg1,' ') - 1
      extension = arg1(1:index_blank)
      for05_file = 'for05.' // trim(extension)
      for05_file_new = for05_file(1:(index_blank+6)) // '.new'
      input_data_file = 'input.' // trim(extension)
      bnorm_file = 'bnorm.' // trim(extension)
      
#if defined(GEOM_ONLY)
#else
      wout_file = 'wout.' // trim(extension)
      scratch_dir = "coilopt" // "_" // trim(extension)
!
!     DO ALL CALCULATIONS IN SCRATCH DIRECTORY
!
      if (myid .eq. master) then                                          !START MPI
         temp = '/bin/rm -Rf ' // scratch_dir
         call system(temp)
         temp = '/bin/mkdir -m 755 ' // scratch_dir
         call system(temp)
         k = chdir(scratch_dir)
         if (k .ne. 0) then
            print *,
     1      'ierr = ',k,': Unable to chdir to working directory!'
            stop
         end if

!
!     COPY INPUT FILE FROM ORIGINAL DIRECTORY
!
         inquire(file="../" // for05_file, exist=isthere)
         if (isthere) then
            temp = "/bin/cp ../" // trim(for05_file) // 
     1                     " ./" // trim(for05_file)
            call system(temp)
         end if

         isfor05there = isthere

         inquire(file="../" // for05_file_new, exist=isthere)
         if (.not.isthere .and. isfor05there) then
            temp = "/bin/cp ../"//trim(for05_file) //
     1                     " ./" // trim(for05_file_new)
            call system(temp)
         else if (isthere) then
            temp = "/bin/cp ../" // trim(for05_file_new) // 
     1                     " ./" // trim(for05_file_new)
            call system(temp)
         end if

         inquire(file="../" // input_data_file, exist=isthere)
         if (isthere) then
            temp = "/bin/cp ../" // trim(input_data_file) // 
     1                     " ./" // trim(input_data_file)
            call system(temp)
         end if

         inquire(file="../" // wout_file, exist=isthere)
         if (isthere) then
            temp = "/bin/ln -s ../" // trim(wout_file) // 
     1                     " ./" // trim(wout_file)
            call system(temp)
         end if

         inquire(file="../" // bnorm_file, exist=isthere)
         if (isthere) then
            temp = "/bin/ln -s ../" // trim(bnorm_file) // 
     1                     " ./" // trim(bnorm_file)
            call system(temp)
         end if
      end if                                                             !END MPI
#endif

#ifdef MPI_OPT
      call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      if (ierr_mpi .ne. 0) stop 'MPI_BARRIER error in PROGRAM COILOPT4'
      if (myid .ne. master) k = chdir(scratch_dir)                       !MPI
#endif
   
      xc = zero
 
!     READ IN PLASMA BOUNDARY. IF INSIDE OPTIMIZER (GEOM_ONLY), WE DO NOT
!     KNOW THE FINAL BOUNDARY SHAPE AT THIS POINT SO READING WOUT FILE IS
!     INAPPROPRIATE. WILL GET THE NFP VARIABLE LATER BY READING INDATA NAMELIST
  
#if defined(GEOM_ONLY)
      lplasma = .false.
      if (numargs .ge. 2) then
!        numchars = getarg(2,arg2)
         call getcarg (2, arg2, numargs)
         if (arg2(1:2).eq.'-P'.or.arg2(1:2).eq.'-p') lplasma = .true.
      endif
#else
      lplasma = .true.
#endif      
  
      if (lplasma) then
         call read_wout(extension)
         nedge = nuv
 
!     COMPUTE SURFACE NORMAL AT A FIXED NUMBER (NEDGE) OF MATCHING POINTS
!     (IN THETA=U, ZETA=V SPACE)
 
         call normal_vector
      else
         nedge = 0
      end if      
  
      nvar = -1
      call initialize_opt(xc, nxc, nvar)
 
!     READ COEFFICIENTS AND EVALUATE BNORM (NOTE: LBNORM ONLY USED IN FOLLOWING BLOCK OF CODE)
!#if defined(GEOM_ONLY)
!#else
      if (lbnorm .and. lplasma) then
         if (myid .eq. master) print *, 'read ', trim(bnorm_file)
         call read_bnorm_coefs
         call evaluate_bnorm
#if defined(GEOM_ONLY)
#else
      else
         if (myid .eq. master) then
            print *, 'LBNORM = FALSE: BNORM FILE WILL BE IGNORED'
            print *, 'IS THIS CORRECT OR SHOULD THIS BE ABORTED?'
            pause
         end if
#endif 
      endif
!#endif 
!     INITIALIZE POINTS FOR ACCESS ZONES
      if (laccess) then
         call initialize_access
      end if
 
!     PERFORM OPTIMIZATION
 
!     ifun needs to be changed when fvec components are added
 
      ifun = nedge + 4*nmod_unique_coils + nmod_coils_per_period
     1     + 9*nsad_unique_coils + n_access + 2*num_vf + 5
 
      call optimize (xc,nvar,numsurf,ifun,nu,nv)
  
      if( myid .eq. master ) then                            !mpi stuff
         iunit = 14
         call safe_open(iunit, istat, 'for05', 'unknown', 
     1        'formatted')
         close(iunit,iostat=istat,status='delete')
 
!     Write b-norm error to file
 
#if defined(GEOM_ONLY)
#else
         iunit = 23
         call safe_open(iunit, istat, 'b_norm.dat', 'unknown', 
     1        'formatted')
         n = 0
         do i=1, nu
            do j=1, nv
               n = n + 1
               write(iunit,1000) nfp*phib(n)/dtwopi, thetab(n)/dtwopi,
     1                        b_error(n)
            end do
            write(iunit,1010)
         end do
 1000 format(1p3e12.4)
 1010 format("")
         close (iunit)
#endif
 
         print 1030, rbphi_avg
 1030 format(/,' <R*Bphi> [T-m] = ',1pe14.4)
 
#if defined(GEOM_ONLY)
#else
 
!        Evaluate max field at coils
 
         if (lmodular) then
            iunit = 29
            call safe_open(iunit, istat, 'surf_norm.dat', 'unknown', 
     1        'formatted')
            call evaluate_max_field (iunit)
            print 1100
            do i=1, nmid
               print 1110, i, b_max(i), curmod(i)
            end do
 1100 format (" Coil",4x,"Bmax (T)",5x,"Imod (Amp)")
 1110 format (i4,1p2e14.4)
            close (iunit)
         end if
#endif
 
!        Write coil and surface data to file
 
         call write_coils (extension)
#if defined(GEOM_ONLY)
#else
         call write_surfaces
         call evaluate_bg_field
#endif         
 
      end if     !if( myid .eq. master )                     !mpi stuff
 
!     FREE MEMORY
 
      if (allocated(xm_bmn)) deallocate (xm_bmn, xn_bmn, bmn)
      if (allocated(rb)) deallocate (rb, zb, phib, thetab, x_p, y_p, 
     1  z_p, rb_ph, rb_th, zb_ph, zb_th, n_r, n_phi, n_z, d_area,
     2  theta_d, phi_d, bnormal_match, b_error, b_mod, bmn_error, 
     3  bsl_error, luv)
      if (allocated(ixm)) deallocate (ixm, ixn, rmnc_b, zmns_b,
     1  xm, xn, gmn_b, lmns_b)
      if (allocated(nbrho_opt)) deallocate (nbrho_opt, mbrho_opt, 
     1  rbc, zbs, rhobc, nrz0_opt, delta_mn)

  
      call second0(finish_time)
      del_time = finish_time - start_time
#ifdef MPI_OPT
      call MPI_Reduce(del_time, del_time_max, 1, MPI_REAL8,
     1  MPI_MAX, master, MPI_COMM_WORLD, ierr)
#else
      del_time_max = del_time      
#endif
      if(myid .eq. master) then
         write(*,'("Elapsed time (from rtc) = ",e15.7)') del_time_max
      endif
#ifdef MPI_OPT
       call MPI_FINALIZE(ierr)      !Close out MPI           !mpi stuff
#endif
      END PROGRAM COILOPT4
