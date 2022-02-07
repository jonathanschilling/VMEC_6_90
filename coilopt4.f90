#!/bin/sh
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
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
EOF
cat > driver.f << "EOF"
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
EOF
cat > coilopt4.f << "EOF"
      subroutine normal_vector
      use cparms
      use boundary
      use bnorm_mod
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, n, ku, kv, mn
      real(rprec) :: cosmu, cosnv, sinmu, sinnv, alu, alv,
     1  cosmn1, sinmn1
!-----------------------------------------------
      allocate (rb(nuv), zb(nuv),
     1  phib(nuv), thetab(nuv), x_p(nuv), y_p(nuv), z_p(nuv),
     2  rb_ph(nuv), rb_th(nuv), zb_ph(nuv), zb_th(nuv),
     3  n_r(nuv), n_phi(nuv), n_z(nuv), d_area(nuv),
     4  theta_d(nuv), phi_d(nuv), bnormal_match(nuv),
     5  b_error(nuv), b_mod(nuv), bmn_error(nuv), bsl_error(nuv),
     6  luv(nuv), stat = i)
  
      if( myid .eq. master ) then                            !mpi stuff
         if (i .ne. 0) stop 'allocation error in normal_vector'
      end if     !if( myid .eq. master )                     !mpi stuff
  
      alu = dtwopi/nu
      alv = dtwopi/(nv*nfp)
  
      z_p = zero
      x_p = zero
      y_p = zero
      d_area = zero
       
!  as part of calculation the normal vector, the array phib calculated for each 
!  of the mnmax boundary points  (nedge  = nuv = nu*nv from vmec or something
!  larger from the assignment statements in the read_wout routine.
!  note that the rb points assume vmec symmetry (m*theta-n*zeta) while
!  the winding surface has NESCOIL symmetry (plus sign).  
 
      if( myid .eq. master ) then                            !mpi stuff
         print *, 'nedge = ', nedge
      end if     !if( myid .eq. master )                     !mpi stuff
      i = 0
      do ku = 1, nu
         do kv = 1, nv
            i = i + 1
            thetab(i) = (ku-1)*alu
            phib(i) = (kv-1)*alv
!  boundary grid coordinates - cylinderical
            rb(i)  = 0
            zb(i)  = 0
            rb_th(i)  = 0
            rb_ph(i)  = 0
            zb_th(i)  = 0
            zb_ph(i)  = 0
            do mn = 1, mnmax
               cosnv = cos(xn(mn)*phib(i))
               sinnv = sin(xn(mn)*phib(i))
               cosmu = cos(xm(mn)*thetab(i))
               sinmu = sin(xm(mn)*thetab(i))
               cosmn1 = cosmu*cosnv + sinmu*sinnv
               sinmn1 = sinmu*cosnv - cosmu*sinnv
               rb(i) = rb (i) + rmnc_b(mn) * cosmn1
               zb(i) = zb (i) + zmns_b(mn) * sinmn1
!  partial derivatives w.r.t. phi, theta
               rb_th(i) = rb_th(i) - xm(mn) * rmnc_b(mn) * sinmn1
               rb_ph(i) = rb_ph(i) + xn(mn) * rmnc_b(mn) * sinmn1
               zb_th(i) = zb_th(i) + xm(mn) * zmns_b(mn) * cosmn1
               zb_ph(i) = zb_ph(i) - xn(mn) * zmns_b(mn) * cosmn1
            end do
         end do
      end do
!  surface normal components
      sum_d_area = 0
      do i=1,nedge
         n_r(i) = rb(i)*zb_th(i)
         n_phi(i) = rb_th(i)*zb_ph(i) - rb_ph(i)*zb_th(i)
         n_z(i) = -rb(i)*rb_th(i)
      end do
      do i=1,nedge
         d_area(i) = sqrt(n_r(i)**2 + n_phi(i)**2 + n_z(i)**2)
         sum_d_area = sum_d_area + d_area(i)
      end do
!     print *, sum_d_area
      do i=1,nedge
         n_r(i) = n_r(i)/d_area(i)
         n_phi(i) = n_phi(i)/d_area(i)
         n_z(i) = n_z(i)/d_area(i)
      end do
!  boundary grid coordinates - cartesian
      do i=1,nedge
         x_p(i)=rb(i)*cos(phib(i))
         y_p(i)=rb(i)*sin(phib(i))
         z_p(i)=zb(i)
      end do
 
      end subroutine normal_vector


      subroutine initialize_opt (xc, nxc, nvar)
      use cparms
      use boundary, only: nfp
      use vmec_input, only: nfp_vmec => nfp, indata
      use modular_coils
      use saddle_coils
      use saddle_surface
      use bcoils_mod
      use tor_field
      use Vcoilpts
      use coilsnamin
      use Vname
      use Vwire   
      use coils   
      use control_mod
      use gade_mod
      use mpi_params                                         !mpi stuff
      use safe_open_mod
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: nxc, nvar
      real(rprec) :: xc(nxc)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr, i, modes, n, nc, k, status
      integer :: iunit=24
      integer :: nvariables
      real(rprec) :: xvariables(nxc)
      logical :: isthere
      character(len=3) :: decision = 'no'
      character*200 :: temp
!-----------------------------------------------
 
!     INITIALIZE CONSTANTS, ARRAYS
 
      call initialize_coilsin

      lbatch = .true.
      mmax_bmn = 10
      nmax_bmn = 18
      nsmid = 0
      nsodd = 0
      m_in=20
      n_in=25

      nwire = nwdim
      nwire1 = nwire+1        
      mbwires = 0
      nvar = 0
 
!     OPEN INPUT FILE
  
      if (lbatch) then
         call safe_open(iunit, ierr, for05_file_new, 'old', 'formatted')
         if (ierr .ne. 0) then            !!for05_file_new is not there!
            call safe_open(iunit, ierr, input_data_file, 'old', 
     1         'formatted')
            if (ierr .eq. 0) then
               read (iunit, nml=indata, iostat=ierr)
               if (ierr .ne. 0) stop 'error reading namelist indata'
            else
               stop 'error opening input data file in xcoilgeom'
            end if
            nfp = nfp_vmec
         end if
      else
 
!     CHECK TO SEE IF A .NEW FILE ALREADY EXISTS
 
         inquire(file=for05_file_new, exist=isthere)
 
         if( myid .eq. master ) then                         !mpi stuff
            if( isthere ) then
               write(*,*)'File ',for05_file_new,' already exists!'
               write(*,*)
     1    'Do you wish to use it (and winding file) as the input? (Y/N)'
              read(*,*) decision
            end if
         end if     !if( myid .eq. master )                  !mpi stuff

#ifdef MPI_OPT
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)              !mpi stuff
#endif
!     OPEN FILE WITH INITIAL NAMELIST NAMIN
 
         if (decision(1:1) .eq. 'N' .or. decision(1:1) .eq. 'n') then
            call safe_open(iunit, ierr, for05_file, 'old', 
     1          'formatted')
        else
            call safe_open(iunit, ierr, for05_file_new, 'old', 
     1          'formatted')
        end if
      end if

      if (ierr .ne. 0) then
        if( myid .eq. master )  write(*,100) for05_file,ierr
        stop
      end if
 100  format(a,' => could not open file! Iostat = ',i4)

!     GENERATE MODE NUMBERS FOR B_ERROR SPECTRUM
 
      mnmax_bmn = (mmax_bmn+1)*(2*nmax_bmn/nfp+1)
      if (mnmax_bmn .ne. 0) then
        allocate (xm_bmn(mnmax_bmn), xn_bmn(mnmax_bmn), bmn(mnmax_bmn),
     1    stat = i)
         if( myid .eq. master ) then                         !mpi stuff
            if (i .ne. 0) stop 'allocation error in initialize_opt'
         end if     !if( myid .eq. master )                  !mpi stuff
        call get_modes (mmax_bmn, nmax_bmn, xm_bmn, xn_bmn, mnmax_bmn)
      end if

 
!     LOAD XC WITH INITIAL GUESS
 
!     IF BOUNDS ARE REQUIRED, USE THE UNCONSTRAINED VARIABLE
!     X AND DEFINE Y(X) = X/SQRT(1+X**2). THEN, VARY X AND LET
!     THE TARGET (R, FOR EXAMPLE) BE WRITTEN:
 
!        R(X) = .5*(RMAX + RMIN) + .5*Y(X)*(RMAX - RMIN)
 
      read (iunit,nml=coilsin,iostat=ierr)
      if (ierr .ne. 0)then
        if (myid .eq. master) then                         !mpi stuff
          print *,' ierr = ',ierr,' in initialize_opt reading coilsin'
        end if     !if( myid .eq. master )                  !mpi stuff
        stop
      endif

#if defined(GEOM_ONLY)
#else
      if (myid .eq. master) then                         !mpi stuff
         inquire(file="../" // bcoil_file, exist=isthere)
         if (isthere) then
            temp = "/bin/cp ../" // trim(bcoil_file) // 
     1                     " ./" // trim(bcoil_file)
            call system(temp)
         end if
      end if

#ifdef MPI_OPT
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                         !MPI
      if (ierr .ne. 0) stop 'MPI_BARRIER error in PROGRAM COILOPT4'
#endif

      if (nopt_alg .gt. 0) then
         call GA_preset
         call DE_preset
         read(iunit,nml=ga_de,iostat=ierr)
         niter_opt = ngen
         if (ierr .ne. 0) then
           if( myid .eq. master ) then                         !mpi stuff
              print *,' ierr = ',ierr,' in initialize_opt reading ga_de'
           end if     !if( myid .eq. master )                  !mpi stuff
           stop
         end if
      end if
#endif
 
!     Initialize and count modular coil variables
 
      if (lmodular) then
         call allocate_modular_coils
 
         call init_modular_coils (nvariables, xvariables)
         xc(1:nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
 
         call init_modular_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
!     Initialize and count saddle coil variables (if lsaddle = true)
 
      if (lsaddle) then
         call allocate_saddle_coils
 
         call init_saddle_coils (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
 
         call init_saddle_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
!     Initialize and count background coil variables (if lbcoil = true)
 
      if (lbcoil) then
         call read_bcoils
         call init_bg_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
!     Initialize and count vf coil variables (if lvf = true)
 
      if (lvf) then
         call allocate_vf_coils
         if (lvfvar) then
            call init_vf_coils (nvariables, xvariables)
            xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
            nvar = nvar + nvariables
         end if
 
         call init_vf_currents (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
!     Initialize and count 1/R coil variables (if ltfc = true)
 
      if (ltfc) then
         call init_tf_coils (nvariables, xvariables)
         if (ltfcv) then
            xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
            nvar = nvar + nvariables
         end if
      end if
 
!     Initialize modular winding surface variables (if lsurfv = true)
 
      if (lsurfv) then
         call init_modular_wsurf (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
!     Initialize saddle winding surface variables (if lsadsfv = true)
 
      if (lsadsfv) then
         call init_saddle_wsurf (nvariables, xvariables)
         xc(nvar+1:nvar+nvariables) = xvariables(1:nvariables)
         nvar = nvar + nvariables
      end if
 
      if (nvar .gt. nxc) stop 'nvar>nxc'
 
      close(iunit)
 
      end subroutine initialize_opt


      subroutine allocate_modular_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use modular_coils
      use Vcoilpts
      use Vwire
      use coils
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nc, i, j, k, status, modes
!-----------------------------------------------
 
!     Allocate modular coil arrays
 
      allocate(modular(nmod_coils_per_period), stat=status)
 
      if( myid .eq. master ) then                            !mpi stuff
        if(status /= 0) stop "Cannot allocate modular"
      end if     !if( myid .eq. master )                     !mpi stuff
 
      nc = nmod_coils_per_period
      if( myid .eq. master ) then                            !mpi stuff
        if (nc .le. 0) stop 'nmod_coils_per_period must > 0!'
      end if     !if( myid .eq. master )                     !mpi stuff
 
      nfper = nfp
      nmod_coils = nc * nfper
      if( myid .eq. master ) then                            !mpi stuff
        if (nmod_coils .gt. ncdim) stop 'nmod_coils > ncdim'
      end if     !if( myid .eq. master )                     !mpi stuff
 
      nmid = (nc+1)/2
      nodd = mod(nc,2)
      if( myid .eq. master ) then                            !mpi stuff
        print *, 'no. modular coils (per period) = ', nc
      end if     !if( myid .eq. master )                     !mpi stuff
 
      ncoils = nmod_coils                !total no. modular coils
      if (ncoils .le. 0) return
 
!     Initialize arrays to values of unique coil parameters
 
!     First consider the case with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).
 
      if ((nodd.eq.0) .and. (.not.lsymm)) then 
        nmid = nmid + 1
 
        if( myid .eq. master ) then                            !mpi stuff
          print *, 'nmid = ', nmid, 'nodd = ', nodd
        end if     !if( myid .eq. master )                     !mpi stuff
       
!     Allocate only those coil componants which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.
 
!     Allocate coil components for symmetry coil at phi = 0.
        i = 1
 
        allocate(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1           modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2           modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3           stat = status)
        if( myid .eq. master ) then                            !mpi stuff
          if(status /= 0) stop "Cannot allocate center coil components"
        end if     !if( myid .eq. master )                     !mpi stuff
          
!     Allocate coil components for coils 2 through nmid-1.
        do i = 2, nmid-1
       
          allocate(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             stat = status)
        if( myid .eq. master ) then                            !mpi stuff
          if(status /= 0) stop "Cannot allocate modular components"
        end if     !if( myid .eq. master )                     !mpi stuff
  
        end do               !! do i = ... loop over half of coils
 
!     Allocate coil components for symmetry coil at phi = pi/nfp.
        i = nmid
 
        allocate(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1           modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2           modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3           stat = status)
        if( myid .eq. master ) then                            !mpi stuff
          if(status /= 0) stop "Cannot allocate center coil components"
        end if     !if( myid .eq. master )                     !mpi stuff
          
!     Allocate coil components phi, rhi for coils nmid + 1 through nc.
        do i = 1, nmid-2
       
          allocate(modular(nmod_coils_per_period-i+1)%phi(nwire1),
     1             modular(nmod_coils_per_period-i+1)%rho(nwire1),
     2             stat = status)
        if( myid .eq. master ) then                            !mpi stuff
          if(status /= 0) stop "Cannot allocate modular components"
        end if     !if( myid .eq. master )                     !mpi stuff
  
        end do               !! do i = ... loop over half of coils
 
!     Symmetry coil at phi = 0.
        i = 1
        modes = 0
        phic(i,modes) = 0
        modular(i)%phic(modes) = 0
        modular(i)%phis(modes) = phis(i,modes)
        do modes = 1, nf_phi
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
        end do
 
        modes = 0
        modular(i)%rhos(modes) = 0
        do modes = 1, nf_rho
          modular(i)%rhos(modes) = rhos(i,modes)
        end do
 
!     Coils 2 through nmid-1.
        do i = 2, nmid-1
          modes = 0
          modular(i)%phic(modes) = phic(i,modes)
          modular(i)%phis(modes) = 0
          do modes = 1,nf_phi
            modular(i)%phic(modes) = phic(i,modes)
            modular(i)%phis(modes) = phis(i,modes)
          end do
 
          modes = 0
          modular(i)%rhos(modes) = 0
          do modes = 1,nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
          end do
        end do
 
!     Symmetry coil at phi = pi/nfp.
        i = nmid
        modes = 0
        phic(i,modes) = 0
        modular(i)%phic(modes) = 0
        modular(i)%phis(modes) = phis(i,modes)
        do modes = 1, nf_phi
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
        end do
 
        modes = 0
        modular(i)%rhos(modes) = 0
        do modes = 1, nf_rho
          modular(i)%rhos(modes) = rhos(i,modes)
        end do
 
!     Set currents from input
 
        do i = 1, nmid
          modular(i)%current = curmod(i)
        end do
 
!     Impose symmetry on coil-coil penalty weights, exponents
 
        do i = 1, nmid-2
          dcc_wgt(nc+1-i) = dcc_wgt(i+1)
          dcc_exp(nc+1-i) = dcc_exp(i+1)
          dcc_tgt(nc+1-i) = dcc_tgt(i+1)
          rc_wgt(nc+1-i) = rc_wgt(i+1)
          rc_exp(nc+1-i) = rc_exp(i+1)
          rc_tgt(nc+1-i) = rc_tgt(i+1)
          r_ext(nc+1-i) = r_ext(i+1)
        end do
 
!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), but not both. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).
 
      else
 
        if( myid .eq. master ) then                            !mpi stuff
          print *, 'nmid = ', nmid
          print *, 'nodd = ', nodd
        end if     !if( myid .eq. master )                     !mpi stuff
 
        do i = 1,nmid-nodd
       
!     Allocate only those coil componants which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.
         
          allocate(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             modular(nmod_coils_per_period-i+1)%phi(nwire1),
     4             modular(nmod_coils_per_period-i+1)%rho(nwire1),
     5             stat = status)
        if( myid .eq. master ) then                            !mpi stuff
          if(status /= 0) stop "Cannot allocate modular components"
        end if     !if( myid .eq. master )                     !mpi stuff
  
        end do               !! do i = ... loop over half of coils
 
        if (nodd .eq. 1) then
          i = nmid
 
          allocate(modular(i)%rhoc(0:nf_rho), modular(i)%rhos(0:nf_rho),
     1             modular(i)%phic(0:nf_phi), modular(i)%phis(0:nf_phi),
     2             modular(i)%phi(nwire1),    modular(i)%rho(nwire1),
     3             stat = status)
          if( myid .eq. master ) then                            !mpi stuff
           if(status /= 0)stop "Cannot allocate center coil components"
          end if     !if( myid .eq. master )                     !mpi stuff
        endif
 
        do i = 1, nmid-nodd
          modes = 0
          modular(i)%phic(modes) = phic(i,modes)
          modular(i)%phis(modes) = 0
          do modes = 1,nf_phi
            modular(i)%phic(modes) = phic(i,modes)
            modular(i)%phis(modes) = phis(i,modes)
          end do
 
          modes = 0
          modular(i)%rhos(modes) = 0
          do modes = 1,nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
          end do
        end do
 
        if (nodd .eq. 1) then
          i = nmid
          modes = 0
          phic(i,modes) = 0
          modular(i)%phic(modes) = 0
          modular(i)%phis(modes) = phis(i,modes)
          do modes = 1, nf_phi
            phic(i,modes) = 0
            modular(i)%phic(modes) = 0
            modular(i)%phis(modes) = phis(i,modes)
          end do
 
          modes = 0
          modular(i)%rhos(modes) = 0
          do modes = 1, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
          end do
        end if
 
!     Set currents from input
 
        do i = 1, nmid
          modular(i)%current = curmod(i)
        end do
 
!     Impose symmetry on coil-coil penalty weights, exponents
 
        do i = 1,nmid-nodd
          dcc_wgt(nc+1-i) = dcc_wgt(i)
          dcc_exp(nc+1-i) = dcc_exp(i)
          dcc_tgt(nc+1-i) = dcc_tgt(i)
          rc_wgt(nc+1-i) = rc_wgt(i)
          rc_exp(nc+1-i) = rc_exp(i)
          rc_tgt(nc+1-i) = rc_tgt(i)
          r_ext(nc+1-i) = r_ext(i)
        end do
 
      end if  ! end if ((nodd .eq. 0) .and. (lsymm .eqv. .false.))
 
      nmod_unique_coils = nmid
 
      end subroutine allocate_modular_coils


      subroutine eval_modular_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use modular_coils
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nc, np, icoil, i, j, k, n, itheta
      real(rprec) ::
     1   x(nwdim), y(nwdim), z(nwdim),
     2   u0, du, crv, di0, r, phi0, phi2, phi1,
     3   phi_tot, zw, rw, dtheta, theta, theta1, theta2,
     4   dssq, dr_ext
!-----------------------------------------------
 
!        LSYMM      NODD       Modular Configuration
!        _____      ____       _____________________
!          F          0        coils on both v=0 and v=1/2 symmetry planes
!          T          0        no coils on either v=0 or v=1/2 planes
!          T          1        coils on v=0 planes
!          F          1        coils on v=1/2 planes
 
      nc = ncoils / nfp                  !coils per field period
      if (ncoils .le. 0) return
 
!     offset in secular part of phi
      if (lsymm) then
!        symmetry coil at phi = 0 for nc odd
         di0 = (nc + 1)*.5_dp
      else
         if (nodd .eq. 1) then
!           symmetry coil at phi = pi/3 for nc odd
            di0 = 0.5_dp
         else
!           symmetry coils at phi = 0 and phi = pi/3 for nc even
            di0 = 1
         end if
      end if
 
      do n = 1, ncoils
        phi_full(n) = dtwopi*real(n - di0,rprec)/ncoils
      end do
 
!     Set coil segment intervals in poloidal angle
 
      dtheta = dtwopi/nwire
 
!     Zero the phi segment points in the modular coil structure
 
      do n = 1,nmid
        do i = 1,nwire1
          modular(n)%phi(i) = zero
        end do
      end do
 
!     First consider the case with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp. This implies that the number of coils
!     per field period must be even (nodd = 0).
 
      if ((nodd .eq. 0) .and. (lsymm .eqv. .false.)) then
 
!     For each modular, compute phi at nwire1 points equally spaced in theta
 
!     Symmetry coil at phi = 0.
        n = 1
        do itheta = 1, nwire1
          theta = dtheta*(itheta-1)
          if (nf_rho .gt. 0) call theta_f (n, theta, theta1, theta2)
          do k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                   + (modular(n)%phis(k)*sin(k*theta))
          end do
        end do
 
!     Coils 2 through nmid - 1
        do n = 2, nmid-1
          do itheta = 1, nwire1
            theta = dtheta*(itheta-1)
            if (nf_rho .gt. 0) call theta_f (n, theta, theta1, theta2)
            modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                             + modular(n)%phic(0)
            do k = 1,nf_phi
               modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                                +(modular(n)%phic(k)*cos(k*theta)
     2                                + modular(n)%phis(k)*sin(k*theta))
            end do
          end do
        end do
 
!     Symmetry coil at phi = pi/nfp.
        n = nmid
        do itheta = 1, nwire1
          theta = dtheta*(itheta-1)
          if (nf_rho .gt. 0) call theta_f (n, theta, theta1, theta2)
          do k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                   + (modular(n)%phis(k)*sin(k*theta))
          end do
        end do
 
!     Compute remaining phi for the field period using stellarator symmetry
 
        do i = 1, nmid
           curcon(i) = modular(i)%current
        end do
 
        do i = 1, nmid-2
        curcon(nc+1-i) = modular(i+1)%current
        do itheta = 1,nwire+1
        modular(nc+1-i)%phi(nwire1-itheta+1) =-modular(i+1)%phi(itheta)
        end do
        end do
 
!     Compute R and Z for first field period
 
        do n = 1, nc
          phi1 = modular(n)%phi(1)
          phi2 = phi1
          do itheta = 1,nwire
             theta = dtheta*(itheta-1)
             phi0 = modular(n)%phi(itheta)
             phi_tot = phi0 + phi_full(n)
             call rz_surf(theta,phi_tot,rw,zw,numsurf,
     1          rmn_sf,zmn_sf,m_num,n_num,nfper)
!            Extend symmetry and adjacent coils for NCSX NBI access
             if (lncsx) then
               call radial_ext (n, theta, dr_ext)
               rw = rw + dr_ext
             end if
             rcoil(n,itheta) = rw
             zcoil(n,itheta) = zw
             rcoil(n,itheta) = abs(rcoil(n,itheta))
             phi0 = modular(n)%phi(itheta)
             phi1 = min(phi1,phi0)
             phi2 = max(phi2,phi0)
          end do
          phimin(n) = phi1
          phimax(n) = phi2
        end do
 
!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), or on neither. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).
 
      else
 
!     For each modular, compute phi at nwire1 points equally spaced in theta
 
      do n = 1,nmid-nodd
        do itheta = 1,nwire1
          theta = dtheta*(itheta-1)
          if (nf_rho .gt. 0) call theta_f (n, theta, theta1, theta2)
          modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                           + modular(n)%phic(0)
          do k = 1,nf_phi
             modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                              +(modular(n)%phic(k)*cos(k*theta)
     2                              + modular(n)%phis(k)*sin(k*theta))
          end do
        end do
      end do
 
      if (nodd.eq.1) then               !central coil for odd no coils
         n = nmid
         do itheta = 1,nwire1
           theta = dtheta*(itheta-1)
           if (nf_rho .gt. 0) call theta_f (n, theta, theta1, theta2)
           do k = 1,nf_phi
              modular(n)%phi(itheta) = modular(n)%phi(itheta)
     1                    + (modular(n)%phis(k)*sin(k*theta))
           end do
         end do
      end if
 
!     Compute remaining phi for the field period using stellarator symmetry
 
      do i = 1, nmid
         curcon(i) = modular(i)%current
      end do
      do i = 1,nmid-nodd
         curcon(nc+1-i) = modular(i)%current
         do itheta = 1,nwire+1
           modular(nc+1-i)%phi(nwire1-itheta+1) =-modular(i)%phi(itheta)
         end do
      end do
 
!     Compute R and Z for first field period
 
      do n = 1, nc
        phi1 = modular(n)%phi(1)
        phi2 = phi1
        do itheta = 1,nwire
           theta = dtheta*(itheta-1)
           phi0 = modular(n)%phi(itheta)
           phi_tot = phi0 + phi_full(n)
           call rz_surf(theta,phi_tot,rw,zw,numsurf,
     1        rmn_sf,zmn_sf,m_num,n_num,nfper)
!          Extend symmetry and adjacent coils for NCSX NBI access
           if (lncsx) then
              call radial_ext (n, theta, dr_ext)
              rw = rw + dr_ext
           end if
 
           rcoil(n,itheta) = rw
           zcoil(n,itheta) = zw
           rcoil(n,itheta) = abs(rcoil(n,itheta))
           phi0 = modular(n)%phi(itheta)
           phi1 = min(phi1,phi0)
           phi2 = max(phi2,phi0)
        end do
        phimin(n) = phi1
        phimax(n) = phi2
      end do
 
      endif   ! end if ((nodd .eq. 0) .and. (lsymm .eqv. .false.))
 
!     Store currents for all field periods
 
      do j = 1,nfp-1
         do i = 1,nc
           curcon(i + j*nc) = curcon(i)
         end do
      end do
 
!     Store x,y,z coordinates for coils as x_mod(i,j,n), ..., where
!     i = segment number, j = filament number (=1 for central filament),
!     and n = coil number
 
      do 100 n = 1, ncoils
         icoil = 1 + mod(n-1,nc)              !coil index, mod nc
         ymin_cls(n) = 1000
         do i = 1,nwire
           r = rcoil(icoil,i)
           z(i) = zcoil(icoil,i)
           x(i) = r*cos(modular(icoil)%phi(i) + phi_full(n))
           y(i) = r*sin(modular(icoil)%phi(i) + phi_full(n))
           x_mod(i,1,n) = x(i)
           y_mod(i,1,n) = y(i)
           z_mod(i,1,n) = z(i)
           if (abs(y(i)) .lt. ymin_cls(n)) ymin_cls(n) = abs(y(i))
         end do
 
         x_mod(nwire1,1,n) = x_mod(1,1,n)
         y_mod(nwire1,1,n) = y_mod(1,1,n)
         z_mod(nwire1,1,n) = z_mod(1,1,n)
 
!     Compute length of coil n
 
         mod_length(n) = zero
         do i = 1,nwire-1
           dssq = (x(i+1) - x(i))**2 + (y(i+1) - y(i))**2
     1          + (z(i+1) - z(i))**2
           mod_length(n) = mod_length(n) + sqrt(dssq)
         end do
         dssq = (x(1) - x(nwire))**2 + (y(1) - y(nwire))**2
     1        + (z(1) - z(nwire))**2
         mod_length(n) = mod_length(n) + sqrt(dssq)
 
 100  continue
 
      end subroutine eval_modular_coils


      subroutine theta_f (icoil, th0, th1, th2)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use modular_coils
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, icoil
      real(rprec) :: theta, th0, th1, th2, t0, t1, t2
      real(rprec) :: rk, sk, ck
 
      i = icoil
      theta = th0
      t0 = theta
      t1 = 1
      t2 = 0
      do k = 1, nf_rho
         rk = modular(i)%rhos(k)
         sk = sin(k*theta)
         ck = cos(k*theta)
         t0 = t0 + rk*sk
         t1 = t1 + k*rk*ck
         t2 = t2 - k**2*rk*sk
      end do
      th0 = t0
      th1 = t1
      th2 = t2
 
      end subroutine theta_f


      subroutine radial_ext (ncoil, theta, dr_ext)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncoil
      real(rprec) :: theta, dr_ext
 
      dr_ext = 0
      if ((theta.le.dpi/2).or.(theta.ge.(3*dpi)/2))
     1  dr_ext = r_ext(ncoil)*cos(theta)
 
      end subroutine radial_ext


      subroutine load_modular_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary, only: nfp
      use modular_coils
      use tor_field
      use Vcoilpts
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nc, i, j, k, modes
!-----------------------------------------------
 
!     load the unique coil parameters with values from variables in
!     optimization
 
!     First consider the case with coils on both symmetry planes at
!     phi = 0 and phi = pi/nfp (lsymm = F). This implies that the number
!     of coils per field period must be even (nodd = 0).
 
      if ((nodd .eq. 0) .and. (lsymm .eqv. .false.)) then
 
!        Symmetry coil at phi = 0.
         i = 1
         do modes = 1, nf_phi
            modular(i)%phis(modes) = phis(i,modes)
         end do
 
         rhos(i,0) = 0
         do modes = 1, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
         end do
 
!        Coils 2 through nmid - 1
         do i = 2, nmid-1
            modes = 0
            modular(i)%phic(modes) = phic(i,modes)
            do modes = 1,nf_phi
               modular(i)%phic(modes) = phic(i,modes)
               modular(i)%phis(modes) = phis(i,modes)
            end do
 
            rhos(i,0) = 0
            do modes = 1,nf_rho
               modular(i)%rhos(modes) = rhos(i,modes)
            end do
         end do
 
!        Symmetry coil at phi = pi/nfp
         i = nmid
         do modes = 1, nf_phi
            modular(i)%phis(modes) = phis(i,modes)
         end do
 
         rhos(i,0) = 0
         do modes = 1, nf_rho
            modular(i)%rhos(modes) = rhos(i,modes)
         end do
 
         do i = 1, nmid
            modular(i)%current = curmod(i)
         end do
 
!     Next consider the cases with a coil on phi = 0 (lsymm = T), or
!     on phi = pi/nfp (lsymm = F), but not both. Then there may be an
!     even number of coils per period (nodd = 0) or an odd number of
!     coils per period (nodd = 1).
 
      else
 
         do i = 1, nmid-nodd
            modes = 0
            modular(i)%phic(modes) = phic(i,modes)
            do modes = 1,nf_phi
               modular(i)%phic(modes) = phic(i,modes)
               modular(i)%phis(modes) = phis(i,modes)
            end do
 
            rhos(i,0) = 0
            do modes = 1,nf_rho
               modular(i)%rhos(modes) = rhos(i,modes)
            end do
         end do
         if (nodd .eq. 1) then
            i = nmid
            do modes = 1, nf_phi
               modular(i)%phis(modes) = phis(i,modes)
            end do
 
            rhos(i,0) = 0
            do modes = 1, nf_rho
               modular(i)%rhos(modes) = rhos(i,modes)
            end do
         end if
 
         do i = 1, nmid
            modular(i)%current = curmod(i)
         end do
 
      endif   ! end if ((nodd .eq. 0) .and. (lsymm .eqv. .false.))
 
!     mod current scale factor for coil spacing, curvature penalties
 
      do i = 1, nmid
         cmod_scl(i) = abs(i_pol)/(one + abs(curmod(i)))
      end do

      end subroutine load_modular_structures


      subroutine allocate_saddle_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use saddle_coils
      use Vcoilpts
      use Vwire
      use coils
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, n, nc, ncoils, modes, status
!-----------------------------------------------
 
!     Allocate vf coil arrays
 
      allocate(saddle(nsad_coils_per_period), stat=status)
      if( myid .eq. master ) then                            !mpi stuff
        if(status /= 0) stop "Cannot allocate saddle"
      end if     !if( myid .eq. master )                     !mpi stuff
 
 
!     Allocate saddle coil x,y,z arrays
 
      allocate(x_sad(nwdim1,ncdim,nfdim), y_sad(nwdim1,ncdim,nfdim),
     1         z_sad(nwdim1,ncdim,nfdim), stat=status)
      if( myid .eq. master ) then                            !mpi stuff
        if(status /= 0) stop "Cannot allocate saddle coils"
      end if     !if( myid .eq. master )                     !mpi stuff
 
      allocate(u_sad(nwdim1,ncdim), v_sad(nwdim1,ncdim),
     1         stat=status)
      if( myid .eq. master ) then                            !mpi stuff
        if(status /= 0) stop "Cannot allocate saddle coils"
      end if     !if( myid .eq. master )                     !mpi stuff
 
      nc = nsad_coils_per_period
      if( myid .eq. master ) then                            !mpi stuff
        if (nc .le. 0) stop 'nsad_coils_per_period must > 0!'
      end if     !if( myid .eq. master )                     !mpi stuff

!     Number of coil types
      nsad_coils = nc * nfp
      if( myid .eq. master ) then                            !mpi stuff
        if (nsad_coils .gt. ncdim) stop 'nsad_coils > ncdim'
      end if     !if( myid .eq. master )                     !mpi stuff

      nsmid = nc/2
      nsodd = mod(nc,2)
      if( myid .eq. master ) then                            !mpi stuff
        if (nsodd .eq. 1) stop 'nsad_coils_per_period should be even'
      end if     !if( myid .eq. master )                     !mpi stuff
 
!     Number of unique coil currents
      num_cursad = 0
      do i = 1, nsmid
         if (nsad_group(i) .gt. num_cursad)
     1      num_cursad = nsad_group(i)
      end do
      if( myid .eq. master ) then                            !mpi stuff
        if (num_cursad .eq. 0) stop 'num_cursad = 0'
        do i = 1, nsmid
           if (nsad_group(i) .le. 0) stop 'nsad_group .le. 0'
        end do
      end if     !if( myid .eq. master )                     !mpi stuff

      do i = 1,nsmid
!     Allocate only those coil componants which are needed. There
!     is no need to store the Fourier coefficients for coils whose
!     locations in real space will be computed from stellarator
!     symmetry.
 
         allocate(saddle(i)%v_c(0:nsad_v),
     1            saddle(i)%v_s(0:nsad_v),
     2            saddle(i)%u_c(0:nsad_u),
     3            saddle(i)%u_s(0:nsad_u),
     4            stat = status)
      if( myid .eq. master ) then                            !mpi stuff
         if(status /= 0) stop "Cannot allocate saddle componants"
      end if     !if( myid .eq. master )                     !mpi stuff
 
      end do               !! do i = ... loop over half of coils
 
      ncoils = nsad_coils                !total no. saddle coils
      nc = ncoils / nfp                  !coils per field period
      if (ncoils .le. 0) return
 
!     initialize the variables to values of unique coil parameters
 
      do n = 1, nsmid
        sad_phi0(n) = twopi*sad_v0(n)/nfp
        sad_theta0(n) = twopi*sad_u0(n)
      end do
 
      do i = 1, nsmid
         modes = 0
         saddle(i)%v_c(modes) = sad_v_c(i,modes)
         saddle(i)%v_s(modes) = 0
         do modes = 1,nsad_v
            saddle(i)%v_c(modes) = sad_v_c(i,modes)
            saddle(i)%v_s(modes) = sad_v_s(i,modes)
         end do
 
         modes = 0
         saddle(i)%u_c(modes) = 0
         saddle(i)%u_s(modes) = 0
         if (nsad_u .gt. 0) then
            saddle(i)%u_c(modes) = sad_u_c(i,modes)
            do modes = 1,nsad_u
               saddle(i)%u_c(modes) = sad_u_c(i,modes)
               saddle(i)%u_s(modes) = sad_u_s(i,modes)
            end do
         end if
      end do
 
!     Set currents from input for each unique current group
 
      do i = 1, nsmid
         saddle(i)%current = cursad(nsad_group(i))*csad_scl(i)
      end do
 
!     Impose symmetry on saddle coil-coil penalty weights, exponents
 
      do i = 1,nsmid
         dsc_wgt(nc+1-i) = dsc_wgt(i)
         dsc_exp(nc+1-i) = dsc_exp(i)
         dsc_tgt(nc+1-i) = dsc_tgt(i)
         rs_wgt(nc+1-i) = rs_wgt(i)
         rs_exp(nc+1-i) = rs_exp(i)
         rs_tgt(nc+1-i) = rs_tgt(i)
      end do
 
      nsad_unique_coils = nsmid
 
      end subroutine allocate_saddle_coils


      subroutine eval_saddle_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary, only: nfp
      use saddle_coils
      use saddle_surface
      use Vcoilpts
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, m, nc, np, status
      real(rprec) :: u,  v,  s, ds, cmax, dlsq, csum, crad
      real(rprec), dimension(2,3,5) :: x,  y,  z
      real(rprec), dimension(2) :: c
      real(rprec) :: r, alf, bet
!-----------------------------------------------
 
      alf = 1
      bet = 0
 
      nc = nsad_coils_per_period
      ds = 1.0_dp/nwire
!     loop over unique coils
      do n=1, nsmid
         m = nc + 1 - n
!        loop over number of poloidal points
         cmax = 0
         csum = 0
         do i=1, nwire + 1
            s = (i - 1)*ds
            call coil_curve(n, s, alf, bet, u, v, r, c, x, y, z)
            cmax = max (cmax, c(1))
            csum = csum + c(1)
!           set coil x,y,z points for all field periods
            do j=1, nfp
               np = (j-1)*nc
               u_sad(i,n + np) = u
               v_sad(i,n + np) = v + 1.0*(j-1)
!              stellarator symmetry
               u_sad(nwire + 2 - i,m + np) = 1.0 - u
               v_sad(nwire + 2 - i,m + np) = -v + 1.0*(j-1)
               do k=1, nfils
                  x_sad(i,n + np,k) = x(1,j,k)
                  y_sad(i,n + np,k) = y(1,j,k)
                  z_sad(i,n + np,k) = z(1,j,k)
!                 stellarator symmetry
                  x_sad(i,m + np,k) = x(2,j,k)
                  y_sad(i,m + np,k) = y(2,j,k)
                  z_sad(i,m + np,k) = z(2,j,k)
               end do
            end do
         end do
         rs_min(n) = 10000.0_dp
         if(cmax .gt. 0.0_dp) rs_min(n) = 1.0_dp/cmax
         cs_sum(n) = csum*ds
      end do
 
!     coil currents
      do n=1, nsmid
         m = nc + 1 - n
         do j=1, nfp
            np = (j-1)*nc
            c_sad(n + np) =  saddle(n)%current
            c_sad(m + np) = -saddle(n)%current
         end do
      end do
 
!     coil lengths
      do n=1, nsmid
         ymin_sad(n) = 1000
         rmax_sad(n) = 0
         sad_length(n) = 0
         do i=1, nwire
            dlsq = (x_sad(i+1,n,1) - x_sad(i,n,1))**2
     1           + (y_sad(i+1,n,1) - y_sad(i,n,1))**2
     2           + (z_sad(i+1,n,1) - z_sad(i,n,1))**2
            sad_length(n) = sad_length(n) + sqrt(dlsq)
            if (abs(y_sad(i,n,1)) .lt. ymin_sad(n))
     1         ymin_sad(n) = abs(y_sad(i,n,1))
            crad = sqrt(x_sad(i,n,1)**2+y_sad(i,n,1)**2)
            if (crad .gt. rmax_sad(n)) rmax_sad(n) = crad
         end do
      end do
 
      end subroutine eval_saddle_coils


      subroutine eval_poloidal_currents
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary
      use modular_coils
      use saddle_coils
      use bcoils_mod
      use tor_field
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n
      real (rprec) :: sum0, csum0
!-----------------------------------------------

!     Evaluate the sum0 of the poloidal currents

      csum0 = 0

      if (lmodular) then
         sum0 = 0
         if ((nodd .eq. 0) .and. (.not. lsymm)) then
            do i = 2, nmid - 1
               sum0 = sum0 + curmod(i)
            end do
            sum0 = 2*sum0 + curmod(1) + curmod(nmid)
         else
            if (nodd .eq. 0) then
               do i = 1, nmid
                  sum0 = sum0 + curmod(i)
               end do
               sum0 = 2*sum0
            end if
            if (lsymm) then
               do i = 1, nmid - 1
                  sum0 = sum0 + curmod(i)
               end do
               sum0 = 2*sum0 + curmod(nmid)
            else
               do i = 2, nmid
                  sum0 = sum0 + curmod(i)
               end do
               sum0 = 2*sum0 + curmod(1)
            end if
         end if
         csum0 = csum0 + nfp*sum0
      end if

      if (lsaddle .and. lsmod) then
         sum0 = 0
         do i = 1, nsmid
            sum0 = sum0 + cursad(nsad_group(i))*csad_scl(i)
         end do
         csum0 = csum0 + 2*nfp*sum0
      end if

      if (lbcoil) then
         sum0 = 0
         do i = 1, mbcoils
            if (lp_bg(i)) sum0 = sum0 + bcoil_cur(i)
         end do
         csum0 = csum0 + sum0
      end if

      if (ltfc) then
         csum0 = csum0 + i_tfc
      end if

      pol_cur = csum0

      end subroutine eval_poloidal_currents


      subroutine coil_curve(n, s, alf, bet, u, v, r, c, x, y, z)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary, only: nfp
      use saddle_coils
      use saddle_surface
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, nf, np, ncl, mk, nk, nu, nv, ifail
      real(rprec), dimension(2,3,5) :: x, y, z
      real(rprec), dimension(2) :: x0, x1, x2
      real(rprec), dimension(2) :: y0, y1, y2
      real(rprec), dimension(2) :: tcx, tcy, tcz
      real(rprec), dimension(2) :: tsx, tsy, tsz
      real(rprec), dimension(2) :: xu, xv, yu, yv
      real(rprec), dimension(2) :: nx, ny, nz, nr, nphi, dn, d1, c
      real(rprec), dimension(2) :: ph0,  ph1,  ph2, phv
      real(rprec), dimension(4) :: fu, fv
      real(rprec), dimension(mfourier) :: vc, vs, uc, us
      real(rprec) :: s, r, u,  v, dph, alf, bet
      real(rprec) :: s0, s1
      real(rprec) :: ck, sk
      real(rprec) :: r0, r1, r2
      real(rprec) :: z0, z1, z2
      real(rprec) :: u0,  u1,  u2
      real(rprec) :: v0,  v1,  v2
      real(rprec) :: ru, rv, zu, zv
!-----------------------------------------------
 
      ncl = n
      dph = twopi/nfp
      s0 = twopi*s
      s1 = twopi
 
!     compute toroidal variable v and derivatives wrt s
      if (lspline) then

         ! b-spline series representation
         nv = nsad_v
         v0 = saddle(ncl)%v_c(0)
         v1 = 0
         v2 = 0
         ! coefficients
         if (nv .ge. 7) then
            do k = 1, nv
               vc(k) = saddle(ncl)%v_c(k)
            end do
            ! breakpoints
            do k = 5, nv
               vs(k) = saddle(ncl)%v_s(k)
            end do
            ! boundary conditions
            do k = 1, 4
               vs(k) = 0
               vs(nv + k) = 1
            end do
            call speval (nv, vs, vc, s, fv, ifail)
            v0 = v0 + fv(1)
            v1 = fv(2)
            v2 = fv(3)
         end if

      else

         ! fourier series representation
         v0 = 0
         v1 = 0
         v2 = 0
         do k = 0,nsad_v
            ck = cos(k*s0)
            sk = sin(k*s0)
            v0 = v0 +  saddle(ncl)%v_c(k)*ck
     1              +  saddle(ncl)%v_s(k)*sk
            v1 = v1 + (saddle(ncl)%v_s(k)*ck
     1              -  saddle(ncl)%v_c(k)*sk)*k
            v2 = v2 + (saddle(ncl)%v_s(k)*sk
     1              +  saddle(ncl)%v_c(k)*ck)*k**2
         end do
         v1 = s1*v1
         v2 = -s1**2*v2

      end if           ! if (lspline)
 
!     compute poloidal variable u and derivatives wrt s
      if (lspline) then

         ! b-spline series representation
         nu = nsad_u
         u0 = 0
         u1 = 0
         u2 = 0
         ! coefficients
         if (nu .ge. 7) then
            do k = 1, nu
               uc(k) = saddle(ncl)%u_c(k)
            end do
            ! breakpoints
            do k = 5, nu
               us(k) = saddle(ncl)%u_s(k)
            end do
            ! boundary conditions
            do k = 1, 4
               us(k) = 0
               us(nu + k) = 1
            end do
            call speval (nu, us, uc, s, fu, ifail)
            u0 = fu(1)
            u1 = fu(2)
            u2 = fu(3)
         end if

      else

         ! fourier series representation
         u0 = 0
         u1 = 0
         u2 = 0
         do k = 0,nsad_u
            ck = cos(k*s0)
            sk = sin(k*s0)
            u0 = u0 +  saddle(ncl)%u_c(k)*ck
     1              +  saddle(ncl)%u_s(k)*sk
            u1 = u1 + (saddle(ncl)%u_s(k)*ck
     1              -  saddle(ncl)%u_c(k)*sk)*k
            u2 = u2 + (saddle(ncl)%u_s(k)*sk
     1              +  saddle(ncl)%u_c(k)*ck)*k**2
         end do
 
         u1 = s1*u1
         u2 = -s1**2*u2

      end if           ! if (lspline)

      if (lsmod) then
!        lsmod = T => modular
         u0 = u0 + s
         u1 = u1 + 1
      end if
 
!     compute R, Z and derivatives wrt s
 
      r0 = 0
      r1 = 0
      r2 = 0
      z0 = 0
      z1 = 0
      z2 = 0
      ru = 0
      rv = 0
      zu = 0
      zv = 0
      do k = 1, numsurf_sad
         mk = m_sad(k)
!        nk = nfp*n_sad(k)
         nk = n_sad(k)
         ck = cos(twopi*(mk*u0 + nk*v0))
         sk = sin(twopi*(mk*u0 + nk*v0))
!     R ...
         r0 = r0 + rmn_sad(k)*ck
         r1 = r1 - rmn_sad(k)*(twopi*(mk*u1 + nk*v1))*sk
         r2 = r2 - rmn_sad(k)*((twopi*(mk*u2 + nk*v2))*sk
     1      + (twopi*(mk*u1 + nk*v1))**2*ck)
         ru = ru - rmn_sad(k)*mk*sk*twopi
         rv = rv - rmn_sad(k)*nk*sk*twopi
!     Z ...
         z0 = z0 + zmn_sad(k)*sk
         z1 = z1 + zmn_sad(k)*(twopi*(mk*u1 + nk*v1))*ck
         z2 = z2 + zmn_sad(k)*((twopi*(mk*u2 + nk*v2))*ck
     1      - (twopi*(mk*u1 + nk*v1))**2*sk)
         zu = zu + zmn_sad(k)*mk*ck*twopi
         zv = zv + zmn_sad(k)*nk*ck*twopi
      end do
 
!     derivatives of toroidal angle wrt s
      ph1(1) =  twopi*v1/nfp
      ph1(2) = -twopi*v1/nfp
      ph2(1) =  twopi*v2/nfp
      ph2(2) = -twopi*v2/nfp
      phv(1) =  twopi/nfp
      phv(2) = -twopi/nfp
 
!     field period invariance
      do np = 1,nfp
 
!     toroidal angle in field period np
      ph0(1) =  twopi*v0/nfp + (np - 1)*dph
      ph0(2) = -twopi*v0/nfp + (np - 1)*dph
!     print 100, u0, v0, ph0(1), ph0(2), r0, z0
!     print 110
! 100 format (1p6e13.4)
! 110 format (" ")
 
!     stellarator symmetry
      do i = 1,2
!     X and derivatives wrt s
      x0(i) = r0*cos(ph0(i))
      x1(i) = r1*cos(ph0(i)) - r0*ph1(i)*sin(ph0(i))
      x2(i) = r2*cos(ph0(i)) - r1*ph1(i)*sin(ph0(i))
     1         - (r0*ph1(i)**2*cos(ph0(i))
     2         + (r0*ph2(i) + r1*ph1(i))*sin(ph0(i)))
 
!     X derivatives wrt u, v
 
      xu(i) = ru*cos(ph0(i))
      xv(i) = rv*cos(ph0(i)) - r0*sin(ph0(i))*phv(i)
 
!     Y and derivatives wrt s
 
      y0(i) = r0*sin(ph0(i))
      y1(i) = r1*sin(ph0(i)) + r0*ph1(i)*cos(ph0(i))
      y2(i) = r2*sin(ph0(i)) + r1*ph1(i)*cos(ph0(i))
     1         + (-r0*ph1(i)**2*sin(ph0(i))
     2         + (r0*ph2(i) + r1*ph1(i))*cos(ph0(i)))
 
!     Y derivatives wrt u, v
 
      yu(i) = ru*sin(ph0(i))
      yv(i) = rv*sin(ph0(i)) + r0*cos(ph0(i))*phv(i)
 
!     compute |dr/ds| = dl/ds
 
      d1(i) = dsqrt (x1(i)*x1(i) + y1(i)*y1(i) + z1*z1)
 
!     end loop for stellarator symmetry
      end do
 
!     tangent vector (to curve) tcx=dx/ds, tcy=dy/ds, tcz=dz/ds
 
      tcx(1) =  x1(1)/d1(1)
      tcy(1) =  y1(1)/d1(1)
      tcz(1) =  z1/d1(1)
!     stellarator symmetry
      tcx(2) =  x1(2)/d1(2)
      tcy(2) =  y1(2)/d1(2)
      tcz(2) = -z1/d1(2)
 
!     components of surface normal vector (cylindrical) nr, nphi, nz
 
      nr(1) =  r0*zu
      nphi(1) =  ru*zv - rv*zu
      nz(1) = -r0*ru
      dn(1) =  dsqrt (nr(1)*nr(1) + nphi(1)*nphi(1) + nz(1)*nz(1))
      nr(1) =  nr(1)/dn(1)
      nphi(1) = nphi(1)/dn(1)
      nz(1) =  nz(1)/dn(1)
!     stellarator symmetry
      nr(2) = -r0*zu
      nphi(2) = -ru*zv + rv*zu
      nz(2) = -r0*ru
      dn(2) =  dsqrt (nr(2)*nr(2) + nphi(2)*nphi(2) + nz(2)*nz(2))
      nr(2) =  nr(2)/dn(2)
      nphi(2) = nphi(2)/dn(2)
      nz(2) =  nz(2)/dn(2)
 
!     components of surface normal vector (cartesian) nx, ny, nz
 
      nx(1) =  yu(1)*zv - zu*yv(1)
      ny(1) =  zu*xv(1) - zv*xu(1)
      nz(1) =  xu(1)*yv(1) - xv(1)*yu(1)
      dn(1) =  dsqrt (nx(1)*nx(1) + ny(1)*ny(1) + nz(1)*nz(1))
      nx(1) =  nx(1)/dn(1)
      ny(1) =  ny(1)/dn(1)
      nz(1) =  nz(1)/dn(1)
!     stellarator symmetry
      nx(2) = -yu(2)*zv + zu*yv(2)
      ny(2) = -zu*xv(2) + zv*xu(2)
      nz(2) =  xu(2)*yv(2) - xv(2)*yu(2)
      dn(2) =  dsqrt (nx(2)*nx(2) + ny(2)*ny(2) + nz(2)*nz(2))
      nx(2) =  nx(2)/dn(2)
      ny(2) =  ny(2)/dn(2)
      nz(2) =  nz(2)/dn(2)
 
!     tangent vector (to surface) tsx, tsy, tsz
 
      tsx(1) = ny(1)*tcz(1) - tcy(1)*nz(1)
      tsy(1) = nz(1)*tcx(1) - tcz(1)*nx(1)
      tsz(1) = nx(1)*tcy(1) - tcx(1)*ny(1)
      dn(1) = dsqrt (tsx(1)*tsx(1) + tsy(1)*tsy(1) + tsz(1)*tsz(1))
      tsx(1) = tsx(1)/dn(1)
      tsy(1) = tsy(1)/dn(1)
      tsz(1) = tsz(1)/dn(1)
!     stellarator symmetry
      tsx(2) = ny(2)*tcz(2) - tcy(2)*nz(2)
      tsy(2) = nz(2)*tcx(2) - tcz(2)*nx(2)
      tsz(2) = nx(2)*tcy(2) - tcx(2)*ny(2)
      dn(2) = dsqrt (tsx(2)*tsx(2) + tsy(2)*tsy(2) + tsz(2)*tsz(2))
      tsx(2) = tsx(2)/dn(2)
      tsy(2) = tsy(2)/dn(2)
      tsz(2) = tsz(2)/dn(2)
 
!     return curvature c = |dr/ds X d2r/ds2|/|dr/ds|^3
 
      c(1) = dsqrt ((y1(1)*z2 - y2(1)*z1)**2
     1         + (z1*x2(1) - z2*x1(1))**2
     2         + (x1(1)*y2(1) - x2(1)*y1(1))**2)/d1(1)**3
!     stellarator symmetry
      c(2) = dsqrt ((-y1(2)*z2 + y2(2)*z1)**2
     1         + (-z1*x2(2) + z2*x1(2))**2
     2         + (x1(2)*y2(2) - x2(2)*y1(2))**2)/d1(2)**3
 
!     return multifilament representation for x, y, z
!     do nf = 1,nfils
!     filament 1
      nf = 1
      x(1,np,nf) =  x0(1)
      y(1,np,nf) =  y0(1)
      z(1,np,nf) =  z0
!     stellarator symmetry
      x(2,np,nf) =  x0(2)
      y(2,np,nf) =  y0(2)
      z(2,np,nf) = -z0
      if (nfils .lt. 4) then
!     Three filament model
!        filament 2
         nf = 2
         x(1,np,nf) =  x0(1) + deln*(alf*nx(1) + bet*tsx(1))
         y(1,np,nf) =  y0(1) + deln*(alf*ny(1) + bet*tsy(1))
         z(1,np,nf) =  z0    + deln*(alf*nz(1) + bet*tsz(1))
!        stellarator symmetry
         x(2,np,nf) =  x0(2) + deln*(alf*nx(2) + bet*tsx(2))
         y(2,np,nf) =  y0(2) + deln*(alf*ny(2) + bet*tsy(2))
         z(2,np,nf) = -z0    + deln*(alf*nz(2) + bet*tsz(2))
!        filament 3
         nf = 3
         x(1,np,nf) =  x0(1) - deln*(alf*nx(1) + bet*tsx(1))
         y(1,np,nf) =  y0(1) - deln*(alf*ny(1) + bet*tsy(1))
         z(1,np,nf) =  z0    - deln*(alf*nz(1) + bet*tsz(1))
!        stellarator symmetry
         x(2,np,nf) =  x0(2) - deln*(alf*nx(2) + bet*tsx(2))
         y(2,np,nf) =  y0(2) - deln*(alf*ny(2) + bet*tsy(2))
         z(2,np,nf) = -z0    - deln*(alf*nz(2) + bet*tsz(2))
      else
!     Four filament model (no current in central filament 1)
!        filament 2
         nf = 2
         x(1,np,nf) =  x0(1) + deln*nx(1) + delt*tsx(1)
         y(1,np,nf) =  y0(1) + deln*ny(1) + delt*tsy(1)
         z(1,np,nf) =  z0    + deln*nz(1) + delt*tsz(1)
!        stellarator symmetry
         x(2,np,nf) =  x0(2) + deln*nx(2) + delt*tsx(2)
         y(2,np,nf) =  y0(2) + deln*ny(2) + delt*tsy(2)
         z(2,np,nf) = -z0    + deln*nz(2) + delt*tsz(2)
!        filament 3
         nf = 3
         x(1,np,nf) =  x0(1) - deln*nx(1) + delt*tsx(1)
         y(1,np,nf) =  y0(1) - deln*ny(1) + delt*tsy(1)
         z(1,np,nf) =  z0    - deln*nz(1) + delt*tsz(1)
!        stellarator symmetry
         x(2,np,nf) =  x0(2) - deln*nx(2) + delt*tsx(2)
         y(2,np,nf) =  y0(2) - deln*ny(2) + delt*tsy(2)
         z(2,np,nf) = -z0    - deln*nz(2) + delt*tsz(2)
!        filament 4
         nf = 4
         x(1,np,nf) =  x0(1) - deln*nx(1) - delt*tsx(1)
         y(1,np,nf) =  y0(1) - deln*ny(1) - delt*tsy(1)
         z(1,np,nf) =  z0    - deln*nz(1) - delt*tsz(1)
!        stellarator symmetry
         x(2,np,nf) =  x0(2) - deln*nx(2) - delt*tsx(2)
         y(2,np,nf) =  y0(2) - deln*ny(2) - delt*tsy(2)
         z(2,np,nf) = -z0    - deln*nz(2) - delt*tsz(2)
!        filament 5
         nf = 5
         x(1,np,nf) =  x0(1) + deln*nx(1) - delt*tsx(1)
         y(1,np,nf) =  y0(1) + deln*ny(1) - delt*tsy(1)
         z(1,np,nf) =  z0    + deln*nz(1) - delt*tsz(1)
!        stellarator symmetry
         x(2,np,nf) =  x0(2) + deln*nx(2) - delt*tsx(2)
         y(2,np,nf) =  y0(2) + deln*ny(2) - delt*tsy(2)
         z(2,np,nf) = -z0    + deln*nz(2) - delt*tsz(2)
      end if
 
!     end loop for multifilament representation
!     end do
 
!     end loop for field period invariance
      end do
 
!     return u, v, r
      u = u0
      v = v0
      r = r0
 
      end subroutine coil_curve


      subroutine speval(n,xk,c,t,f,ifail)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: n, ifail
      real(rprec), dimension(*) :: xk, c, f
      real(rprec) :: t
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l, ip1
      real(rprec) :: a10, a11, a12, a20, a21
      real(rprec) :: b02, b12, b03, b13, b23, b04, b14, b24, b34
      real(rprec) :: d10, d11, d12, d20, d21, d30
      real(rprec) :: tm0, tm1, tm2, tp1, tp2, tp3
 
!     evaluate the cubic spline with knots xk and b-spline
!     coefficients c at the point t. return function value
!     in f(1) and derivatives in f(2)-f(4).
 
      ifail=1
      do j=1,4
         f(j)=0
      end do
      if((t.lt.xk(4)) .or. (t.gt.xk(n+1))) return
      i=4
      ip1=n+1
   20 l=(i+ip1)/2
      if(ip1-i.le.1) go to 40
      if(t.lt.xk(l)) go to 30
      i=l
      go to 20
   30 ip1=l
      go to 20
   40 tm2=t-xk(i-2)
      tm1=t-xk(i-1)
      tm0=t-xk(i)
      tp1=xk(i+1)-t
      tp2=xk(i+2)-t
      tp3=xk(i+3)-t
      d10=tp1+tm0
      d11=tp1+tm1
      d20=tp2+tm0
      d12=tp1+tm2
      d21=tp2+tm1
      d30=tp3+tm0
      b12=tp1/d10
      b02=tm0/d10
      b23=tp1*b12/d11
      b13=tm1*b12/d11+tp2*b02/d20
      b03=tm0*b02/d20
      b34=tp1*b23/d12
      b24=tm2*b23/d12+tp2*b13/d21
      b14=tm1*b13/d21+tp3*b03/d30
      b04=tm0*b03/d30
      f(1)=c(i-3)*b34+c(i-2)*b24+c(i-1)*b14+c(i)*b04
      a12=(c(i-2)-c(i-3))/d12
      a11=(c(i-1)-c(i-2))/d21
      a10=(c(i)-c(i-1))/d30
      f(2)=3*(a12*b23+a11*b13+a10*b03)
      a21=(a11-a12)/d11
      a20=(a10-a11)/d20
      f(3)=6*(a21*b12+a20*b02)
      f(4)=6*(a20-a21)/d10
      ifail=0

      end subroutine speval


      subroutine load_saddle_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary, only: nfp
      use saddle_coils
      use tor_field
      use Vcoilpts
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, modes, nv, nu, ifail
      real(rprec), dimension(4) :: fv, fu
      real(rprec), dimension(mfourier) :: vc0, vc, vs, uc0, uc, us
      real(rprec) :: b1_0_0, b1_1_0, b1_2_0, b2_1_0, b2_2_0, b3_2_0
      real(rprec) :: bnm0_0_1, bnm0_1_1, bnm0_2_1
      real(rprec) :: bnm1_1_1, bnm1_2_1, bnm2_2_1
      real(rprec) :: s0, s1
!-----------------------------------------------
 
!     load the coil structures with values from variables in
!     optimization
 
      do i = 1, nsmid
         modes = 0
         saddle(i)%v_c(modes) = sad_v_c(i,modes)
         sad_v_s(i,modes) = 0
         saddle(i)%v_s(modes) = 0
         do modes = 1,nsad_v
            saddle(i)%v_c(modes) = sad_v_c(i,modes)
            saddle(i)%v_s(modes) = sad_v_s(i,modes)
         end do
 
         if (nsad_u .gt. 0) then
            modes = 0
            saddle(i)%u_c(modes) = sad_u_c(i,modes)
            sad_u_s(i,modes) = 0
            saddle(i)%u_s(modes) = 0
            do modes = 1,nsad_u
               saddle(i)%u_c(modes) = sad_u_c(i,modes)
               saddle(i)%u_s(modes) = sad_u_s(i,modes)
            end do
         end if
      end do
 
!     saddle coil currents for each unique current group
      do i=1, nsmid
         saddle(i)%current = cursad(nsad_group(i))*csad_scl(i)
      end do

      if (lspline .and. (nsad_v .ge. 7)) then
!     set boundary conditions for spline representation of v

         nv = nsad_v

         do i = 1, nsmid
            ! breakpoints
            do k = 5, nv
               vs(k) = saddle(i)%v_s(k)
            end do
            ! boundary conditions on breakpoints
            do k = 1, 4
               vs(k) = 0
               vs(nv + k) = 1
            end do

            ! boundary conditions on coefficients
            s0 = 0
            s1 = 1
            vc0 = 0

            vc0(1) = 1
            call speval (nv, vs, vc0, s0, fv, ifail)
            b1_0_0 = fv(1)
            b1_1_0 = fv(2)
            b1_2_0 = fv(3)
            vc0(1) = 0

            vc0(2) = 1
            call speval (nv, vs, vc0, s0, fv, ifail)
            b2_1_0 = fv(2)
            b2_2_0 = fv(3)
            vc0(2) = 0

            vc0(3) = 1
            call speval (nv, vs, vc0, s0, fv, ifail)
            b3_2_0 = fv(3)
            vc0(3) = 0

            vc0(nv) = 1
            call speval (nv, vs, vc0, s1, fv, ifail)
            bnm0_0_1 = fv(1)
            bnm0_1_1 = fv(2)
            bnm0_2_1 = fv(3)
            vc0(nv) = 0

            vc0(nv-1) = 1
            call speval (nv, vs, vc0, s1, fv, ifail)
            bnm1_1_1 = fv(2)
            bnm1_2_1 = fv(3)
            vc0(nv-1) = 0

            vc0(nv-2) = 1
            call speval (nv, vs, vc0, s1, fv, ifail)
            bnm2_2_1 = fv(3)
            vc0(nv-2) = 0

            vc = 0
            do k = 1, 3
               vc(k) = saddle(i)%v_c(k)
            end do
            vc(nv) = vc(1)*(b1_0_0/bnm0_0_1)
            vc(nv-1) = (vc(1)*b1_1_0 + vc(2)*b2_1_0 - vc(nv)*bnm0_1_1)
     1                 /bnm1_1_1
            vc(nv-2) = (vc(1)*b1_2_0 + vc(2)*b2_2_0 + vc(3)*b3_2_0
     1                - vc(nv-1)*bnm1_2_1 - vc(nv)*bnm0_2_1)/bnm2_2_1
            do k = nv-2, nv
               saddle(i)%v_c(k) = vc(k)
               sad_v_c(i,k) = vc(k)
            end do

         end do     ! i=1,nsmid

      end if        ! if (lspline .and. (nsad_v .ge. 7))
 
      if (lspline .and. (nsad_u .ge. 7)) then
!     set boundary conditions for spline representation of u

         nu = nsad_u

         do i = 1, nsmid
            ! breakpoints
            do k = 5, nu
               us(k) = saddle(i)%u_s(k)
            end do
            ! boundary conditions on breakpoints
            do k = 1, 4
               us(k) = 0
               us(nu + k) = 1
            end do

            ! boundary conditions on coefficients
            s0 = 0
            s1 = 1
            uc0 = 0

            uc0(1) = 1
            call speval (nu, us, uc0, s0, fu, ifail)
            b1_0_0 = fu(1)
            b1_1_0 = fu(2)
            b1_2_0 = fu(3)
            uc0(1) = 0

            uc0(2) = 1
            call speval (nu, us, uc0, s0, fu, ifail)
            b2_1_0 = fu(2)
            b2_2_0 = fu(3)
            uc0(2) = 0

            uc0(3) = 1
            call speval (nu, us, uc0, s0, fu, ifail)
            b3_2_0 = fu(3)
            uc0(3) = 0

            uc0(nu) = 1
            call speval (nu, us, uc0, s1, fu, ifail)
            bnm0_0_1 = fu(1)
            bnm0_1_1 = fu(2)
            bnm0_2_1 = fu(3)
            uc0(nu) = 0

            uc0(nu-1) = 1
            call speval (nu, us, uc0, s1, fu, ifail)
            bnm1_1_1 = fu(2)
            bnm1_2_1 = fu(3)
            uc0(nu-1) = 0

            uc0(nu-2) = 1
            call speval (nu, us, uc0, s1, fu, ifail)
            bnm2_2_1 = fu(3)
            uc0(nu-2) = 0

            uc = 0
            do k = 1, 3
               uc(k) = saddle(i)%u_c(k)
            end do
            uc(nu) = uc(1)*(b1_0_0/bnm0_0_1)
            uc(nu-1) = (uc(1)*b1_1_0 + uc(2)*b2_1_0 - uc(nu)*bnm0_1_1)
     1                 /bnm1_1_1
            uc(nu-2) = (uc(1)*b1_2_0 + uc(2)*b2_2_0 + uc(3)*b3_2_0
     1                - uc(nu-1)*bnm1_2_1 - uc(nu)*bnm0_2_1)/bnm2_2_1
            do k = nu-2, nu
               saddle(i)%u_c(k) = uc(k)
               sad_u_c(i,k) = uc(k)
            end do

         end do     ! i=1,nsmid

      end if        ! if (lspline .and. (nsad_u .ge. 7))

      end subroutine load_saddle_structures


      subroutine allocate_vf_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      use coils
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: status
!-----------------------------------------------
 
!     Allocate vf coil arrays
 
      allocate(vertical(num_vf), stat=status)
      if( myid .eq. master ) then                            !mpi stuff
        if(status /= 0) stop "Cannot allocate vertical"
      end if     !if( myid .eq. master )                     !mpi stuff

      end subroutine allocate_vf_coils


      subroutine eval_vf_coils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary
      use vf_coils
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, j1, j2
      real(rprec) :: x(nwdim1), y(nwdim1), z(nwdim1)
      real(rprec) :: phi, dphi, rvfc, zvfc, argt, rvfmax
!-----------------------------------------------
 
      dphi = dtwopi/nwire
      nvf = 2*num_vf
 
      do n = 1, num_vf
         j1 = 2*(n - 1) + 1
         j2 = j1 + 1
         zvfc = zc_vf(n)
         rvfc = rc_vf(n)
         rvfmax = rc_vf(n)
         do i = 1, nwire + 1
            phi = (i-1)*dphi
            if (nrvf_c .gt. 0) then
               do k = 1, nrvf_c
                  argt = k*nfp*phi
                  rvfc = rvfc + rcfc_vf(n,k)*cos(argt)
     1                        + rcfs_vf(n,k)*sin(argt)
               end do
               if (rvfc .gt. rvfmax) rvfmax = rvfc
            end if

            x_vf(i,1,j1) = rvfc*cos(phi)
            y_vf(i,1,j1) = rvfc*sin(phi)
            z_vf(i,1,j1) = zvfc

            x_vf(i,1,j2) = rvfc*cos(-phi)
            y_vf(i,1,j2) = rvfc*sin(-phi)
            z_vf(i,1,j2) = -zvfc
         end do
         rvf_max(n) = rvfmax
      end do

!     vf coil currents
      k = 0
      do i = 1, nvf, 2
         k = k + 1
         cvf(i)   =  cc_vf(k)
         cvf(i+1) = -cc_vf(k)
      end do
!     norm of vf current vector
      cvf_ssq = 0
      do i=1, nvf
         cvf_ssq = cvf_ssq + cc_vf(i)**2
      end do
 
      end subroutine eval_vf_coils


      subroutine load_vf_structures
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vf_coils
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
 
!     load the coil parameters with values from variables
 
      do i=1, num_vf
         vertical(i)%current = cc_vf(i)
      end do
 
      end subroutine load_vf_structures


      subroutine read_bnorm_coefs
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use boundary
      use bnorm_mod
      use mpi_params                                         !mpi stuff
      use safe_open_mod
      implicit none
!-----------------------------------------------
!   L o c a l  P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: imnmax = 1000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iostat, ierr, ibnorm = 22
      integer :: i, im, in
      real (rprec) :: bn
!-----------------------------------------------
      call safe_open(ibnorm, ierr, bnorm_file, 'old', 'formatted')
      if (ierr .ne. 0) then
         if( myid .eq. master ) then                         !mpi stuff
           print *, 'Error opening ' // trim(bnorm_file)
     1                               // ': ierr =', iostat
         end if                                              !mpi stuff
         stop
      end if
 
!     Read and load coefficients into boundary
 
      mnbn_max = 0
      do i=1,imnmax
          read (ibnorm, *, end=120, iostat=ierr) im, in, bn
          mnbn_max = i
          xbn_m(i) = im
          xbn_n(i) = in*nfp                                  ! note nfp
          bn_coef(i) = bn
      end do
  120 continue

      close(ibnorm)
      if( myid .eq. master ) 
     1   print *, 'number of bnorm coefficients read = ', mnbn_max

      end subroutine read_bnorm_coefs


      subroutine evaluate_bnorm
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use bnorm_mod
      use tor_field
      use safe_open_mod
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ku, kv, n, nmatch, iunit=29, ierr
      real (rprec) :: pscale, theta, phi, bnrml
!-----------------------------------------------
!     Write field normal at plasma edge vs u and v for plot
      if (myid .eq. master) call safe_open(iunit,
     1    ierr, 'b_norm_eq.dat', 'unknown', 'formatted')
      pscale = dtwopi/nfp
 
      bnormal_match = zero
      n = 0
      do ku = 1, nu
          do kv = 1, nv
              n = n + 1
              phi = phib(n)
              theta = thetab(n)
              bnrml = zero
              do i=1,mnbn_max
                  bnrml = bnrml + bn_coef(i)*
     1                    sin(xbn_m(i)*theta + xbn_n(i)*phi)
              end do
!     bnrml (Tesla) = bnrml*i_pol*mu, where i_pol = Ipol/nfp is the
!     total poloidal current per period obtained from R*Bt = mu*Ipol/dtwopi
!     (R*BT = R-BTOR(s=1) is output from vmec run producing bnrml)
!     For m3.b15, i_pol = 3.0159*mu (whete R*Bt = 1.44 T-m)
              bnormal_match(n) = bnrml*i_pol*dmu0
              if( myid .eq. master ) write(iunit,1000) phi/pscale, 
     1            theta/dtwopi, bnrml*i_pol*dmu0
 1000         format(1p3e16.7)
          end do
          if( myid .eq. master ) write(iunit,1010)
 1010     format("")
      end do
      close(iunit)
      nmatch = n

      end subroutine evaluate_bnorm


      subroutine initialize_access
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use bcoils_mod
      use Vcoilpts
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n, status
      real (rprec) :: dx_access, dy_access, dz_access
!-----------------------------------------------
 
      n_access_pts = nwdim1
 
      allocate(x_access(ncdim, nwdim1), y_access(ncdim, nwdim1),
     1         z_access(ncdim, nwdim1), stat=status)
      if(status /= 0) stop "Cannot allocate access points"
 
      do n = 1, n_access
         dx_access = (x1_access(n) - x0_access(n))/(n_access_pts - 1)
         dy_access = (y1_access(n) - y0_access(n))/(n_access_pts - 1)
         dz_access = (z1_access(n) - z0_access(n))/(n_access_pts - 1)
         do i = 1, n_access_pts
            x_access (n,i) = x0_access (n) + (i-1)*dx_access
            y_access (n,i) = y0_access (n) + (i-1)*dy_access
            z_access (n,i) = z0_access (n) + (i-1)*dz_access
         end do
      end do

      end subroutine initialize_access


      subroutine read_bcoils
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use bcoils_mod
      use Vcoilpts
      use safe_open_mod
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iostat, ierr, icoils=21
      integer :: i, n, mbw
      real (rprec) :: bn
!-----------------------------------------------
      call safe_open (icoils, ierr, trim(bcoil_file), 'old','formatted')
      if (ierr .ne. 0) then
         print *, 'Error opening bcoil_file: ierr = ', ierr,
     1   ' for processor ', myid         
         stop
      end if
 
!     Read coil filament coordinates
 
      read (icoils, *) mbcoils
      if (mbcoils .gt. ncdim) then
         if( myid .eq. master ) then                         !mpi stuff
           print *, 'Error: no. of bg coils > ncdim ', mbcoils
         end if     !if( myid .eq. master )                  !mpi stuff
         stop
      end if
 
      do n=1,mbcoils
         read (icoils, *) mbwires(n)
         mbw = mbwires(n) + 1
         if (mbw .gt. nwdim1) then
            if( myid .eq. master ) then                         !mpi stuff
               print *, 'Error: no. of bg wires > nwdim ', mbw
            end if     !if( myid .eq. master )                  !mpi stuff
            stop
         end if
         do i=1,mbw
            read (icoils, *) bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i)
         end do
      end do
 
      close (icoils)
 
      end subroutine read_bcoils


      subroutine bfield (nw, xw, yw, zw, cur, xp, yp, zp, bx, by, bz)
 
      use Vcoilpts
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, nw
      real (rprec) :: cur, xp, yp, zp, bx, by, bz
      real (rprec) :: ax, ay, az, r12, fac
      real (rprec), dimension (nwdim1) :: xw, yw, zw
      real (rprec), dimension (nwdim1) :: dxp, dyp, dzp, rwp
      real (rprec), dimension (nwdim1) :: fap, dxw, dyw, dzw
      real (rprec), dimension (nwdim1) :: vxw, vyw, vzw
!-----------------------------------------------
 
      fac = 1.e-7_dp
 
      do i=1, nw
        dxp(i) = xp - xw(i)
        dyp(i) = yp - yw(i)
        dzp(i) = zp - zw(i)
        rwp(i) = sqrt (dxp(i)**2 + dyp(i)**2 + dzp(i)**2)
      end do
 
      do i=1, nw-1
        r12 = rwp(i+1)*rwp(i)
        fap(i) = (rwp(i+1) + rwp(i))/
     1    (r12*(r12 + dxp(i+1)*dxp(i)
     2              + dyp(i+1)*dyp(i)
     3              + dzp(i+1)*dzp(i)))
      end do
 
      do i=1, nw-1
        dxw(i) = (xw(i+1) - xw(i))*cur*fac
        dyw(i) = (yw(i+1) - yw(i))*cur*fac
        dzw(i) = (zw(i+1) - zw(i))*cur*fac
        vxw(i) = yw(i)*dzw(i) - zw(i)*dyw(i)
        vyw(i) = zw(i)*dxw(i) - xw(i)*dzw(i)
        vzw(i) = xw(i)*dyw(i) - yw(i)*dxw(i)
      end do
 
      ax = sum (fap(:nw-1)*dxw(:nw-1))
      ay = sum (fap(:nw-1)*dyw(:nw-1))
      az = sum (fap(:nw-1)*dzw(:nw-1))
 
      bx = sum (fap(:nw-1)*vxw(:nw-1)) - yp*az +zp*ay
      by = sum (fap(:nw-1)*vyw(:nw-1)) - zp*ax +xp*az
      bz = sum (fap(:nw-1)*vzw(:nw-1)) - xp*ay +yp*ax
 
      end subroutine bfield

#if defined(GEOM_ONLY)
#else
      subroutine evaluate_bg_field
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use modular_coils
      use saddle_coils
      use tor_field
      use vf_coils
      use bcoils_mod
      use Vcoilpts
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, ks, n, nw, mbw
      real (rprec), dimension (nwdim1) :: xw, yw, zw
      real (rprec) :: cur, xp, yp, zp, bx, by, bz
      real (rprec) :: bxp, byp, bzp
!-----------------------------------------------
      xp = 1.4_dp
      yp = 0
      zp = 0
      bx = 0
      by = 0
      bz = 0
 
!     Field due to background coils
 
      do n=1, mbcoils
         mbw = mbwires(n) + 1
         cur = bcoil_cur(n)
!        print *, 'bcoil_cur = ', bcoil_cur(n)
         do i=1, mbw
            xw(i) = bcoil_x(n,i)
            yw(i) = bcoil_y(n,i)
            zw(i) = bcoil_z(n,i)
         end do
         call bfield (mbw, xw, yw, zw, cur, xp, yp, zp,
     1                bxp, byp, bzp)
         bx = bxp + bx
         by = byp + by
         bz = bzp + bz
      end do
 
!     Field due to other coils
 
      nw = nwire1
      if (lmodular) then
!        compute field due to modulars
         do j=1, ncoils
            cur = curcon(j)
            do i=1, nw
               xw(i) = x_mod(i,1,j)
               yw(i) = y_mod(i,1,j)
               zw(i) = z_mod(i,1,j)
            end do
            call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                   bxp, byp, bzp)
            bx = bx + bxp
            by = by + byp
            bz = bz + bzp
         end do
      end if
      if (lsaddle) then
!        add field due to saddle coils
         do j=1, nsad_coils
            if (nfils .lt. 4) then
!           One or three filament model
               ks = 1
               cur = c_sad(j)/nfils
            else
!           Five filament model (no current in central filament 1)
               ks = 2
               cur = c_sad(j)/(nfils - 1)
            end if
            do k=ks, nfils
               do i=1, nw
                  xw(i) = x_sad(i,j,k)
                  yw(i) = y_sad(i,j,k)
                  zw(i) = z_sad(i,j,k)
               end do
               call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                      bxp, byp, bzp)
               bx = bx + bxp
               by = by + byp
               bz = bz + bzp
            end do
         end do
      end if
      if (lvf) then
!        add field due to vf coils
         do j=1, nvf
            cur = cvf(j)
            do i=1, nw
               xw(i) = x_vf(i,1,j)
               yw(i) = y_vf(i,1,j)
               zw(i) = z_vf(i,1,j)
            end do
            call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                   bxp, byp, bzp)
            bx = bx + bxp
            by = by + byp
            bz = bz + bzp
         end do
      end if
      if (ltfc) then
!        add field due to tf coil (1/R)
         do j=1, mtfcoil
            cur = tfc_cur(j)
            do i=1, mtfwire
               xw(i) = tfc_x(j,i)
               yw(i) = tfc_y(j,i)
               zw(i) = tfc_z(j,i)
            end do
            call bfield (mtfwire, xw, yw, zw, cur, xp, yp, zp,
     1                   bxp, byp, bzp)
            bx = bx + bxp
            by = by + byp
            bz = bz + bzp
         end do
      end if
      print 1010, xp, yp, zp
      print 1020, bx, by, bz
 1010 format('          r [m] = ',1p3e14.4)
 1020 format('       B(r) [T] = ',1p3e14.4)

      end subroutine evaluate_bg_field
#endif

      subroutine evaluate_field_error
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use modular_coils
      use saddle_coils
      use tor_field
      use vf_coils
      use bcoils_mod
      use Vcoilpts
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: iunit=23
      integer :: i, j, k, ks, n, nw, mbw
      real (rprec), dimension (nwdim1) :: xw, yw, zw
      real (rprec) :: cur, xp, yp, zp, bx, by, bz
      real (rprec) :: bxp, byp, bzp
      real (rprec) :: sinphi, cosphi, br, bphi
      real (rprec) :: b_error_rms
!-----------------------------------------------
      nw = nwire1
      rbphi_avg = 0
      do n=1, nedge
         xp = x_p(n)
         yp = y_p(n)
         zp = z_p(n)
         sinphi = sin(phib(n))
         cosphi = cos(phib(n))
         bx = 0
         by = 0
         bz = 0
         if (lmodular) then
!        compute field due to modulars
            do j=1, ncoils
               cur = curcon(j)
               do i=1, nw
                  xw(i) = x_mod(i,1,j)
                  yw(i) = y_mod(i,1,j)
                  zw(i) = z_mod(i,1,j)
               end do
               call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                      bxp, byp, bzp)
               bx = bx + bxp
               by = by + byp
               bz = bz + bzp
            end do
         end if
         if (lsaddle) then
!        add field due to saddle coils
            do j=1, nsad_coils
               if (nfils .lt. 4) then
!              One - three filament model
                  ks = 1
                  cur = c_sad(j)/nfils
               else
!              Five filament model (no current in central filament 1)
                  ks = 2
                  cur = c_sad(j)/(nfils - 1)
               end if
               do k=ks, nfils
                  do i=1, nw
                     xw(i) = x_sad(i,j,k)
                     yw(i) = y_sad(i,j,k)
                     zw(i) = z_sad(i,j,k)
                  end do
                  call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                         bxp, byp, bzp)
                  bx = bx + bxp
                  by = by + byp
                  bz = bz + bzp
               end do
            end do
         end if
         if (lvf) then
!        add field due to vf coils
            do j=1, nvf
               cur = cvf(j)
               do i=1, nw
                  xw(i) = x_vf(i,1,j)
                  yw(i) = y_vf(i,1,j)
                  zw(i) = z_vf(i,1,j)
               end do
               call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1                      bxp, byp, bzp)
               bx = bx + bxp
               by = by + byp
               bz = bz + bzp
            end do
         end if
         if (ltfc) then
!        add field due to tf coil (1/R)
            do j=1, mtfcoil
               cur = tfc_cur(j)
               do i=1, mtfwire
                  xw(i) = tfc_x(j,i)
                  yw(i) = tfc_y(j,i)
                  zw(i) = tfc_z(j,i)
               end do
               call bfield (mtfwire, xw, yw, zw, cur, xp, yp, zp,
     1                      bxp, byp, bzp)
               bx = bx + bxp
               by = by + byp
               bz = bz + bzp
            end do
         end if
         if (lbcoil) then
!        add field due to background coils
            do j=1, mbcoils
               mbw = mbwires(j) + 1
               cur = bcoil_cur(j)
               do i=1, mbw
                  xw(i) = bcoil_x(j,i)
                  yw(i) = bcoil_y(j,i)
                  zw(i) = bcoil_z(j,i)
               end do
               call bfield (mbw, xw, yw, zw, cur, xp, yp, zp,
     1                      bxp, byp, bzp)
               bx = bx + bxp
               by = by + byp
               bz = bz + bzp
            end do
         end if
         br =    cosphi*bx + sinphi*by
         bphi = -sinphi*bx + cosphi*by
         b_mod(n) = sqrt (br**2 + bphi**2 + bz**2)
         b_error(n) = (br*n_r(n) + bphi*n_phi(n) + bz*n_z(n)
     1              +  bnormal_match(n)) / b_mod(n)
         rbphi_avg = rbphi_avg + d_area(n)*rb(n)*bphi
      end do
      b_error_rms = sqrt(sum(b_error(1:nedge)**2)/nedge)
      rbphi_avg = rbphi_avg/sum_d_area
 
      end subroutine evaluate_field_error


      subroutine read_wout(extension)
      use read_wout_mod, nfp_w => nfp, ns_w => ns, mpol_w => mpol,
     1   ntor_w => ntor, mnmax_w => mnmax, xm_w => xm, xn_w => xn
      use boundary
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character*(*) :: extension
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat
!-----------------------------------------------
      call read_wout_file('wout.' // trim(extension), istat)      
      if (istat .ne. 0) then
         if( myid .eq. master ) then                         !mpi stuff
           print *,' Error reading wout file: ierr = ', istat
         end if     !if( myid .eq. master )                  !mpi stuff
         stop
      end if   
  
!     Load values from read_wout_mod module into boundary
 
      nfp = nfp_w
      ns  = ns_w
      mpol= mpol_w
      ntor= ntor_w
      mnmax = mnmax_w
 
!     ALLOCATE ARRAYS FIRST TIME THRU
 
      allocate (ixm(mnmax), ixn(mnmax), rmnc_b(mnmax), zmns_b(mnmax),
     1  xm(mnmax), xn(mnmax), gmn_b(mnmax), lmns_b(mnmax),
     2  stat = istat)
        if( myid .eq. master ) then                         !mpi stuff
          if (istat .ne. 0) stop 'allocation error in read_wout'
        end if     !if( myid .eq. master )                  !mpi stuff
 
!     ONLY SAVE BOUNDARY COEFFICIENTS
 
      xm = zero
      xn = zero
      ixm = zero
      ixn = zero
      rmnc_b = zero
      zmns_b = zero
      lmns_b = zero
      gmn_b = zero
 
      xm (:) = xm_w(:)
      xn (:) = xn_w(:)
      ixm(:) = nint(xm(:))
      ixn(:) = nint(xn(:))
      rmnc_b(:) = rmnc(:,ns)
      zmns_b(:) = zmns(:,ns)
      lmns_b(:) = 1.5_dp*lmns(:,ns) - 0.5_dp*lmns(:,ns-1)
      gmn_b(:) = 1.5_dp*gmn(:,ns) - 0.5_dp*gmn(:,ns-1)
      iota_b = 1.5_dp*iotas(ns) - 0.5_dp*iotas(ns-1)
      phip_b = 1.5_dp*phip(ns) - 0.5_dp*phip(ns-1)
 
      if( myid .eq. master ) then                         !mpi stuff
        print 100, iota_b, phip_b
      end if     !if( myid .eq. master )                  !mpi stuff
  100 format (" iota = ",f7.5,", phip = ",f7.5)
  
      nu = max(2*mpol + 6, 32)
      nv = max(2*ntor + 4, 32)
       
!      nu = 36
!      nv = 36
       
      nuv = nu * nv
 
      if( myid .eq. master ) then                         !mpi stuff
        print '(a,i6,a,i6)', 'nu =',nu,'  nv =',nv
      end if     !if( myid .eq. master )                  !mpi stuff

!     Deallocate memory
 
 1000 call read_wout_deallocate
  
      end subroutine read_wout


      subroutine loadparams (xc, nvar, mfun)
      use modular_coils
      use saddle_coils
      use bcoils_mod
      use saddle_surface
      use tor_field
      use vf_coils
      use Vcoilpts
      use coils
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: nvar, mfun
      real(rprec), target :: xc(nvar)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, nvariables
      real(rprec), pointer :: xvariables(:)
!-----------------------------------------------
 
!     Load modular coefficients
 
      n = 0
      if (lmodular) then
         xvariables => xc(n+1:)
         call load_modular_coils (nvariables, xvariables)
         n = n + nvariables
 
         
         if (lmodcur) xvariables => xc(n+1:)
!        Need to call load_modular_currents even if lmodcur = F, in case
!        TF current has changed
         call load_modular_currents (nvariables, xvariables)
         n = n + nvariables
 
         call load_modular_structures
      end if
 
!     Load saddle coil coefficients (if lsaddle = true)
 
      if (lsaddle) then
         xvariables => xc(n+1:)
         call load_saddle_coils (nvariables, xvariables)
         n = n + nvariables
 
!        For now, can not vary TF current when lsaddle = T, unless
!        saddle currents also vary (need to fix)
         if (lsadcur) then
            xvariables => xc(n+1:)
            call load_saddle_currents (nvariables, xvariables)
            n = n + nvariables
         end if
 
         call load_saddle_structures
      end if
 
!     Load background coil currents (if lbcoil = true)
 
      if (lbcoil) then
         xvariables => xc(n+1:)
         call load_bg_currents (nvariables, xvariables)
         n = n + nvariables
      end if
 
!     Load vf coil currents (if lvf = true)
 
      if (lvf) then
         if (lvfvar) then
            xvariables => xc(n+1:)
            call load_vf_coils (nvariables, xvariables)
            n = n + nvariables
         end if

         xvariables => xc(n+1:)
         call load_vf_currents (nvariables, xvariables)
         n = n + nvariables
!
         call load_vf_structures
      end if
 
!     Load tf coil currents (if ltfc = true)
 
      if (ltfcv) then
         xvariables => xc(n+1:)
         call load_tf_coils (nvariables, xvariables)
         n = n + nvariables
      end if
 
!     Load modular winding surface coefficients (if lsurfv = true)
 
      if (lsurfv) then
         xvariables => xc(n+1:)
         call load_modular_wsurf (nvariables, xvariables)
         n = n + nvariables
      end if
 
!     Load saddle winding surface coefficients (if lsadsfv = true)
 
      if (lsadsfv) then
         xvariables => xc(n+1:)
         call load_saddle_wsurf (nvariables, xvariables)
         n = n + nvariables
      end if
 
      if( myid .eq. master ) then                            !mpi stuff
        if (n .ne. nvar) stop 'loadparams: n .ne. nvar'
      end if     !if( myid .eq. master )                     !mpi stuff
 
      end subroutine loadparams


      subroutine save_for05(for05file)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use coilsnamin
      use safe_open_mod
      use gade_mod
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character for05file*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iunit=15, ierr
!-----------------------------------------------
 
      call safe_open(iunit, ierr, for05file, 'replace', 'formatted')
      call write_coilsin (iunit, ierr)
      if (nopt_alg .gt. 0) call write_gade_nml(iunit)
      close(iunit)
  
!     print *,'minimum chi-sq state saved in file ',
!    1    trim(for05file)

      end subroutine save_for05


      subroutine plasma_mod_coil_distance
      use boundary
      use modular_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat, i, n, j
      real(rprec) :: xpl, ypl, zpl, xcl, ycl, zcl, d_cp, cpd
!-----------------------------------------------
 
      p_d_min=1000
 
      do i=1,ncoils
        do j = 1,nwire
           xcl=x_mod(j,1,i)
           ycl=y_mod(j,1,i)
           zcl=z_mod(j,1,i)
           cpd = 1000
           do n = 1, nedge
              xpl=rb(n)*cos(phib(n))
              ypl=rb(n)*sin(phib(n))
              zpl=zb(n)
              d_cp=sqrt((xcl-xpl)**2+(ycl-ypl)**2+(zcl-zpl)**2)
              cpd = min(cpd, d_cp)
              p_d_min=min(p_d_min,d_cp)
           end do
         end do
      end do
 
      end subroutine plasma_mod_coil_distance


      subroutine plasma_sad_coil_distance
      use boundary
      use saddle_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat, i, n, j
      real(rprec) :: xpl, ypl, zpl, xcl, ycl, zcl, d_cp, cpd
!-----------------------------------------------
 
      p_s_min=1000
 
      do i=1,nsad_coils
        do j = 1,nwire
           xcl=x_sad(j,i,1)
           ycl=y_sad(j,i,1)
           zcl=z_sad(j,i,1)
           cpd = 1000
           do n = 1, nedge
              xpl=rb(n)*cos(phib(n))
              ypl=rb(n)*sin(phib(n))
              zpl=zb(n)
              d_cp=sqrt((xcl-xpl)**2+(ycl-ypl)**2+(zcl-zpl)**2)
              cpd = min(cpd, d_cp)
              p_s_min=min(p_s_min,d_cp)
           end do
         end do
      end do
 
      end subroutine plasma_sad_coil_distance


      subroutine spline_bkp_distance
      use saddle_coils
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n
      real(rprec), parameter :: vs1 = 1, us1 = 1
      real(rprec) :: bkp_tst
!-----------------------------------------------
 
      bkp_min=1000

      do i = 1, nsmid

         if (nsad_v .gt. 5) then
            do n = 4, nsad_v - 1 
               bkp_tst = sad_v_s(i,n+1) - sad_v_s(i,n)
               bkp_min = min(bkp_min, bkp_tst)
            end do
            bkp_tst = vs1 - sad_v_s(i,nsad_v)
            bkp_min = min(bkp_min, bkp_tst)
         end if

         if (nsad_u .gt. 5) then
            do n = 4, nsad_u - 1
               bkp_tst = sad_u_s(i,n+1) - sad_u_s(i,n)
               bkp_min = min(bkp_min, bkp_tst)
            end do
            bkp_tst = us1 - sad_u_s(i,nsad_u)
            bkp_min = min(bkp_min, bkp_tst)
         end if

      end do

      end subroutine spline_bkp_distance


      subroutine evaluate_access
      use boundary, only: nfp
      use bcoils_mod
      use modular_coils
      use saddle_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat, i, j, k, n, nc, ns
      real(rprec) :: x_c0, y_c0, z_c0, x_c1, y_c1, z_c1
      real(rprec) :: d_cc, ci_min
!-----------------------------------------------
 
      nc = nmod_coils_per_period*nfp
      ns = nsad_coils_per_period*nfp
      do n=1,n_access
         ci_min = 1000
         do i = 1,n_access_pts
            x_c0=x_access(n,i)
            y_c0=y_access(n,i)
            z_c0=z_access(n,i)
!           Distance to all modular coils
            if (lmodular) then
               do j = 1,nc
                  do k = 1,nwire
                     x_c1=x_mod(k,1,j)
                     y_c1=y_mod(k,1,j)
                     z_c1=z_mod(k,1,j)
                     d_cc = (x_c0-x_c1)**2
     1                    + (y_c0-y_c1)**2
     2                    + (z_c0-z_c1)**2
                     ci_min = min(ci_min, d_cc)
                  end do
               end do
            end if
!           Distance to all saddle coils
            if (lsaddle) then
               do j = 1,ns
                  do k = 1,nwire
                     x_c1=x_sad(k,j,1)
                     y_c1=y_sad(k,j,1)
                     z_c1=z_sad(k,j,1)
                     d_cc = (x_c0-x_c1)**2
     1                    + (y_c0-y_c1)**2
     2                    + (z_c0-z_c1)**2
                     ci_min = min(ci_min, d_cc)
                  end do
               end do
            end if
         end do
         acc_min(n) = sqrt(ci_min)
      end do
 
      end subroutine evaluate_access


      subroutine mod_coil_distance
      use boundary, only: nfp
      use modular_coils
      use saddle_coils
      use vf_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat, i, j, k, n, nc, ns
      real(rprec) :: x_c0, y_c0, z_c0, x_c1, y_c1, z_c1
      real(rprec) :: d_cc, ci_min
!-----------------------------------------------
 
      nc = nmod_coils_per_period*nfp
      ns = nsad_coils_per_period*nfp
      do i=1,nmid
        ci_min = 1000
        do n = 1,nc
!       Distance to other modular coils
           if (i .ne. n) then
              do j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 do k = 1,nwire
                    x_c1=x_mod(k,1,n)
                    y_c1=y_mod(k,1,n)
                    z_c1=z_mod(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = min(ci_min, d_cc)
                 end do
              end do
           end if
        end do
!       Distance to all saddle coils
        if (lsaddle) then
           do n = 1,ns
              do j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 do k = 1,nwire
                    x_c1=x_sad(k,n,1)
                    y_c1=y_sad(k,n,1)
                    z_c1=z_sad(k,n,1)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = min(ci_min, d_cc)
                 end do
              end do
           end do
        end if
!       Distance to all vf coils
        if (lvf) then
           do n = 1,nvf
              do j = 1,nwire
                 x_c0=x_mod(j,1,i)
                 y_c0=y_mod(j,1,i)
                 z_c0=z_mod(j,1,i)
                 do k = 1,nwire
                    x_c1=x_vf(k,1,n)
                    y_c1=y_vf(k,1,n)
                    z_c1=z_vf(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = min(ci_min, d_cc)
                 end do
              end do
           end do
        end if
        cc_min(i) = sqrt(ci_min)
      end do
 
      end subroutine mod_coil_distance


      subroutine sad_coil_distance
      use cparms
      use boundary
      use saddle_coils
      use modular_coils
      use vf_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istat, i, j, k, l, m, n, nmc, nsc
      integer, dimension(ncdim) :: nsv, jsv, ksv
      real(rprec) :: x_c0, y_c0, z_c0, x_c1, y_c1, z_c1
      real(rprec) :: xpl, ypl, zpl, xpls, ypls, zpls, cpd
      real(rprec) :: dphib, d_cc, ci_min, d_cp, d_cxp
!-----------------------------------------------
 
      dphib = dtwopi/nfp
      nsc = nsad_coils_per_period*nfp
      nmc = nmod_coils_per_period*nfp
!     For each unique coil i ...
      do i=1,nsmid
        ci_min = 1000
        do n = 1,nsc
!          Find the min distance to all other coils n .ne. i
           if (i .ne. n) then
              do j = 1,nwire
                 do l = 1, nfils
                    x_c0 = x_sad(j,i,1) + 2*(x_sad(j,i,l)-x_sad(j,i,1))
                    y_c0 = y_sad(j,i,1) + 2*(y_sad(j,i,l)-y_sad(j,i,1))
                    z_c0 = z_sad(j,i,1) + 2*(z_sad(j,i,l)-z_sad(j,i,1))
!                Find the nearest point k, coil n, to point j on coil i
                    do k = 1,nwire
                       do m = 1, nfils
                          x_c1 = x_sad(k,n,1)
     1                         + 2*(x_sad(k,n,m)-x_sad(k,n,1))
                          y_c1 = y_sad(k,n,1)
     1                         + 2*(y_sad(k,n,m)-y_sad(k,n,1))
                          z_c1 = z_sad(k,n,1)
     1                         + 2*(z_sad(k,n,m)-z_sad(k,n,1))
                          d_cc = (x_c0-x_c1)**2
     1                         + (y_c0-y_c1)**2
     2                         + (z_c0-z_c1)**2
                          if (d_cc .lt. ci_min) then
                             ci_min = d_cc
!                            Save indices of points j, k and coil n
                             ns = n
                             js = j
                             ks = k
                          end if
                       end do
                    end do
                 end do
              end do
           end if
        end do
!       Distance to all modular coils
        if (lmodular) then
           do n = 1,nmc
              do j = 1,nwire
                 x_c0=x_sad(j,i,1)
                 y_c0=y_sad(j,i,1)
                 z_c0=z_sad(j,i,1)
                 do k = 1,nwire
                    x_c1=x_mod(k,1,n)
                    y_c1=y_mod(k,1,n)
                    z_c1=z_mod(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = min(ci_min, d_cc)
                 end do
              end do
           end do
        end if
!       Distance to all vf coils
        if (lvf) then
           do n = 1,nvf
              do j = 1,nwire
                 x_c0=x_sad(j,i,1)
                 y_c0=y_sad(j,i,1)
                 z_c0=z_sad(j,i,1)
                 do k = 1,nwire
                    x_c1=x_vf(k,1,n)
                    y_c1=y_vf(k,1,n)
                    z_c1=z_vf(k,1,n)
                    d_cc = (x_c0-x_c1)**2
     1                   + (y_c0-y_c1)**2
     2                   + (z_c0-z_c1)**2
                    ci_min = min(ci_min, d_cc)
                 end do
              end do
           end do
        end if
        sc_min(i) = sqrt(ci_min)
        nsv(i) = ns
        jsv(i) = js
        ksv(i) = ks
      end do

      if (nedge .gt. 0) then
!     Find size of min coil-coil vector crossed with min coil-plasma vector
         do i=1,nsmid
            ns = nsv(i)
            js = jsv(i)
            ks = ksv(i)
            x_c0=x_sad(js,i,1)
            y_c0=y_sad(js,i,1)
            z_c0=z_sad(js,i,1)
            x_c1=x_sad(ks,ns,1)
            y_c1=y_sad(ks,ns,1)
            z_c1=z_sad(ks,ns,1)
            cpd = 1000
!           Min distance to from coil i, point js to plasma
            do m = 1, nedge
               zpl=zb(m)
               do k = 1,nfp
                  xpl=rb(m)*cos(phib(m) + (k-1)*dphib)
                  ypl=rb(m)*sin(phib(m) + (k-1)*dphib)
                  d_cp = (x_c0-xpl)**2
     1                 + (y_c0-ypl)**2
     2                 + (z_c0-zpl)**2
                  if (d_cp .lt. cpd) then
                     cpd = d_cp
                     xpls = xpl
                     ypls = ypl
                     zpls = zpl
                  end if
               end do
            end do
            d_cxp = ((y_c0 - y_c1)*(zpls - z_c0)
     1            -  (ypls - y_c0)*(z_c0 - z_c1))**2
     2            + ((z_c0 - z_c1)*(xpls - x_c0)
     3            -  (zpls - z_c0)*(x_c0 - x_c1))**2
     4            + ((x_c0 - x_c1)*(ypls - y_c0)
     5            -  (xpls - x_c0)*(y_c0 - y_c1))**2
            scxp_min(i) = sqrt(d_cxp/((xpls - x_c0)**2
     1                              + (ypls - y_c0)**2
     2                              + (zpls - z_c0)**2))
         end do
      end if
 
      end subroutine sad_coil_distance


      subroutine coil_curvature
      use boundary
      use modular_coils
      use Vwire
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, ncpts
      real(rprec) :: crv_max, crv, csum
      real(rprec) :: u0, v0, du
      real(rprec) :: tx, ty, tz
!-----------------------------------------------
 
      ncpts = 4*nwire
      du = 1.0_dp/ncpts
      do i=1,nmid
        crv_max = 0
        csum = 0
        do j = 1,ncpts
          u0 = (j-1)*du
          call modular_curve (i, u0, v0, tx, ty, tz, crv)
          crv_max = max(crv_max, crv)
          csum = csum + crv
        end do
        rc_min(i) = 10000.0_dp
        if (crv_max .gt. 0.0_dp) rc_min(i) = 1.0_dp/crv_max
        cu_sum(i) = csum*du
      end do
 
      end subroutine coil_curvature


      subroutine modular_curve (n, u, v, tx, ty, tz, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use modular_coils
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, ncl
      real(rprec) :: tx, ty, tz
      real(rprec) :: ck, sk, u, v,
     1   ph0, ph1, ph2, r0, r1, r2, z0, z1, z2,
     2   x0, x1, x2, y0, y1, y2, s1, s2, c
      real(rprec) :: th0, th1, th2
      real(rprec) :: fth0, fth1, fth2
!-----------------------------------------------
 
!     compute theta and derivative wrt u
 
      ncl = n
      th0 = dtwopi*u
      th1 = dtwopi
 
!     compute phi and derivatives wrt u
 
      if (nf_rho .gt. 0) then
!     case with possible windbacks - call theta_f to modify arguement th0
!     put theta and derivatives in temporary variables fth0, fth1, fth2 so
!     the winding law calculation will not interfere with the winding 
!     surface evaluation below
         fth0 = th0
         call theta_f (ncl, fth0, fth1, fth2)
         ph0 = phi_full(ncl)
         ph1 = 0
         ph2 = 0
         do k = 0,nf_phi
            ck = cos(k*fth0)
            sk = sin(k*fth0)
            ph0 = ph0 +  modular(ncl)%phic(k)*ck
     1                +  modular(ncl)%phis(k)*sk
            ph1 = ph1 + (modular(ncl)%phis(k)*ck
     1                -  modular(ncl)%phic(k)*sk)*k
            ph2 = ph2 + (modular(ncl)%phis(k)*sk
     1                +  modular(ncl)%phic(k)*ck)*k**2
         end do
!     note the order of the next two lines is important
         ph2 = dtwopi**2*(ph2*fth1**2 + ph1*fth2)
         ph1 = dtwopi*fth1*ph1
      else
!     case without windbacks
         ph0 = phi_full(ncl)
         ph1 = 0
         ph2 = 0
         do k = 0,nf_phi
            ck = cos(k*th0)
            sk = sin(k*th0)
            ph0 = ph0 +  modular(ncl)%phic(k)*ck
     1                +  modular(ncl)%phis(k)*sk
            ph1 = ph1 + (modular(ncl)%phis(k)*ck
     1                -  modular(ncl)%phic(k)*sk)*k
            ph2 = ph2 + (modular(ncl)%phis(k)*sk
     1                +  modular(ncl)%phic(k)*ck)*k**2
         end do
         ph1 = th1*ph1
         ph2 = -th1**2*ph2
      end if
 
!     compute R, Z and derivatives wrt u
 
      r0 = 0
      r1 = 0
      r2 = 0
      z0 = 0
      z1 = 0
      z2 = 0
      do k = 1, numsurf
         ck = cos(m_num(k)*th0 + nfper*n_num(k)*ph0)
         sk = sin(m_num(k)*th0 + nfper*n_num(k)*ph0)
!     R ...
         r0 = r0 + rmn_sf(k)*ck
         r1 = r1 - rmn_sf(k)*(m_num(k)*th1 + nfper*n_num(k)*ph1)*sk
         r2 = r2 - rmn_sf(k)*(nfper*n_num(k)*ph2*sk
     1      + (m_num(k)*th1 + nfper*n_num(k)*ph1)**2*ck)
!     Z ...
         z0 = z0 + zmn_sf(k)*sk
         z1 = z1 + zmn_sf(k)*(m_num(k)*th1 + nfper*n_num(k)*ph1)*ck
         z2 = z2 + zmn_sf(k)*(nfper*n_num(k)*ph2*ck
     1      - (m_num(k)*th1 + nfper*n_num(k)*ph1)**2*sk)
      end do
 
!     compute X and derivatives wrt u
 
      x0 = r0*cos(ph0)
      x1 = r1*cos(ph0) - r0*ph1*sin(ph0)
      x2 = r2*cos(ph0) - r1*ph1*sin(ph0)
     1   - (r0*ph1**2*cos(ph0) + (r0*ph2 + r1*ph1)*sin(ph0))
 
!     compute Y and derivatives wrt u
 
      y0 = r0*sin(ph0)
      y1 = r1*sin(ph0) + r0*ph1*cos(ph0)
      y2 = r2*sin(ph0) + r1*ph1*cos(ph0)
     1   + (-r0*ph1**2*sin(ph0) + (r0*ph2 + r1*ph1)*cos(ph0))
 
!     compute ds/du
 
      s1 = dsqrt (x1*x1 + y1*y1 + z1*z1)
 
!     tangent vector dx/ds, dy/ds, dz/ds
 
      tx = x1/s1
      ty = y1/s1
      tz = z1/s1
 
!     compute curvature
 
      c = dsqrt ((y1*z2 - y2*z1)**2
     1         + (z1*x2 - z2*x1)**2
     2         + (x1*y2 - x2*y1)**2)/s1**3
 
!     return toroidal v
      v = nfper*ph0/dtwopi
 
      end subroutine modular_curve


      subroutine modular_surface (u, v, x, y, z, nx, ny, nz)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use modular_coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k
      real(rprec) :: u, v
      real(rprec) :: x, y
      real(rprec) :: nx, ny, nz, nmag
      real(rprec) :: th0, th1, ph0, ph1
      real(rprec) :: r, ru, rv
      real(rprec) :: z, zu, zv
      real(rprec) :: ck, sk
!-----------------------------------------------
 
!     theta and derivative wrt u
 
      th0 = dtwopi*u
      th1 = dtwopi
 
!     phi and derivative wrt v
 
      ph0 = dtwopi*v/nfper
      ph1 = dtwopi/nfper
 
!     compute R, Z and derivatives wrt u and v
 
      r   = 0
      ru  = 0
      rv  = 0
      z   = 0
      zu  = 0
      zv  = 0
      do k = 1, numsurf
         ck = cos(m_num(k)*th0 + nfper*n_num(k)*ph0)
         sk = sin(m_num(k)*th0 + nfper*n_num(k)*ph0)
!     R ...
         r  = r  + rmn_sf(k)*ck
         ru = ru - rmn_sf(k)*m_num(k)*sk
         rv = rv - rmn_sf(k)*n_num(k)*sk
!     Z ...
         z  = z  + zmn_sf(k)*sk
         zu = zu + zmn_sf(k)*m_num(k)*ck
         zv = zv + zmn_sf(k)*n_num(k)*ck
      end do
      ru = ru*th1
      zu = zu*th1
      rv = rv*ph1*nfper
      zv = zv*ph1*nfper
 
      sk = sin(ph0)
      ck = cos(ph0)
 
!     return X and Y
 
      x = r*ck
      y = r*sk
 
!     compute Nx, Ny, Nz
 
      nx = ru*zv*sk - zu*(rv*sk + r*ck)
      ny = zu*(rv*ck - r*sk) - ru*zv*ck
      nz = r*ru
      nmag = dsqrt (nx**2 + ny**2 + nz**2)
      nx = -nx/nmag
      ny = -ny/nmag
      nz = -nz/nmag
 
      end subroutine modular_surface


      subroutine evaluate_max_field (iunit)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use bcoils_mod
      use modular_coils
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iunit, i, j, k, n, nw, ncpts, icl
      real (rprec), dimension (nwdim1) :: xw, yw, zw
      real (rprec) :: cur, xp, yp, zp, bx, by, bz
      real (rprec) :: bxp, byp, bzp, bpmax, bp
      real (rprec) :: u0, v0, du, delta
      real (rprec) :: tx, ty, tz, crv
      real (rprec) :: nx, ny, nz, x, y, z
!-----------------------------------------------
      nw = nwire1                         ! PROBLEM ?
      ncpts = 2*nwire
      delta = 0.06_dp
      du = 1.0_dp/ncpts
 
! loop over unique modular coils
 
      do i=1, nmid
         icl = i
         bpmax = 0
         do j=1, ncpts
            u0 = (j-1)*du
            call modular_curve (icl, u0, v0, tx, ty, tz, crv)
            call modular_surface (u0, v0, x, y, z, nx, ny, nz)
 
! determine xp, yp, zp
 
            xp = x - delta*nx
            yp = y - delta*ny
            zp = z - delta*nz
 
! compute field at xp, yp, zp
 
            bxp = 0
            byp = 0
            bzp = 0
            do n = 1, ncoils
               cur = curcon(n)
               do k = 1, nw
                  xw(k) = x_mod(k,1,n)
                  yw(k) = y_mod(k,1,n)
                  zw(k) = z_mod(k,1,n)
               end do                            ! over k
               call bfield (nw, xw, yw, zw, cur, xp, yp, zp,
     1            bx, by, bz)
               bxp = bxp + bx
               byp = byp + by
               bzp = bzp + bz
            end do                               ! over n
            bp = sqrt (bxp**2 + byp**2 + bzp**2)
            bpmax = max (bpmax, bp)
            write (iunit, 1010) x, y, z, xp, yp, zp, bp
 
! add field due to vf''s
 
         end do                                  ! over j
         write (iunit, 1000)
         b_max(i) = bpmax
      end do                                     ! over i
 1000 format ("")
 1010 format (7f11.4)
      end subroutine evaluate_max_field


      subroutine write_coils (extension)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use modular_coils
      use saddle_coils
      use tor_field
      use vf_coils
      use coiltypes
      use coils
      use bcoils_mod
      use Vwire
      use safe_open_mod
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character*(*) :: extension
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: cunit=15, munit, gunit, sunit, ierr
      integer :: n, i, j, js, nc, ncm, ncp, mbw, nmax, ncoil_count
      integer :: nz
      real(rprec) :: ph_w, th_w, zz, sfils
      character*7 :: status      
      character :: ch_sad*10
!-----------------------------------------------
#if defined(GEOM_ONLY)
      status = 'scratch'
#else      
      status = 'replace'
#endif

      call safe_open(cunit, ierr, 'coilxyz.dat', status, 'formatted')
      if (ierr .ne. 0) stop 'Error opening coilxyz.dat in write_coils'

      munit = cunit+1
      call safe_open(munit, ierr, 'coils.' // trim(extension), 
     1    'replace', 'formatted')
      if (ierr .ne. 0) stop 'Error opening coils file in write_coils'
     
      gunit = munit+1
      call safe_open(gunit, ierr, 'coilgnu.dat', status, 'formatted')
      if (ierr .ne. 0) stop 'Error opening coilgnu.dat in write_coils'
      
      sunit = gunit+1
      call safe_open(sunit, ierr, 'coilsad.dat', status, 'formatted')
      if (ierr .ne. 0) stop 'Error opening coilsad.dat in write_coils'
 
      zz = 0
      nz = 0

      write (munit,300) nfp
      write (munit,310)
      write (munit,320)
  300 format("periods ",i2)
  310 format("begin filament")
  320 format("mirror NUL")
 
      nmax = 0
      ncoil_count = nmax
      ncp = 0
 
      if (lmodular) then
!     Write modular coils, currents
        do n = 1, ncoils
           write(cunit,'(1pe22.14,2i6)') curcon(n), nwire, nz
           write(gunit,'(/)')
           ncm = 1 + mod(n-1,nmod_coils_per_period)
           ncp = ncp + 1
           if (ncp .eq. nmid) ncp = 0
!          compute coil group number
           if (lmodcur) then
!             only good for lsymm = t, nodd = 1 for now
              nc = ncoil_count + ncm
              if (ncm .ge. nmid) nc = ncoil_count + nmid - ncp
              if (ncm .eq. nmod_coils_per_period) ncp = 0
              nmax = max(nmax, nc)
           else
              nc = ncoil_count + 1
              nmax = nc
           end if
 
           do i = 1,nwire1
             write(cunit,'(1p3e22.14)')
     1          x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n)
             if (i .eq. nwire1) then
                write(munit,'(1p4e22.14,i4,a9)')
     1             x_mod(1,1,n), y_mod(1,1,n), z_mod(1,1,n), zero, nc,
     2             "  Modular"
             else
                write(munit,'(1p4e22.14)')
     1             x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n), curcon(n)
             end if
             ph_w = modular(nc)%phi(i) + phi_full(n)
             th_w = dtwopi*(i-1)/nwire
             write(gunit,'(3f14.6,i6,2e14.6)')
     1         x_mod(i,1,n), y_mod(i,1,n), z_mod(i,1,n),i,
     2         ph_w/(dtwopi/nfper), th_w/dtwopi
           end do
        end do
      end if
 
      if (lsaddle) then
        ncoil_count = nmax
        ncp = 0
        if (lsmod) then
           ch_sad = ' Modular'
        else
           ch_sad = ' Saddle'
        end if
 
        js = 1
        if (nfils .gt. 4) then
          js = 2
          sfils = 4.0_dp
        else if (nfils .gt. 1) then
          sfils = 3.0_dp
        else
          sfils = 1.0_dp
        end if

!       Write saddle coils, currents
        do n = 1, nsad_coils
           ncm = 1 + mod(n-1,nsad_coils_per_period)
           ncp = 1 + mod(n-1,nsmid)
!          compute coil group number
           nc = ncoil_count + nsad_group(ncp)
           if (ncm .gt. nsmid) nc = ncoil_count
     1        + nsad_group(nsmid + 1 - ncp)
           nmax = max(nmax, nc)
 
           do j = js, nfils
              write(cunit,'(1pe22.14,2i6)') c_sad(n)/sfils, nwire, nz
              write(sunit,'(/)')
              do i = 1,nwire1
                 write(cunit,'(1p3e22.14)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j)
                 if (i .eq. nwire1) then
                    write(munit, '(1p4e22.14,i4,a8)')
     1              x_sad(1,n,j), y_sad(1,n,j), z_sad(1,n,j), zz, nc,
     2              ch_sad
                 else
                    write(munit,'(1p4e22.14)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j),
     2              c_sad(n)/sfils
                 end if
                 write(sunit,'(3f14.6,2e14.6)')
     1              x_sad(i,n,j), y_sad(i,n,j), z_sad(i,n,j),
     2              v_sad(i,n), u_sad(i,n)
              end do
           end do
        end do
      end if
 
      if (lvf) then
      ncoil_count = nmax
 
!     Write vf coils, currents
        do n = 1, nvf
         write(gunit,'(/)')
!        compute coil group number
         nc = ncoil_count + (n+1)/2
         nmax = max(nmax, nc)
 
         write(cunit,'(1pe22.14,2i6)') cvf(n), nwire, nz
         do i = 1,nwire1
           write(cunit,'(1p3e22.14)')
     1        x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n)
           if (i .eq. nwire1) then
              write(munit, '(1p4e22.14,i4,a4)')
     1          x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), zz, nc, " VF"
           else
              write(munit,'(1p4e22.14)')
     1          x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), cvf(n)
           end if
           write(gunit,'(3f14.6,i6,2e14.6)')
     1      x_vf(i,1,n), y_vf(i,1,n), z_vf(i,1,n), i, zz, zz
          end do
        end do
      end if
 
      if (ltfc) then
        ncoil_count = nmax
 
!       compute coil group number
        nc = ncoil_count + 1
        nmax = max(nmax, nc)
 
!       Write tf coils, currents
        do n = 1, mtfcoil
          write(gunit,'(/)')
          write(cunit,'(1pe22.14,2i6)') tfc_cur(n), mtfwire-1, nz
          do i = 1,mtfwire
            write(cunit,'(1p3e22.14)')
     1        tfc_x(n,i), tfc_y(n,i), tfc_z(n,i)
            if (i .eq. mtfwire) then
              write(munit, '(1p4e22.14,i4,a4)')
     1          tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), zz, nc, " TF"
            else
              write(munit,'(1p4e22.14)')
     1          tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), tfc_cur(n)
            end if
            write(gunit,'(3f14.6,i6,2e14.6)')
     1       tfc_x(n,i), tfc_y(n,i), tfc_z(n,i), i, zz, zz
          end do
        end do
      end if
 
      if (lbcoil) then
      ncoil_count = nmax
      
!     compute coil group number
      nc = ncoil_count + 1
      nmax = max(nmax, nc)
!     Write background coils, currents
        do n = 1, mbcoils
         mbw = mbwires(n) + 1
         write(cunit,'(1pe22.14,2i6)') bcoil_cur(n), mbwires(n),
     1      mc_bg(n)
         write(gunit,'(/)')
         do i = 1,mbw
           write(cunit,'(1p3e22.14)')
     1        bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i)
           if (i .eq. mbw) then
              write(munit, '(1p4e22.14,i4,a4)')
     1          bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), zz,
     2          nc + mc_bg(n), " BG"
           else
              write(munit,'(1p4e22.14)')
     1          bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), bcoil_cur(n)
           end if
           write(gunit,'(3f14.6,i6,2e14.6)')
     1       bcoil_x(n,i), bcoil_y(n,i), bcoil_z(n,i), i, zz, zz
         end do
        end do
      end if
 
      if (laccess) then
!     Write access zones to coils file
        do n = 1, n_access
         write(gunit,'(/)')
         do i = 1, n_access_pts
          write(gunit,'(3f14.6,i6,2e14.6)')
     1      x_access(n,i), y_access(n,i), z_access(n,i), i, zz, zz
         end do
        end do
      end if
       
      write (munit, '(a3)') "end"
 
      close(munit)
      close(cunit)
      close(gunit)
      close(sunit)
 
      end subroutine write_coils

#if defined(GEOM_ONLY)
#else
      subroutine write_surfaces
      use cparms
      use boundary
      use modular_coils
      use saddle_coils
      use saddle_surface
      use Vwire
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, parameter :: nphi = 9
      integer :: i, k, m, n, nc, mn, itheta, iphi, nscmax
      integer :: npt
      integer, dimension(ncdim) :: npts
      real(rprec) :: v0, phi0, theta, dtheta, dphi, rws, zws
      real(rprec) :: ck, sk, cosmu, sinmu, cosnv, sinnv, rps, zps
      real(rprec), dimension (nwdim1,nphi) :: rw, zw, rs, zs, rp, zp
      real(rprec), dimension (nwdim1,2) :: xw, yw, xp, yp
      real(rprec), dimension (ncdim, 5*ncdim) :: rpts, zpts
      real(rprec), dimension (5*ncdim, nphi) :: rplt, zplt
!-----------------------------------------------
      open(unit=22,file='surfaces.dat',status='unknown')
      open(unit=25,file='topview.dat',status='unknown')
      open(unit=26,file='rzpts.dat',status='unknown')
 
      nscmax = ncdim/2
      dtheta = dtwopi/(nwire-1)
      dphi = dtwopi/(16*nfp)
      theta = zero

      do n = 1, nphi
        phi0 = -dpi/nfp + (n-1)*dphi

        if (lsaddle) then
           v0 = nfp*phi0/dtwopi
           call eval_saddle_rz (v0, npts, rpts, zpts)
           m = 0
           do k = 1, nsmid
              npt = npts(k)
              do i = 1, nscmax
                 m = m + 1
                 if (i .le. npt) then
                    rplt(m,n) = rpts(k,i)
                    zplt(m,n) = zpts(k,i)
                 else
                    rplt(m,n) = 1000
                    zplt(m,n) = 1000
                 end if
              end do
           end do
        end if

        do itheta = 1,nwire
           theta = dtheta*real(itheta-1)
 
!     Modular winding surface at phi    
           call rz_surf(theta,phi0,rws,zws,numsurf,
     1                   rmn_sf,zmn_sf,m_num,n_num,nfp)
           rw(itheta,n) = rws
           zw(itheta,n) = zws
 
!     Saddle winding surface edge at phi
           rws = 0.0
           zws = 0.0
           do k = 1, numsurf_sad
              ck = cos(m_sad(k)*theta + nfp*n_sad(k)*phi0)
              sk = sin(m_sad(k)*theta + nfp*n_sad(k)*phi0)
              rws = rws + rmn_sad(k)*ck
              zws = zws + zmn_sad(k)*sk
           end do
           rs(itheta,n) = rws
           zs(itheta,n) = zws
 
!     Plasma surface at phi
           rps = 0.0
           zps = 0.0
           do mn = 1, mnmax
              cosmu = cos(xm(mn)*theta)
              sinmu = sin(xm(mn)*theta)
              cosnv = cos(xn(mn)*phi0)
              sinnv = sin(xn(mn)*phi0)
              rps = rps + rmnc_b(mn)*(cosmu*cosnv + sinmu*sinnv)
              zps = zps + zmns_b(mn)*(sinmu*cosnv - cosmu*sinnv)
           end do
           rp(itheta,n) = rps
           zp(itheta,n) = zps
 
        end do
 
      end do
      do i = 1,nwire
        write(22,1000)(rw(i,mn),zw(i,mn), mn=1,9)
      end do
      write(22,'(/)')
      do i = 1,nwire
        write(22,1000)(rs(i,mn),zs(i,mn), mn=1,9)
      end do
      write(22,'(/)')
      do i = 1,nwire
        write(22,1000)(rp(i,mn),zp(i,mn), mn=1,9)
      end do
      if (lsaddle) then
         do m = 1,nsmid*nscmax
            write(26,1000) (rplt(m,n), zplt(m,n), n=1,9)
         end do
      end if
 1000 format(1p18e12.4) 
 
!     Now write data to plot topview
 
      dphi = dtwopi/(nwire-1)
      do n = 1,2
        theta = 0.5*(n-1)*dtwopi
        do iphi = 1,nwire
           phi0 = dphi*(iphi-1)
 
!     Modular winding surface edge at phi    
           call rz_surf(theta,phi0,rws,zws,numsurf,
     1                   rmn_sf,zmn_sf,m_num,n_num,nfp)
           xw(iphi,n) = rws*cos(phi0)
           yw(iphi,n) = rws*sin(phi0)
 
!     Plasma surface edge at phi
           rps = 0.0
           zps = 0.0
           do mn = 1, mnmax
              cosmu = cos(xm(mn)*theta)
              sinmu = sin(xm(mn)*theta)
              cosnv = cos(xn(mn)*phi0)
              sinnv = sin(xn(mn)*phi0)
              rps = rps + rmnc_b(mn)*(cosmu*cosnv + sinmu*sinnv)
              zps = zps + zmns_b(mn)*(sinmu*cosnv - cosmu*sinnv)
           end do
           xp(iphi,n) = rps*cos(phi0)
           yp(iphi,n) = rps*sin(phi0)
 
        end do
      end do
 
!     do i = 1,nwire
!       write(25,1000) xw(i,1), yw(i,1)
!     end do
!     write(25,'(/)')
      do i = 1,nwire
        write(25,1000) xp(i,1), yp(i,1)
      end do
      write(25,'(/)')
!     do i = 1,nwire
!       write(25,1000) xw(i,2), yw(i,2)
!     end do
!     write(25,'(/)')
      do i = 1,nwire
        write(25,1000) xp(i,2), yp(i,2)
      end do
 
      close(22)
      close(25)
      close(26)
 
!     write winding surface coefficients
      open (unit=27,file='rz_coeff.dat',status='unknown')
      open (unit=28,file='rz_coeff2.dat',status='unknown')
      write (27,115) numsurf
      write (27,120) (m_num(i), n_num(i), rmn_sf(i), zmn_sf(i),
     1                i=1, numsurf)
      write (28,115) numsurf_sad
      write (28,120) (m_sad(i), n_sad(i), rmn_sad(i), zmn_sad(i),
     1                i=1, numsurf_sad)
  115 format(i4)
  120 format(2i4,1p2e16.8)
      close (27)
      close (28)
 
!     write plasma surface coefficients
      open (unit=31,file='rz_plasma.dat',status='unknown')
      open (unit=32,file='rz_plasma2.dat',status='unknown')
      write (31,135) mnmax
      write (31,140) (xm(i), -xn(i)/nfp, rmnc_b(i), zmns_b(i),
     1                i=1, mnmax)
      write (32,135) mnmax
      write (32,145) (xm(i), -xn(i)/nfp, rmnc_b(i), zmns_b(i),
     1                i=1, mnmax)
  135 format(i4)
  140 format(2f8.0,2e22.12)
  145 format(2f10.0,2f20.12)
      close (31)
      close (32)
 
      end subroutine write_surfaces


      subroutine eval_saddle_rz (v0, npts, rpts, zpts)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary, only:nfp
      use saddle_coils
      use saddle_surface
      use Vcoilpts
      use Vwire
      use coils
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n, m, nc, np, npt, npoints, status
      integer, dimension(ncdim) :: npts
      real(rprec) :: u, v0, v1, v2, ds, cmax, dlsq, csum
      real(rprec) :: r1, s1, d1
      real(rprec) :: r2, s2, d2
      real(rprec), dimension(2,3,5) :: x1,  x2, y1,  y2, z1, z2
      real(rprec), dimension(2) :: c
      real(rprec), dimension(ncdim, 5*ncdim) :: rpts, zpts
      real(rprec) :: r, alf, bet
!-----------------------------------------------
 
      alf = 1
      bet = 0
      npoints = 2*nwire
 
      ds = 1.0/npoints
!     loop over unique coils
      do n=1, nsmid
         npt = 0
!        loop over number of poloidal points
         do i=1, npoints + 1
            s1 = (i - 1)*ds
            call coil_curve(n, s1, alf, bet, u, v1, r1, c, x1, y1, z1)
            s2 = i*ds
            call coil_curve(n, s2, alf, bet, u, v2, r2, c, x2, y2, z2)
            d1 = v1 - v0
            d2 = v2 - v0
            if (d1*d2 .lt. 0) then
               npt = npt + 1
               rpts(n,npt) = (r1 + r2)/2
               zpts(n,npt) = (z1(1,1,1) + z2(1,1,1))/2
               if (nfils .ge. 3) then
                  npt = npt + 1
                  r1 = sqrt(x1(1,1,2)**2 + y1(1,1,2)**2)
                  r2 = sqrt(x2(1,1,2)**2 + y2(1,1,2)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,2) + z2(1,1,2))/2
                  npt = npt + 1
                  r1 = sqrt(x1(1,1,3)**2 + y1(1,1,3)**2)
                  r2 = sqrt(x2(1,1,3)**2 + y2(1,1,3)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,3) + z2(1,1,3))/2
               end if
               if (nfils .eq. 5) then
                  npt = npt + 1
                  r1 = sqrt(x1(1,1,4)**2 + y1(1,1,4)**2)
                  r2 = sqrt(x2(1,1,4)**2 + y2(1,1,4)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,4) + z2(1,1,4))/2
                  npt = npt + 1
                  r1 = sqrt(x1(1,1,5)**2 + y1(1,1,5)**2)
                  r2 = sqrt(x2(1,1,5)**2 + y2(1,1,5)**2)
                  rpts(n,npt) = (r1 + r2)/2
                  zpts(n,npt) = (z1(1,1,5) + z2(1,1,5))/2
               end if
            end if
         end do
         npts(n) = npt
      end do
 
      end subroutine eval_saddle_rz
#endif

      subroutine optimize(xc,nopt,numsurf,mopt,nu,nv)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use cparms
      use control_mod
      use Vname, ONLY: extension, num_processors, num_levmar_params
      use mpi_params                                         !mpi stuff
      use modular_coils, ONLY: epsfcn, niter_opt, epsfcn_array, 
     1                         niter_array
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nopt, mopt, numsurf, nu, nv
!      integer :: nfev_end                                       !mpi stuff
      real(rprec), dimension(*) :: xc
!      real(rprec) :: epsfcn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: info, lwa, ierr_mpi, nfev_end, iter
      integer :: ierr, mode = 1                              !mpi stuff
      real(rprec), dimension(:), allocatable :: fvec
      real(rprec), dimension(:), allocatable :: diag    !mpi stuff
      real(rprec) :: tol
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
      external lsfun1
!-----------------------------------------------
#ifdef MPI_OPT
      call MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi)                         !MPI
      if (ierr_mpi .ne. 0) stop 'MPI_BARRIER error in OPTIMIZE'
#endif
 
      if (myid .eq. master) then                            !mpi stuff
        print *, 'no. independent variables = ', nopt
!       print *, 'no. poloidal match points = ', nu
!       print *, 'no. toroidal match points = ', nv
      end if     !if( myid .eq. master )                     !mpi stuff
 
 
      allocate(diag(nopt), fvec(mopt), stat = ierr)      
      lwa = mopt*(nopt + 1) + 5*nopt

      if( myid .eq. master ) then                               !mpi stuff
         if (ierr .ne. 0) stop 'allocation error in optimize'   !mpi stuff
      end if     !if( myid .eq. master )                        !mpi stuff
 
      tol = 1.e-6_dp
      mode = 1                                                               

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     BY PASS FOR COILGEOM ONLY (NO OPTIMIZATION, ONLY WRITE OUT COILS FILE
!
#if defined(GEOM_ONLY)
      call lsfun1(mopt, nopt, xc, fvec, info, nfev_end)
#else 
      if( niter_array(1) > 0 .and. epsfcn_array(1) .ne. 0) then
         nfev_end = niter_array(1)
         epsfcn = epsfcn_array(1)
      else
         nfev_end = niter_opt
      endif

      do iter = 1, 100
         if(NOPT_ALG .eq. 0) then

#if defined(MPI_OPT)

            call lmdif1 (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,     
     1            nfev_end, diag, mode, info, lwa)
#else
            call lmdif1 (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
     1             nfev_end, diag, mode, info, lwa, 
     2             num_processors, num_levmar_params)
#endif

         else if(NOPT_ALG .eq. 1) then
            call ga_driver (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
     1       nfev_end, num_processors, extension, info, lwa, lrestart )


         else if(NOPT_ALG .eq. 2) then
            call de_driver (lsfun1, mopt, nopt, xc, fvec, tol, epsfcn,
     1       nfev_end, num_processors, extension, info, lwa, lrestart )

         else
            write(6,*) "NOPT_ALG undefined, unable to proceed"
            stop
         endif

         if( iter==1 .and. (niter_array(1)<=0 .or. 
     1       epsfcn_array(1) == 0)) exit

         if( niter_array(iter+1)<=0 .or. epsfcn_array(iter+1)==0) then
            exit
         else
            nfev_end = niter_array(iter+1)
            epsfcn = epsfcn_array(iter+1)
         endif
      enddo
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if( myid .eq. master ) then                            !mpi stuff
        write (*, 67) info
   67   format(' info=',i5)
      end if     !if( myid .eq. master )                     !mpi stuff
      deallocate(diag, fvec)
 
      end subroutine optimize


      subroutine rz_surf(theta, phi, rw, zw, kmodes, rmn, zmn, m_num, 
     1   n_num, nfp)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use cparms
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer kmodes, nfp
      integer, dimension(*) :: m_num, n_num
      real(rprec) :: theta, phi, rw, zw
      real(rprec), dimension(*) :: rmn, zmn
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
 
      rw = zero
      zw = zero
 
! 2/16/98 WHM The NESCOIL sign convention is used here
 
      do i=1,kmodes
         rw = rw + rmn(i)*cos(m_num(i)*theta+nfp*n_num(i)*phi)
         zw = zw + zmn(i)*sin(m_num(i)*theta+nfp*n_num(i)*phi)
      end do
  
      end subroutine rz_surf


      subroutine get_modes (mmax, nmax, xm, xn, mnmax)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
      integer :: m, n, mn, mmax, nmax, mnmax
      real(rprec) :: xm(*), xn(*)
!-----------------------------------------------
 
      m = 0
      mn = 0
      do n = 0, nmax, 3
         mn = mn + 1
         xm(mn) = m
         xn(mn) = n
      end do
 
      do m = 1, mmax
         do n = -nmax, nmax, 3
            mn = mn + 1
            xm(mn) = m
            xn(mn) = n
         end do
      end do
 
      mnmax = mn
 
      end subroutine get_modes


#if defined(GEOM_ONLY)
#else
      subroutine b_error_spectrum
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use bnorm_mod
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
      integer :: n, mn
      real(rprec) :: sinmn
!-----------------------------------------------
 
!     Compute lamda  over u,v
      do n = 1, nedge
         luv(n) = zero
         do mn = 1, mnmax
            sinmn = sin (xm(mn)*thetab(n) - xn(mn)*phib(n))
            luv(n) = luv(n) + lmns_b(mn)*sinmn
         end do
      end do
 
!     Spectrum of b-error
      do mn = 1, mnmax_bmn
         bmn(mn) = zero
         do n = 1, nedge
            sinmn = sin (xm_bmn(mn)*thetab(n) - xn_bmn(mn)*phib(n))
            bmn(mn) = bmn(mn) + b_error(n)*sinmn
         end do
         bmn(mn) = 2.0_dp*bmn(mn)/nedge
!        if( myid .eq. master ) then                             !mpi stuff
!          print 100, nint(xm_bmn(mn)), nint(xn_bmn(mn))/3, bmn(mn)
!        end if     !if( myid .eq. master )                      !mpi stuff
      end do
 
!     Evaluate b-error in straight line coordinates 
      do n = 1, nedge
         bsl_error(n) = zero
         do mn = 1, mnmax
            sinmn = sin (xm(mn)*(thetab(n) + luv(n)) - xn(mn)*phib(n))
            bsl_error(n) = bsl_error(n) + bmn(mn)*sinmn
         end do
      end do
 
!     Spectrum of b-error in straight line coordinates
      do mn = 1, mnmax_bmn
         bmn(mn) = zero
         do n = 1, nedge
            sinmn = sin (xm_bmn(mn)*(thetab(n) + luv(n))
     1            - xn_bmn(mn)*phib(n))
            bmn(mn) = bmn(mn) + bsl_error(n)*sinmn*(1.0_dp + luv(n))
         end do
         bmn(mn) = 2.0_dp*bmn(mn)/nedge
         if( myid .eq. master ) then                             !mpi stuff
           print 100, nint(xm_bmn(mn)), nint(xn_bmn(mn))/3, bmn(mn)
         end if     !if( myid .eq. master )                     !mpi stuff
      end do
  100 format(2i5,1pe15.5)
 
      end subroutine b_error_spectrum


      subroutine b_error_mode (m)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cparms
      use boundary
      use bnorm_mod
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
 
      integer :: ku, kv, n, m, mn
      real(rprec) :: sinmn
!-----------------------------------------------
 
      open (unit=33,file='bmnerr.dat',status='unknown')
 
      mn = m
 
      bmn_error(:) = zero
      n = 0
      do ku = 1, nu
         do kv = 1, nv
            n = n + 1
            sinmn = sin (xm_bmn(mn)*(thetab(n) + luv(n))
     1            - xn_bmn(mn)*phib(n))
            bmn_error(n) = bmn(mn)*sinmn
            write (33,100) nfp*phib(n)/twopi, thetab(n)/twopi,
     1         bmn_error(n)
         end do
         write (33,110)
      end do
 
  100 format (1p3e15.6)
  110 format ("")
      close (33)
 
      end subroutine b_error_mode
#endif

      subroutine lsfun1 (mfun, nvar, xc, fvec, iflag, niter)
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use cparms
      use boundary
      use bnorm_mod
      use bcoils_mod
      use modular_coils
      use saddle_coils
      use vf_coils
      use tor_field
      use Vwire
      use Vname
      use times_stuff
      use safe_open_mod
      use control_mod, ONLY: nopt_alg
      use gade_mod, ONLY: save_space
      use mpi_params                                         !mpi stuff
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: mfun, nvar, iflag, niter
      real(rprec), dimension(nvar) :: xc
      real(rprec), dimension(mfun) :: fvec
      integer, parameter :: cleanup_flag = -100
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n, nfun, ntotal, ierr, n_coil, ncper
      integer, save :: erunit = 70, icnt = 0
      real (rprec) ::  chisq_total = min_chisq, chisq_max = min_chisq,
     1  rms_pd_err = min_chisq, rms_field_err = min_chisq,
     2  avg_field_err = min_chisq, max_field_err = min_chisq,
     3  chisq_field_err = min_chisq, chisq_cp_err = min_chisq
      real (rprec) :: cp_penalty, chisq_sc_err
      real (rprec) :: chisq_lc_err, chisq_cc_err, chisq_rc_err
      real (rprec) :: chisq_rs_err, chisq_ac_err, chisq_ymin_err
      real (rprec) :: chisq_cvf_err, chisq_cs_err, chisq_cu_err
      real (rprec) :: sp_penalty,   chisq_sp_err, chisq_scxp_err
      real (rprec) :: pc_penalty, chisq_pc_err, chisq_rvf_err
      real (rprec) :: chisq_csc_err, chisq_scd_err
      real (rprec) :: bkp_penalty, chisq_bkp_err, chisq_rmax_err
      real (rprec) :: mxb_penalty, chisq_mxb_err
      real (rprec), dimension(ncdim) :: lc_penalty, cc_penalty
      real (rprec), dimension(ncdim) :: rc_penalty, sc_penalty
      real (rprec), dimension(ncdim) :: rs_penalty, ac_penalty
      real (rprec), dimension(ncdim) :: cs_penalty, cu_penalty
      real (rprec), dimension(ncdim) :: ymin_penalty, cvf_penalty
      real (rprec), dimension(ncdim) :: scxp_penalty, csc_penalty,
     1                                  scd, scd_penalty
      real (rprec), dimension(ncdim) :: rmax_penalty, rvf_penalty
      real (rprec) :: csphi, snphi                        
      real(rprec), save :: total_time=0
      character*200 :: temp_input, input_file, output_file, opt_ext
      integer :: istat, iunit
      logical :: isthere, lscreen, lprint
!-----------------------------------------------
!
!         calculates the target functions at xc and
!         return this vector in fvec.
!
!       lsfun1 must be declared EXTERNAL in call to lmdif1
!
!       mfun is a positive integer input variable set to the number
!         of functions.
!
!       nvar is a positive integer input variable set to the number
!         of variables. nvar must not exceed mfun.
!
!       xc is an array of length nvar. on input xc must contain
!         an initial estimate of the solution vector. on output xc
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length mfun which contains
!         the functions evaluated at the output xc.
! 
!****************************************************************************
!****************************************************************************
!                                                                           *
!   **** NOTE CAREFULLY ****                                                *
!                                                                           *
!  if the sequence of values printed by this routine is changed,            *
!  the subroutine CHISQ_COILGEOM in Stellopt (which reads this output)      *
!   must also be changed to match                                           *
!                                                                           *
!****************************************************************************
!****************************************************************************


      if (iflag == cleanup_flag) then
#if defined(GEOM_ONLY) 
         iflag = 0
         return
#else
         if(nopt_alg .gt. 0 .and. save_space) then
         iflag = 0
         return
         else
         call clean_up (nvar, iflag)
         iflag = 0
         return
         end if
#endif
      end if
      
      fvec(:mfun) = zero
      nfun = 0
      call second0(strt_time)
       
!
!     compute unique extension for parallelized operation
!
#if defined(GEOM_ONLY) 
      icnt = icnt +1
      if (icnt .ge. niter_opt) iflag = -1
      if ((icnt .eq. 1) .or. (iflag .eq. 1) .or.
     1    (iflag .eq. -1)) lscreen=.true.
#else
      if (nopt_alg .gt. 0 .and. save_space) then
         icnt = icnt + 1
         if ((icnt .eq. 1) .or. (iflag .eq. 1) .or.
     1       (iflag .eq. -1)) lscreen=.true.
      else
         istat = iflag
         if (iflag .eq. -1) istat = 0
         write (temp_input,'(i5)') istat
         opt_ext = trim(extension)//'_opt'//trim(adjustl(temp_input))

         input_file = 'for05.' // trim(opt_ext)
         output_file = 'coil_targets.'// trim(opt_ext)

         if (iflag .ge. 0) then
            icnt = niter + iflag - 1
            lscreen = .false.
         else
            icnt = niter
            lscreen = .true.
         end if
      end if
#endif
!     Load model parameters with latest values of xc
 
      call loadparams (xc, nvar, mfun)

!     Evaluate coil filaments
       
      if (lmodular) call eval_modular_coils
      if (lsaddle) call eval_saddle_coils
      if (lvf) call eval_vf_coils

!
!     IF IN COILGEOM MODE (CALL FROM STELLOPT), RETURN UNLESS -P OPTION PASSED
!     FROM COMMAND LINE TO EVALUATE PLASMA-DEPENDENT PENALTIES BELOW
!     (FROM STELLOPT, THIS IS THE SECOND CALL TO COILGEOM)
!
      if (nedge .le. 0) return
 
!#if defined(GEOM_ONLY)
!      print *, 'coilgeom penalties'
!#else

!     Residuals due to normal magnetic field error

      call evaluate_field_error
      rms_field_err = sqrt (sum (b_error(1:nedge)**2)/nedge)
      avg_field_err = sum (d_area(1:nedge)*abs(b_error(1:nedge)))
     1               / sum (d_area(1:nedge))
      max_field_err = sqrt(maxval(b_error(1:nedge)**2))
      chisq_field_err = sum (b_error(1:nedge)**2)
      fvec (nfun+1:nfun+nedge) = b_error(1:nedge)
      nfun = nfun + nedge

#if defined(GEOM_ONLY)
      print *, 'coilgeom penalties'
      print *, 2, 'Berr penalties'
      print *, rms_field_err, max_field_err

#else 
!     Find b-error spectrum
 
!        if (mnmax_bmn .ne. 0) then
!           call b_error_spectrum
!           call b_error_mode (64)
!        end if

#endif

!     Residuals due to modular coil length penalties
 
      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_lc_err = zero
      if (lmodular .and. (.not. lsaddle)) then
         lc_penalty(1:ncper) =
     1      lmod_wgt(1:ncper)*(mod_length(1:ncper) - lmod_tgt(1:ncper))
         chisq_lc_err = sum (lc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = lc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'modular length penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper

!     Residuals due to saddle coil length penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_lc_err = zero
      if (lsaddle) then
!        For inequality constraint, try this
         do i = 1, ncper
           if (lsad_tgt(i) .lt. sad_length(i)) then
             lc_penalty(i) = lsad_wgt(i)*(sad_length(i) - lsad_tgt(i))
           else
             lc_penalty(i) = 0
           end if
         end do
!        For equality constraint, try this
!        lc_penalty(1:ncper) =
!    1      lsad_wgt(1:ncper)*(sad_length(1:ncper) - lsad_tgt(1:ncper))
         chisq_lc_err = sum (lc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = lc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle length penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to modular coil ymin penalties
 
      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_ymin_err = zero
      if ((lmodular .eqv. .true.) .and. (lsaddle .eqv. .false.)) then
         do i = 1, ncper
           if (ymin_tgt(i) .gt. ymin_cls(i)) then
             ymin_penalty(i) = ymin_wgt(i)*(ymin_tgt(i) - ymin_cls(i))
           else
             ymin_penalty(i) = 0
           end if
         end do
         chisq_ymin_err = sum (ymin_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ymin_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'modular ymin penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to saddle coil ymin penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_ymin_err = zero
      if (lsaddle) then
         do i = 1, ncper
           if (ymin_tgt(i) .gt. ymin_sad(i)) then
             ymin_penalty(i) = ymin_wgt(i)*(ymin_tgt(i) - ymin_sad(i))
           else
             ymin_penalty(i) = 0
           end if
         end do
         chisq_ymin_err = sum (ymin_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ymin_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle ymin penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to minimum modular coil-coil distance penalties
 
      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cc_err = zero
      if (lmodular) then
         call mod_coil_distance
         cc_penalty(1:ncper) =
     1      dcc_wgt(1:ncper)*exp(-dcc_exp(1:ncper)*(cc_min(1:ncper)
     2         - dcc_tgt(1:ncper))*cmod_scl(1:ncper))
         chisq_cc_err = sum (cc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'modular coil-coil distance penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to minimum saddle coil-coil distance penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_sc_err = zero
      if (lsaddle) then
         call sad_coil_distance
         do i = 1,ncper
            if (dsc_exp(i) .ge. 0) then    ! use exp penalty
               sc_penalty(i) =
     1            dsc_wgt(i)*exp(-dsc_exp(i)
     2               *(sc_min(i) - dsc_tgt(i)))
            else                              ! use abs penalty
               if (sc_min(i) .lt. dsc_tgt(i)) then
                  sc_penalty(i) = dsc_wgt(i)*
     1               (dsc_tgt(i) - sc_min(i))
               else
                  sc_penalty(i) = 0
               end if
            end if
         end do
!        sc_penalty(1:ncper) = 
!    1      dsc_wgt(1:ncper)*exp(-dsc_exp(1:ncper)*(sc_min(1:ncper)
!    2         - dsc_tgt(1:ncper)))
         chisq_sc_err = sum (sc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = sc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle coil-coil distance penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to minimum saddle c-c dist. X c-p dist. penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_scxp_err = zero
      if (lsaddle) then
!     Note: scxp_min calculated in sad_coil_distance (called above)
         do i = 1,ncper
            if (dscxp_exp(i) .ge. 0) then    ! use exp penalty
               scxp_penalty(i) = 
     1            dscxp_wgt(i)*exp(-dscxp_exp(i)
     2               *(scxp_min(i) - dscxp_tgt(i)))
            else                              ! use abs penalty
               if (scxp_min(i) .lt. dscxp_tgt(i)) then
                  scxp_penalty(i) = dscxp_wgt(i)*
     1               (dscxp_tgt(i) - scxp_min(i))
               else
                  scxp_penalty(i) = 0
               end if
            end if
         end do

         chisq_scxp_err = sum (scxp_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = scxp_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle c-c dist. X c-p dist. penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to minimum modular coil radius of curvature penalties
 
      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rc_err = zero
      if (lmodular) then
         call coil_curvature
         rc_penalty(1:ncper) =
     1      rc_wgt(1:ncper)*exp(-rc_exp(1:ncper)*(rc_min(1:ncper)
     2         - rc_tgt(1:ncper))*cmod_scl(1:ncper))
         chisq_rc_err = sum (rc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'modular coil radius of curvature penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to minimum saddle radius of curvature penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rs_err = zero
      if (lsaddle) then
         do i = 1,ncper
            if (rs_exp(i) .ge. 0) then    ! use exp penalty
               rs_penalty(i) =
     1            rs_wgt(i)*exp(-rs_exp(i)
     2               *(rs_min(i) - rs_tgt(i)))
            else                              ! use abs penalty
               if (rs_min(i) .lt. rs_tgt(i)) then
                  rs_penalty(i) = rs_wgt(i)*
     1               (rs_tgt(i) - rs_min(i))
               else
                  rs_penalty(i) = 0
               end if
            end if
         end do
!        rs_penalty(1:ncper) =
!    1      rs_wgt(1:ncper)*exp(-rs_exp(1:ncper)*(rs_min(1:ncper)
!    2         - rs_tgt(1:ncper)))
         chisq_rs_err = sum (rs_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rs_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle coil radius of curvature penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to modular coil integrated curvature penalties
 
      ncper = nmod_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cu_err = zero
      if (lmodular) then
         cu_penalty(1:ncper) =
     1      cu_wgt(1:ncper)*(cu_sum(1:ncper) - cu_tgt(1:ncper))
         chisq_cu_err = sum (cu_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cu_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'modular coil integrated curvature penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to saddle coil integrated curvature penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cs_err = zero
      if (lsaddle) then
         cs_penalty(1:ncper) =
     1      cs_wgt(1:ncper)*(cs_sum(1:ncper) - cs_tgt(1:ncper))
         chisq_cs_err = sum (cs_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cs_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle coil integrated curvature penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
      
!     Residuals due to saddle coil current regularization penalty
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_csc_err = zero
      if (lsaddle) then
         csc_penalty(1:ncper) = csc_wgt(1:ncper)
     1                         *(cursad(1:ncper) - csc_tgt(1:ncper))
         chisq_csc_err = sum (csc_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = csc_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'Saddle current regularization penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper

!     Residuals due to saddle maximum current density penalty, using
!     saddle c-c dist. X c-p dist. calculation
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_scd_err = zero
      if (lsaddle) then
!     Note: scxp_min calculated in sad_coil_distance (called above)
         do i = 1,ncper

            if( scxp_min(i) .lt. 1.e-10) then
               scd_penalty(i) = 1.e10
               scd(i) = 1.e10

            else
               scd(i) = abs(cursad(i)/scxp_min(i))

               if (scd(i) .gt. scd_tgt(i)) then
                  scd_penalty(i) = scd_wgt(i)*
     1                  (scd(i) - scd_tgt(i))
               else
                  scd_penalty(i) = 0
               end if
            end if
         end do

         chisq_scd_err = sum (scd_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = scd_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle lin. current density penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
         end if
         nfun = nfun + ncper
 
!     Residuals due to saddle coil rmax penalties
 
      ncper = nsad_unique_coils
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rmax_err = zero
      if (lsaddle) then
         do i = 1, ncper
           if (rmax_tgt(i) .lt. rmax_sad(i)) then
             rmax_penalty(i) = rmax_wgt(i)*(rmax_sad(i) - rmax_tgt(i))
           else
             rmax_penalty(i) = 0
           end if
         end do
         chisq_rmax_err = sum (rmax_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rmax_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'saddle rmax penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper

!     Residual due to minimum plasma-modular coil distance penalty
 
      fvec (nfun+1) = zero
      chisq_cp_err = zero
      if (lmodular) then
         call plasma_mod_coil_distance
         cp_penalty = dcp_wgt*exp(-dcp_exp*(p_d_min-dcp_tgt))
         chisq_cp_err = cp_penalty**2
         fvec (nfun+1) = cp_penalty
#if defined(GEOM_ONLY)
         print *, '1 plasma-modular coil distance penalty'
         print *, fvec (nfun+1)
#else
#endif
      end if
      nfun = nfun + 1
 
!     Residual due to minimum plasma-saddle coil distance penalty
 
      fvec (nfun+1) = zero
      chisq_sp_err = zero
      if (lsaddle) then
         call plasma_sad_coil_distance
         if (dcp_exp .ge. 0) then
            sp_penalty = dcp_wgt*exp(-dcp_exp*(p_s_min-dcp_tgt))
         else
            if (p_s_min .lt. dcp_tgt) then
               sp_penalty = dcp_wgt*(dcp_tgt - p_s_min)
            else
               sp_penalty = 0
            end if
         end if
         chisq_sp_err = sp_penalty**2
         fvec (nfun+1) = sp_penalty
#if defined(GEOM_ONLY)
         print *, '1 plasma-saddle coil distance penalty'
         print *, fvec (nfun+1)
#else
#endif
      end if
      nfun = nfun + 1
 
!     Residuals due to minimum access penalties
 
      ncper = n_access
      fvec (nfun+1:nfun+ncper) = zero
      chisq_ac_err = zero
      if (laccess) then
         call evaluate_access
         do i = 1,ncper
            if (dac_exp(i) .ge. 0) then    ! use exp penalty
               ac_penalty(i) =
     1            dac_wgt(i)*exp(-dac_exp(i)
     2               *(acc_min(i) - dac_tgt(i)))
            else                              ! use abs penalty
               if (acc_min(i) .lt. dac_tgt(i)) then
                  ac_penalty(i) = dac_wgt(i)*
     1               (dac_tgt(i) - acc_min(i))
               else
                  ac_penalty(i) = 0
               end if
            end if
         end do
!        ac_penalty(1:ncper) =
!    1      dac_wgt(1:ncper)*exp(-dac_exp(1:ncper)*(acc_min(1:ncper)
!    2         - dac_tgt(1:ncper)))
         chisq_ac_err = sum (ac_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = ac_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'minimum access penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper
 
!     Residuals due to vf coil current regularization penalty
 
      ncper = num_vf
      fvec (nfun+1:nfun+ncper) = zero
      chisq_cvf_err = zero
      if (lvf) then
         cvf_penalty(1:ncper) = cvf_wgt(1:ncper)*rc_vf(1:ncper)
     1                         *(cc_vf(1:ncper) - cvf_tgt(1:ncper))
         chisq_cvf_err = sum (cvf_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = cvf_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'VF regularization penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper

!     Residuals due to vf coil radius penalty
 
      ncper = num_vf
      fvec (nfun+1:nfun+ncper) = zero
      chisq_rvf_err = zero
      if (lvfvar) then
         do i = 1,ncper
            if (rvf_max(i) .gt. rvf_tgt(i)) then
               rvf_penalty(i) = rvf_wgt(i)*
     1            (rvf_max(i) - rvf_tgt(i))
            else
               rvf_penalty(i) = 0
            end if
         end do
         chisq_rvf_err = sum (rvf_penalty(1:ncper)**2)
         fvec (nfun+1:nfun+ncper) = rvf_penalty(1:ncper)
#if defined(GEOM_ONLY)
         print *, ncper, 'VF radius penalties'
         print *, (fvec (i), i=nfun+1, nfun+ncper)
#else
#endif
      end if
      nfun = nfun + ncper

!     Residual due to total poloidal current penalty

      fvec (nfun+1) = zero
      chisq_pc_err = zero
      if (lpolcur) then
         call eval_poloidal_currents
         pc_penalty = dpc_wgt*(pol_cur/nfp - i_pol)/i_pol
         chisq_pc_err = pc_penalty**2
         fvec (nfun+1) = pc_penalty
#if defined(GEOM_ONLY)
         print *, '1  Total poloidal current penalty'
         print *, fvec (nfun+1)
#else
#endif
      end if
      nfun = nfun + 1

!     Residual due to spline breakpoint separation penalty

      fvec (nfun+1) = zero
      chisq_bkp_err = zero
      if (lspline .and. lsplbkp) then
         call spline_bkp_distance
         if (bkp_min .lt. bkp_tgt) then
            bkp_penalty = bkp_wgt*(bkp_min - bkp_tgt)
         else
            bkp_penalty = 0
         end if
         chisq_bkp_err = bkp_penalty**2
         fvec (nfun+1) = bkp_penalty
#if defined(GEOM_ONLY)
!        print *, '1  Spline breakpoint separation penalty'
!        print *, fvec (nfun+1)
#else
#endif
      end if
      nfun = nfun + 1

!     Residual due to maximum field error penalty

      mxb_penalty = mxb_wgt*max_field_err
      chisq_mxb_err = mxb_penalty**2
      fvec (nfun+1) = mxb_penalty
#if defined(GEOM_ONLY)
!     print *, '1  Maximum field error penalty'
!     print *, fvec (nfun+1)
#else
#endif
      nfun = nfun+1
      if (nfun .gt. mfun) stop 'NFUN > MFUN in LSFUN1'

!     Total error
 
      chisq_total = chisq_field_err + chisq_lc_err   + chisq_cc_err
     1            + chisq_rc_err    + chisq_cp_err   + chisq_sc_err
     2            + chisq_rs_err    + chisq_ac_err   + chisq_ymin_err
     3            + chisq_cvf_err   + chisq_cs_err   + chisq_cu_err
     4            + chisq_sp_err    + chisq_pc_err   + chisq_scxp_err
     5            + chisq_bkp_err   + chisq_mxb_err  + chisq_rvf_err
     6            + chisq_rmax_err  + chisq_csc_err  + chisq_scd_err
 
#if defined(GEOM_ONLY)
#else
      if (nopt_alg.gt.0 .and. save_space) then
         if (chisq_total .lt. chisq_min) then
            chisq_min = chisq_total
            inquire(file=for05_file_new,exist=isthere)
            if (myid.eq.master .and. isthere)
     1         call save_for05 (for05_file_new)
         end if
      else
         call save_for05(input_file)
      end if
#endif
!     Output section
 
      if (myid .eq. master ) then
         call second0(fin_time)
         total_time=total_time+(fin_time-strt_time)      
         lprint = (((icnt .eq. 1) .or. (iflag .eq. 1) .or.
     1              (iflag .eq. -1)) .and. lscreen) 
#if defined(GEOM_ONLY)
         do iunit = 6, 6
#else
         do iunit = 6, 15, 9
#endif      
            if (iunit.eq.6 .and. .not.lprint) cycle
            if (iunit.eq.15 .and. 
     1         (nopt_alg.ne.0 .and. save_space)) cycle
            if (iunit.eq.15) call safe_open(iunit,
     1           istat, output_file, 'replace', 'formatted')
            write(iunit,160) icnt
            write(iunit,100)
            write(iunit,110)rms_field_err, avg_field_err, max_field_err,
     1                       chisq_total
            write(iunit,130)
            write(iunit,110)chisq_field_err, chisq_mxb_err,chisq_cc_err,
     1                       chisq_rc_err
            write(iunit,135)
            write(iunit,110) chisq_cp_err, chisq_ac_err, chisq_ymin_err,
     1                       chisq_cvf_err
            write(iunit,310)
            write(iunit,110) chisq_pc_err, chisq_rvf_err, chisq_csc_err, 
     1                       chisq_scd_err
            if (lmodular) then
               write(iunit,340)
               write(iunit,110) (curmod(i), i=1, nmid)
            end if
            if (lsaddle) then
               if (lsmod) then
                  write(iunit,340)
               else
                  write(iunit,350)
               end if
               write(iunit,110) (cursad(nsad_group(i))
     1                          *csad_scl(i), i=1, nsmid)
            end if
            if (lvf) then
               write(iunit,360)
               write(iunit,110) (cc_vf(i), i=1, num_vf)
            end if
            if (lbcoil) then
               write(iunit,370)
               write(iunit,110) (bcoil_cur(i), i=1, mbcoils)
            end if
            write(iunit,320)
            write(iunit,110) pol_cur/nfp, i_pol
            if (lmodular) then
               write(iunit,145)
               write(iunit,110) chisq_cu_err, chisq_lc_err
               write(iunit,250)
               write(iunit,170)
               write(iunit,110) (mod_length(i), i=1, nmod_unique_coils)
               write(iunit,180)
               write(iunit,110) (lmod_tgt(i), i=1, nmod_unique_coils)
               write(iunit,140)
               write(iunit,110) (cc_min(i), i=1, nmod_unique_coils)
               write(iunit,190)
               write(iunit,110) (rc_min(i), i=1, nmod_unique_coils)
               write(iunit,280)
               write(iunit,110) (cu_sum(i), i=1, nmod_unique_coils)
               write(iunit,120)
               write(iunit,110) p_d_min
            end if
            if (lsaddle) then
               if (lsmod) then
                  write(iunit,250)
               else
                  write(iunit,220)
               end if
               write(iunit,230)
               write(iunit,110)chisq_sc_err, chisq_rs_err, chisq_cs_err,
     1                         chisq_sp_err
               write(iunit,235)
               write(iunit,110) chisq_scxp_err, chisq_lc_err,
     1                          chisq_rmax_err
               write(iunit,200)
               write(iunit,110) (sc_min(i), i=1, nsad_unique_coils)
               write(iunit,205)
               write(iunit,110) (scxp_min(i), i=1, nsad_unique_coils)
               write(iunit,210)
               write(iunit,110) (rs_min(i), i=1, nsad_unique_coils)
               write(iunit,280)
               write(iunit,110) (cs_sum(i), i=1, nsad_unique_coils)
               write(iunit,240)
               write(iunit,110) (sad_length(i), i=1, nsad_unique_coils)
               write(iunit,290)
               write(iunit,110) (ymin_sad(i), i=1, nsad_unique_coils)
               write(iunit,380)
               write(iunit,110) (rmax_sad(i), i=1, nsad_unique_coils)
               write(iunit,'(6x,a)') 
     1                     "max saddle linear current density (A/m)"
               write(iunit,110) (scd(i), i=1, nsad_unique_coils)
               write(iunit,120)
               write(iunit,110) p_s_min
            end if
            if (laccess) then
               write(iunit,260)
               write(iunit,270)
               write(iunit,110) (acc_min(i), i=1, n_access)
            end if
            if (lspline .and. lsplbkp) then
               write (iunit, 125)
               write (iunit, 110) bkp_min, chisq_bkp_err
            end if
            write(iunit,150)
            if (iunit .ne. 6) close(iunit)
         end do
      end if

!     Reset iflag to no-error condition
      iflag = 0

  100 format (7x,'rms error',7x,'avg error',7x,'max error',
     1 7x,'chisq_tot')
  110 format (5e16.4)
  120 format (6x,' min plasma-coil distance (m)')
  125 format (6x,' bkp min',6x,'chisq_bkp')
  130 format (6x,'chisq_field',7x,'chisq_mxb',7x,'chisq_cc',
     1 8x,'chisq_rc')
  135 format (8x,'chisq_cp',7x,'chisq_acc',7x,'chisq_ymin',
     1 6x,'chisq_cvf')
  145 format (8x,'chisq_cu',8x,'chisq_lc')
  140 format (6x,' min modular coil-coil distance (m)')
  150 format (18x,'__________________________________________')
  160 format (/,' N=',i6)
  170 format (6x,' mod coil lengths (m)')
  180 format (6x,' mod coil length targets (m)')
  190 format (6x,' mod coil min radius of curvature (m)')
  200 format (6x,' min saddle coil-coil distance (m)')
  205 format (6x,' |min. del(c-c) X min. del(c-p)|/|min. del(c-p)| (m)')
  210 format (6x,' min saddle radius of curvature (m)')
  220 format (18x,'_______________Saddle coils_______________')
  230 format (6x,'chisq_sc',8x,'chisq_rs',8x,'chisq_cs',8x,'chisq_sp')
  235 format (6x,'chisq_scxp',8x,'chisq_lc',9x,'chisq_rmax')
  240 format (6x,' sad coil lengths (m)')
  250 format (18x,'_______________Modular coils______________')
  260 format (18x,'_______________Access zones_______________')
  270 format (6x,' min distance to coils (m)')
  280 format (6x,' integrated coil curvature (1/m)')
  290 format (6x,' minimum y - coils (m)')
  310 format (6x,'chisq_pc',6x,'chisq_rvf',7x,'chisq_csc',6x,
     1           'chisq_scd')
  320 format (6x,' total poloidal current per fp, required i-pol (amp)')
  340 format (6x,' modular currents (amp)')
  350 format (6x,' saddle currents (amp)')
  360 format (6x,' vertical field currents (amp)')
  370 format (6x,' background currents (amp)')
  380 format (6x,' maximum r - coils (m)')
 
      end subroutine lsfun1
#if defined(GEOM_ONLY)
#else
      subroutine clean_up (nvar, iflag)
      use control_mod
      use Vname
      use safe_open_mod
      use system_mod
      use mpi_params
      use gade_mod, only: npopsiz
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                                   !MPI
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nvar, iflag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iunit, icount, istat, istat2, itemp, icmax, iteration
      real(rprec) ::  chisq_tot
      logical :: lexist, lnew_min
      character*200 :: temp_input, save_min_ext, testfile, opt_ext
      real(rprec) :: rms_err, avg_err, max_err
C-----------------------------------------------

      iflag = 0
      save_min_ext = " "
      lnew_min = .false.

      if (myid .ne. master) go to 2000

      iunit = 14
      chisq_tot = chisq_min
      icmax = max(nvar, num_levmar_params)
      if (nopt_alg .gt. 0) icmax = max(icmax, npopsiz)

      findmin: do icount = 0, icmax
         write (temp_input,'(i5)') icount
         opt_ext =
     1      trim(extension)//'_opt'//trim(adjustl(temp_input))
         temp_input = 'coil_targets.' // trim(opt_ext)

         call safe_open(iunit, istat, temp_input, 'old',
     1             'formatted')

         do while (istat .eq. 0)
            read  (iunit,'(a)', iostat=istat) temp_input
            if (index(temp_input,'chisq_tot').ne.0 ) exit
         end do

         if (istat .eq. 0) then
            read (iunit,'(5e16.4)',iostat=istat)
     1             rms_err, avg_err, max_err, chisq_tot            


            if (chisq_tot .lt. chisq_min) then
               chisq_min = chisq_tot
               save_min_ext = trim(opt_ext)
               lnew_min = .true.
            end if

         end if

         close (iunit)

      end do findmin

!
!     STORE MIN STATE FILES, BUT ONLY IF ONE WAS FOUND
!
      if (lnew_min) then
 
         temp_input = "/bin/ls -1 *." // trim(save_min_ext)
     1                // " > min_file"
         call system(temp_input, istat)
         itemp = iunit+1
         call safe_open(itemp, istat, "min_file", "old", "formatted")
         do while (.true.)
            read (itemp, '(a)', end=50) testfile
            icount = index(testfile, trim(save_min_ext)) - 1
            temp_input = "/bin/mv -f " // trim(testfile)
     1                // " " // testfile(1:icount) // trim(extension)
     2                // ".min"
            call system(temp_input, istat)
         end do
 50      continue
         close (itemp, status='delete')

      end if              !lnew_min

      temp_input = "/bin/ls -1 > all_files"
      call system(temp_input)
      itemp = iunit + 1
      call safe_open(itemp, istat, "all_files", "old", "formatted")
      do while (.true.)
         read (itemp, '(a)', end=60) testfile
!      Avoid wiping out necessary files for continuing run
         if ((nopt_alg.eq.1 .and. index(testfile,"ga_restart").gt.0)
     1     .or. (nopt_alg.eq.2 .and. index(testfile,"de_restart").gt.0)
     2     .or. (trim(testfile) .eq. "all_files")
     3     .or. (index(testfile,".dat") .gt. 0)
     4     .or. (index(testfile,"wout.") .gt. 0)
     5     .or. (index(testfile,"bnorm.") .gt. 0)
     6     .or. (index(testfile,".new") .gt. 0)) cycle

         icount = index(testfile, ".min", back=.true.)
         if (icount .eq. (len_trim(testfile) - 3)) cycle
         temp_input = "rm -f " // trim(testfile)
         call system(temp_input, istat)
      end do

 60   close (itemp, status='delete', iostat=istat)

 2000 continue

#ifdef MPI_OPT
      call MPI_BCAST(chisq_min, 1, MPI_REAL8, master, MPI_COMM_WORLD,
     1     istat)
      if (istat .ne. 0) stop 'MPI_BCAST error in STELLOPT CLEAN_UP'
      call MPI_BCAST(iflag, 1, MPI_INTEGER, master, MPI_COMM_WORLD,
     1     istat)
      if (istat .ne. 0) stop 'MPI_BCAST error in STELLOPT CLEAN_UP'
#endif
      end subroutine clean_up
#endif
EOF
EOC
