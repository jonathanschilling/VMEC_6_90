       program cobra
c___________________________________________________________________________________
c                                                                                   |
c              COBRA (COde for Ballooning Rapid Analysis)                           |
c                     ==       =          =     =                                   |
c                                                                                   |
c      VERSION: 3.0  (FORTRAN 90: fixed format; standard I/O; VMEC coordinates)     |
c      Last update: 01/26/00                                                        | 
c                                                                                   |
c      AUTHORS: R. Sanchez(*) and S.P. Hirshman(**)                                 |
c      (*)    Universidad Carlos III de Madrid, Madrid 28911, SPAIN                 |
c      (**)   Oak Ridge National Laboratory, Oak Ridge, TN 37831-9701, USA          |
c                                                                                   |
c      REFERENCES:                                                                  |
c                                                                                   |
c        1."COBRA: and optimized code for fast analysis of ideal ballooning         |
c        stability of 3-D magnetic equilibria", R. Sanchez, S.P. Hirshman, J.C.     |
c        Whitson and A.S. Ware, submitted to Journal of Computational Physics (1999)|
c                                                                                   |
c        2."Improved magnetic coordinate representation for ideal ballooning        |
c        stability calculations with the COBRA code", R. Sanchez, S.P. Hirshman     |
c        H.V. Wong, to be submitted to Journal of Computational physics (2000)      |
c                                                                                   |
c      DISCLAIMER:  This code is under development by R.Sanchez at the Departamento |
c        de Fisica, Universidad Carlos III de Madrid, SPAIN. As a BETA version, the |
c        code is supplied on "as it is" basis, and the non-existence of "bugs" of   |
c        any kind is NOT guaranteed. Any problem or comment should be reported to   |
c        R. Sanchez at rsanchez@fis.uc3m.es.                                        |
c___________________________________________________________________________________|
c
       use kind_spec
       use normalize_data
       use ballooning_data
       use readin_data
       use safe_open_mod
       use date_and_computer
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: unit_cobra = 12, input_cobra = 14
      character*(*), parameter :: banner =
     1   ' THIS IS THE (VMEC-BASED) COBRA BALLOONING CODE Version '
      character*(5), parameter :: cobra_version = '3.00'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec), dimension(:), allocatable :: grate 
       real(rprec) :: t1,t2
       integer :: j, nlis, i, iunit_out, numargs, iunit_in, ij, jk
       integer :: istat, ierr, ns_surf, nini_zeta, nini_theta
       integer,dimension(:),allocatable:: in_surf, in_idx, bsurf
       integer :: imon

       character*120 :: extension  
       character*10 :: date0, time0, zone0
       character*120 :: dateloc, arg1, arg2
!-----------------------------------------------
!    
!     Read command line argument to get input file or sequence file
!
      lscreen = .true.
      
      call getcarg(1, arg1, numargs)
      if (numargs .gt. 1) call getcarg(2, arg2, numargs)

      if (numargs .lt. 1) then
         stop 'Invalid command line'
      else if (arg1 .eq. '-h' .or. arg1 .eq. '/h') then
         print *,
     1   ' ENTER INPUT FILE NAME ON COMMAND LINE'
         print *,' For example: xcobra in_cobra.tftr'
         print *
         print *,' Optional command line argument:'
         print *,' xcobra <input> (T or F)'
         print *
         print *,' where F will suppress screen output'
         stop
      else if (numargs .gt. 1) then
         if (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.   
      endif


!      INPUT FILE
!      1st line:   nini_tot(number starting angle pairs), extension of WOUT file
!      2nd line:   init_zeta_v, init_theta_v (initial toroidal and poloidal angle, in degrees) vectors
!      3rd line:   ns_surf (number surfaces to compute growth rate on)
!      4rd line:   surfaces to be computed (read in ORDER_INPUT)

      iunit_in = input_cobra
      call safe_open(iunit_in, istat, trim(arg1), 'old', 'formatted')

      if (istat .ne. 0) stop ' Error opening input file in COBRA'

!     read number initial points, extension of WOUT VMEC file
      read (iunit_in, *, iostat=istat) nini_zeta, nini_theta, extension      
       
      allocate (init_zeta_v(nini_zeta), init_theta_v(nini_theta), 
     1          stat=istat) 

!     read vector of ballooning initial angle points     
      read (iunit_in, *, iostat=istat) (init_zeta_v(j), j=1,nini_zeta)
      read (iunit_in, *, iostat=istat) (init_theta_v(j),j=1,nini_theta)

      read (iunit_in, *, iostat=istat) ns_surf                  ! read no surfaces

      allocate (in_surf(ns_surf), in_idx(ns_surf))

      read (iunit_in, *, iostat=istat) in_surf(1:ns_surf)       ! read surfaces on the full grid to be computed
      close (iunit_in)
       
      iunit_out = unit_cobra
      call safe_open(iunit_out, istat, 'cobra_grate.'//extension, 
     1     'replace', 'formatted')

      if (istat .ne. 0) then
          write (iunit_out, *) ' Error opening COBRA_GRATE'
          close (iunit_out)
          stop ' Error opening COBRA_GRATE'
      end if    

      if (lscreen) then
         print 48
         print *
         call date_and_time(date0,time0,zone0)
         read (date0(5:6),'(i2)') imon
         write (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1      time0(1:2),time0(3:4),time0(5:6)
         write (*,'(1x,2a,/,1x,2a)') banner, cobra_version, 
     1      computer, dateloc
         print *
      end if
 100  format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      nlis = 0
      do i = 1, ns_surf
         if (in_surf(i) .gt. 1) then
            nlis = nlis+1
            in_idx(nlis) = in_surf(i)
         endif
      enddo

      allocate (bsurf(nlis))
      bsurf(:nlis) = in_idx(:nlis)
       
!...  read equilibrium data and form surface list 

      call order_input(extension, nlis, bsurf, istat)
      if (istat .ne. 0) then
         extension = ' Error reading WOUT.' // trim(extension)
     1         // ' in COBRA'
         write (iunit_out, *) trim(extension)
         close (iunit_out)
         print *, trim(extension)
         stop 
      end if   
       
!..   solve ballooning equation

      allocate (grate(ns_cob))               ! allocate growth rate vector     
      call second0(t1)            

      theta: do ij = 1, nini_theta
         zeta: do jk = 1, nini_zeta

            init_zeta = init_zeta_v(jk)
            init_theta = init_theta_v(ij)
            call get_ballooning_grate(grate)       ! get growth rates
            call second0(t2)

!...   generate output
      
            if (lscreen) then
               write(*,120) init_zeta, init_theta, t2-t1
               write(*,48)
            end if
 120  format(/,'ZETA0 = ', 1pe10.3,' THETA0 = ', 1pe10.3,
     1        ' TIME IN COBRA CODE:',1pe10.2,' SEC')
48    format('====================================================')

            write (iunit_out, 12) init_zeta, init_theta, nlis
            write (iunit_out, 14) (bsurf(j), grate(bsurf(j)),j=1,nlis)        ! write output file
       
         end do zeta
      end do theta

12    format(1p2e10.3,i5)
14    format(i4, 1pe16.8)
      close (unit=iunit_out)
      deallocate (hiota, hpres, hphip, rmncf, zmnsf, list,
     1   lmnsh, bmnch, bsupumnh, bsupvmnh, mercierf, xn_v, xm_v)
      deallocate (grate, bsurf, in_surf, in_idx, init_theta_v, 
     1       init_zeta_v) ! deallocate growth rate vector     

      end program cobra
       

       subroutine order_input(extension, nl, bsurf, ierr)
       use kind_spec
       use read_wout_mod
       use normalize_data
       use readin_data
       use fmesh_quantities, ONLY: rmncf, zmnsf, mercierf
       use vparams, only: dmu0
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       character*(*), intent(IN) :: extension
       integer, intent(IN) :: nl
       integer, intent(OUT) :: ierr
       integer, dimension(nl), intent(IN):: bsurf
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer:: mn, lk, i, j, k
       integer, dimension(:), allocatable :: plist
       real(rprec) :: aspect_v, r0max_v, r0min_v, betaxis_v,
     1   rdum
       logical, allocatable, dimension(:):: lballoon
!-----------------------------------------------
!
!      Variables allocated in this routine (deallocated elsewhere)
!         hiota, hpres, hphip, xm_v, xn_v, rmncf, zmnsf, mercierf,
!         lmnsh, bmnch, list, bsupumnh, bsupvmnh
!
!      All half mesh quantities have first surface (js =1) equal to zero
!
!...   read surface quantities and normalization data

       call read_wout_file('wout.'//trim(extension), ierr)
       if (ierr .ne. 0) then
         print *,' ierr = ', ierr,
     1   'error in read_wout_file called from COBRA'
       end if
       if (ierr_vmec.ne.0 .or. ierr.ne.0) goto 1000
       
       nfp_v = nfp
       ns_cob = ns
       aspect_v = aspect
       r0max_v = rmax_surf 
       r0min_v = rmin_surf
       betaxis_v = betaxis
       mnmax_v = mnmax

!...   IOTA, PRES, PHIP are ALL on HALF-MESH!!!   
    
       allocate (hiota(ns_cob), hpres(ns_cob), hphip(ns_cob), 
     1   mercierf(ns_cob), stat=k)
       if (k .ne. 0) stop 'Allocation error 1 in cobra order_input'

       hiota = iotas(1:ns_cob)
       if (version_ .gt. 6.0) then
          hpres = dmu0*pres(1:ns_cob)  ! in VMEC versions > 6.00 pressure is given in pascals
       else
         hpres = pres(1:ns_cob)        ! in VMEC versions  .LE. 6.00 pressure is p= mu_o * p(pascals)
       endif
       hphip = phip(1:ns_cob) 

       mercierf = Dmerc(1:ns_cob)

       r0 = (r0max_v+r0min_v)/2
       amin = r0/aspect_v
       beta0 = betaxis_v
       if (beta0 .le. zero) then
          beta0 = epsilon(beta0)
          hpres = beta0*( /(1 - (j - 1)/real(ns - 1,rprec), j=1,ns)/ )
       end if 
       b0_v = sqrt((2._dp/beta0)*(1.5_dp*hpres(2)-.5_dp*hpres(3)))

       allocate (xm_v(mnmax_v), xn_v(mnmax_v), plist(ns_cob), stat=k)
       if (k .ne. 0) stop 'Allocation error 2 in cobra order_input'
       xm_v = xm(1:mnmax_v)
       xn_v = xn(1:mnmax_v)              

       mndim = ns_cob*mnmax_v

!...   RMN, ZMN are on FULL-MESH; LMNS, BSUPUMN (c), BSUPCMN (c), 
!      and BMNC are on HALF-MESH!!!

       allocate (rmncf(mndim), zmnsf(mndim), lmnsh(mndim), 
     1    bmnch(mndim),  bsupumnh(mndim), bsupvmnh(mndim), 
     2    stat=k)
       if (k .ne. 0) stop 'Allocation error 3 in cobra order_input' 

       zmnsf=0; rmncf=0; lmnsh=0; bmnch=0; bsupvmnh=0; bsupumnh=0

       do k=1, ns_cob
         lk = (k-1)*mnmax_v
         rmncf(lk+1:lk+mnmax_v) = rmnc(1:mnmax_v,k)
         zmnsf(lk+1:lk+mnmax_v) = zmns(1:mnmax_v,k)
         bmnch(lk+1:lk+mnmax_v) = bmn(1:mnmax_v,k)
         lmnsh(lk+1:lk+mnmax_v) = lmns(1:mnmax_v,k)
         bsupvmnh(lk+1:lk+mnmax_v) = bsupvmn(1:mnmax_v,k)
         bsupumnh(lk+1:lk+mnmax_v) = bsupumn(1:mnmax_v,k)
       enddo

!...   identify wanted surfaces

       allocate (lballoon(ns_cob), stat=k)
       lballoon=.false.
       do i=1,nl
         lballoon(bsurf(i))=.true.
       enddo

!...   check for feasibility of ballooning calculation
 
       nlist=0
       do i = 2, ns_cob - 1                 ! exclude Boundary (i=ns)
          if (lballoon(i)) then
             nlist=nlist+1
             plist(nlist)=i              ! i-th surface on full is between i-th and (i+1)-th on half
          endif
       enddo

!...   store NLIST surfaces where calculation can be done in LIST(1:NLIST) vector
    
       allocate (list(nlist), stat=k)
       list(1:nlist) = plist(1:nlist)
       deallocate (plist,lballoon)
1000   call read_wout_deallocate
       if (ierr .eq. 0) then
          ierr = -ierr_vmec
       else
          print *,' ierr = ', ierr, ' writing in READ_WOUT_BOOZ'
       end if  

       end subroutine order_input
 

       subroutine get_ballooning_grate(grate)
       use kind_spec
       use ballooning_data
       use normalize_data
       use readin_data
       use general_dimensions
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       real(rprec), intent(OUT), dimension(ns_cob) :: grate 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec), dimension(ns_cob) ::  radii_flux, h0, xmax, xmin
       real(rprec) :: fullp, semip, h, eigf, err, rerr, ohs2, ohs, pi
       integer :: np, npm1, npm2, npm3, inc
       integer :: nturns, j, kl, i, k, l
       integer :: dinf, dinf1, dinfp, dinfp1, dinfpp, dinfm1
       real(rprec),allocatable, dimension(:) :: xf, xh, pf, ph,
     1    qf, qh, rf, rh, feigfun, eigfun, i_xf, i_xh, i_pf, i_ph,
     2    i_qf, i_qh, i_rf, i_rh, hv, eigv, a1f, a2f, a3f
       real(rprec) :: feigenv
       logical :: l_success
!-----------------------------------------------
 
      grate = 0
      do i=1,ns_cob
!        radii(i)=sqrt(real(i-1,rprec)/(ns_cob-1)     !  r/a= (radii_flux)**1/2
         radii_flux(i) = real(i-1,rprec)/(ns_cob-1)   !  s = radii_flux
      enddo

!=================
!      BEGIN MIGRATION TO FULL MESH
!=================

!...  Store surface quantites on RADIAL full mesh

      allocate(iotaf(ns_cob), presf(ns_cob), phipf(ns_cob), stat=k)
      if (k .ne. 0) stop 'Allocation error 1 in get_ballooning_grate'

      iotaf = 0 ; presf = 0; phipf = 0

      iotaf(2:ns_cob-1) = 0.5_dp*(hiota(2:ns_cob-1) + hiota(3:ns_cob))
      presf(2:ns_cob-1) = 0.5_dp*(hpres(2:ns_cob-1) + hpres(3:ns_cob))
      phipf(2:ns_cob-1) = 0.5_dp*(hphip(2:ns_cob-1) + hphip(3:ns_cob))

!...  Evaluate and store surface quantities derivatives on RADIAL full mesh
       
      allocate(iotapf(ns_cob), prespf(ns_cob), stat=k)
      if (k .ne. 0) stop 'Allocation error 2 in get_ballooning_grate'

      iotapf = 0; prespf = 0

      ohs = ns_cob-1                            ! ds to differentiate on radial HALF mesh
      ohs2 = ohs/2                           ! 2 comes from differentiating on radial FULL mesh
      
      iotapf(2:ns_cob-1) = ohs*(hiota(3:ns_cob) - hiota(2:ns_cob-1))
      prespf(2:ns_cob-1) = ohs*(hpres(3:ns_cob) - hpres(2:ns_cob-1))


!...  EVALUATE AND STORE VMEC Fourier coefficients and
!           their derivatives on RADIAL full mesh

      allocate (lmnsf(mndim), bmncf(mndim), rmncpf(mndim), 
     1          zmnspf(mndim), lmnspf(mndim), bmncpf(mndim),
     2          bsupvmnf(mndim), bsupumnf(mndim), stat=k) 

      if (k .ne. 0) stop 'Allocation error 3 in get_ballooning_grate'
 
      lmnsf=0; bmncf=0; bsupvmnf=0; bsupumnf=0
      rmncpf=0; zmnspf=0; lmnspf=0; bmncpf=0
       
      do kl = 1, nlist
         dinf  = (list(kl)-1)*mnmax_v
         dinf1 = dinf+1
         dinfp = dinf+mnmax_v
         dinfp1= dinfp+1
         dinfpp= dinf+2*mnmax_v
         dinfm1= dinf-mnmax_v+1

!...   VMEC Fourier coefficients on RADIAL full mesh

         lmnsf(dinf1:dinfp) = 0.5_dp*(lmnsh(dinfp1:dinfpp)
     1    +lmnsh(dinf1:dinfp))
         bmncf(dinf1:dinfp) = 0.5_dp*(bmnch(dinfp1:dinfpp)
     1    +bmnch(dinf1:dinfp))
         bsupvmnf(dinf1:dinfp) = 0.5_dp*(bsupvmnh(dinfp1:dinfpp)
     1    +bsupvmnh(dinf1:dinfp))
         bsupumnf(dinf1:dinfp) = 0.5_dp*(bsupumnh(dinfp1:dinfpp)
     1    +bsupumnh(dinf1:dinfp))

!...   VMEC Fourier coefficients radial derivatives on RADIAL full mesh      
       
         rmncpf(dinf1:dinfp) = ohs2*(rmncf(dinfp1:dinfpp)
     1    -rmncf(dinfm1:dinf))
         zmnspf(dinf1:dinfp) = ohs2*(zmnsf(dinfp1:dinfpp)
     1    -zmnsf(dinfm1:dinf))
         lmnspf(dinf1:dinfp) = ohs*(lmnsh(dinfp1:dinfpp)
     1    -lmnsh(dinf1:dinfp))
         bmncpf(dinf1:dinfp) = ohs*(bmnch(dinfp1:dinfpp)
     1    -bmnch(dinf1:dinfp))
       
       enddo

!...   Check for stellarator symmetry
    
       pi = 4*atan(one)
       xmax(1) = 0
       xmax(ns_cob) = 0
       xmax(2:ns_cob-1)=pi*(2*k_w-1)/(2*nfp_v*iotaf(2:ns_cob-1))   ! account for K_w potential wells!!!
       fullp = 360._dp/nfp_v
       semip = fullp/2
       nturns = int(init_zeta/fullp)
       init_zeta = init_zeta-fullp*nturns
       if ((init_zeta.eq.zero .or. init_zeta.eq.semip .or.
     1   init_zeta.eq.fullp) .and. init_theta.eq.zero) tsymm = 0
       if (tsymm .eq. 0) then
         xmin = 0
         np0  = 6*np0_in/12+1       ! always odd + 3-multiple  number of points; (1/2) in symmetric case
         h0   = xmax/(np0-1)
         inc  = 1                ! solution vector includes F(0)
       else
         np0  = 6*(np0_in/6)+1      ! always odd + 3-multiple number of points
         xmin = -xmax
         h0   = 2*(xmax/(np0-1))
         inc  = 0                ! solution vector does not include F(0) since F(0)=0.
       endif

!================
!      END MIGRATION TO FULL MESH
!================
       if (lscreen) then
          write(*,48)
          write(*,49)
          write(*,48)
       end if

!=============
!     BEGIN BALLOONING LOOP
!=============

       ballooning: do kl=1,nlist
               
         allocate (hv(lmax+1), eigv(lmax), stat=k)
         if (k .ne. 0) stop 'Allocation 4 error in get_ballooning_grate'

         hv(1) = 1
         h = h0(list(kl))             ! initialize step size
         np = np0                     ! initialize number of points
         npm1 = np-1
         npm2 = np-2
         npm3 = np-3
         l_success =.false. 

!-----------------------------
!      BEGIN RICHARDSON'S LOOP
!-----------------------------

         richardson: do j=1,lmax    ! lmax= max. #refinements in Richardson's scheme

           if (j .gt. 1) then       ! reduce step ALWAYS by 2, so
             h = h/2                ! new full mesh  =  previous half+full meshes
             np = np*2-1            ! np =  total number of points in current iteration
             npm1 = np-1
             npm2 = np-2
             npm3 = np-3
           endif

           allocate(a1f(npm2+inc), a2f(npm3+inc), a3f(npm3+inc),     !a1= diag; a2=sup-diag; a3=sub-diag 
     1      xf(np), xh(npm1), pf(np), ph(npm1), qf(np), qh(npm1),    !xf=LINE full mesh; xh=LINE half mesh
     2      rf(np), rh(npm1), eigfun(np),feigfun(npm2+inc), stat=k)
           if (k .ne. 0)
     1      stop 'Allocation error 5 in get_ballooning_grate'
     
           if(j.gt.1)then                   
              xf(1:np:2) = i_xf             
              pf(1:np:2) = i_pf
              qf(1:np:2) = i_qf
              rf(1:np:2) = i_rf                             ! if refinement, use previous....
              xf(2:npm1:2) = i_xh                           ! ... evaluations to form new LINE full mesh
              pf(2:npm1:2) = i_ph
              qf(2:npm1:2) = i_qh
              rf(2:npm1:2) = i_rh
              deallocate(i_xf, i_pf, i_qf, i_rf, i_xh, i_ph, 
     1          i_qh, i_rh, stat=k)                         ! deallocate passing vectors
           endif
      
           call getmatrix(list(kl), a1f, a2f, a3f, h, np, j, xf,
     1          xh, pf, ph, qf, qh, rf, rh, inc, xmin(list(kl)))
           call geteigm(a1f, a2f, a3f, npm2+inc, feigenv, feigfun) ! get eigenvalue and eigenvector: 2-nd order
          
           if (tsymm .eq. 0) then
             eigfun(1:np-1) = feigfun(1:np-1)
           else
             eigfun(1) = 0
             eigfun(2:np-1) = feigfun(1:np-2)
           endif
           eigfun(np) = 0
           call variat_eig_full(np, h, eigfun, pf, qf, rf, eigv(j)) !variational eigenvalue: 4th-order accurate

           deallocate (a1f, a2f, a3f, feigfun, stat=k)                   ! deallocate matrix memory space
         
           if (j .ge. krich) then
              call polint(hv(j-km), eigv(j-km), krich, zero, eigf, 
     1          err)       ! carry out extrapolation
              rerr = abs(err/eigf)
              if(rerr .lt. tole) then
                l_success =.true. 
                if (eigf.lt.0) grate(list(kl)) = sqrt(-eigf)             ! if conveged, exit richardson loop
                if (eigf.ge.0) grate(list(kl)) =-sqrt(eigf)
                exit richardson
              endif
           endif

           hv(j+1) = 0.0625_dp*hv(j)  ! 2nd. order discretization =>  variat.eigenv. error is also 4th. order

           allocate (i_xf(np), i_xh(npm1), i_pf(np), i_ph(npm1),
     1        i_qf(np), i_qh(npm1), i_rf(np), i_rh(npm1), stat=k)        ! allocate passing vectors
           if (k .ne. 0)
     1       stop 'Allocation error 6 in get_ballooning_grate'

           i_xf = xf 
           i_xh = xh 
           i_pf = pf 
           i_ph = ph                                                     ! store LINE mesh info in passing vectors
           i_qf = qf 
           i_qh = qh 
           i_rf = rf
           i_rh = rh

           deallocate(xf, xh, pf, ph, qf, qh, rf, rh, eigfun, stat=k)    ! deallocate old LINE mesh vectors
         
         end do richardson

!--------------------------
!     END RICHARDSON'S LOOP
!--------------------------

         deallocate (hv, eigv, stat=k)                                  ! deallocate interpolation variables
         if (allocated(i_xf))
     1      deallocate(i_xf,i_xh,i_pf,i_ph,i_qf,i_qh,i_rf,i_rh, stat=k)
         if (allocated(xf)) 
     1      deallocate (xf, xh, pf, ph, qf, qh, rf, rh, eigfun, stat=k)
        
         if (lscreen) write(*, 50) list(kl), radii_flux(list(kl)), 
     1    init_zeta, init_theta, grate(list(kl)), j, np, xmax(list(kl)),
     2    l_success, tsymm, presf(list(kl)), mercierf(list(kl))


       end do ballooning 

!============
!     END BALLOONING LOOP
!============

       if (lscreen) write(*, 48)

49     format(3x,'NS',5x,'FLUX-s',5x,'ZT_0',8x,'TH_0',7x,
     1  'GR. RATE',4x,'IT',3x,'POINTS',6x,'XMAX',3x,'OK?',
     2  2x,'SYMM',4x,'PRES',8x,'MERC')
48     format(112('-'))
50     format(i5, 2x, 1pe10.2, 2(1pe10.2,2x), 1pe12.4, 3x, i2,
     1   2x, i6, 3x, 1pe10.2, 2x, l1, 3x, i2, 2(2x,1pe10.2))

!...  deallocate all variables

      deallocate (lmnsf, bmncf, rmncpf, zmnspf, 
     1  lmnspf, bmncpf, iotapf, prespf, iotaf,
     2  presf, phipf, bsupvmnf, bsupumnf, stat=k)  

      end subroutine get_ballooning_grate
      
      
       subroutine getmatrix(nsurf, ad, asup, asub, h, n, ik, xf, 
     1  xh, pf, ph, qf, qh, rf, rh, inc, xmin)
       use kind_spec
       use ballooning_data
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer, intent(in) :: nsurf, n, ik, inc
       real(rprec), intent(in):: h, xmin
       real(rprec), intent(out):: ad(n-2+inc), asup(n-3+inc),
     1  asub(n-2+inc)
       real(rprec),intent(inout):: xf(n), xh(n-1), pf(n), ph(n-1),
     1  qf(n), qh(n-1), rf(n), rh(n-1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec) :: h2
       integer :: k,i
!------------------------------------------------------------------------

       h2=h*h
       if(ik.eq.1)then                              !   if first call then.....
         xf=(/(xmin+(k-1)*h,k=1,n)/)                ! .... build full mesh and ...
         call coeffs(nsurf, n, xf, pf, qf, rf)      ! ...evaluate coefficients on full mesh
       endif
 
       xh(1:n-1)=xf(2:n)-h/2                        ! build half mesh and ...
       call coeffs(nsurf, n-1, xh, ph, qh, rh)      ! ...evaluate coefficients on half mesh
      
         
       if(tsymm .eq. 0)then

         ad(1)=(-2*ph(1)+qf(1)*h2)/(rf(1)*h2)       !  B.C.: F(1)=F(-1) , i.e. F'(0)=0
         asup(1)=2*ph(1)/(rf(1)*h2)                 !  size of matrix is then NPTS-1  

         ad(2:n-1)=(-ph(2:n-1)-ph(1:n-2)+qf(2:n-1)*h2)
     1            /(rf(2:n-1)*h2)                                  ! build A-diagonal
         asup(2:n-2)=ph(2:n-2)/(rf(2:n-2)*h2)                      ! build A-supdiagonal
         asub(1:n-2)=ph(1:n-2)/(rf(2:n-1)*h2)                      ! build A-subdiagonal

       else                                                        !  size of matrix is NPTS-2
 
         ad(1:n-2)=(-ph(2:n-1)-ph(1:n-2)+qf(2:n-1)*h2)
     1            /(rf(2:n-1)*h2)                                  ! build A-diagonal
         asup(1:n-3)=ph(2:n-2)/(rf(2:n-2)*h2)                      ! build A-supdiagonal
         asub(1:n-3)=ph(2:n-2)/(rf(3:n-1)*h2)                      ! build A-subdiagonal

       endif
        
       end subroutine getmatrix          


       subroutine coeffs(nsurf, npt, x, p, q, r)
       use kind_spec
       use normalize_data
       use ballooning_data
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer,intent(in) :: nsurf,npt
       real(rprec),intent(in), dimension (npt) :: x
       real(rprec),intent(out), dimension(npt) :: p, q, r
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec), dimension(npt) :: k_s, c2, b2, bsuppar,
     2   c2n, b2n, bsupparn
       real(rprec) :: presn, prespn, phipn2, eps2,
     1  factor, phipc, iotac, prespc, iotapc 
!------------------------------------------------------------------------
         integer::j

!...   values of surface quantities at current surface
   
       phipc=phipf(nsurf)                            ! toroidal magnetic flux
       iotac=iotaf(nsurf)                            ! iota
       iotapc=iotapf(nsurf)                          ! radial iota derivative
       prespc=prespf(nsurf)                          ! radial pressure gradient
      
       call summodosd(nsurf, npt, x, c2, k_s, bsuppar, b2,
     1   iotac, iotapc, prespc)
      
       b2n = b2/(b0_v**2)                              ! magnetic field normalized to B_0
       bsupparn = r0*bsuppar/b0_v
       c2n = c2*(amin**2)                              ! perpend. wave vector squared
       prespn = 2*prespc/(beta0*b0_v**2)               ! pressure normalized to p_0
       phipn2=(phipc/(b0_v*amin**2))**2                ! perpend. lengths normalized to a
       eps2=(r0/amin)**2
       factor=eps2*beta0*prespn/phipn2

!...   P, Q and R:  ballooning coefficients

       p=bsupparn*c2n/b2n                                          
       q=factor*k_s/bsupparn
       r=c2n/(b2n*bsupparn)                         

       end subroutine coeffs  
      
      
       subroutine summodosd(nsurf, npt, x, c2, k_s, bsuppar, 
     1            b2, iotac, iotapc, prespc)
       use kind_spec
       use ballooning_data
       use general_dimensions
       use fmesh_quantities
       use summod
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer,intent(in) :: nsurf, npt
       real(rprec),intent(in), dimension (npt) :: x
       real(rprec),intent(in):: iotac, prespc, iotapc 
       real(rprec),intent(out), dimension(npt) :: c2,
     1   k_s, b2, bsuppar
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: j, lj, istat
       real(rprec) :: zetacn, thetacn, twopi, alpha
!------------------------------------------------------------------------
       call alloc_summod(npt)

       twopi   = 8*atan(one)
       zetacn  = (twopi*init_zeta)/360
       thetacn = (twopi*init_theta)/360
       call obtain_field_line(zetacn, thetacn, iotac, nsurf,
     1   alpha)
       zetang  = zetacn+x
       call obtain_theta(alpha, iotac, thetacn, nsurf, npt, 
     1   zetang, thetang, lambdath)

!=============
!   BEGIN FOURIER INVERSION
!=============

!...   initialize for Fourier inversion
      
       bfield = 0; bfieldze = 0; bfields = 0; bfieldth = 0
       rboo = 0; rs = 0; rze = 0; rth = 0; zs = 0; zze = 0
       zth = 0; lambdaze = 0; lambdas = 0; bsupze = 0; bsupth = 0

       fourier: do j = 1,mnmax_v             ! Fourier invert 

         arg = xm_v(j)*thetang-zetang*xn_v(j)
         ccosi = cos(arg)
         ssine = sin(arg)
         lj = mnmax_v*(nsurf-1)+j
         bfield = bfield+bmncf(lj)*ccosi                !magnetic field magnitude
         bfields = bfields+bmncpf(lj)*ccosi             ! ..... radial derivative
         bfieldze = bfieldze+xn_v(j)*bmncf(lj)*ssine    ! ..... zeta derivative
         bfieldth = bfieldth-xm_v(j)*bmncf(lj)*ssine    ! ..... theta derivative
         rboo = rboo+rmncf(lj)*ccosi                    ! cylindrical R
         rth = rth-rmncf(lj)*xm_v(j)*ssine              ! ..... radial derivative
         rze = rze+rmncf(lj)*xn_v(j)*ssine              ! ..... zeta derivative
         rs = rs+rmncpf(lj)*ccosi                       ! ..... theta derivative
         zth = zth+zmnsf(lj)*xm_v(j)*ccosi              ! cylindrical Z: theta derivative
         zze = zze-zmnsf(lj)*xn_v(j)*ccosi              ! ..... zeta derivative
         zs = zs+zmnspf(lj)*ssine                       ! ..... radial derivative
         lambdas = lambdas+lmnspf(lj)*ssine             ! lambda radial derivative
         lambdaze = lambdaze-xn_v(j)*lmnsf(lj)*ccosi    ! ..... zeta derivative
         bsupth= bsupth+ bsupumnf(lj)*ccosi             ! contravariant theta-comp. magnetic field
         bsupze= bsupze+ bsupvmnf(lj)*ccosi             ! contravariant zeta-comp. magnetic field

       enddo fourier
       
!============
!   END FOURIER INVERSION
!============

!...   auxiliar quantities
     
       rboo2 = rboo**2
       bfieldi = one/bfield
       b2 = bfield**2
       bfield2i = bfieldi**2
 
!...   VMEC lower metric elements

       gtssub = rth*rs+zth*zs
       gstsub = gtssub
       gzssub = rze*rs+zze*zs
       gszsub = gzssub                                   
       gtzsub = rth*rze+zth*zze
       gztsub = gtzsub
       gttsub = (rth**2)+(zth**2)
       gsssub = (rs**2)+(zs**2)
       gzzsub = (rze**2)+(zze**2)+rboo2

!...    VMEC jacobian

       jacob2 = gsssub*gttsub*gzzsub +
     1    2*gstsub*gtzsub*gszsub -
     2    gttsub*gszsub**2 -gzzsub*gstsub**2-
     3    gsssub*gztsub**2
       rjac2i = one/jacob2

!...   VMEC upper metric elements
 
       gsssup = (gttsub*gzzsub-gtzsub*gztsub)*rjac2i
       gttsup = (gsssub*gzzsub-gszsub*gzssub)*rjac2i
       gzzsup = (gsssub*gttsub-gstsub*gtssub)*rjac2i  
       gstsup = (gztsub*gzssub-gstsub*gzzsub)*rjac2i       
       gtssup = gstsup
       gszsup = (gtssub*gtzsub-gttsub*gzssub)*rjac2i
       gzssup = gszsup
       gtzsup = (gstsub*gzssub-gsssub*gztsub)*rjac2i
       gztsup = gtzsup
       
!...   covariant B-components

       bsubs=  bsupze*gszsub+bsupth*gstsub
       bsubth=  bsupze*gtzsub+bsupth*gttsub
       bsubze=  bsupze*gzzsub+bsupth*gtzsub
       
!...   VMEC covariant curvature components
       
       aux = (bsupze*bfieldze+bsupth*bfieldth)*bfield2i
       cks_vmec = bfieldi*(bfields + bfieldi*prespc - bsubs*aux)
       ckth_vmec = bfieldi*(bfieldth - bsubth*aux)


!...   normal curvature in (s,alpha,phi)-coordinates!

       lam1 = 1 + lambdath
       lam2 = - iotac + lambdaze      
       lam3 = - iotapc*(x-x0)+ lambdas
       k_s = cks_vmec - ckth_vmec*lam3/lam1

!...   normal vector squared

       c2 = gzzsup*lam2**2 + gttsup*lam1**2 + 
     1   gsssup*lam3**2 + 2*lam2*lam3*gzssup +
     1   2*lam3*lam1*gtssup + 2*lam1*lam2*gztsup

!...   Contravariant zeta-component of field going to B.grad!!
 
       bsuppar = bsupze 

       call free_summod

       end subroutine summodosd


       subroutine obtain_field_line(zetacn, thetacn, iotac, 
     1    nsurf, alpha)
       use kind_spec
       use ballooning_data
       use general_dimensions
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer, intent(IN) :: nsurf
       real(rprec), intent(IN) :: zetacn, thetacn, iotac
       real(rprec), intent(OUT) :: alpha
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: j, lj
       real(rprec):: lambda0, ssine, arg
!-----------------------------------------------
       
       lambda0 = 0
       fourier: do j = 1, mnmax_v                       ! Fourier invert lambda at initial point (thetacn, zetacn)
         arg = xm_v(j)*thetacn-xn_v(j)*zetacn
         ssine = sin(arg)
         lj = mnmax_v*(nsurf-1)+j
         lambda0 = lambda0 + lmnsf(lj)*ssine     
       enddo fourier
       
       alpha = thetacn + lambda0 - iotac*zetacn         ! obtain field line label value
       
       end subroutine obtain_field_line


       subroutine obtain_theta(alpha, iotac, thetacn, nsurf, npt, 
     1   zetang, thetang, lambdath)
       use kind_spec
       use ballooning_data
       use general_dimensions
       use fmesh_quantities
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer, intent(IN):: nsurf, npt
       real(rprec), intent(IN) :: iotac, alpha, thetacn
       real(rprec), intent(IN), dimension(npt):: zetang
       real(rprec), intent(OUT), dimension(npt):: thetang,
     1   lambdath
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec), parameter :: guard_band = 0.95_dp
       integer :: j, lj, l, kj       
       real(rprec):: lambda0, dlambda0, ccosi, ssine, arg, 
     1   aux1, theta0, dtheta, fun0, dfun0
!-----------------------------------------------
     
       theta_loop: do l = 1, npt
         aux1 = iotac*zetang(l) + alpha         
         if(l == 1)  then
           theta0 = thetacn
         else
          theta0 = thetang(l-1)
         endif

! 
!        Use Newton-Raphson to find angle theta that is the zero of 
!        F(theta) = theta + lambda(theta) - (alpha + iota*zeta)
!
         newton_raphson: do kj = 1, jnewton_max             
           lambda0 = 0
           dlambda0 = 0

           fourier: do j = 1, mnmax_v                        

             arg = xm_v(j)*theta0-xn_v(j)*zetang(l)
             ssine = sin(arg)
             ccosi = cos(arg)
             lj = mnmax_v*(nsurf-1)+j
             lambda0 = lambda0 + lmnsf(lj)*ssine            ! Fourier invert lambda and d(lambda)/dtheta
             dlambda0 = dlambda0 + xm_v(j)*lmnsf(lj)*ccosi

           enddo fourier
       
           fun0 = theta0 - aux1 + lambda0
           dfun0 = 1 + dlambda0

           if (kj .ge. jnewton_max/2) dfun0 = max(dfun0, guard_band)    ! Added by SPH (10/00)

           dtheta = fun0/dfun0
           theta0 = theta0 - dtheta

           if (abs(dtheta) < newton_tol) then
             thetang(l) = theta0
             lambdath(l) = dlambda0
             exit
           else if (kj == jnewton_max) then
             stop 'COBRA: NEWTON fails!'
           endif

         enddo newton_raphson       

       enddo theta_loop

       end subroutine obtain_theta
       
             
       subroutine variat_eig_full(n, h, eigfun, pf, qf, rf,
     1   eigenv)
       use kind_spec
       use ballooning_data
       implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer,intent(IN) :: n
       real(rprec), intent(IN) :: h
       real(rprec), dimension(n), intent(IN):: eigfun, qf, 
     1  rf, pf
       real(rprec), intent(OUT) :: eigenv       
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: istat
       real(rprec), dimension(:), allocatable :: intr, intq,
     1   intp, coef, dery
       real(rprec) :: frac,frac1
!-----------------------------------------------
       allocate (intr(n), intq(n), intp(n), coef(n), dery(n), 
     1           stat=istat)
       if (istat .ne. 0) stop 'Allocation error in variat_eig_full'
   
       frac = 2._dp/(3*h)
       frac1= 1._dp/(12*h)

       if(tsymm.eq.0)then
         dery(1) = 0
         dery(2)=frac*(eigfun(3)-eigfun(1))
     1                 -frac1*(eigfun(4)-eigfun(2))
       else
         dery(1)=frac*eigfun(2)-frac1*eigfun(3)
         dery(2)=frac*eigfun(3)-frac1*eigfun(4)
       endif
       dery(3:n-2)=frac*(eigfun(4:n-1)-eigfun(2:n-3))
     1  -frac1*(eigfun(5:n)-eigfun(1:n-4))
       dery(n-1)=-frac*eigfun(n-2)-frac1*eigfun(n-3)
       dery(n)=-frac*eigfun(n-1)-frac1*eigfun(n-2)

       intp=pf*(dery**2)
       intr=rf*(eigfun**2)
       intq=qf*(eigfun**2)

       coef=1
       coef(1)=0.354166666667_dp       ! use 4-th oder alternative extended Simpson rule
       coef(2)=1.229166666667_dp       ! Numerical Recipes, Chapter 4, page 108 
       coef(3)=0.895833333333_dp
       coef(4)=1.020833333333_dp
       coef(n:n-3:-1)=coef(1:4)
!      coef(1)=.5_dp                   ! 2-nd order extended trapezoidal rule can be chosen
!      coef(n)=.5_dp                   ! commenting the previous lines and uncommenting these

       intp=coef*intp
       intr=coef*intr
       intq=coef*intq

       eigenv=(SUM(intp)-(SUM(intq)))/SUM(intr)

       deallocate (intr, intq, intp, coef, dery)
       
       end subroutine variat_eig_full
       

       subroutine polint(xa, ya, n, x, y, dy)
       use kind_spec
       implicit none                              
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer,intent(IN) :: n
       real(rprec),intent(IN), dimension(n) :: xa, ya
       real(rprec),intent(IN) :: x
       real(rprec),intent(OUT) :: y, dy
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
       integer,parameter :: nmax=10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(rprec), dimension(nmax) :: c, d  
       real(rprec) ::  dif, dift, h0, hp, w, den
       integer :: ns, i, j, m
!-----------------------------------------------
!      Given arrays XA and YA of length N, and given X,
!      returns Y, and an error estimate DY.
!      If P(x) is the polynomial of degree N-1 such that 
!      P(AXi)=AYi, i=1,...,N, then the returned Y=P(X).
!      Interpolation is done using Neville's algorithm.
!      (NUMERICAL RECIPES, p. 82).
!
       ns=1
       dif=abs(x-xa(1))
       do i=1,n
       dift=abs(x-xa(i))
         if(dift.lt.dif)then
           ns=i
           dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
       enddo
       y=ya(ns)
       ns=ns-1
       do m=1,n-1
         do i=1,n-m
           h0=xa(i)-x
           hp=xa(i+m)-x
           w=c(i+1)-d(i)
           den=h0-hp
           if (den .eq. 0._dp)then
             write(*,*)'Two identical xa in input'
             return
           endif
           den=w/den
           d(i)=hp*den
           c(i)=h0*den
         enddo
         if(2*ns.lt.n-m)then
           dy=c(ns+1)
         else
           dy=d(ns)
           ns=ns-1
         endif
         y=y+dy
       enddo
       
       end subroutine polint
      

       subroutine geteigm(ad, asup, asub, n, eigm, eigf)
       use kind_spec
       use ballooning_data
       implicit none   
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       integer, intent(IN) :: n
       real(rprec), intent(IN), dimension(n) :: ad
       real(rprec), intent(IN), dimension(n-1) ::asup, asub
       real(rprec), intent(OUT) :: eigm
       real(rprec), intent(OUT),dimension(n) :: eigf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: i, j, nm
       real(rprec), dimension(:), allocatable :: w, w2
       real(rprec):: maxm, minm, norm
!-------------------------------------------------------------------------------------
       allocate (w(n), w2(4*n), stat = j)
       if (j .ne. 0) stop 'Allocation error in COBRA geteigm'
       
       call TVAL(eigm, -kth, asub, ad, asup, n, w)
       call TVECT(eigm, eigf, asub, ad, asup, n, w2)

       deallocate (w, w2)
       
       maxm=maxval(eigf)
       minm=minval(eigf)
       norm=max(abs(maxm),abs(minm))
      
       if(norm.eq.abs(maxm))norm=maxm
       if(norm.eq.abs(minm))norm=minm
       eigf=eigf/norm                    ! normalize eigenfunction to unity at maximum
       eigm=-eigm                        ! convert to Sturm-Liouville eigenvalue
       
       end subroutine geteigm
