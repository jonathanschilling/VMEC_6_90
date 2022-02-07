cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
       module normalize_data
       use kind_spec
       implicit none
       integer :: nfp_v
       real(rprec), parameter :: mu_0=1.2566368e-6_dp
       real(rprec) :: r0,b0_v,amin,beta0, twopi
       end module normalize_data
 
       module general_dimensions
       integer :: ns_cob, mnmax_v, mndim, nlist 
       end module general_dimensions

       module readin_data
       use kind_spec
       use general_dimensions
       integer,dimension(:),allocatable :: list   
       real(rprec), dimension(:), allocatable :: hiota, hphip,
     1   hpres
       real(rprec), dimension(:), allocatable :: xn_v, xm_v
       real(rprec), dimension(:), allocatable :: lmnsh,bmnch,
     1   bsupvmnh, bsupumnh
       logical :: lscreen
       end module readin_data
       
       module ballooning_data
       use kind_spec
       implicit none
       integer :: tsymm = 1
       real(rprec), allocatable, dimension(:) :: 
     1    init_zeta_v, init_theta_v
       real(rprec) :: init_zeta, init_theta
       real(rprec), parameter :: tole = 1.e-3_dp, x0 = 0
       integer, parameter :: np0_in = 101, kth = 1, k_w = 7
       integer,parameter :: lmax=8, krich=3, km=krich-1
       integer, parameter :: jnewton_max=200
       integer :: np0
       real(rprec), parameter :: newton_tol = tole
       
!      lmax=#max iterations in Richardson's; krich= min.#samples to interpolate
!      x0=radial wave number (default=0.)
!      k_w: number of potential wells included in integration domain!
!      tole= tolerance for eigenvalue fit error
!      np0_in=coarsest grid used in Richardson's; kth-largest eigenvalue obtained
!      jnewton_max= maximum number of iterations in Newton-Raphson to get theta
!      newton_tol= tolerance in Newton_Raphson

       end module ballooning_data
       
       module fmesh_quantities
       use kind_spec
       use general_dimensions
       use readin_data, ONLY: xn_v, xm_v
       implicit none
       real(rprec), dimension(:), allocatable:: iotaf, phipf,
     1   presf, mercierf
       real(rprec), dimension(:), allocatable:: iotapf, prespf
       real(rprec), dimension(:), allocatable:: rmncf, zmnsf,
     1   lmnsf, bmncf, bsupvmnf, bsupumnf
       real(rprec), dimension(:), allocatable:: rmncpf, zmnspf,
     1   lmnspf, bmncpf
       end module fmesh_quantities

       module summod
       use kind_spec
       implicit none 
              
       real(rprec), dimension(:), allocatable :: bfields, bfieldze,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, aux, 
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, lam2, cks_vmec, ckth_vmec, bfieldi, bfield2i,
     7  lam3

       contains

       subroutine alloc_summod(npt)
       integer :: npt, istat

       allocate(bfields(npt), bfieldze(npt), bfieldth(npt), rboo(npt),
     1  rs(npt), rze(npt), rth(npt), zs(npt), zze(npt), zth(npt), 
     2  ccosi(npt), ssine(npt), lam1(npt), zetang(npt), thetang(npt), 
     3  arg(npt), lambdaze(npt), lambdath(npt), bsubze(npt), bsubs(npt),
     4  bsubth(npt), lambdas(npt), rboo2(npt),jacob2(npt),bfield2i(npt), 
     5  rjac2i(npt), bfield(npt), bsupth(npt), bsupze(npt), gtssub(npt),
     6  gstsub(npt), gzssub(npt), gszsub(npt), gtzsub(npt), gztsub(npt),
     7  gttsub(npt), gzzsub(npt), gsssub(npt), gttsup(npt), gzzsup(npt),
     8  gtssup(npt), gstsup(npt), gzssup(npt), gszsup(npt), gtzsup(npt),
     9  gztsup(npt), gsssup(npt), lam2(npt),  bfieldi(npt), aux(npt),
     A  lam3(npt), cks_vmec(npt), ckth_vmec(npt), stat= istat)

       if (istat .ne. 0) stop 'Allocation error in COBRA summod...'

       end subroutine alloc_summod

       subroutine free_summod
       integer :: istat
       
       deallocate( bfields, bfieldze, bfieldi, bfield2i, aux, lam3,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, lam2,
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, cks_vmec, ckth_vmec, stat=istat)

       end subroutine free_summod

       end module summod
EOF

cat > cobra_source.f << "EOF"
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
EOF

cat > cobra_lib.f <<"EOF"
!
! ===================================
! NIST Guide to Available Math Software.
! Fullsource for module TVAL from package NAPACK.
! Retrieved from NETLIB on Fri Aug 7 17:41:59 1998.
! ===================================
!
      SUBROUTINE TVAL(E, K, L, D, U, N, W) 

!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE (REPLACE K |
!     |                 BY -K TO GET K-TH LARGEST EIGENVALUE)  |
!     |                                                        |
!     |         L     --SUBDIAGONAL (CAN BE IDENTIFIED WITH U) |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST N)         |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX                        |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP,STM                 |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: K, N 
      REAL(rprec):: E 
      REAL(rprec), DIMENSION(N) :: L, D, U, W 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, M 
      REAL(rprec) :: S, T, Y, Z 
!-----------------------------------------------
      IF (N > 1) THEN 
!     -------------------------------------------------------
!     |*** BOUND EIGENVALUES USING GERSCHGORIN'S THEOREM ***|
!     -------------------------------------------------------
         M = N - 1 
         Y = D(1) 
         Z = D(1) 
         S = 0 
         DO I = 1, M 
            W(I) = L(I)*U(I) 
            IF (W(I) < ZERO) GO TO 50 
            S = S + ABS(U(I)) 
            T = D(I) + S 
            Z = MAX(T,Z) 
            T = D(I) - S 
            Y = MIN(T,Y) 
            S = ABS(L(I)) 
         END DO 
         T = D(N) + S 
         Z = MAX(T,Z) 
         T = D(N) - S 
         Y = MIN(T,Y) 
         M = K 
         IF (K .LE. 0) M = N + K + 1 
         T = ONE 
         T = T/2 
         S = ONE + T 
 1003    CONTINUE 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         IF (S .LE. ONE) GO TO 1002 
         T = T/2 
         S = ONE + T 
         GO TO 1003 
 1002    CONTINUE 
         T = 8*N*T*MAX(ABS(Y),ABS(Z)) 
         IF (T .NE. ZERO) THEN 
            CALL STM (E, Y, Z, M, T, D, W, N) 
            RETURN  
         ENDIF 
      ENDIF 
      E = D(1) 
      RETURN  

   50 CONTINUE 
      WRITE (6, *) 'ERROR: SUBROUTINE TVAL CAN ONLY BE APPLIED' 
      WRITE (6, *) 'WHEN THE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE' 
      STOP  

      END SUBROUTINE TVAL 


      SUBROUTINE STM(E, Y, Z, K, T, D, P, N) 
!      ________________________________________________________
!     |                                                        |
!     |   FIND THE K-TH SMALLEST EIGENVALUE OF A TRIDIAGONAL   |
!     |  MATRIX WHOSE CROSS-DIAGONAL PRODUCTS ARE NONNEGATIVE  |
!     |USING BOTH THE BISECTION METHOD AND A NEWTON-LIKE METHOD|
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         Y,Z   --ENDS OF AN INTERVAL THAT CONTAINS THE  |
!     |                 DESIRED EIGENVALUE                     |
!     |                                                        |
!     |         K     --INDEX OF DESIRED EIGENVALUE            |
!     |                                                        |
!     |         T     --TOLERANCE (ITERATIONS CONTINUE UNTIL   |
!     |                 THE ERROR IN THE EIGENVALUE .LE. T)    |
!     |                                                        |
!     |         D     --DIAGONAL OF THE COEFFICIENT MATRIX A   |
!     |                                                        |
!     |         P     --CROSS-DIAGONAL PRODUCTS A SUB I+1,I    |
!     |                 TIMES A SUB I,I+1                      |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,MAX,SIGN,SQRT              |
!     |    PACKAGE SUBROUTINES: CP,EQL,INP                     |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer:: K, N 
      REAL(rprec):: E, Y, Z, T 
      REAL(rprec), DIMENSION(N) :: D, P 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, J, M 
      REAL(rprec) :: A,B,FL,FR,L,Q,R,S,U,V,
     1   W,D2,D3,D4,P2,P3,P4,D34,D42,D23 
      REAL(rprec), DIMENSION(4) :: F, H, X 
      REAL(rprec), DIMENSION(4) :: G 
      REAL(rprec) :: C, GL, GR 
!-----------------------------------------------
      L = Y
      R = Z
      IF ( L .LT. R ) GOTO 10
      L = Z
      R = Y
10    IF ( N .EQ. 1 ) GOTO 210
      S = T + T
      U = ONE
20    U = U/2
      A = ONE + U
      IF (A .GT. ONE) GOTO 20
      U = 5*U
      V = U/2
      CALL CP(L,FL,GL,I,D,P,N)
      CALL CP(R,FR,GR,J,D,P,N)
      IF ( I .GE. K ) GOTO 190
      IF ( J .LT. K ) GOTO 200
C     --------------------------------
C     |*** ISOLATE THE EIGENVALUE ***|
C     --------------------------------
30    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
      IF ( J .EQ. I+1 ) GOTO 70
40    A = .5*(L+R)
      CALL CP(A,B,C,M,D,P,N)
      IF ( K .LE. M ) GOTO 50
      L = A
      I = M
      FL = B
      GL = C
      GOTO 30
50    R = A
      J = M
      FR = B
      GR = C
      GOTO 30
60    E = R - L
      IF ( E .LE. S ) GOTO 180
      IF ( E .LE. U*(ABS(L)+ABS(R)) ) GOTO 180
70    X(1) = L
      F(1) = FL
      G(1) = GL
      X(2) = R
      F(2) = FR
      G(2) = GR
      CALL EQL(X,F,G,H,2)
      IF ( H(1) .EQ. H(2) ) GOTO 160
C     ---------------------
C     |*** SECANT STEP ***|
C     ---------------------
      A = X(1) - H(1)*(X(1)-X(2))/(H(1)-H(2))
      Q = A
      W = max(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) A = L + W
      IF ( ABS(A-R) .LT. W ) A = R - W
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 80
      R = A
      FR = B
      GR = C
      GOTO 90
80    L = A
      FL = B
      GL = C
90    X(3) = A
      F(3) = B
      G(3) = C
      W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      CALL EQL(X,F,G,H,3)
C     --------------------------------------
C     |*** QUADRATIC INTERPOLATION STEP ***|
C     --------------------------------------
      CALL INP(A,X(1),X(2),X(3),H(1),H(2),H(3),L,R)
      B = L
      IF ( ABS(A-L) .GT. ABS(A-R) ) B = R
C     ------------------------------------
C     |*** APPLY PSEUDO-NEWTON METHOD ***|
C     ------------------------------------
100   Q = A
      W = max(T,V*(ABS(L)+ABS(R)))
      IF ( ABS(A-L) .LT. W ) GOTO 110
      IF ( ABS(A-R) .GT. W ) GOTO 130
110   IF ( A+A .GT. L+R ) GOTO 120
      A = L + W
      GOTO 130
120   A = R - W
130   IF ( A .LE. L ) GOTO 160
      IF ( A .GE. R ) GOTO 160
      E = .5*E
      IF ( E .LT. ABS(B-A) ) GOTO 160
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 140
      R = A
      FR = B
      GR = C
      GOTO 150
140   L = A
      FL = B
      GL = C
150   W = R - L
      IF ( W .LE. S ) GOTO 220
      IF ( W .LE. U*(ABS(L)+ABS(R)) ) GOTO 220
      X(4) = A
      F(4) = B
      G(4) = C
      CALL EQL(X,F,G,H,4)
      IF ( X(1) .LT. L ) GOTO 160
      IF ( X(1) .GT. R ) GOTO 160
      B = X(1)
      D4 = X(4) - B
      D3 = X(3) - B
      D2 = X(2) - B
      D34 = X(3) - X(4)
      D42 = X(4) - X(2)
      D23 = X(2) - X(3)
      P2 = D2*(ONE+((D2/D3)*D42+(D2/D4)*D23)/D34)
      P3 = D3*(ONE+((D3/D2)*D34+(D3/D4)*D23)/D42)
      P4 = D4*(ONE+((D4/D2)*D34+(D4/D3)*D42)/D23)
      IF (P2 .NE. ZERO) P2 = (H(2)-H(1))/P2
      IF (P3 .NE. ZERO) P3 = (H(3)-H(1))/P3
      IF (P4 .NE. ZERO) P4 = (H(4)-H(1))/P4
      P2 = P2 + P3 + P4
      IF (P2 .EQ. ZERO) GOTO 160
      A = B - H(1)/P2
      GOTO 100
C     --------------------------
C     |*** BISECTION METHOD ***|
C     --------------------------
160   A = (L+R)/2
      CALL CP(A,B,C,J,D,P,N)
      IF ( I .GE. J ) GOTO 170
      R = A
      FR = B
      GR = C
      GOTO 60
170   L = A
      FL = B
      GL = C
      GOTO 60
180   E = (L+R)/2
      RETURN
190   E = L
      RETURN
200   E = R
      RETURN
210   E = D(1)
      IF ( L .GT. E ) E = L
      IF ( R .LT. E ) E = R
      RETURN
220   E = Q
      RETURN
      END SUBROUTINE STM


      SUBROUTINE INP(A, X, Y, Z, U, V, W, L, R) 
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec) :: A, X, Y, Z, U, V, W, L, R 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: B, C, P, Q, S, T 
!-----------------------------------------------
      S = Z - X 
      T = (Y - X)/S 
      A = (U - V)/T + (W - V)/(1. - T) 
      IF (A .NE. ZERO) THEN 
         B = .5_dp*(W - U)/A - .5_dp 
         C = U/A 
         T = SQRT(ABS(C)) 
         IF (ABS(B) < SIGN(T,C)) GO TO 60 
         T = MAX(T,ABS(B)) 
         IF (T .EQ. ZERO) GO TO 50 
         Q = ONE/T 
         P = SQRT((Q*B)**2 - Q*C*Q) 
         P = T*P 
         IF (ABS(P + B) .LE. ABS(P - B)) THEN 
            Q = P - B 
         ELSE 
            Q = -(B + P) 
         ENDIF 
         P = C/Q 
         Q = X + S*Q 
         P = X + S*P 
         IF (Q .GE. L) THEN 
            IF (Q .LE. R) THEN 
               A = Q 
               RETURN  
            ENDIF 
         ENDIF 
         A = P 
         RETURN  
      ENDIF 
      IF (U .NE. W) THEN 
         A = X + S*U/(U - W) 
         RETURN  
      ENDIF 
   50 CONTINUE 
      A = L 
      RETURN  
   60 CONTINUE 
      A = X - S*B 
      RETURN  

      END SUBROUTINE INP 


      SUBROUTINE CP(T, B, V, L, D, P, N) 
!------------------------------------------------
!*** EVALUATE CHARACTERISTIC POLYNOMIAL AND ***
!              SIGN ALTERNATIONS               
!------------------------------------------------
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: L, N 
      REAL(rprec) ::  T, B 
      REAL(rprec) :: V 
      REAL(rprec), DIMENSION(N) :: D, P 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I, J, M 
      REAL(rprec) :: A, C, E, F, U, Z 
!-----------------------------------------------
      L = 0 
      Z = 0 
      V = Z 
      IF (N .LE. 1) THEN 
         B = D(1) - T 
         IF (B .LE. Z) L = 1 
         RETURN  
      ENDIF 
      F = 65536
      F = F**4
      E = ONE/F 
      U = 16 
      M = 0 
      I = 1 
   20 CONTINUE 
      C = ONE 
      B = D(I) - T 
   30 CONTINUE 
      IF (B .LE. Z) GO TO 70 
      IF (I .GE. N) GO TO 130 
   40 CONTINUE 
      J = I 
      I = I + 1 
      A = (D(I)-T)*B - P(J)*C 
      C = B 
      B = A 
   50 CONTINUE 
      A = ABS(B) 
      IF (A > F) GO TO 60 
      IF (A > E) GO TO 30 
      IF (A .EQ. Z) GO TO 70 
      C = C*F 
      B = B*F 
      V = V - U 
      GO TO 50 
   60 CONTINUE 
      C = C*E 
      B = B*E 
      V = V + U 
      GO TO 50 
   70 CONTINUE 
      L = L + 1 
      IF (I .GE. N) GO TO 130 
      IF (B < Z) GO TO 90 
      IF (P(I) > ZERO) GO TO 90 
      I = I + 1 
      M = 1 
      V = Z 
      GO TO 20 
   80 CONTINUE 
      IF (B .GE. Z) GO TO 120 
      IF (I .GE. N) GO TO 130 
   90 CONTINUE 
      J = I 
      I = I + 1 
      A = (D(I)-T)*B - P(J)*C 
      C = B 
      B = A 
  100 CONTINUE 
      A = ABS(B) 
      IF (A > F) GO TO 110 
      IF (A > E) GO TO 80 
      IF (A .EQ. Z) GO TO 120 
      C = C*F 
      B = B*F 
      V = V - U 
      GO TO 100 
  110 CONTINUE 
      C = C*E 
      B = B*E 
      V = V + U 
      GO TO 100 
  120 CONTINUE 
      L = L + 1 
      IF (I .GE. N) GO TO 130 
      IF (B > Z) GO TO 40 
      IF (P(I) > ZERO) GO TO 40 
      I = I + 1 
      M = 1 
      V = Z 
      GO TO 20 
  130 CONTINUE 
      IF (M .NE. 1) THEN 
         IF (B .NE. ZERO) THEN 
            A = ONE/U 
            IF (ABS(B) .GE. A) THEN 
  140          CONTINUE 
               IF (ABS(B) < ONE) RETURN  
               B = B*A 
               V = V + ONE 
               GO TO 140 
            ENDIF 
            B = B*U 
            V = V - ONE
 1003       CONTINUE 
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            IF (ABS(B) .GE. A) GO TO 1002 
            B = B*U 
            V = V - ONE
            GO TO 1003 
 1002       CONTINUE 
            RETURN  
         ENDIF 
      ENDIF 
      V = Z 
      B = Z 
      RETURN  

      END SUBROUTINE CP 


      SUBROUTINE EQL(X, F, G, H, N) 
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer N 
      REAL(rprec), DIMENSION(N) :: X, F, H 
      REAL(rprec), DIMENSION(N) :: G 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: I1, J2, I, J, K, L, M, NM1 
      REAL(rprec) :: R, S, Z, U 
!-----------------------------------------------
      Z = 0 
      I = 1 
      R = X(N) 
      S = F(N) 
      U = G(N) 
      IF (S .EQ. Z) GO TO 30 
      I1 = I 
      J2 = MAX(N - 1,I1) 
      DO I = I1, J2 
         IF (F(I) .EQ. Z) CYCLE  
         IF (U < G(I)) GO TO 30 
         IF (U > G(I)) CYCLE  
         IF (ABS(F(I)) .GE. ABS(S)) GO TO 30 
      END DO 
      GO TO 50 
   30 CONTINUE 
      M = N + I 
      NM1 = N - 1 
      DO J = I, NM1 
         K = M - J 
         L = K - 1 
         X(K) = X(L) 
         F(K) = F(L) 
         G(K) = G(L) 
      END DO 
      X(I) = R 
      F(I) = S 
      G(I) = U 
   50 CONTINUE 
      U = G(N) 
      WHERE ((G(:N)-U).LE.(-99._dp) .OR. F(:N).EQ.Z) H(:N) = Z 
      WHERE (F(:N).NE.Z .AND. (G(:N)-U)>(-99._dp))
     1    H(:N) = F(:N)*16._dp**(G(:N)-U) 

      RETURN  
      END SUBROUTINE EQL 


      SUBROUTINE TVECT(E, X, L, D, U, N, W) 
!      ________________________________________________________
!     |                                                        |
!     |     COMPUTE EIGENVECTOR CORRESPONDING TO GIVEN REAL    |
!     |        EIGENVALUE FOR A REAL TRIDIAGONAL MATRIX        |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |         E     --EIGENVALUE                             |
!     |                                                        |
!     |         L     --SUBDIAGONAL                            |
!     |                                                        |
!     |         D     --DIAGONAL                               |
!     |                                                        |
!     |         U     --SUPERDIAGONAL                          |
!     |                                                        |
!     |         N     --MATRIX DIMENSION                       |
!     |                                                        |
!     |         W     --WORK ARRAY (LENGTH AT LEAST 4N)        |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |         E     --IMPROVED ESTIMATE FOR EIGENVALUE       |
!     |                                                        |
!     |         X     --EIGENVECTOR                            |
!     |                                                        |
!     |    BUILTIN FUNCTIONS: ABS,SQRT                         |
!     |________________________________________________________|
!
      use kind_spec
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer:: N 
      REAL(rprec):: E 
      REAL(rprec), DIMENSION(N) :: X, L, D, U
      REAL(rprec), DIMENSION(4*N) :: W
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: G1, J1, J2, J4, J3, F, G, H, I, J, K, M 
      REAL(rprec) :: O, P, Q, R, S, T, V, Y, Z 
!-----------------------------------------------
      IF (N .LE. 1) THEN 
         E = D(1) 
         X(1) = ONE 
         RETURN  
      ENDIF 
      M = N - 1 
      J = 2 
!     ---------------------------------------------------------
!     |*** STORE MATRIX IN W ARRAY AND SUBTRACT EIGENVALUE ***|
!     ---------------------------------------------------------
      X(1:M) = 0
      DO I = 1,M
         W(J) = U(I)
         W(J+1) = D(I) - E
         W(J+2) = L(I)
         J = J + 4
      END DO   
      W(J) = 0 
      W(J+1) = D(N) - E 
      O = 65536
      O = O**(-4) 
      T = O/2
      S = T 
      T = T/2
      P = S 
      S = S + T 
      IF (S .GE. O) GO TO 40 
      DO WHILE(S + T > S) 
         T = T/2
         P = S 
         S = S + T 
         IF (S .GE. O) GO TO 40 
      END DO 
   40 CONTINUE 
      R = W(3) 
      V = ABS(R) + ABS(W(4)) 
      F = 4 
      M = 4*N - 3 
      G = -2 
!     ---------------------------
!     |*** FACTOR THE MATRIX ***|
!     ---------------------------
      G = G + 4 
      G1 = G 
      DO G = G1, M, 4 
         H = G - 1 
         I = G + 2 
         J = G + 5 
!     --------------------
!     |*** FIND PIVOT ***|
!     --------------------
         Q = W(I) 
         Y = ABS(Q) 
         Z = ABS(R) 
         IF (Z < Y) THEN 
!     -------------------
!     |*** SWAP ROWS ***|
!     -------------------
            IF (V < Y) THEN 
               V = Y 
               F = I 
            ENDIF 
            T = W(G) 
            W(G) = W(J) 
            W(J) = T 
            T = R/Q 
            K = G + 1 
            W(K) = Q 
            K = J - 1 
            S = W(K) 
            IF (S .NE. ZERO) THEN 
               IF (S .EQ. O) S = P 
               W(K) = -S*T 
               W(H) = S 
               GO TO 100 
            ENDIF 
            W(K) = S 
            W(H) = O 
            GO TO 100 
         ENDIF 
         W(H) = 0 
         IF (V < Z) THEN 
            V = Z 
            F = I 
         ENDIF 
         IF (R .EQ. ZERO) GO TO 120 
         T = Q/R 
!     -------------------
!     |*** ELIMINATE ***|
!     -------------------
  100    CONTINUE 
         R = W(J) - T*W(G) 
         W(J) = R 
         W(I) = T 
      END DO 
      IF (ABS(R) < V) THEN 
         V = R 
         F = J + 1 
!     ---------------------------------------------------
!     |*** COMPUTE INITIAL EIGENVECTOR APPROXIMATION ***|
!     ---------------------------------------------------
      ENDIF 
  120 CONTINUE 
      J = F/4 
      X(J) = ONE 
      IF (J .NE. 1) THEN 
         K = F - 5 
         J = J - 1 
         X(J) = (X(J)-W(K-1)*X(J+1))/W(K) 
         J1 = J 
         IF (J1 - 1 > 0) THEN 
            DO J = J1, 2, -1 
               K = K - 4 
               T = W(K-2) 
               IF (T .EQ. O) T = 0 
               X(J-1) = (X(J-1)-W(K-1)*X(J)-T*X(J+1))/W(K) 
            END DO 
         ENDIF 
         GO TO 140 
      ENDIF 
  140 CONTINUE 
      IF (V .NE. ZERO) THEN 
         S = MAXVAL(ABS(X(:N)))
         S = ONE/S 
         X(:N) = S*X(:N)
         R = SUM(X(:N)*X(:N))
         K = 0 
         J = 1 
         Y = X(1) 
!     -----------------------------------------------------
!     |*** APPLY ONE ITERATION OF INVERSE POWER METHOD ***|
!     -----------------------------------------------------
         J2 = J 
         J4 = MAX(N - 1,J2) 
         DO J = J2, J4 
            K = K + 4 
            I = J 
            S = W(K) 
            W(K) = Y 
            Y = X(J+1) 
            IF (W(K-3) .NE. ZERO) THEN 
               T = X(J+1) 
               X(J+1) = X(I) 
               X(I) = T 
            ENDIF 
            X(J+1) = X(J+1) - S*X(I) 
         END DO 
!     ---------------------------
!     |*** BACK SUBSTITUTION ***|
!     ---------------------------
         S = X(J)/W(K+3) 
         X(J) = S 
         T = ABS(S) 
         V = S*Y 
         J = J - 1 
         K = K - 1 
         S = (X(J)-W(K-1)*S)/W(K) 
         X(J) = S 
         T = MAX(ABS(S),T) 
         V = V + S*W(K+1) 
         J3 = J 
         IF (J3 - 1 > 0) THEN 
            DO J = J3, 2, -1 
               K = K - 4 
               Z = W(K-2) 
               IF (Z .EQ. O) Z = 0 
               S = (X(J-1)-W(K-1)*S-Z*X(J+1))/W(K) 
               X(J-1) = S 
               T = MAX(ABS(S),T) 
               V = V + S*W(K+1) 
            END DO 
         ENDIF 
         IF (V .NE. ZERO) V = R/V 
         T = ONE/T 
         S = SUM((X(:N)*V)**2) 
         Z = SUM((T*X(:N))**2) 
         T = T/SQRT(Z) 
         X(:N) = T*X(:N) 
!     --------------------------------------------------------------
!     |*** USE RAYLEIGH QUOTIENT TO IMPROVE EIGENVALUE ESTIMATE ***|
!     --------------------------------------------------------------
         IF (R + R .GE. S) E = E + V 
         RETURN  
      ENDIF 
      T = MAXVAL(ABS(X(:N)))
      T = ONE/T 
      Z = SUM((T*X(:N))**2) 
      T = T/SQRT(Z) 
      X(:N) = T*X(:N) 

      END SUBROUTINE TVECT 
EOF
EOC
