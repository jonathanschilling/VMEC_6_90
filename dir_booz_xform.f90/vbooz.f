      program boozer_xform
      use booz_params
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, jrad, i, iunit, iread, numargs
      integer, dimension(200) :: jlist
      real(rprec) :: t1, t2
      character*38, parameter :: version = 
     1   "Boozer Transformation Code Version 1.0"
      character*(120) :: arg1, arg2
      character :: extension*120
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the boozer coordinates are required
!     writes the boozer coordinates to a file, boozmn.extension
!           
!     call this as follows:   
!
!          xbooz_xform input.boz [T or F]
!
!     where input.boz contains the mboz, nboz, wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     mboz   nboz
!     FILE_EXTENSION
!     1  3   5   10  12
!
!     The optional command line argument, (T) or (F), allows the user to turn off screen 
!     output if set to F.
!
!     call xbooz_xform -h brings up a help screen
!
      lscreen = .true.          !!Default, write to screen
      jlist = 0
!    
!     Read command line argument to get input file name
!      
      call getcarg(1, arg1, numargs)
      if (numargs .gt. 1) call getcarg(2, arg2, numargs)

      if (numargs .lt. 1) then
         stop 'Invalid command line in calling xbooz_xform'
      else if (arg1 .eq. '-h' .or. arg1 .eq. '/h') then
         print *,' ENTER INPUT FILE NAME ON COMMAND LINE'
         print *,' For example: xbooz_xform in_booz.qos'
         print *
         print *,' where in_booz.qos is the input file'
         print *
         print *,' Optional command line argument'
         print *,' xbooz_xform <infile> (T or F)'
         print *
         print *,' where F suppresses output to the screen'
         stop
      else if (numargs .gt. 1) then
         if (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.   
      endif

      iread = unit_booz-1
      call safe_open (iread, istat, trim(arg1), 'old', 'formatted')
      if (istat .ne. 0) stop 'Error opening input file in booz_xform'

      read (iread, *, iostat = istat) mboz, nboz
      read (iread, *, iostat = istat) extension
      if (istat .ne. 0) stop 'Error reading input file in booz_xform'

      iunit = iread+1
      call safe_open (iunit, istat, 'boozmn.' // extension, 'replace',
     1     'unformatted')
      if (istat .ne. 0) then
         print *,' istat = ', istat
         stop 'Error opening boozmn file in XBOOZ_XFORM!'
      end if  

!
!     READ IN PARAMETERS, DATA FROM WOUT FILE 
!
      call read_wout_booz(extension, iunit, istat)
      if (istat .ne. 0) then
         print *,' ierr_vmec !=0 in booz_xform read_wout_booz'
         goto 1010
      end if   

      write (iunit, iostat=istat, err=1010) mboz, nboz, mnboz
      write (iunit, iostat=istat, err=1010) version
      
      read (iread, *, iostat=istat) jlist
      close (unit=iread)
      call second0(t1)
            
      do i = 1, 200
        jrad = jlist(i)
        if (jrad.le.0 .or. jrad.gt.ns) cycle
        write (iunit, iostat=istat, err=1010) jrad
        call boozer_coords(jrad, iunit)
      end do

 1010 continue
      
!
!     FREE MEMORY : USER MUST CALL BOOZER_COORDS AT END WITH LDEALLOC = TRUE
!     OTHERWISE MEMORY ALLOCATED WILL NOT BE AVAILABLE
!
 
      call free_mem_boozer
      close(unit=iunit, iostat=istat)
      
      call second0(t2)
      
      if (lscreen) print 120, t2-t1
 120  format(/,' TIME IN BOOZER TRANSFORM CODE:',1pe12.2,' SEC')
       
      end program boozer_xform
      
      
      subroutine boozer_coords(jrad, iunit)
      use booz_params
      use booz_extern
      use booz_persistent
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer jrad, iunit
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nparity, mn, istat1, nv2_b, mb, nb
      integer, save :: icount_save = 0
      real(rprec) ::  b1, bmodv(4), bmodb(4), err(4), jacfac
      real(rprec), dimension(:), allocatable ::
     1   r1, rt, rz, z1, zt, zz, rodd, rtodd, rzodd, zodd, ztodd, zzodd, 
     2   gsqrt, r12, z12, p1, q1, xjac, lt, lz, lam, wt, wz, wp
      real(rprec), dimension(:,:), allocatable ::
     1   cosmm, sinmm, cosnn, sinnn
C-----------------------------------------------
c       jrad         radial point where Boozer coords. are needed
c       ns           number of vmec radial grid points
c       nu_boz       number of theta points in integration
c       nv_boz       number of zeta points in integration
c       mpol         number of theta harmonics from vmec
c       ntor         number of zeta harmonics from vmec (no. zeta modes = 2*ntor+1)
c       mboz         number of boozer theta harmonics desired
c       nboz         number of boozer zeta harmonics desired
c

 
      if (icount_save .eq. 0) then
         icount_save = icount_save + 1
!
!        ALLOCATE GLOBAL ARRAYS
!

         call setup_booz (ixm, ntorsum, nunv, ns, mnmax, ohs, nfp, xmb, 
     1      xnb, sfull, scl, mboz, nboz, mnboz, nu2_b, nv_boz)

!
!        SET UP FIXED ANGLE ARRAYS
!

         nunv = nu2_b*nv_boz               !!only need on top half of theta mesh
         call foranl (cosm_b, cosn_b, sinm_b, sinn_b, thgrd, ztgrd,
     1      nu2_b, nv_boz, mpol1, ntor, nfp, nunv)

         if (lscreen) write(*,50) mboz-1, -nboz, nboz, nu_boz, nv_boz
  50     format('  0 <= mboz <= ',i4,3x,i4,' <= nboz <= ',i4,/,
     1          '  nu_boz = ',i5,' nv_boz = ',i5,//,
     1         13x,'OUTBOARD (u=0)',11x,'JS',10x,'INBOARD (u=pi)'
     2         /,70('-')/,'  v     |B|vmec   |B|booz   Error',12x,
     3         '|B|vmec   |B|booz   Error'/)

      endif

!
!      ALLOCATE LOCAL MEMORY 
! 
      allocate( r12(nunv), z12(nunv), r1(nunv*2), rt(nunv*2), 
     1   rodd(nunv*2), rtodd(nunv*2), rzodd(nunv*2), rz(nunv*2), 
     2   z1(nunv*2), zt(nunv*2), zodd(nunv*2), ztodd(nunv*2), 
     3   zzodd(nunv*2), zz(nunv*2), gsqrt(nunv), lt(nunv), lz(nunv),
     4   p1(nunv), q1(nunv), xjac(nunv), stat=istat1 )
      if (istat1 .ne. 0) stop 'memory allocation error in boozer_coords'

!
!     COMPUTE FOURIER COEFFICIENTS (pmn) OF THE
!     BOOZER-TO-VMEC TRANSFORMATION FUNCTION P:
!
!     Theta-Booz = Theta-VMEC + Lambda + Iota*p
!     Zeta-Booz  = Zeta-VMEC  + p
!
      call transpmn (pmn, xm, xn, bsubtmn(1,jrad), 
     1 bsubzmn(1,jrad), gpsi, ipsi, mnmax, ixm, ixn, hiota, jrad)

!
!     BEGIN CALCULATION OF BOOZER QUANTITIES AT HALF-RADIAL
!     MESH POINT JRAD
!     (ALL TRANSFORMED QUANTITIES MUST BE ON HALF-MESH FOR ACCURACY)
!
!
!     COMPUTE EVEN (in poloidal mode number) AND ODD COMPONENTS
!     OF R,Z, LAMDA AND P(transformation function) IN REAL SPACE
!
      nparity = 0
 
      allocate (lam(nunv), wt(nunv), wz(nunv), wp(nunv), stat=istat1)
      if (istat1 .ne. 0) stop 'allocation error in boozer_coord'
      
      call vcoords (rmnc, zmns, lmns, pmn, xm, xn, ixm, 
     1   ixn, ntorsum, jrad, ns, mnmax, r1, rt, rz, z1, zt, zz,
     2   lt, lz, lam, wt, wz, wp, cosm_b, sinm_b, cosn_b, sinn_b, 
     3   sfull, nparity, nunv, mpol1, ntor, nfp)
 
      nparity = 1
 
      call vcoords (rmnc, zmns, lmns, pmn, xm, xn, ixm, 
     1   ixn, ntorsum, jrad, ns, mnmax, rodd, rtodd, rzodd, zodd, 
     2   ztodd, zzodd, lt, lz, lam, wt, wz, wp, cosm_b, sinm_b, 
     3   cosn_b, sinn_b, sfull, nparity, nunv, mpol1, ntor, nfp)
 

!
!     COMPUTE MAPPING FUNCTION P AND MAPPING JACOBIAN (XJAC)
!     FOR DOING FOURIER INTEGRALS IN VMEC COORDINATES
!
      call harfun (jacfac, hiota, gpsi, ipsi, jrad, nunv, 
     1   lt, lz, lam, wt, wz, wp, p1, q1, xjac)
      deallocate (lam, wt, wz, wp, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in boozer_coord'

!
!     COMPUTE JACOBIAN (gsqrt) BETWEEN VMEC FLUX COORDINATES
!     AND CYLIDRICAL COORDINATES (NEED TO COMPUTE CO- AND
!     CONTRA-VARIANT COMPONENTS OF MAGNETIC FIELD)
!
      call booz_jac (r1, z1, rt, zt, rodd, zodd, rtodd, ztodd, gsqrt, 
     1   r12, z12, ohs, jrad, nunv)

!
!     COMPUTE MOD-B IN VMEC COORDINATES
!
      call bmetrc (hiota, phip, lz, lt, gsqrt, bmod_b, ohs, r1, rt, rz,
     1   zt, zz, rodd, rtodd, rzodd, ztodd, zzodd, nunv, jrad)
!     Store VMEC-Space fixed point values for |B| for checking accuracy later
      nv2_b = nv_boz/2 + 1       !Index of v=pi
      bmodv(1) = bmod_b(1,1)           !(0,0)
      bmodv(2) = bmod_b(1,nu2_b)       !(0,pi)
      bmodv(3) = bmod_b(nv2_b,1)       !(pi,0)
      bmodv(4) = bmod_b(nv2_b,nu2_b)   !(pi,pi)
 
! 
!     COMPUTE BOOZER-SPACE FOURIER COEFFICIENTS FOR R,Z,P, AND |B|
!
      deallocate (r1, rt, rodd, rtodd, rzodd, rz, 
     1   z1, zt, zodd, ztodd, zzodd, zz, gsqrt, lt, lz,
     2   stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in boozer_coord'

      allocate (cosmm(nunv,0:mboz), sinmm(nunv,0:mboz),
     1   cosnn(nunv,0:nboz), sinnn(nunv,0:nboz), stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in boozer_coord'

      call boozer (thgrd, ztgrd, bmod_b, r12, z12, xmb, xnb, 
     1   bmnb, rmnb, zmnb, pmnb, gmnb, scl, p1, q1, xjac,
     2   cosmm, sinmm, cosnn, sinnn, mnboz, nunv, mboz, nboz, 
     3   nfp, nu2_b, nv_boz, jacfac)
 
      deallocate (r12, z12, p1, q1, xjac, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in boozer_coord'

!
!     COMPUTE BOOZER-SPACE MOD-B
!
!     call modbooz(bmnb, rmnb, zmnb, pmnb,  bmod_b, xmb, xnb, thgrd, 
!    1   ztgrd, cosmm, sinmm, cosnn, sinnn, 
!    2   mnboz, nunv, mboz, nboz, nfp)

      deallocate (cosmm, sinmm, cosnn, sinnn, stat=istat1)

!     NOTE: THIS ONLY IS CORRECT FOR STELLARATOR-SYMMETRIC FIELDS
!     Summation rules: |B|(vmec) = |B|(boozer) at fixed points of map
!       SUM (Bmn)(vmec) = SUM (Bmn)(boozer) (at u=0,v=0 fixed point), etc.

      bmodb = zero
      if (icount_save .eq. 1) then
         do mn = 1,mnboz
            mb = nint(xmb(mn))
            nb = nint(xnb(mn))/nfp
            write (iunit, iostat=istat1) nfp*nb, mb
            icount_save = icount_save + 1
         end do   
      end if  

!
!     WRITE rmnb, zmnb, bmnb, pmnb, gmnb (ALL on half grid)
!
      do mn = 1, mnboz
         write(iunit, iostat=istat1) bmnb(mn), rmnb(mn), 
     1      zmnb(mn), pmnb(mn), gmnb(mn)
         if (istat1 .ne. 0) stop ' Error writing out boozmn file'
      end do     
      
      do mn = 1,mnboz
        mb = nint(xmb(mn))
        nb = nint(xnb(mn))/nfp
        b1 = bmnb(mn)
        bmodb(1) = bmodb(1) + b1
        bmodb(2) = bmodb(2) + b1 * (-1)**mb
        bmodb(3) = bmodb(3) + b1 * (-1)**nb
        bmodb(4) = bmodb(4) + b1 * (-1)**(nb - mb)
      end do
 
      err = abs(bmodb - bmodv)/max(bmodb,bmodv)
      if (lscreen) write(*,100)
     1          bmodv(1),bmodb(1),err(1),jrad,bmodv(2),bmodb(2),err(2),
     2          bmodv(3),bmodb(3),err(3),bmodv(4),bmodb(4),err(4)      
 100  format('  0  ',1p3e10.2,i5,2x,1p3e10.2,/' pi  ',2(1p3e10.2,7x) )
 
      end subroutine boozer_coords


      subroutine allocate_boozer
      use booz_params
      use booz_extern
      use booz_persistent
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1=0, istat2=0
c-----------------------------------------------
      if (.not.allocated(bsubtmn)) allocate(
     1    bsubtmn(mnmax,ns), bsubzmn(mnmax,ns), 
     2    rmnc(mnmax,ns), zmns(mnmax,ns), lmns(mnmax,ns), 
     3    xm(mnmax), xn(mnmax), ixm(mnmax), ixn(mnmax), 
     4    hiota(ns), phip(ns), gpsi(ns), ipsi(ns), pmn(mnmax), 
     5    rmnb(mnboz), zmnb(mnboz), pmnb(mnboz), gmnb(mnboz),
     6    bmnb(mnboz), bmod_b(nv_boz,nu_boz), stat=istat1 )
     
      if (.not.allocated(cosm_b)) allocate(
     1    cosm_b(nunv,0:mpol1), sinm_b(nunv,0:mpol1),
     2    cosn_b(nunv,0:ntor),  sinn_b(nunv,0:ntor),
     3    sfull(ns), scl(mnboz), xmb(mnboz), xnb(mnboz), 
     4    thgrd(nunv), ztgrd(nunv), stat=istat2 ) 
     
      if (istat1.ne.0 .or. istat2.ne.0) then
          print *,' problem allocating boozer memory'
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          stop
      endif  

      end subroutine allocate_boozer


      subroutine free_mem_boozer
      use booz_params
      use booz_extern
      use booz_persistent
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0
C-----------------------------------------------
      if (allocated(bsubtmn)) deallocate(
     1    bsubtmn, bsubzmn, rmnc, zmns, lmns, 
     2    xm, xn, ixm, ixn, hiota, phip, gpsi, ipsi, pmn, 
     3    rmnb, zmnb, pmnb, gmnb, bmnb, bmod_b, stat=istat1 )
     
      if (allocated(cosm_b)) deallocate(
     1    cosm_b, sinm_b, cosn_b, sinn_b,
     2    sfull, scl, xmb, xnb, thgrd, ztgrd, stat=istat2) 

      if (istat1 .ne.0 .or. istat2 .ne. 0)
     1  print *,' Deallocation error in Free_mem_boozer'       
 
      end subroutine free_mem_boozer


      subroutine read_wout_booz(extension, iunit, ierr)
      use read_wout_mod
      use booz_params, mpol_in => mpol, nfp_in => nfp, ns_in => ns,
     1   ntor_in => ntor, mnmax_in => mnmax, phip_in => phip      
      use booz_extern
      use booz_persistent, rmnc_in => rmnc, zmns_in => zmns, 
     1   lmns_in => lmns, xm_in => xm, xn_in => xn 
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ierr, iunit
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js
C----------------------------------------------- 

      call read_wout_file('wout.' // trim(extension), ierr)
      
      if (ierr .ne. 0) then
         print *,' ierr = ', ierr,
     1   ' error in read_wout_file called from xbooz_xform'
      end if
      if (ierr_vmec.ne.0 .or. ierr.ne.0) goto 1000
 
      mpol_in = mpol
      mpol1   = mpol - 1 
      ntor_in = ntor
      mnmax_in= mnmax
      nfp_in  = nfp
      ns_in   = ns
*
*     COMPUTE ACTUAL NO. THETA, PHI POINTS FOR INTEGRATIONS
*     NEEDED FOR DYNAMIC MEMORY ALLOCATION
*
      mboz = max(6*mpol, 2, mboz)                 !USER COULD INPUT MORE HERE
      nboz = max(2*ntor - 1, 1, nboz)             !USER COULD INPUT MORE HERE
      nu_boz   = 6*mboz+2
      nv_boz   = 4*nboz+2
      nu_boz   = nu_boz + mod(nu_boz,2)                !nu_boz, nv_boz MUST be even
      nv_boz   = nv_boz + mod(nv_boz,2)
      nunv = nu_boz*nv_boz
      mnboz = nboz + 1 + (mboz - 1)*(1 + 2*nboz)
      nu2_b = nu_boz/2 + 1
 
*
*     ALLOCATE ARRAYS FIRST TIME THRU
*
      call allocate_boozer
 
      xm_in(:mnmax) = xm(:mnmax)
      xn_in(:mnmax) = xn(:mnmax)
      ixm(:mnmax) = nint(xm(:mnmax))
      ixn(:mnmax) = nint(xn(:mnmax))

      rmnc_in(:mnmax,:ns) = rmnc(:mnmax,:ns)
      zmns_in(:mnmax,:ns) = zmns(:mnmax,:ns)
      lmns_in(:mnmax,:ns) = lmns(:mnmax,:ns)
      bsubtmn(:mnmax,:ns) = bsubumn(:mnmax,:ns)
      bsubzmn(:mnmax,:ns) = bsubvmn(:mnmax,:ns)    
      hiota(2:ns) = iotas(2:ns)
      phip_in(2:ns) = phip(2:ns)
      
!    
!     write out file surface quantities needed by ballooning code
!
      write(iunit, iostat=ierr, err=1000) 
     1   nfp, ns, aspect, rmax_surf, rmin_surf, betaxis

      do js = 2, ns
         write(iunit, iostat=ierr, err=1000) hiota(js), pres(js), 
     1   beta_vol(js), phip(js),  phi(js), bvco(js),  buco(js)
      end do

!
!     Deallocate memory in read_wout module
!
 1000 call read_wout_deallocate
      if (ierr .eq. 0) then
         ierr = -ierr_vmec
      else 
         print *,' ierr = ', ierr, ' writing in READ_WOUT_BOOZ'
      end if   
       
      end subroutine read_wout_booz


      subroutine boozer(thgrd, ztgrd, bmod, rad, zee, xm, xn,
     1   bmnb, rmnb, zmnb, pmnb, gmnb, scl, uboz, vboz, xjac, cosmm, 
     2   sinmm, cosnn, sinnn, mnmax, nznt, mboz, nboz, nfp, nu2, nv,
     3   jacfac)
C...MODIFIED 6/98 by A. WARE to speed up by factor of 8
C
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: mnmax, nznt, mboz, nboz, nfp, nu2, nv
      real(rprec), dimension(nznt), intent(in) ::
     1   thgrd, ztgrd, bmod, rad, zee, uboz, vboz, xjac
      real(rprec), dimension(mnmax), intent(in) ::
     1   xm, xn, scl
      real(rprec), dimension(nznt,0:mboz), intent(out) ::
     1   cosmm, sinmm 
      real(rprec), dimension(nznt,0:nboz), intent(out) ::
     1   cosnn, sinnn 
      real(rprec), dimension(mnmax), intent(out) ::
     1   bmnb, rmnb, zmnb, pmnb, gmnb
      real(rprec), intent(in) :: jacfac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, m, n, imax, i
      real(rprec), dimension(nznt) :: cost, sint, bbjac
      real(rprec) :: sgn
C-----------------------------------------------

      cosmm(:,0) = one
      sinmm(:,0) = zero
      cosmm(:,1) = cos(thgrd + uboz)
      sinmm(:,1) = sin(thgrd + uboz)

      cosnn(:,0) = one
      sinnn(:,0) = zero
      cosnn(:,1) = cos((ztgrd + vboz)*nfp)
      sinnn(:,1) = sin((ztgrd + vboz)*nfp)

      do m = 2,mboz
         cosmm(:,m) = cosmm(:,m-1)*cosmm(:,1)
     1              - sinmm(:,m-1)*sinmm(:,1)
         sinmm(:,m) = sinmm(:,m-1)*cosmm(:,1)
     1              + cosmm(:,m-1)*sinmm(:,1)
      end do

!     ONLY INTEGRATE IN U HALF WAY AROUND

      i = nv*(nu2 - 1)+1
      imax = i-1 + nv
      do m = 0,mboz
         cosmm(1:nv,m)    = p5*cosmm(1:nv,m)
         cosmm(i:imax,m)  = p5*cosmm(i:imax,m)
         sinmm(1:nv,m)    = p5*sinmm(1:nv,m)
         sinmm(i:imax,m)  = p5*sinmm(i:imax,m)
      end do

      do n = 2,nboz
         cosnn(:,n) = cosnn(:,n-1)*cosnn(:,1)
     1              - sinnn(:,n-1)*sinnn(:,1)
         sinnn(:,n) = sinnn(:,n-1)*cosnn(:,1)
     1              + cosnn(:,n-1)*sinnn(:,1)
      end do

!     jacobian from Euclidean to Boozer coords, with SPECIAL
!     radial variable s = (toroidal flux)/twopi (phip = 1)
      bbjac = jacfac/(bmod*bmod)   

      do mn = 1,mnmax
        m = nint(xm(mn))
        n = nint(abs(xn(mn)/nfp))
        sgn = sign(one,xn(mn))
        cost = (cosmm(:,m)*cosnn(:,n)
     1       +  sinmm(:,m)*sinnn(:,n)*sgn)*xjac
        sint = (sinmm(:,m)*cosnn(:,n)
     1       -  cosmm(:,m)*sinnn(:,n)*sgn)*xjac
        bmnb(mn) = dot_product(bmod,cost)
        rmnb(mn) = dot_product(rad ,cost)
        zmnb(mn) = dot_product(zee, sint)
        pmnb(mn) =-dot_product(vboz,sint)
        gmnb(mn) = dot_product(bbjac, cost)
      end do

      bmnb = scl*bmnb
      rmnb = scl*rmnb
      zmnb = scl*zmnb
      pmnb = scl*pmnb
      gmnb = scl*gmnb
      
      end subroutine boozer
!
!      subroutine modbooz(bmnb, rmnb, zmnb, pmnb, bmod, xm, xn,
!     1   thgrd, ztgrd, cosmm, sinmm, cosnn, sinnn,
!     2   mnmax, nznt, mboz, nboz, nfp)
!      use kind_spec
!      implicit none
!C-----------------------------------------------
!C   D u m m y   A r g u m e n t s
!C-----------------------------------------------
!      integer :: mnmax, nznt, mboz, nboz, nfp
!      real(rprec), dimension(mnmax), intent(in) ::
!     1   bmnb, rmnb, zmnb, pmnb
!      real(rprec), dimension(nznt), intent(out) :: bmod
!      real(rprec), dimension(mnmax), intent(in) :: xm, xn
!      real(rprec), dimension(nznt), intent(in) :: thgrd, ztgrd
!      real(rprec), dimension(nznt,0:mboz), intent(out) ::
!     1   cosmm, sinmm
!      real(rprec), dimension(nznt,0:nboz), intent(out) ::
!     1   cosnn, sinnn
!C-----------------------------------------------
!C   L o c a l   P a r a m e t e r s
!C-----------------------------------------------
!      real(rprec), parameter :: one = 1.0_dp, zero = 0.0_dp
!C-----------------------------------------------
!C   L o c a l   V a r i a b l e s
!C-----------------------------------------------
!      integer :: mn, m, n, i
!      real(rprec), dimension(nznt) :: cost, sint
!      real(rprec) :: sgn
!C-----------------------------------------------
!
!      bmod = zero
!
!      cosmm(:,0) = one
!      sinmm(:,0) = zero
!      cosmm(:,1) = cos(thgrd)
!      sinmm(:,1) = sin(thgrd)
!
!      cosnn(:,0) = one
!      sinnn(:,0) = zero
!      cosnn(:,1) = cos(ztgrd*nfp)
!      sinnn(:,1) = sin(ztgrd*nfp)
!
!      do m = 2,mboz
!         cosmm(:,m) = cosmm(:,m-1)*cosmm(:,1)
!     1              - sinmm(:,m-1)*sinmm(:,1)
!         sinmm(:,m) = sinmm(:,m-1)*cosmm(:,1)
!     1              + cosmm(:,m-1)*sinmm(:,1)
!      end do
!
!      do n = 2,nboz
!         cosnn(:,n) = cosnn(:,n-1)*cosnn(:,1)
!     1              - sinnn(:,n-1)*sinnn(:,1)
!         sinnn(:,n) = sinnn(:,n-1)*cosnn(:,1)
!     1              + cosnn(:,n-1)*sinnn(:,1)
!      end do
!
!      do mn=1,mnmax
!        m = nint(xm(mn))
!        n = nint(abs(xn(mn)/nfp))
!        sgn = sign(one,xn(mn))
!        cost = cosmm(:,m)*cosnn(:,n)
!     1       + sinmm(:,m)*sinnn(:,n)*sgn
!c       sint = sinmm(:,m)*cosnn(:,n)
!c    1       - cosmm(:,m)*sinnn(:,n)*sgn
!        bmod = bmod + bmnb(mn)*cost
!c       r12  = r12  + rmnb(mn)*cost
!c       z12  = z12  + zmnb(mn)*sint
!c       p12  = p12  + pmnb(mn)*sint
!      end do
!
!      end subroutine modbooz
!      end
!
 
      subroutine harfun(jacfac, hiota, gpsi, ipsi, js, nznt, xlt, 
     1   xlz, xl, wt, wz, w, uboz, vboz, xjac)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: js, nznt
      real(rprec) :: jacfac
      real(rprec), dimension(*) :: hiota, gpsi, ipsi
      real(rprec), dimension(nznt) ::
     1   xl, xlt, xlz, w, wt, wz, uboz, vboz, xjac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(nznt) ::
     1   uboz1u, vboz1u, uboz2u, vboz2u
      real(rprec) :: dem, gpsi1, hiota1, ipsi1
C-----------------------------------------------
!
!     HERE, W IS THE PART OF THE TRANSFORMATION FUNCTION P
!     WHICH DEPENDS ON THE COVARIANT B COMPONENTS. THE FULL
!     TRANSFORMATION FUNCTION IS THEN
!
!     P = ( w - Ipsi*Lambda ) / (g + iota*Ipsi)
!
!     ALSO,
!
!     uboz = lambda + iota*p   (non-secular piece of boozer theta in vmec conates)
!
!            = (gpsi*lambda + iota*w)/(gpsi + iota*Ipsi)
!
!     vboz = p                 (non-secular piece of boozer phi in vmec coortes)
!
!          = (-Ipsi*lambda + w)/(gpsi + iota*Ipsi)
!
!     FINALLY, XJAC IS THE JACOBIAN BETWEEN BOOZER, VMEC COORDINATES
!
!
!     NOTE THAT LAMBDA = UBOZ - IOTA*VBOZ
!
      jacfac = gpsi(js) + hiota(js)*ipsi(js)
      dem = one/jacfac
      gpsi1 = gpsi(js)*dem
      hiota1 = hiota(js)*dem
      ipsi1 = ipsi(js)*dem


      uboz = gpsi1*xl + hiota1*w
      vboz = dem*w - ipsi1*xl
      uboz1u(:nznt) = gpsi1*xlt + hiota1*wt
      vboz1u(:nznt) = dem*wt - ipsi1*xlt
      uboz2u(:nznt) = gpsi1*xlz + hiota1*wz
      vboz2u(:nznt) = dem*wz - ipsi1*xlz
      xjac = one + uboz1u(:nznt)*(one + vboz2u(:nznt)) 
     1     + vboz2u(:nznt) - uboz2u(:nznt)*vboz1u(:nznt)
     
      dem = minval(xjac)
      if (dem .lt. zero) print *,
     1   ' Jacobian xjac changed sign in harfun in xbooz_xform'

      end subroutine harfun

      subroutine transpmn(pmn, xm, xn, bsubtmn, bsubzmn, gpsi, ipsi, 
     1   mnmax, ixm, ixn, hiota, js)
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mnmax, js
      integer, dimension(*) :: ixm, ixn
      real(rprec), dimension(mnmax) :: pmn, bsubtmn, bsubzmn
      real(rprec), dimension(*) :: xm, xn
      real(rprec), dimension(*) :: gpsi, ipsi, hiota
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn
C-----------------------------------------------
!
!     COMPUTE THE PART OF Pmn WHICH IS INDEPENDENT OF Lmn
!       (not including the 1/i phase, where i=sqrt(-1))
!     THE NET P IS GIVEN AS FOLLOWS:
!
!     P(FINAL) = [ Pmn * sin(m*theta - n*zeta) - Ipsi * Lambda ] / D
!
!     where D = gpsi + iota*Ipsi
!
      do mn=1,mnmax
        if (ixn(mn) .ne. 0) then
          pmn(mn) = -bsubzmn(mn)/xn(mn)
        else if (ixm(mn) .ne. 0) then
          pmn(mn) = bsubtmn(mn)/xm(mn)
        else
          pmn(mn) = zero
          gpsi(js) = bsubzmn(mn)
          Ipsi(js) = bsubtmn(mn)
        endif
      end do

      end subroutine transpmn
      

      subroutine vcoords(rmnc, zmns, xlmn, pmn, xm, xn, 
     1   ixm, ixn, ntorsum, jrad, ns, mnmax, r, rt, rz, z, zt, zz, 
     2   lt, lz, lam, wt, wz, w, cosm_b, sinm_b, cosn_b, sinn_b, 
     3   sfull, nparity, nznt, mpol1, ntor, nfp)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer jrad, ns, mnmax, nparity, nznt, mpol1, ntor, nfp
      integer, dimension(mnmax) :: ixm, ixn
      integer, dimension(0:1) :: ntorsum
      real(rprec), dimension(mnmax,*) :: rmnc, zmns, xlmn
      real(rprec), dimension(mnmax) :: xm, xn, pmn
      real(rprec), dimension(nznt,2) :: r, rt, rz, z, zt, zz
      real(rprec), dimension(nznt) :: lam, lt, lz, w, wt, wz
      real(rprec), dimension(ns) :: sfull
      real(rprec), dimension(nznt,0:mpol1) :: cosm_b, sinm_b
      real(rprec), dimension(nznt,0:ntor) :: cosn_b, sinn_b
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: jlo = 1, jhi = 2
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, js1, mn, m, n
      real(rprec) :: t1, t2, rc, zs, sgn
      real(rprec), dimension(nznt) :: tsin, tcos
C-----------------------------------------------
    
      js = jrad
      js1= js-1
      if (js .le. 1) stop 'js must be > 1!'
 
      r (:,jlo:jhi) = zero; rt(:,jlo:jhi) = zero; rz(:,jlo:jhi) = zero
      z (:,jlo:jhi) = zero; zt(:,jlo:jhi) = zero; zz(:,jlo:jhi) = zero
 
 
*
*       Compute Reven, Rodd and Zeven, Zodd in Real Space
*       on FULL radial grid
*       Lambda, W (p transformation) on half grid
*       (even, nparity = 0; odd, nparity = 1)
*
      if (nparity .eq. 0) then
         t1 = one
         t2 = one
         lt  = zero
         lz  = zero
         lam = zero
         w    = zero
         wt   = zero
         wz   = zero
      else if (js .gt. 2) then
         t1 = one/sfull(js)
         t2 = one/sfull(js1)
      else
         t1 = one/sfull(2)
         t2 = one
         rmnc(1+ntorsum(0):ntorsum(1),1) = 2*rmnc(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - rmnc(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
         zmns(1+ntorsum(0):ntorsum(1),1) = 2*zmns(1+ntorsum(0):
     1      ntorsum(1),2)/sfull(2) - zmns(1+ntorsum(0):ntorsum(1),3)/
     2      sfull(3)
      endif

      do mn = 1, mnmax
         m = ixm(mn)
         n = abs(ixn(mn))/nfp
         sgn = sign(one,xn(mn))
         if (mod(m,2) .eq. nparity) then
            tcos(:nznt) = cosm_b(:,m)*cosn_b(:,n) 
     1                  + sinm_b(:,m)*sinn_b(:,n) * sgn
            tsin(:nznt) = sinm_b(:,m)*cosn_b(:,n)
     1                  - cosm_b(:,m)*sinn_b(:,n) * sgn
            rc = t1*rmnc(mn,js)
            zs = t1*zmns(mn,js)
            rt(:,jhi) = rt(:,jhi) - tsin(:nznt)*rc*xm(mn)
            rz(:,jhi) = rz(:,jhi) + tsin(:nznt)*rc*xn(mn)
            r (:,jhi) = r(:,jhi)  + tcos(:nznt)*rc
            zt(:,jhi) = zt(:,jhi) + tcos(:nznt)*zs*xm(mn)
            zz(:,jhi) = zz(:,jhi) - tcos(:nznt)*zs*xn(mn)
            z (:,jhi) = z(:,jhi)  + tsin(:nznt)*zs
            lt  = lt  + tcos(:nznt)*xlmn(mn,js)*xm(mn)
            lz  = lz  - tcos(:nznt)*xlmn(mn,js)*xn(mn)
            lam = lam + tsin(:nznt)*xlmn(mn,js)
            w   = w   + tsin(:nznt)*pmn(mn)
            wt  = wt  + tcos(:nznt)*pmn(mn)*xm(mn)
            wz  = wz  - tcos(:nznt)*pmn(mn)*xn(mn)
            rc = t2*rmnc(mn,js1)
            zs = t2*zmns(mn,js1)
            rt(:,jlo) = rt(:,jlo) - tsin(:nznt)*rc*xm(mn)
            rz(:,jlo) = rz(:,jlo) + tsin(:nznt)*rc*xn(mn)
            r (:,jlo) = r(:,jlo)  + tcos(:nznt)*rc
            zt(:,jlo) = zt(:,jlo) + tcos(:nznt)*zs*xm(mn)
            zz(:,jlo) = zz(:,jlo) - tcos(:nznt)*zs*xn(mn)
            z (:,jlo) = z(:,jlo)  + tsin(:nznt)*zs
          endif
      end do

      end subroutine vcoords

      subroutine booz_jac(r, z, rt, zt, rodd, zodd, rtodd, ztodd, 
     1   gsqrt, r12, z12, ohs, js, nznt)
      use kind_spec
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer js, nznt
      real(rprec) ohs
      real(rprec), dimension(nznt,2) :: r, z, rt, zt,
     1   rodd, zodd, rtodd, ztodd
      real(rprec), dimension(nznt) :: gsqrt, r12, z12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: hs, shalf
      real(rprec) , dimension(nznt) ::
     1   rt121u, zt121u, rs1u, zs1u
C-----------------------------------------------
c
c       (RS, ZS)=(R, Z) SUB S, (RT12, ZT12)=(R, Z) SUB THETA,
c       AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
c       IN THE SECOND ARGUMENT (J) OF FULL MESH QUANTITIES, X(L,J),
c       2 = JS'th FULL MESH POINT, 1 = JS-1'th FULL MESH POINT
c
      hs = one/ohs
      shalf = sqrt(hs*abs(js - 1 - p5))

      rt121u = p5*(rt(:,2)+rt(:,1)+shalf*(rtodd(:,2)+rtodd(:,1)))
      zt121u = p5*(zt(:,2)+zt(:,1)+shalf*(ztodd(:,2)+ztodd(:,1)))
      rs1u = ohs*(r(:,2)-r(:,1)+shalf*(rodd(:,2)-rodd(:,1)))
      zs1u = ohs*(z(:,2)-z(:,1)+shalf*(zodd(:,2)-zodd(:,1)))
      r12 = p5*(r(:,2)+r(:,1)+shalf*(rodd(:,2)+rodd(:,1)))
      z12 = p5*(z(:,2)+z(:,1)+shalf*(zodd(:,2)+zodd(:,1)))
      gsqrt = r12*(rt121u*zs1u - rs1u*zt121u + p25*(rtodd(:,2)*zodd(:,2)
     1   +rtodd(:,1)*zodd(:,1)-ztodd(:,2)*rodd(:,2)-ztodd(:,1)*rodd(:,1)
     2   +(rt(:,2)*zodd(:,2)+rt(:,1)*zodd(:,1)-zt(:,2)*rodd(:,2)-zt(:,1)
     3   *rodd(:,1))/shalf))

      end subroutine booz_jac

      subroutine bmetrc(hiota, phip, xlz, xlt, gsqrt, bmod, ohs, r, rt, 
     1   rz, zt, zz, rodd, rtodd, rzodd, ztodd, zzodd, nznt, js)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nznt, js
      real(rprec) ohs
      real(rprec), dimension(*) :: hiota, phip
      real(rprec), dimension(nznt,*) :: r, rt, rz, zt, zz, rodd,
     1   rtodd, rzodd, ztodd, zzodd
      real(rprec), dimension(nznt) :: xlz, xlt, gsqrt, bmod
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: is
      real(rprec) :: hs, shalf, s2
      real(rprec) , dimension(:), allocatable :: phip1u, bst1u,
     1   bsz1u, bsub1u, bsub2u, guu, guv, gvv
C-----------------------------------------------
c
      allocate(phip1u(nznt), bst1u(nznt), bsz1u(nznt), bsub1u(nznt), 
     1 bsub2u(nznt), guu(nznt), guv(nznt), gvv(nznt), stat=is)
      if (is .ne. 0) stop 'allocation error in bmetrc'

      hs = one/ohs
 
C     COMPUTE METRIC COEFFICIENTS ON HALF MESH
 
      guu = zero
      guv = zero
      gvv = zero
 
      shalf = sqrt(hs*(js - 1 - p5))
 

      do is = 1,2
         s2 = p5*hs*(js-3+is)
         guu = guu + p5*(rt(:,is)**2 + zt(:,is)**2)
     1         + s2     *(  rtodd(:,is)**2 + ztodd(:,is)**2)
     2         + shalf   *(rt(:,is)*rtodd(:,is) + zt(:,is)*ztodd(:,is))
         guv = guv + p5*(rt(:,is)*rz(:,is)+zt(:,is)*zz(:,is))
     1         + s2*(rtodd(:,is)*rzodd(:,is) + ztodd(:,is)*zzodd(:,is))
     2         + p5*shalf*(rt(:,is)*rzodd(:,is) + rz(:,is)*rtodd(:,is)
     3         +            zt(:,is)*zzodd(:,is) + zz(:,is)*ztodd(:,is))
         gvv = gvv + p5*(rz(:,is)**2 + zz(:,is)**2)
     1         + s2      *(rzodd(:,is)**2 + zzodd(:,is)**2)
     2         + shalf*(rz(:,is)*rzodd(:,is) + zz(:,is)*zzodd(:,is))
     3         + p5*r(:,is)**2 + s2*rodd(:,is)**2
     4         + shalf*r(:,is)*rodd(:,is)
      enddo
 
      
      phip1u = phip(js)/gsqrt
      bst1u = phip1u*(hiota(js)-xlz)
      bsz1u = phip1u*(one + xlt)
      bsub1u = guu*bst1u + guv*bsz1u
      bsub2u = guv*bst1u + gvv*bsz1u
      bmod = sqrt(bst1u*bsub1u + bsz1u*bsub2u)
       
      deallocate(phip1u, bst1u, bsz1u, bsub1u, 
     1   bsub2u, guu, guv, gvv, stat=is)

      end subroutine bmetrc


      subroutine foranl(cosm_b, cosn_b, sinm_b, sinn_b, thgrd, ztgrd,
     1   nu, nv, mpol1, ntor, nfp, nunv)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nu, nv, nfp, nunv, mpol1, ntor
      real(rprec), dimension(nunv) :: thgrd, ztgrd
      real(rprec), dimension(nunv,0:mpol1) :: cosm_b, sinm_b
      real(rprec), dimension(nunv,0:ntor) :: cosn_b, sinn_b
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lk, lt, lz, m, n
      real(rprec) :: dth, dzt, twopi
C-----------------------------------------------
 
      twopi = 8*atan(1.0_dp)
!
!     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
!
!     dth = twopi/nu                  !USE THIS FOR FULL 2-pi
      dth = twopi/(2*(nu-1))          !Half-around in theta
      dzt = twopi/(nv*nfp)
      lk = 0
 
      do lt = 1, nu
         do lz=1, nv
           lk = lk + 1
           thgrd(lk) = (lt-1)*dth
           ztgrd(lk) = (lz-1)*dzt
          end do
      end do
      
      do m = 0, mpol1
        cosm_b(:,m) = cos(m*thgrd)
        sinm_b(:,m) = sin(m*thgrd)
      end do  

      do n = 0, ntor
        cosn_b(:,n) = cos(n*nfp*ztgrd)
        sinn_b(:,n) = sin(n*nfp*ztgrd)
      end do  
          
 
      end subroutine foranl
 

      subroutine setup_booz(ixm, ntorsum, nunv, ns, mnmax, ohs, nfp, 
     1   xmb, xnb, sfull, scl, mboz, nboz, mnboz, nu2_b, nv_boz)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nunv, ns, mnmax, nfp, mboz, nboz, mnboz, nu2_b, nv_boz
      real(rprec) ohs
      integer, dimension(mnmax) :: ixm
      integer, dimension(0:1) :: ntorsum
      real(rprec), dimension(mnboz) :: xmb, xnb, scl
      real(rprec), dimension(ns) :: sfull
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mnboz0, n2, m, n1, n, i
      real(rprec) :: fac, hs
C-----------------------------------------------
!
!     SETUP BOOZER M,N ARRAYS
!
      mnboz0 = 0
      n2 = nboz
      do m = 0, mboz - 1
         n1 = -nboz
         if (m .eq. 0) n1 = 0
         do n = n1, n2
            mnboz0 = mnboz0 + 1
            if (mnboz0 .gt. mnboz) then
               stop 'mnboz exceeds limit in booz xform'
            endif
            xnb(mnboz0) = n*nfp
            xmb(mnboz0) = m
         end do
      end do
 
      if (mnboz0 .ne. mnboz) mnboz = mnboz0
 
!     fac = 2.0_dp/nunv
      fac = 2.0_dp/((nu2_b-1)*nv_boz)
      scl(:mnboz) = fac
      where (nint(abs(xnb(:mnboz))+abs(xmb(:mnboz))) .eq. 0) 
     1  scl(:mnboz) = 0.5_dp*fac

      ntorsum(0) = 0
      ntorsum(1) = 0

      do i=1,mnmax
        if( ixm(i).eq.0 ) ntorsum(0) = ntorsum(0) + 1
        if( ixm(i).le.1 ) ntorsum(1) = ntorsum(1) + 1
      end do

      ohs = (ns-1)
      hs = one/ohs
      do i=2,ns
        sfull(i) = sqrt(hs*(i-1))
      end do
 
      end subroutine setup_booz
