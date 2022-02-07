      program plotout
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use Vpname1
      use Vpname2
      use Vpname3
      use Vindat2
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: lastplot = 23
      integer, parameter :: maxlist = 100
      logical, parameter :: lwstbf = .TRUE.
      character*(*), parameter :: version = 
     1   ' VMEC Plotter Version 6.10  Jan-1999 '
#ifndef CRAY
      character*(*), parameter :: tempfile = 'QqZzXLftemp'
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: numargs, numchars, lenlist, i, 
     1   js, lj, mn, l, j, n, ii, k, 
     2   ntor0, lp, ierr
      real(rprec), dimension(:), allocatable :: szc, szb, dummy
      real(rprec) :: r0
      character*100 :: input_id, ictrans_sel, explist
      character, dimension(1:lastplot) :: pagedir*50
#ifdef CRAY
      character*8 cday,ctime
#else
      character*30 timeloc
#endif
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      integer , EXTERNAL :: iargc
C-----------------------------------------------
*
*                 THIS PROGRAM - PROUT - ACCEPTS THE OUTPUT
*                 FROM THE EQUILIBRIUM PROGRAM VMEC (WOUT FILE)
*                 AND PRINTS/PLOTS THE APPROPRIATE DATA
*
*                GET FILE ID EXTENSION FROM COMMAND LINE
*
*                 Modified (Aug, 1993) ... by M. Vachharajan, R Wieland
*                 Removed DISSPLA Graphics, replaced by NCAR Graphics (V3.0),
*                 using a modified version (ncarplts1.)of the PPPL "ncarplts"
*                 graphics package by Derek Fox.
*                 Graphics calls in this program are a mix of AUTOGRAPH
*                 calls (for 2d plots) and NCARPLTS calls for contours
*
*                To generate all plots %xprout fileid
*                To generate selected plots %xprout fileid "1,5 7 12 14,$"
*                To generate no plots %xprout fileid 0
*
*                 Modified Jan 1994 by R. Wieland to generate fileid.EQI
*                 file for stability code analysis. (Cf lwstbf vbl in WRFCN)
*************
*               Bsubs,Bsubu,Bsubv,Gmn  come in on the half-mesh
*               Vp,Mass,Pres(*),phip   come in on the half-mesh
*               Iotas(*), Lmns         come in on the half-mesh
*               Rmnc,Zmns              come in on the full-mesh
*               Jcuru,Jcurv            come in on the full-mesh
*               Gmn ==> Gsqrt                  on the half mesh
*               sknots,ystark,y2stark  spline knots come in on the sqrt(s) grid
*               pknots,ythom,y2thom    spline knots come in on the sqrt(s) grid
*************
*               (*)  Use Iota Spline Knots instead
*               (*)  Use Pressure Spline Knots instead
 
      pagedir(1) = '3D Outer magnetic flux surface'
      pagedir(2) = 'Flux Contours for 4 Toroidal Cuts'
      pagedir(3) = 'Flux, Mod-B, Jacobian Contours'
      pagedir(4) = 'J Contours'
      pagedir(5) = 'VMEC Convergence Criterion'
      pagedir(6) = 'Mercier criterion'
      pagedir(7) = 'Ballooning Growth Rates'
      pagedir(8) = 'Bootstrap vs VMEC <J*B>'
      pagedir(9) = 'Fluxes, Iota, Well'
      pagedir(10) = 'More Soln Profiles'
      pagedir(11) = 'More Soln Profiles'
      pagedir(12) = 'Magnetic Data Comparison'
      pagedir(13) = 'MSE Data Comparison'
      pagedir(14) = 'Iota Profile Comparison'
      pagedir(15) = 'Pressure Profile Comparison'
      pagedir(16) = 'Midplane Q Profile'
      pagedir(17) = 'Midplane J Profile'
      pagedir(18) = 'Midplane Shear Profile'
      pagedir(19) = 'Midplane Alpha-Prime Profile'
      pagedir(20) = 'Midplane Enclosed ITOR Profile'
      pagedir(21) = '{R,Z,L} set - the DEFAULT is to skip these'
      pagedir(22) = 'Limiter & Coil Positions'
      pagedir(23) = 'Poloidal Flux'
 
      print *, version
#ifndef CRAY
      timeloc = 'date >> '//tempfile
      call system (timeloc)
      open(unit=99, file=tempfile, status='old', err=100)
      read (99, 5) cdate
    5 format(a)
  100 continue
      close(unit=99, status='delete')
#else
      call date(cday)
      call clock(ctime)
      cdate = cday // ' ' // ctime
#endif
 
      numchars = 0
 
!       Get the Command Line Arguments (Case & Plot Selection)
#ifdef VAX
      call lib$get_foreign(input_id,,numchars)
#else
      numargs = iargc()
      if (numargs >= 1) then
         call getarg (1, input_id)
         ictrans_sel = '1,15 17,$'     ! the default is ALL PLOTS PLEASE
                                   ! except for the {R,Y,L} set
      endif
      if (input_id .eq. '-h' .or. numargs .eq. 0) then
         write (*, 104)
  104    format(/' ** To generate all plots %xprout fileid'/
     1' ** To generate selected plots %xprout fileid "1,5 12 14,$"'/
     2      ' ** To generate no plots %xprout fileid 0')
         write (*, 105) (i,pagedir(i),i=1,lastplot)
  105    format(/,' *** Plot Directory ***',2/,(i5,3x,a))
         stop 
      endif
      if (numargs .eq. 2) call getarg (2, ictrans_sel)
 
!     Parse the Plot Selection Part, If Present
      call icsel1 (ictrans_sel, lastplot, explist, maxlist, lenlist)
#endif
!
!     Does the Input File Exist ?
!
      numchars = len_trim(input_id)
      if (numchars .eq. 0) then
         print *, ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         stop 
      endif
      if (index(input_id,'wout.') .eq. 1) then
         input_id = input_id(6:)
         numchars = numchars - 6
      endif
      runlabel = input_id(1:numchars+1)//cdate(1:len_trim(cdate)+1)//
     1   '$'
      threed2_file = 'threed2.' // trim(input_id)
      spreadsheet  = 'spreadsheet.' // trim(input_id)
      answers_file = 'answers.' // trim(input_id)
      bal_grate_file = 'cobra_grate.' // trim(input_id)
      bal_input_file = 'in_cobra.' // trim(input_id)
      gmeta_file = 'gmeta.' // trim(input_id)

      call read_wout_file('wout.'//trim(input_id),ierr)

      if (ierr .eq. 1) stop 'could not read wout file: check extension'
      if (ierr.ne.0 .and. ierr.le.10) stop 'error in plotter read_wout'
      if (niter .le. 0) stop 'VMEC CODE DID NOT RUN PROPERLY!'
      
      
      if (imse .ne. (-1)) then
         i = index(mgrid_file,'.')
         if (i .eq. 0) then
            write (*, *) 
     1      'MGRID_FILE (in WOUT & INPUT) has incorrect format!'
             stop 
         endif
         tokid = mgrid_file(i+1:)
      else
         tokid = ' '
         mgrid_file = ' '
      endif

      if (ireconstruct > 0) then
        nbrbz = nbfld(nbrbz_id)/2          ! count in pairs
        if (nbrbz.gt.0 .and. ireconstruct .gt. 0) then
          allocate (brcoil(nbrbz), plbrfld(nbrbz), brbc(nbrbz),
     1    bzcoil(nbrbz), plbzfld(nbrbz), bzbc(nbrbz))
          ii = 0
          do i = 1, 2*nbrbz, 2
             ii = ii +1
             brcoil(ii)  = bcoil(i,1)
             bzcoil(ii)  = bcoil(i+1,1)
             brbc(ii)    = bbc(i,1)
             bzbc(ii)    = bbc(i+1,1)
             plbrfld(ii) = plbfld(i,1)
             plbzfld(ii) = plbfld(i+1,1)
          end do
        end if
      end if

      nmirnovset = 0
      nmirnov = 0
      do n = 1, nbsets
         if (n .eq. nmirnov_id) then
            nmirnovset = n
            nmirnov = nbfld(nmirnovset)
            exit 
         endif
      end do
 
      allocate (sqrt_phimod(mnmax*ns), phimod(mnmax*ns), dbzcods(ns), 
     1   ffp(ns), ub(ns), iotaf(ns), darea(ns), iotazb(ns), psi(ns),
     2   itors(ns), ipols(ns), szc(ns), szb(ns), dummy(ns),
     3   ixm(mnmax), ixn(mnmax), stat=mn)
      if (mn .ne. 0) stop 'Allocation error in plotout'        

      ixm(:mnmax) = nint(xm(:mnmax))
      ixn(:mnmax) = nint(xn(:mnmax))

      do mn = 1, mnmax
         if (ixm(mn).eq.0 .and. ixn(mn).eq.0) mn0 = mn
      end do

      ntor0 = 1 + ntor
      nrt = ns*nthpts
 
*****************************
*        COMPUTE FOURIER COEFFICIENTS OF d Bu/ds and d Bs/du
*        ON ZONE BNDRY GRID
*        Bsubu,v on 1/2 grid; dBsubu,v on full grid
*****************************
      hs = 1./(ns - 1)
      phip(1) = phip(2)
      
 
      do mn = 1, ntor0
         bmn(mn,1) = 1.5*bmn(mn,2) - 0.5*bmn(mn,3)
      end do

      do mn = 1, ntor
         gmn(mn,1) = 1.5*gmn(mn,2) - 0.5*gmn(mn,3)
      end do

      currvmn(1+ntor0:mnmax,1) = 0.
      bmn(1+ntor0:mnmax,1) = 0.
      gmn(1+ntor0:mnmax,1) = 0.
      darea(2:ns) = vp(2:ns)*overr(2:ns)
      overr(2:ns-1) = .5*(overr(2:ns-1)+overr(3:ns))
      iotaf(2:ns-1) = .5*(iotas(2:ns-1)+iotas(3:ns))
      jcuru(1) = 2.*jcuru(2) - jcuru(3)
      jcurv(1) = 2.*jcurv(2) - jcurv(3)
      overr(1) = 2.*overr(2) - overr(3)
      jcuru(ns) = 2.*jcuru(ns-1) - jcuru(ns-2)
      jcurv(ns) = 2.*jcurv(ns-1) - jcurv(ns-2)
      overr(ns) = 2.*overr(ns-1) - overr(ns-2)
      iotas(1) = 1.5*iotas(2) - 0.5*iotas(3)
      iotaf(1) = iotas(1)
      iotaf(ns) = 1.5*iotas(ns) - 0.5*iotas(ns-1)
 
*     Spline Knots are on sqrt(s) grid ...
      if (ireconstruct .gt. 0) then
         do i = 2, ns
            szc(i) = sqrt(hs*(i - 1.5))
            szb(i) = sqrt(hs*(i - 1))
         end do
         szc(1) = 0.
         szb(1) = 0.
         call splint (sknots, ystark, y2stark, isnodes, szc, iotas, 
     1      dummy, ns)
         call splint (sknots, ystark, y2stark, isnodes, szb, iotazb, 
     1      dummy, ns)
         iotaf(:ns) = iotazb(:ns)
         call splint (pknots, ythom, y2thom, ipnodes, szc, pres, 
     1      dummy, ns)
      endif
*     Forming <J-dot-GradPhi> / < R**-1 > ...
      where (overr .ne. 0.) jcurv = jcurv/overr
      hs = 1.0/(ns - 1)
      ohs = 1.0/hs
      r0 = rmnc(mn0,ns)

      call wrfcn (input_id)
      call plotter (r0, explist, lenlist, input_id, pagedir)
      call wrfcn2 (input_id)
*
*                   CALL EQI STABILITY FILE GENERATOR
*
*     Poloidal Flux / TwoPi
      if (lwstbf .and. ireconstruct.gt.0) then
         psi(1) = 0.
         do j = 2, ns
            psi(j) = psi(j-1) + hs*phip(j)*iotas(j)
         end do
*
         call wrstab (ns, psi, pres, iotaf, input_id, mnmax, 
     1      rmnc, zmns, bvco(ns), itor, ffp, jdotb, betatot, betapol,
     2      betator, betaxis)
      endif
!
!     Deallocate memory
!
      call read_wout_deallocate
      deallocate (sqrt_phimod, phimod, dbzcods, ffp, ub, iotaf, 
     1   darea, iotazb, psi, itors, ipols, szc, szb, dummy,
     2   ixm, ixn)
      if (allocated(brcoil)) deallocate (brcoil, plbrfld, brbc,
     1    bzcoil, plbzfld, bzbc)

      end program plotout
 

      subroutine wrfcn(input_id)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use Vpname1
      use Vpname2
      use Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character input_id*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: nboot = 33
      integer :: iloop, nmax1, mn, n1, j, jp, i, k, jcount
      real(rprec) :: es, tb1, tm1, tv1, tp1, ti1, ub1, phi_es,
     1    init_zeta, init_theta, bigno
      character, dimension(3) :: ichar1*5, ichar2*5
      character :: line*200, form_line*2000
      logical :: lexist
C-----------------------------------------------
 
      data ichar1/'RmnC(', 'ZmnS(', 'LmnS('/
      data ichar2/'RmnS(', 'ZmnC(', 'LmnC('/
      twopi = 8*atan(1.0)

!     Parse answers file if it exists to read in bootstrap current      
      open (unit=nboot, file=answers_file, status='old', iostat=i)
      if (i .eq. 0) then
         allocate (jboots(ns))
         jboots = huge(jboots)
         do while (i .eq. 0)
            read (nboot, '(a)', iostat=i, end=10) line
            j = index (line, 'Beta')
            if (j .gt. 0) exit
         end do

         do while (i .eq. 0)
            read(nboot, '(a)', iostat=i, end=10) line
            read(line, *) j, es, es, es, es, es, ub1
            if (j.gt.1 .and. j .le. ns) jboots(j) = ub1
         end do
      end if

 10   continue         
      close (unit=nboot)

!     Parse ballooning growth rate file if it exists (read first 10 values only)
      inquire (file=trim(bal_grate_file), exist = lexist)
      if (.not.lexist) then
         inquire(file=trim(bal_input_file), exist = lexist)
         if (lexist) then
            print *,' Running xcobravmec code ...'
            call system('xcobravmec ' // trim(bal_input_file)
     1        // ' >> /dev/null', k)
            if (k .ne. 0) lexist = .false.
         end if
      end if

      i = 0
      if (.not. lexist) goto 20
      open(unit=nboot, file=trim(bal_grate_file), status='old',iostat=k)
      if (k .eq. 0) allocate (balloon_grate(ns, max_no_angles_p),
     1    init_zeta_v(max_no_angles_p), init_theta_v(max_no_angles_p),
     2    ns_ball(ns), stat=k)
      balloon_grate = huge(balloon_grate)
      do while (k .eq. 0)
         read (nboot, *, iostat=k, end=20) init_zeta, init_theta, jcount
         if (jcount .gt. ns) k = -1
         if (k.eq.0 .and. i.lt.max_no_angles_p) then
            i = i + 1
            init_zeta_v(i) = init_zeta;   init_theta_v(i) = init_theta
            read (nboot, *, iostat=k) (ns_ball(jp), 
     1          balloon_grate(ns_ball(jp), i), jp = 1, jcount)
         else
            k = -1
         end if     
      end do
      
 20   continue

      max_no_angles = i
      if (lexist) close(unit=nboot)

!
!     PRINT SYMMETRIC TERMS -- RmnC and ZmnS
!     INTERPOLATE LAMBDA ONTO FULL MESH FIRST
!
      open(unit=nthreed2, file=threed2_file, status='replace')

      do mn = 1, mnmax
         if (ixm(mn) .ne. 0) then
            lmns(mn,1) = 0.
         else
            lmns(mn,1) = 1.5*lmns(mn,2) - 0.5*lmns(mn,3)
         end if
      end do
               
      do j = 2, ns-1
         do mn = 1,mnmax
            lmns(mn,j) = 0.5*(lmns(mn,j) + lmns(mn,j+1))
         end do         
      end do   
      
      do mn = 1,mnmax
         lmns(mn,ns) = 2.0*lmns(mn,ns-1) - lmns(mn,ns-2)
      end do      

      write (nthreed2, 200)
      do iloop = 1, 3
         nmax1 = 5
         do mn = 1, mnmax, 6
            if (mn > mnmax - 6) nmax1 = mnmax - mn
            write (nthreed2, 210) (ichar1(iloop),ixm(mn+n1),ixn(mn+n1),
     1         n1=0, nmax1)
            do j = 1, ns
!              es = (j - 1)*hs
               es = phi(j)/phi(ns)
               select case (iloop) 
               case default
                  write (nthreed2, 220) es, (rmnc(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (2) 
                  write (nthreed2, 220) es, (zmns(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (3) 
                  write (nthreed2, 220) es, (lmns(mn+n1,j),n1=0,nmax1)
               end select
            end do
         end do
      end do
 
!
!                 PRINT ASYMMETRIC TERMS -- RmnS and ZmnC
!
      if (iasym .eq. 1) then
         write (nthreed2, 200)
         do iloop = 1, 3
            nmax1 = 5
            do mn = 1, mnmax, 6
               if (mn > mnmax - 6) nmax1 = mnmax - mn
               write (nthreed2, 210) (ichar2(iloop),ixm(mn+n1),
     1            ixn(mn+n1), n1=0,nmax1)
               do j = 1, ns
!                 es = (j - 1)*hs
                  es = phi(j)/phi(ns)
                  select case (iloop) 
                  case default
                     write (nthreed2, 220) es,(rmns(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (2) 
                     write (nthreed2, 220) es,(zmnc(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (3) 
                     write (nthreed2, 220) es,(lmnc(mn+n1,j),n1=0,nmax1)
                  end select
               end do
            end do
         end do
      endif
!
!     DETERMINE RADIAL BETA PROFILE (SURFACE AVERAGED)
!
      phi(1) = 0.
      sqrt_phimod = 0.
      phimod = 0.
      do j = 1, ns
         sqrt_phimod(mn0+mnmax*(j-1)) = sqrt(abs(phi(j)))
         phimod(mn0+mnmax*(j-1)) = phi(j)
      end do
*     d(BVCO)/ds
      dbzcods(2:ns-1) = (bvco(3:ns)-bvco(2:ns-1))/hs
      dbzcods(ns) = 2*dbzcods(ns-1) - dbzcods(ns-2)
      dbzcods(1) = 2*dbzcods(2) - dbzcods(3)
      beta_vol(1) = 1.5*beta_vol(2) - 0.5*beta_vol(3)
      tb1 = 1.5*beta_vol(ns) - 0.5*beta_vol(ns-1)
      bvco(1) = 1.5*bvco(2) - 0.5*bvco(3)
      buco(1) = 0.
      ub(2:ns) = vp(2:ns)/phip(2:ns)
      mass(2:ns) = mass(2:ns)/abs(phip(2:ns))**gamma
      mass(1) = 1.5*mass(2) - 0.5*mass(3)
      pres(1) = 1.5*pres(2) - 0.5*pres(3)
      vp(1)   = 1.5*vp(2) - 0.5*vp(3)
      tm1 = 1.5*mass(ns) - .5*mass(ns-1)
      tv1 = 1.5*vp(ns) - .5*vp(ns-1)
      tp1 = 1.5*pres(ns) - .5*pres(ns-1)
      ti1 = 1.5*iotas(ns) - .5*iotas(ns-1)
      ub1 = 1.5*ub(ns) - .5*ub(ns-1)
      ub(1) = ub(2)*twopi
CDIR$   IVDEP
      do j = 2, ns - 1
         jp = j + 1
         buco(j) = .5*(buco(jp)+buco(j))
         ub(j) = .5*(ub(jp)+ub(j))*twopi
         vp(j) = .5*(vp(jp)+vp(j))
         mass(j) = .5*(mass(jp)+mass(j))
         pres(j) = 0.5*(pres(jp) + pres(j))
         beta_vol(j) = .5*(beta_vol(jp)+beta_vol(j))
         bvco(j) = .5*(bvco(jp)+bvco(j))
      end do
      ub(ns) = twopi*ub1
      vp(ns) = tv1
      mass(ns) = tm1
      pres(ns) = tp1
      beta_vol(ns) = tb1
      buco(ns) = 2.*buco(ns) - buco(ns-1)
      bvco(ns) = 2.*bvco(ns) - bvco(ns-1)
      write (nthreed2, 230)
!     Normalized FF-Prime : ffp(j)
      do j = 1, ns
!        es = (j - 1)*hs
         es = phi(j)/phi(ns)
         itor = -twopi*buco(j)
         itors(j) = itor/dmu0
         ipol = twopi*(bvco(j)-bvco(ns))
         ipols(j) = ipol/dmu0
         ffp(j) = bvco(j)*dbzcods(j)/phip(j)/iotaf(j)/
     1      bvco(ns)**2
         write (nthreed2, 220) es, bvco(j), buco(j), itor, ipol, ffp(j)
      end do
      write (nthreed2, 240)
      do i = 1, ns
!        es = (i - 1)*hs
         es = phi(i)/phi(ns)
         write (nthreed2, 220) es, twopi**2*vp(i), ub(i), mass(i), 
     1      pres(i), iotaf(i), beta_vol(i)
      end do
 
      return 
 
!  Do the above before the call to PLOTTER; do the below after
 
      entry wrfcn2 (input_id)
 
      write (nthreed2, 260)
      do i = 1, ns
!        es = (i - 1)*hs
         es = phi(i)/phi(ns)
         write (nthreed2, 220) es, jcurv(i), jdotb(i)/bdotgradv(i)
      end do
 
      write (nthreed2, 250) beta_vol(1)

      close(unit=nthreed2)
      
!
!     WRITE SPREADSHEET-FORMAT FILE (SPREAD....)
!
      open(unit=nthreed2, file=spreadsheet, status='replace')

      write (form_line, *)
     1   "(9x,' ES',9x,'PHI',5x,'<BZETA>',4x,'<BTHETA>',8x,'ITOR',8x,
     1   'IPOL',9x,'FF''', 10x,'V''',5x,'dV/dPHI',8x,'MASS',11x,'P',
     2   8x,'IOTA',8x,'BETA',6x,'<JTOR>',6x,'<JPOL>',7x,'<J*B>',
     3   5x,'MERCIER', 5x,'SHEAR_M', 6x,'WELL_M',6x, 
     4   'CURR_M', 6x,'GEOD_M'"

      if (allocated(jboots)) then
         write (line, *) ",4x,'<JBOOTS>'"
         form_line = trim(form_line) // trim(line)
      end if

      if (allocated(balloon_grate)) then
         write (line, *) "(", max_no_angles, 
     1      "(a6,i3,1x,a2,i3,a1))"
         write (line, trim(line)) (",'  U=",
     1          init_theta_v(i),"V=", init_zeta_v(i),"'",
     2          i = 1, max_no_angles)
         form_line = trim(form_line) // trim(line)
      end if

      form_line = trim(form_line) // "/)"      

      write (nthreed2, form_line)

      do j = 1, ns
         es = (j - 1)*hs
         phi_es = phi(j)/phi(ns)
         itor = -twopi*buco(j)
         ipol = twopi*(bvco(j)-bvco(ns))
         write (nthreed2, '(1p100e12.4)', advance='no')
     1   es, phi_es, bvco(j), buco(j), itor, ipol, ffp(j),
     2   twopi**2*vp(j), ub(j), mass(j), pres(j), iotaf(j), beta_vol(j),
     3   jcurv(j), jcuru(j), jdotb(j), dmerc(j), 
     4   dshear(j), dwell(j), dcurr(j), dgeod(j)

         if (allocated(jboots)) 
     1      write (nthreed2, '(1pe12.4)', advance='no') jboots(j)

         if (allocated(balloon_grate)) then
            bigno = huge(balloon_grate(1,1))
            where (balloon_grate .eq. bigno) balloon_grate = 0
            write (nthreed2, '(1p100e13.4)', advance='no') 
     2        (balloon_grate(j,i), i = 1, max_no_angles)
         end if
     
         write (nthreed2, *)          !Advance to next record
      end do

      if (allocated(jboots)) deallocate (jboots, stat=k)
      if (allocated(balloon_grate)) deallocate (balloon_grate, 
     1    init_zeta_v, init_theta_v, ns_ball, stat=k)

      close (nthreed2)

 
  200 format(//,40x,'FOURIER COEFFICIENTS X(m,n)',/)
  210 format(//,9x,'PHI',4x,6(4x,a5,i1,',',i3,')'),/)
  220 format(1p8e15.3)
  230 format(//,25x,'COVARIANT COMPONENTS OF B',
     1   ' AND INTEGRATED CURRENTS',2/,9x,'PHI',10x,'<BZETA>',8x,
     2   '<BTHETA>',8x,'ITOR',11x,'IPOL',12x,'FF''',/)
  240 format(//,9x,'PHI',11x,'VP',12x,'dV/dPHI',10x,'MASS',13x,'P',10x,
     1   'IOTA',12x,'BETA',/)
  250 format(//,'  BETA ON AXIS (SUM OVER MODES) = ',1pe10.3)
  260 format(//,9x,'PHI',11x,'<JTOR>',2x,'<JdotB>/<bdotgradv>',/,9x,3x
     1   ,11x,'[A/M2]',10x,'[A/M]',/)
      
 
      end subroutine wrfcn
      

      subroutine plotter(r0, explist, lenlist, input_id, pagedir)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod
      use Vpname2
      USE Vpltcn2
      USE Vpltcn6
      USE Vrzarray
      USE Vplotdata
      USE Vmagaxis
      USE Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: lenlist
      real(rprec) :: r0
      character :: explist*(*), input_id*(*)
      character, dimension(*) :: pagedir*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ntdu = 31, ntdv = 21
      integer, parameter :: nplot1 = 25, ntheta1 = 2, 
     1  nox = 1, noy = 3, lx = -1, ly = 0, iend = 3,
     2  ilog = 1, icart = 0, oneppg = 2, fourppg = 1
      integer, dimension(4) :: ivar = (/0, 0, 2, 0/),
     1  ncon = (/30, 50, 30, 30/), ndeg_param = 
     2    (/ 0, 45, 90, 180 /)     
      real(rprec), parameter :: grsize1 = 3.5, grsize2 = 2.5, 
     1  grfull = 7.5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, nplots, les, kz, k, kt, j, l, js, 
     1   ndeg, ns1, mn, mp, nt, ichar, istat, nrange,
     2   page_count
      real(rprec), dimension(ns) :: br, bz, phin
      real(rprec), dimension(2*ns) :: oqmid
      real(rprec), dimension(100) :: time
      real(rprec), dimension(nplot1) :: raxis_v, zaxis_v
      real(rprec), dimension(nfloops + 1) :: delflm
      real(rprec), dimension(nbloops) :: delbc
      real(rprec), dimension(nfloops + 1) :: indflm
      real(rprec), dimension(nbloops) :: indbc
      real(rprec), dimension(:,:), allocatable :: dummy1
      real(rprec), dimension(:), allocatable :: modb, sqrt_phim, 
     1  phim, gsqrt, torcur, dummy2, r12, z12, ru12, zu12
      real(rprec) :: dth, denom, phiangle, offset,
     1   phinorm, t1, a0, an, cosphi, sinphi, bigno
      character :: page*10, pagedesc*50, mchar*100, nchar*100, 
     1   pchar*100, ititlet*100, fname*80, fdum*80
      logical :: lneed
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      logical , EXTERNAL :: makeplot
C-----------------------------------------------

      newframe = .FALSE.

      mn = mnmax*ns; i = max(ns, nobd+1, nbrbz, nmirnov, imse, itse)
      allocate (dummy1(mnmax,ns), dummy2(nrt), r(nrt), z(nrt), 
     1   modb(nrt), sqrt_phim(nrt), phim(nrt), gsqrt(nrt),torcur(nrt), 
     2   dbsubudstz(nrt), bsubutz(nrt), bsubvtz(nrt), r12(nrt),
     3   z12(nrt), rdata(i), pldata(i), stat = istat)
     
      if (istat .ne. 0) stop 'allocation error in plotter'


      dummy1 = 0
      dummy2 = 0
 
      noxc = 0
      noyc = noy - 1
      lxc = lx
      lyc = ly
 
!
!                 SETUP NCAR GRAPHICS
!
      call grafinit (1, input_id)

      call setclr ('yellow')                      !default color for AUTOGRAPH-only plots
      call pcsetc ('FC - function code delimiter', '!')
      write (pchar, 150) '!PGU!U$'
 
!
!                 COMPUTE R,Z SCALES
!

      mchar = '$'
      nplots = 1
      nrange = max(2*ntor, 1)
      page_count = 0
      
      if (ntor .ne. 0) nplots = nrange/2 + 1
      les = 1 + nthpts*(ns - 1)
      do kz = 1, nplots
         call totz(nthpts,ns,nrange,kz,r,z,rmnc,zmns,rmns,zmnc)
         if (kz .eq. 1) then
            rmax = r(les)
            rmin = r(les)
            zmax = z(les)
            zmin = z(les)
         endif
         rmax = max(rmax,maxval(r(les:nrt)))
         rmin = min(rmin,minval(r(les:nrt)))
         zmax = max(zmax,maxval(z(les:nrt)))
         zmin = min(zmin,minval(z(les:nrt)))
      end do
 
!
!     3-D Surface
!      
      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      if (makeplot(explist,lenlist,page,pagedesc)) 
#ifndef RISC
     1   call tdplotter(rmnc(1,ns), zmns(1,ns), rmns, zmnc)
#else
     1   print *,' 3-D plotting not available on RISC (yet)!'
#endif      

      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      if (makeplot(explist,lenlist,page,pagedesc)) then
         nt = 4          !!No. plots on this page
         if (ntor.le.1) nt = 1
         noxc = 0
         do kz = 1, nt
            ndeg = ndeg_param(kz)
            k = 1 + ndeg
            write (mchar, 30) 'N!B!f!NPGL!u!PRU! = ', 
     1       ndeg, '!S!o$'
            call totz (nthpts, ns, 360, k ,r, z,rmnc, zmns, rmns, zmnc)
            call totz (nthpts, ns, 360, k , sqrt_phim, dummy2, 
     1      sqrt_phimod, zmns, dummy1, dummy1)

            rmagaxis = r(1)
            zmagaxis = z(1)

            call contour (r(nthpts+1), z(nthpts+1), sqrt_phim(nthpts+1),
     1         mchar, 'CONTOURS OF !PRL!$!PGU!U$', ncon(1), ivar(1), 
     2         grsize1, ns - 1, runlabel)
         end do     
      end if

      lneed = .false.
      do i = 1, 2
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) lneed = .true.
      end do

      if (.not. lneed) goto 8050
      page_count = page_count - 2
      
!       ---Transform from (m,n) Space to (u,v) Space---
      cont_plots: do kz = 1, nplots
         call totz (nthpts,ns,nrange,kz,r,z,rmnc,zmns,rmns,zmnc)
         call totz (nthpts, ns, nrange, kz, modb, dummy2, bmn, 
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, sqrt_phim, dummy2, 
     1      sqrt_phimod, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, phim, dummy2, phimod, 
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, gsqrt, dummy2, gmn,
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, torcur, dummy2, 
     1      currvmn, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, bsubutz, dummy2, 
     1      bsubumn, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nrange, kz, bsubvtz, dummy2, 
     1      bsubvmn, zmns, dummy1, dummy1)
 
         do kt = 1, nthpts
            r12(kt) = r(kt)
            z12(kt) = z(kt)
c        do this for r1,z1 because they go like sqrt(s) near origin
            r12(kt+nthpts) = .5*(sqrt(2.) - 1.)*(r(kt+nthpts)-r(kt))
            z12(kt+nthpts) = .5*(sqrt(2.) - 1.)*(z(kt+nthpts)-z(kt))
            do j = 2, ns
               l = kt + nthpts*(j - 1)
               r12(l) = 0.5*(r(l) + r(l-nthpts))
               z12(l) = 0.5*(z(l) + z(l-nthpts))
            end do
         end do

         rmagaxis = r12(1)
         zmagaxis = z12(1)
 
         ndeg = nint(360.*(kz - 1)/nrange)
         if (ntor .ne. 0) write (mchar, 30) 'N!B!f!NPGL!u!PRU! = ', 
     1    ndeg, '!S!o$'
   30    format(a,i3,a)
         noxc = 0
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            call contour (r(nthpts+1), z(nthpts+1), sqrt_phim(nthpts+1),
     1         'CONTOURS OF !PRL!$!PGU!U$', mchar, ncon(1), ivar(1), 
     2         grsize1, ns - 1, runlabel)
            call contour (r12(nthpts+1), z12(nthpts+1), modb(nthpts+1),
     1         'MOD-B CONTOURS$', mchar, 
     2         ncon(2), ivar(2), grsize1, ns-1, runlabel)
            call contour (r, z, phim, 
     1         '!PGU!U!PRL! and !PGL!H!PRU! CONTOURS$', mchar, ncon(3),
     2         ivar(3), grsize1, ns, runlabel)
            call contour (r12(nthpts+1), z12(nthpts+1), gsqrt(nthpts+1),
     1         'JACOBIAN CONTOURS$', mchar, 
     2         ncon(4), ivar(4), grsize1, ns-1, runlabel)
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
            newframe = .TRUE.
            call contour (r, z, torcur, 
     1       'CONTOURS OF RJ!S!.!NIRL!/!PGL!U$',
     1        mchar, ncon(1), ivar(2), grsize1, ns, runlabel)
         end if
         
         if (kz .lt. nplots) page_count = page_count - 2
      end do cont_plots
 
 8050 continue
 
      phinorm = maxval(abs(phi(:ns)))
      phin(:ns) = abs(phi(:ns)/phinorm)

      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      itfsq = min(itfsq, 100)
      if (itfsq.gt.1. and. makeplot(explist,lenlist,page,pagedesc)) then
         wdot(1) = wdot(2)
         fsqt(1) = fsqt(2)
      
         do j = 1, itfsq
            time(j) = real(niter*(j - 1))/real(itfsq - 1)
         end do

         write (nchar, 150) '!PRL!=!PRU!F!S!2!N! dV$'
         noxc = 1
         ndata = 0
         call tableau ('RESIDUAL FORCE$', 'ITERATIONS$', nchar, 
     1     fsqt, fsqt, 1, itfsq, time, grsize1, lxc, lyc, 
     2     noxc, noyc, fourppg, ilog, runlabel)
         write (nchar, 150) '-dW/dt$'
         call tableau ('ENERGY MINIMIZATION$', 'ITERATIONS$', nchar, 
     1     wdot, wdot, 1, itfsq, time, grsize1, lxc, lyc, 
     2     noxc, noyc, fourppg, ilog, runlabel)
*
*                 COMPUTE R,Z MAGNETIC AXES VS. TOROIDAL ANGLE
*
         if (ntor .ne. 0) then
            ns1 = 1
            do kz = 1, nplot1
               call totz (ntheta1, ns1, nplot1-1, kz, r, z, rmnc,
     1         zmns, rmns, zmnc)
               time(kz) = real(kz - 1)/real(nplot1-1)
               raxis_v(kz) = r(ntheta1)
               zaxis_v(kz) = z(ntheta1)
            end do
            write (nchar, 150) 'R!B!M!N!(!PGL!u!PRU!)$'
            call tableau ('MAGNETIC AXIS (R)$', 
     1      'N!B!F!NPGL!u!PRU!/2!PGL!P$', nchar, raxis_v, 
     2      raxis_v, 1, nplot1, time, grsize1, lxc, lyc, noxc, 
     3      noyc, fourppg, icart, runlabel)
            write (nchar, 150) 'Z!B!M!N!(!PGL!u!PRU!)$'
            call tableau ('MAGNETIC AXIS (Z)$', 
     1      'N!B!F!NPGL!u!PRU!/2!PGL!P$', nchar, zaxis_v, 
     2      zaxis_v, 1, nplot1, time, grsize1, lxc, lyc, noxc, 
     3      noyc, fourppg, icart, runlabel)
         endif

         specw(1) = 2*specw(2) - specw(3)
         call tableau ('SPECTRAL WIDTH$', pchar, '<M>$', specw, 
     1   specw, 1, ns, phin, grsize1, lxc, lyc, noxc, noyc, 
     2   fourppg, icart, runlabel)

      end if

!
!                 COMPUTE EQUILIBRIUM PROFILES
!
      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      if (makeplot(explist,lenlist,page,pagedesc)) then
         call tableau ('MERCIER$', pchar, 'DM$', dmerc(3:), 
     1      dmerc(3:), 1, ns - 4, phin(3), grsize1, 
     2      lxc, lyc, noxc, noyc, fourppg, icart, runlabel)
         call tableau ('SHEAR and WELL$', pchar, 'DS, DW$', 
     1      dshear(3:), dwell(3:), 2, ns - 4, phin(3),
     2      grsize1, lxc, lyc, noxc, noyc, fourppg, icart, runlabel)
         call tableau ('CURRENT and GEODESIC$', pchar, 'Dcur, Dgeo$', 
     1      dcurr(3:), dgeod(3:), 2, ns - 4, phin(3), 
     2      grsize1, lxc, lyc, noxc, noyc, fourppg, icart, runlabel)
      endif
 
!
!     Ballooning growth rates
!
      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)

      if (makeplot(explist,lenlist,page,pagedesc)) then
         if (.not.allocated(balloon_grate)) then
            print *,' Can not locate ballooning growth rate data'
         else   
            angles: do kz = 1, max_no_angles
               bigno = huge(balloon_grate(1,kz))
               ndata = 0
               do j = 1, ns
                  if (balloon_grate(j,kz) .eq. bigno) cycle
                  ndata = ndata + 1
                  rdata(ndata) = phin(j)
                  pldata(ndata) = balloon_grate(j, kz)
               end do     
               write (nchar, '(a,i4,a,i4,a)') 
     1           'BALLOONING GROWTH RATE, U=', nint(init_theta_v(kz)),
     2           ' V=', nint(init_zeta_v(kz)), '$'
               noxc = 0
               noyc = 0
               br(1:ns) = balloon_grate(1:ns, kz)
               bz(1:ns) = phin(1:ns)
               where (br .eq. bigno) br = 0
               call tableau (nchar, pchar, '$', br, br,
     1         0, ns, bz, grfull, lxc, lyc, noxc, noyc, 
     2         oneppg, icart, runlabel)
            end do angles
         end if
      end if
      
!     Bootstrap data vs VMEC current
      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
 
      if (makeplot(explist,lenlist,page,pagedesc)) then
         if (.not. allocated(jboots)) then
            print *,' Unable to located bootstrap data'
         else
            bigno = huge(jboots(1))
            ndata = 0
            do j = 2, ns
               if (jboots(j) .eq. bigno) cycle
               ndata = ndata+1
               rdata(ndata) = 0.5*(phin(j) + phin(j-1))               !Half-mesh from bs code
               pldata(ndata) = jboots(j)
            end do     
            write (nchar, 150) '<J*B>$'
            noxc = 0
            noyc = 0
            call tableau ('BOOTSTRAP CONSISTENCY$', pchar, nchar, jdotb,
     1      jdotb, 1, ns, phin, grfull, lxc, lyc, noxc, noyc, oneppg, 
     2      icart, runlabel)
         end if
      end if


      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)

      if (makeplot(explist,lenlist,page,pagedesc)) then
         do j = 1, ns
            br(j) = hs*(j - 1)
         end do
         ndata = 0
         call tableau ('TOROIDAL FLUX$', ' s$', pchar, phi, phi, 1, 
     1      ns, br, grsize1, lxc, lyc, noxc, noyc, fourppg, 
     2      icart, runlabel)
         t1 = 0.5*twopi*hs
         do j = 2, ns
            bz(j) = iotas(j)
            br(j) = br(j-1) + t1*phip(j)*(iotas(j)+iotas(j-1))
         end do
         call tableau ('POLOIDAL FLUX$', pchar, '!PGL!V$', br, br, 1, 
     1      ns, phin, grsize1, lxc, lyc, noxc, noyc, fourppg, 
     2      icart, runlabel)
         a0 = 1.5*bz(2) - .5*bz(3)
         an = 1.5*bz(ns) - .5*bz(ns-1)

         bz(2:ns-1) = .5*(bz(2:ns-1)+bz(3:ns))
         br(2:ns-1) = 1./bz(2:ns-1)
         bz(1) = a0
         if (a0 .ne. 0.) br(1) = 1./a0
         bz(ns) = an
         br(ns) = 1./an
         call tableau ('IOTA$', pchar, 'q,!PGL!I$', br, bz, 2, ns, 
     1      phin, grsize1, lxc, lyc, noxc, noyc, fourppg, 
     2      icart, runlabel)
         bz(:ns) = (ub(:ns)-ub(1))/ub(1)
         write (nchar, 150) '(V''(!PGU!U!PRU!)-V''(0))/V''(0)$'
         call tableau ('MAGNETIC WELL$', pchar, nchar, bz, bz, 1, ns, 
     1      phin, grsize1, lxc, lyc, noxc, noyc, fourppg, icart, 
     2      runlabel)
 
      endif
 

      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      if (makeplot(explist,lenlist,page,pagedesc)) then
 
         ndata = 0

         call tableau ('MASS$', pchar, 'M(!PGU!U!PRU!)$', mass, 
     1       mass, 1, ns, phin, grsize1, lxc, lyc, noxc, 
     2       noyc, fourppg, icart, runlabel)
         call tableau ('PRESSURE$', pchar, 'P(!PGU!U!PRU!)$',pres,
     1       pres, 1, ns, phin, grsize1, lxc, lyc, noxc, 
     2       noyc, fourppg, icart, runlabel)
         call tableau ('BETA$', pchar, '<!PGL!B!PRU!>$', beta_vol, 
     1      beta_vol, 1, ns, phin, grsize1, lxc, lyc, noxc, 
     2      noyc, fourppg, icart, runlabel)
         call tableau ('POLOIDAL CURRENT DENSITY$', pchar, 
     1      '<J!S!.!NIRL!/!PGL!h!PRU!>$', jcuru, jcuru, 1, ns, phin, 
     2      grsize1, lxc, lyc, noxc, noyc, fourppg, icart, runlabel)
 
      endif
 

      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
**      Finally ... <J-dot-B> / <B-dot-gradPhi>
      where (bdotgradv .ne. 0.) br = jdotb/bdotgradv
      
      if (makeplot(explist,lenlist,page,pagedesc)) then
 
         ndata = 0
         call tableau ('FLUX-AV. JTOR (A/m!S!2!N!)$', pchar, 
     1      '<J!S!.!NIRL!/!PGL!U!PRU!>/<R!S!-1!N!>$', jcurv, 
     2      jcurv, 1, ns, phin, grsize1, lxc, lyc, 
     3      noxc, noyc, fourppg, icart, runlabel)
         call tableau ('<J!S!.!N!B>/<B!S!.!NIRL!/!PGL!U!PRU!>$', pchar, 
     1      'Amps/m$', br, br, 1, ns, phin, 
     2      grsize1, lxc, lyc, noxc, noyc, fourppg, icart, runlabel)
         call tableau ('TOROIDAL CURRENT$', pchar, 'ITOR$', itors, 
     1      itors, 1, ns, phin, grsize1, lxc, lyc, noxc, noyc, 
     2      iend, icart, runlabel)
 
      endif
 

      page_count = page_count+1
      if (ireconstruct .gt. 0) then
         write (page, *) page_count
         pagedesc = pagedir(page_count)

         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = 0
            if (ireconstruct > 0 .and. nobd > 0) then
 
*        Index 0 = Diamagnetic Flux Comparison
c                                                !mWb
               delflm(1) = (phidiam - delphid)*1.E3
               indflm(1) = 0
               rdata(1) = indflm(1)
               pldata(1) = delflm(1)
               ndata = 1
 
               do i = 1, nobd
                  ndata = ndata + 1
                  delflm(ndata) = (dsiobt(i)-(dsiext(i)+plflux(i)))*1.E3
                  indflm(ndata) = i
                  rdata(ndata) = indflm(i)
                  pldata(ndata) = delflm(i)
               end do
 
               write (ititlet, 8055) flmwgt*100
 8055          format('Flux Loops : RMS Error= ',f5.2,' %$')
               call tableau (ititlet, 'Index (0=DMG)$', 'Delta (mWb)$', 
     1            delflm, delflm, 1, ndata, indflm, grsize1, lxc, lyc, 
     2            noxc, noyc, fourppg, icart, runlabel)
 
            endif
 
 
            if (nbrbz > 0) then
               do i = 1, nbrbz
c                                                !mT
                  delbc(i) = (brbc(i)-(brcoil(i)+plbrfld(i)))*1.E3
                  indbc(i) = i
               end do
               ndata = nbrbz
               rdata(:ndata) = indbc(:ndata)
               pldata(:ndata) = delbc(:ndata)
               write (ititlet, 8056) bcwgt*100
 8056          format('Br Loops : Br/Bz RMS Error= ',f5.2,' %$')
               call tableau (ititlet, 'Index$', 'Delta (mTesla)$', 
     1            delbc, delbc, 1, nbrbz, indbc, grsize1, lxc, lyc, 
     2            noxc, noyc, fourppg, icart, runlabel)
 
               do i = 1, nbrbz
c                                                !mT
                  delbc(i) = (bzbc(i)-(bzcoil(i)+plbzfld(i)))*1.E3
                  indbc(i) = i
               end do
               ndata = nbrbz
               rdata(:ndata) = indbc(:ndata)
               pldata(:ndata) = delbc(:ndata)
               write (ititlet, 8057) bcwgt*100
 8057          format('Bz Loops : Br/Bz RMS Error= ',f5.2,' %$')
               call tableau (ititlet, 'Index$', 'Delta (mTesla)$', 
     1            delbc, delbc, 1, nbrbz, indbc, grsize1, lxc, lyc, 
     2            noxc, noyc, iend, icart, runlabel)
            endif
 
            if (nmirnov > 0) then
               do i = 1, nmirnov
                  delbc(i) = (bbc(i,nmirnovset)-(bcoil(i,nmirnovset)+
     1               plbfld(i,nmirnovset)))*1.E3 !mT
                  indbc(i) = i
               end do
               ndata = nmirnov
               rdata(:ndata) = indbc(:ndata)
               pldata(:ndata) = delbc(:ndata)
               write (ititlet, 8058) bcwgt*100
 8058          format('Mirnov Loops : RMS Error= ',f5.2,' %$')
               call tableau (ititlet, 'Index$', 'Delta (mTesla)$', 
     1            delbc, delbc, 1, nmirnov, indbc, grsize1, lxc, lyc, 
     2            noxc, noyc, iend, icart, runlabel)
            endif
 
            ndata = 0
         
         endif
      endif
 
 
      page_count = page_count+1
      if (ireconstruct > 0 .and. (nobd.ne.0 .or. nbfldn.ne.0)) then
         write (page, *) page_count
         pagedesc = pagedir(page_count)
 
         if (makeplot(explist,lenlist,page,pagedesc) 
     1     .and. allocated(brcoil)) call bgraph (
     1      phidiam, delphid, dsiobt, dsiext, plflux, nobd, brbc, 
     2      brcoil, plbrfld, bzbc, bzcoil, plbzfld, nbrbz, 
     3      bbc(1,nmirnovset), bcoil(1,nmirnovset), 
     4      plbfld(1,nmirnovset), nmirnov, flmwgt, bcwgt)
      endif
 
      if (imse>2 .or. itse>0) then
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = imse
            rdata(:ndata) = rstark(:ndata)
            pldata(:ndata) = datastark(:ndata)
            write (ititlet, 8071) msewgt*100
 8071       format('MSE PITCH : RMS Error= ',f5.2,' %$')
            call tableau (ititlet, 'RMID$', 'ANGLE (DEG)$', anglemse, 
     1         anglemse, 1, 2*ns - 1, rmid, grfull, lxc, lyc, noxc, 
     2         noyc, oneppg, icart, runlabel)
 
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = imse
            rdata(:ndata) = rstark(:ndata)
            pldata(:ndata) = qmeas(:ndata)
            oqmid(:2*ns-1) = 1.0/qmid(:2*ns-1)
            call tableau ('IOTA PROFILE$', 'RMID$', 'IOTA$', oqmid, 
     1         oqmid, 1, 2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, 
     2         oneppg, icart, runlabel)
 
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = itse
            rdata(:ndata) = rthom(:ndata)
            pldata(:ndata) = 0.001*datathom(:ndata)*pfac
            presmid(:2*ns-1) = 0.001*presmid(:2*ns-1)
            write (mchar, 95) 'P-MID [KPa] $'
c       write(mchar,95)'P-MID [KPa] (data scaled by ',pfac,')$'
   95       format(a,f6.2,a)
            write (ititlet, 8073) tswgt*100
 8073       format('MIDPLANE PRESSURE : RMS Error= ',f5.2,' %$')
            call tableau (ititlet, 'RMID$', mchar, presmid, presmid, 1, 
     1         2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, oneppg, 
     2         icart, runlabel)
 
         endif
 
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = 0
            call tableau ('MIDPLANE Q$', 'RMID$', 'Q-MID$', qmid, qmid, 
     1         1, 2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, oneppg, 
     2         icart, runlabel)
 
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = 0
            curmid(:2*ns-1) = 1.E-6*curmid(:2*ns-1)
            call tableau ('TOROIDAL CURRENT DENSITY$', 'RMID$', 
     1         'JTOR-MID [MA/m!S!2!N!]$', curmid, curmid, 1, 2*ns - 1, 
     2         rmid, grfull, lxc, lyc, noxc, noyc, oneppg, icart, 
     3         runlabel)
 
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = 0
            call tableau ('SHEAR$', 'RMID$', 'S(RMID)$', shear, shear, 1
     1         , 2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, oneppg, 
     2         icart, runlabel)
 
         endif
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
 
            ndata = 0
            call tableau ('ALPHA(P'')$', 'RMID$', 'ALPHA(RMID)$', alfa, 
     1         alfa, 1, 2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, 
     2         oneppg, icart, runlabel)
         endif
 
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (.not.makeplot(explist,lenlist,page,pagedesc)) go to 8130
 
         br(:ns) = 0.
         do i = 2, ns
            br(i) = .5*twopi*hs*darea(i)*(jcurv(i)+jcurv(i-1)) + br(i-1)
         end do
         do i = 1, 2*ns - 1
            if (i <= ns) then
               curmid(i) = br(ns+1-i)
            else
               curmid(i) = br(i-ns+1)
            endif
         end do
         ndata = 0
         call tableau ('ENCLOSED ITOR (A)$', 'RMID$', 'ITOR$', curmid, 
     1      curmid, 1, 2*ns - 1, rmid, grfull, lxc, lyc, noxc, noyc, 
     2      oneppg, icart, runlabel)
 
      else
      
         page_count = page_count+8

      endif             ! end of if( (imse.gt.2) .or. (itse.gt.0) ) then
*
*                 PLOT R,Z, LAM AND MAGNETIC FIELD FOURIER PROFILES
*
 8130 continue
 
      page_count = page_count+1
      write (page, *) page_count
      pagedesc = pagedir(page_count)
      if (makeplot(explist,lenlist,page,pagedesc)) then
 
         do mn = 1, mnmax
            lxc = lx
            lxc = ly
            mp = nint(xm(mn))
            nt = nint(xn(mn))
            write (mchar, 120) ' m =', mp, ' n =', nt
  120       format(a,i3,a,i3,'$')
            if (mn .ne. mn0) then
               write (nchar, 150) 'RmC,ZmS$'
            else
               write (nchar, 150) 'R-R!B!0!N!,Z$'
            endif
            br(:ns) = rmnc(mn,:ns)
            if (mn .eq. mn0) 
     1        br(:ns) = br(:ns) - r0
            bz(:ns) = zmns(mn,:ns)
            call tableau (mchar, pchar, nchar, br, bz, 2, ns, phin, 
     1         grsize2, lxc, lyc, nox, noy, fourppg, icart, runlabel)
            br(:ns) = lmns(mn,:ns)
            if (mn .eq. mn0) then
               lxc = lx + 1
            else
               call tableau (mchar, pchar, 'LAMBDA-Sin$', br, br, 1, 
     1            ns, phin, grsize2, lxc, lyc, nox, noy, fourppg, 
     2            icart, runlabel)
            endif
         end do
 
         if (iasym .eq. 1) then
            do mn = 1, mnmax
               lxc = lx
               mp = nint(xm(mn))
               nt = nint(xn(mn))
               write (mchar, 120) ' m =', mp, ' n =', nt
               write (nchar, 150) 'RmS,ZmC$'
               br(:ns) = rmns(mn,:ns)
               bz(:ns) = zmnc(mn,:ns)
               call tableau (mchar, pchar, nchar, br, bz, 2, ns, phin, 
     1            grsize2, lxc, lyc, nox, noy, fourppg, icart, runlabel)
               br(:ns) = lmnc(mn,:ns)
               if (mn .eq. mn0) then
                  lxc = lx + 1
               else
                  call tableau (mchar, pchar, 'LAMBDA-Cos$', br, br, 1, 
     1               ns, phin, grsize2, lxc, lyc, nox, noy, fourppg, 
     2               icart, runlabel)
               endif
            end do
         endif
 
!
!       PLOT LIMITER AND COILS IF AVAILABLE
!
      endif
 
      if (imse .gt. 1) then
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) call limplot
 
         page_count = page_count+1
         write (page, *) page_count
         pagedesc = pagedir(page_count)
         if (makeplot(explist,lenlist,page,pagedesc)) then
            allocate (ru(nrt), zu(nrt), ru12(nrt), zu12(nrt))        
            kz = 1    !Could loop over toroidal planes, eventually
            call totzu(nthpts,ns,nrange,kz,ru,zu,rmnc,zmns,rmns,zmnc)
            do kt = 1, nthpts
               do j = 2, ns
                  l = kt + nthpts*(j - 1)
                  ru12(l) = 0.5*(ru(l)+ru(l-nthpts))
                  zu12(l) = 0.5*(zu(l)+zu(l-nthpts))
               end do
            end do
            call pfluxeng (gsqrt, r12, z12, ru12, zu12)
            deallocate (ru, zu, ru12, zu12)
         end if
 
 150     format(a)
 
      else
         page_count = page_count+2
      endif

      call grafinit (2, input_id)

      deallocate (dummy1, dummy2, r, z, modb, sqrt_phim, phim, 
     1   gsqrt, torcur, dbsubudstz, bsubutz, bsubvtz, r12, z12,
     2   rdata, pldata)


      end subroutine plotter

#ifndef RISC
      subroutine tdplotter(rmnb, zmnb, rdum, zdum)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nfp, rprec
      use Vpltcn2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*), intent(in) :: 
     1   rmnb, zmnb, rdum, zdum
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
c
c define a parameter specifying the maximum size of the array required
c for the list of triangles.
c
        integer, parameter :: mtri0=5000
c
c define values of pi, two times pi, and pi over 180.
c
        real, parameter :: pi=3.14159265358979323846,
     1        twopi=2.*pi,dtor=pi/180.
c
c set the values determining the resolution of the grid over which the
c surface is generated.
c
        integer, parameter :: idim = 51
c set the desired minimum and maximum values of u and v (for the grid
c over which the surface is generated).
c
        real, parameter :: umin = 0., umax = twopi
        real, parameter :: vmin = 0., vmax = twopi
c
c set the desired values of parameters determining the eye position.
c ang1 is a bearing angle, ang2 is an elevation angle, and rmul is a
c multiplier of the length of the diagonal of the data box, specifying
c the distance from the center of the box to the eye.
c
        real, parameter :: ang1 = 215., ang2 = 35., rmul = 2.9
c
c iste is a flag that says whether to do a simple image (iste=0),
c a one-frame stereo image (iste=-1), or a two-frame stereo image
c (iste=+1). igrids is a flag for drawing grids and labels (igrids=1) or
c suppressing them (igrids=0).
c
        integer, parameter :: iste = 0, igrids = 0
c
c aste is the desired angle (in degrees) between the lines of sight for
c a pair of stereo views.
c
        real, parameter :: aste = 4.
c
c wosw is the width of the stereo windows to be used in one-frame stereo
c images; the width is stated as a fraction of the width of the plotter
c frame.  (the windows are centered vertically; horizontally, they are
c placed as far apart as possible in the plotter frame.)  the value used
c must be positive and non-zero; it may be slightly greater than .5, if
c it is desired that the stereo windows should overlap slightly.
c
        real, parameter :: wosw = .5
c
c set the desired value of the flag that says whether the basic color
c scheme will be white on black (ibow=0) or black on white (ibow=1).
c
        integer, parameter :: ibow = 1
c
c set the desired value of the flag that says whether shading of the
c surfaces will be done using gray scales (iclr=0) or colors (iclr=1).
c
        integer, parameter :: iclr = 1
c
c set the desired values of the shading parameters.  values of shde
c near 0 give brighter colors and values near 1 give pastel shades.
c values of shdr near 0 give a narrow range of shades and values near
c 1 give a wide range of shades.
c
        real, parameter :: shde = 0.1, shdr = 0.8
c
c set the desired value of the rendering-style index for the surface.
c the ordering index iord = 0 gives fastest rendering; use iord=1 if
c overlapping, intersecting surfaces possible
c
        integer, parameter :: isrs = 3, iord = 0
c
c declare variables to hold labels.
c
        character*64 :: xnlb,ynlb,znlb,xilb,yilb,zilb
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
        integer :: jdim, mtri, ntri, jplt, j, i, jvvmi, jvvma,
     1     iuvmi, iuvma, lnlg                
        integer, dimension(:), allocatable :: itwk
        real, dimension(:,:), allocatable :: rval, hval
c
c declare a local array to hold the triangle list and a couple of
c temporary variables to be used in sorting the list.
c
        real, dimension(:,:), allocatable :: rtri, rtwk
        real :: xmin, xmax, ymin, ymax, xmid, ymid, zmid,
     1     xeye, yeye, zeye, r, p, h, otep,
     2     xvpl, xvpr, yvpb, yvpt, xwdl, xwdr, ywdb, ywdt,
     3     rv00, pv00, hv00, rv01, pv01, hv01, rv10, pv10, hv10, 
     4     rv11, pv11, hv11, pval, xval, yval, zval
C-----------------------------------------------
c
c define labels for the edges of the box.
c
        data xnlb / ' -4 -3 -2 -1 0 1 2 3 4 ' /
        data ynlb / ' -4 -3 -2 -1 0 1 2 3 4 ' /
        data znlb / ' -4 -3 -2 -1 0 1 2 3 4 ' /
c
        data xilb / 'x coordinate values' /
        data yilb / 'y coordinate values' /
        data zilb / 'z coordinate values' /

c
c define arithmetic statement functions for r, phi, and h as functions
c of u and v.
c
c       rval(u,v)=cos(u)+ a*cos(u-v)+b
        pval(i,j)=vmin+(real(j-1)/real(jdim-1))*(vmax-vmin)
c       hval(u,v)=sin(u)+ a*sin(u-v)
c
c define arithmetic statement functions to transform cylindrical
c coordinates into cartesian coordinates.
c
        xval(r,p,h)=r*cos(p)
        yval(r,p,h)=r*sin(p)
        zval(r,p,h)=h

        jplt = 51/nfp                           !!v points per period
        if (jplt .lt. 10) jplt = 10
        jdim = 1 + jplt*nfp
c
c
c set the desired minimum and maximum values of x, y, and z.
c
        xmin =-rmax
        xmax = rmax
        ymin =-rmax
        ymax = rmax
        
        mtri = mtri0 * nfp
        allocate (rtri(10,mtri),rtwk(mtri,2),itwk(mtri), 
     1      rval(idim,2), hval(idim,2), stat=i)
        if (i .ne. 0) stop 'allocation error in 3dplot'

c
c turn clipping off.
c
        call gsclip (0)
c
c double the line width.
c
        call gslwsc (2.)
c
c define colors to use.
c
        call tdclrs (1,ibow,shde,shdr,11,42,4)
c
c select font number 25, turn on the outlining of filled fonts, set the
c line width to 1, and turn off the setting of the outline color.
c
        call pcseti ('fn - font number',25)
        call pcseti ('of - outline flag',1)
        call pcsetr ('ol - outline line width',1.)
        call pcsetr ('oc - outline line color',-1.)
c
c make tdpack characters a bit bigger.
c
        call tdsetr ('cs1',1.25)
c
c define tdpack rendering styles 1 through 7, using black-and-white
c shading or colored shading, whichever is selected.  the indices
c 1-7 can then be used as final arguments in calls to tditri, tdstri,
c and tdmtri.
c
        if (iclr.eq.0) then
          call tdstrs (1,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (2,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (3,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (4,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (5,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (6,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (7,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
        else
          call tdstrs (1,43,74, 43, 74,-1,-1,1,0.,0.,0.) ! gray/gray
          call tdstrs (2,43,74, 75,106,-1,-1,1,0.,0.,0.) ! gray/red
          call tdstrs (3,43,74,107,138,-1,-1,1,0.,0.,0.) ! gray/green
          call tdstrs (4,43,74,139,170,-1,-1,1,0.,0.,0.) ! gray/blue
          call tdstrs (5,43,74,171,202,-1,-1,1,0.,0.,0.) ! gray/cyan
          call tdstrs (6,43,74,203,234,-1,-1,1,0.,0.,0.) ! gray/magenta
          call tdstrs (7,43,74,235,266,-1,-1,1,0.,0.,0.) ! gray/yellow
        end if
c
c initialize the count of triangles in the triangle list.
c
        ntri=0
c
c for each box on a rectangular grid in the uv plane, generate two
c triangles and add them to the triangle list.  each triangle is
c transformed from cylindrical coordinates to cartesian coordinates.
c
        j = 1
        call totz(idim,1,jplt,j,rval(1,1),hval(1,1),
     1              rmnb,zmnb,rdum,zdum)
        do 102 j=1,jdim-1
          call totz(idim,1,jplt,j+1,rval(1,2),hval(1,2),
     1              rmnb,zmnb,rdum,zdum)
          jvvmi=1
          jvvma=2
          do 101 i=1,idim-1
            iuvmi=i
            iuvma=i+1
            rv00=rval(iuvmi,jvvmi)
            pv00=pval(iuvmi,j)
            hv00=hval(iuvmi,jvvmi)
            rv01=rval(iuvmi,jvvma)
            pv01=pval(iuvmi,j+1)
            hv01=hval(iuvmi,jvvma)
            rv10=rval(iuvma,jvvmi)
            pv10=pval(iuvma,j)
            hv10=hval(iuvma,jvvmi)
            rv11=rval(iuvma,jvvma)
            pv11=pval(iuvma,j+1)
            hv11=hval(iuvma,jvvma)
            if (ntri.lt.mtri) then
              ntri=ntri+1
              rtri(1,ntri)=xval(rv10,pv10,hv10)
              rtri(2,ntri)=yval(rv10,pv10,hv10)
              rtri(3,ntri)=zval(rv10,pv10,hv10)
              rtri(4,ntri)=xval(rv00,pv00,hv00)
              rtri(5,ntri)=yval(rv00,pv00,hv00)
              rtri(6,ntri)=zval(rv00,pv00,hv00)
              rtri(7,ntri)=xval(rv01,pv01,hv01)
              rtri(8,ntri)=yval(rv01,pv01,hv01)
              rtri(9,ntri)=zval(rv01,pv01,hv01)
              rtri(10,ntri)=real(isrs)
            end if
            if (ntri.lt.mtri) then
              ntri=ntri+1
              rtri(1,ntri)=xval(rv01,pv01,hv01)
              rtri(2,ntri)=yval(rv01,pv01,hv01)
              rtri(3,ntri)=zval(rv01,pv01,hv01)
              rtri(4,ntri)=xval(rv11,pv11,hv11)
              rtri(5,ntri)=yval(rv11,pv11,hv11)
              rtri(6,ntri)=zval(rv11,pv11,hv11)
              rtri(7,ntri)=xval(rv10,pv10,hv10)
              rtri(8,ntri)=yval(rv10,pv10,hv10)
              rtri(9,ntri)=zval(rv10,pv10,hv10)
              rtri(10,ntri)=real(isrs)
            end if
  101     continue
          rval(:,1) = rval(:,2)
          hval(:,1) = hval(:,2)
  102   continue
c
c find the midpoint of the data box (to be used as the point looked at).
c
        xmid=.5*(xmin+xmax)
        ymid=.5*(ymin+ymax)
        zmid=.5*(zmin+zmax)
c
c determine the distance (r) from which the data box will be viewed and,
c given that, the eye position.
c
        r=rmul*sqrt((xmax-xmin)**2+(ymax-ymin)**2+(zmax-zmin)**2)
c
        xeye=xmid+r*cos(dtor*ang1)*cos(dtor*ang2)
        yeye=ymid+r*sin(dtor*ang1)*cos(dtor*ang2)
        zeye=zmid+r*sin(dtor*ang2)
c
c initialize the stereo offset argument to do either a single view or
c a left-eye view (whichever is selected by the value of iste).
c
        if (iste.eq.0) then
          otep=0.                    !  (single view)
        else
          otep=-r*tan(dtor*aste/2.)  !  (left-eye view)
        end if
c
c initialize tdpack.
c
  109   call tdinit (xeye,yeye,zeye,xmid,ymid,zmid,
     +                              xmid,ymid,zmid+r,otep)
c
c if stereo views are being done, do the requested thing, either by
c redoing the set call to put them side by side on the same frame,
c or by calling frame to put them on separate frames.
c
        if (otep.ne.0.) then
          if (iste.lt.0) then
            call getset (xvpl,xvpr,yvpb,yvpt,xwdl,xwdr,ywdb,ywdt,lnlg)
            if (otep.lt.0.) then
              call set  (1.-wosw,1.,.5-.5*wosw,.5+.5*wosw,
     +                           xwdl,xwdr,ywdb,ywdt,lnlg)
            else
              call set  (  0., wosw,.5-.5*wosw,.5+.5*wosw,
     +                           xwdl,xwdr,ywdb,ywdt,lnlg)
            end if
          else
            if (otep.gt.0.) call frame
          end if
        end if
c
c order the triangles in the triangle list.
c
        call tdotri (rtri,mtri,ntri,rtwk,itwk,iord)
c
        if (ntri.eq.mtri) then
          print * , 'triangle list overflow in tdotri'
          stop
        end if
c
c draw labels for the axes.
c
        if (igrids .eq. 1) then
        call tdlbls (xmin,ymin,zmin,xmax,ymax,zmax,
     +               xnlb,ynlb,znlb,xilb,yilb,zilb,1)
c
c draw the sides of the box that could be hidden.
c
        call tdgrds (xmin,ymin,zmin,xmax,ymax,zmax,
     +               .1*(xmax-xmin),.1*(ymax-ymin),.1*(zmax-zmin),
     +                                                       12,1)
        end if
c
c draw the triangles in the triangle list.
c
        call tddtri (rtri,mtri,ntri,itwk)
c
c draw the sides of the box that could not be hidden.
c
        if (igrids .eq. 1) 
     +  call tdgrds (xmin,ymin,zmin,xmax,ymax,zmax,
     +               .1*(xmax-xmin),.1*(ymax-ymin),.1*(zmax-zmin),
     +                                                       12,0)
c
c if a left-eye view has just been done, loop back for a right-eye view.
c
        if (otep.lt.0.) then
          otep=-otep
          go to 109
        end if
c
c advance the frame.
c
        call frame

        deallocate (rtri,rtwk,itwk, rval, hval, stat=i)

c
c done.
c
c
      end subroutine tdplotter
#endif


      subroutine contour(xplot, yplot, func, itop, mtitle, ncon, ivar, 
     1   grsz, ns1, runlabel1)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Vpname1
      USE Vpltcn2
      USE Vpltcn6
      USE Vrzarray
      USE Vmagaxis
      USE Vplotdata
      USE Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ncon, ivar, ns1
      real(rprec) :: grsz
      character itop*(*), mtitle*(*), runlabel1*(*)
      real(rprec), dimension(nthpts*ns1) :: xplot, yplot, func
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nstep, nrav, nzav, n, k, l, kt, j
      real, dimension(ns1) :: rthet, zthet
      real :: step, exstep, fac, xmin, xmax, ymin, ymax, 
     1  delmax, delmin, delf
      real, dimension(:), allocatable :: f_real, x_real, y_real
      character :: lblx*5, lbly*5
C-----------------------------------------------
*
*     LOCATE PHYSICAL ORIGIN
*
*     ROUND-OFF AXIS
*
      step = max(rmax - rmin,zmax - zmin)
      exstep = 2. - nint(log10(step))
      fac = 10.**exstep
      nstep = nint(.25*(fac*step)*1.05)
      nrav = nint(.50*(rmax + rmin)*fac)
      nzav = nint(.50*(zmax + zmin)*fac)
      step = nstep/fac
      xmin = nrav/fac - 2.*step
      xmax = nrav/fac + 2.*step
      ymin = nzav/fac - 2.*step
      ymax = nzav/fac + 2.*step
      delmax = maxval(func(:nthpts*ns1))
      delmin = minval(func(:nthpts*ns1))
      delf = (delmax - delmin)/ncon
      lblx = 'R$'
      lbly = 'Z$'
      
      allocate (f_real(nthpts*ns1), x_real(nthpts*ns1), 
     1          y_real(nthpts*ns1))
 
      f_real = real(func)
      x_real = real(xplot)
      y_real = real(yplot)
      
      call contgraf1 (f_real, nthpts, ns1, xmin, xmax, ymin, ymax, 
     1    1, -2, delf, 4, x_real, y_real, 2, 0, itop, lblx, 
     2    lbly, runlabel1, mtitle, 6, 6, gwnd(1,noxc+1))
 
c         Plot the Boundary -- Last Closed Flux Surface
      l = nthpts*(ns1 - 1) + 1
      call agcurv (x_real(l), 1, y_real(l), 1, nthpts, 1)
 

      deallocate (f_real, x_real, y_real)
      
c         Plot the Magnetic Axis ?
      if (ivar .eq. 0) then
c        cf pg 104 in NCAR Fundamentals manual
         call gsmk (3)
         call gsmksc (1.)
         call gpm (1, real(rmagaxis), real(zmagaxis))
      endif
 
c         Plot Constant THETA Contours?
      if (ivar .eq. 2) then
         do kt = 1, nthpts, 2
            do j=1,ns1
              rthet(j) = r(kt+nthpts*(j-1))
              zthet(j) = z(kt+nthpts*(j-1))
            end do
            call agcurv (rthet, 1, zthet, 1, ns1, 1)
         end do
      endif
 
      if (noxc.eq.3 .or. newframe) then
         call frame
         newframe = .FALSE.
         noxc = 0
      else
         noxc = noxc + 1
      endif

      end subroutine contour
      

      subroutine totz(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nfp, mnmax, xm, xn, rprec, dp
      use Vpname1, only: nrt
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ntheta, ns1, nplots, kz
      real(rprec), dimension(*), intent(out) :: r, z
      real(rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1._dp/(ntheta - 1)
      piz = 1._dp/(nfp*nplots)
      nrth = ns1*ntheta
      if (nrth .gt. nrt) stop 'nrth > nrt in totz'

      r(:nrth) = 0
      z(:nrth) = 0

      do j = 1, ns1
         jes = ntheta*(j - 1)
         mes = mnmax*(j - 1)
         do mn = 1, mnmax
            xm0 = xm(mn)*pit
            xn0 = xn(mn)*piz
            do kt = 1, ntheta
               l = kt + jes
               m = mn + mes
               arg = xm0*(kt - 1) - xn0*(kz - 1)
               iarg = arg
               arg = twopi*(arg - iarg)
               r(l) = r(l) + rmnc(m)*cos(arg) + rmns(m)*sin(arg)
               z(l) = z(l) + zmns(m)*sin(arg) + zmnc(m)*cos(arg)
            end do
         end do
      end do

      end subroutine totz


      subroutine totzu(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: dp, rprec, nfp, mnmax, xm, xn
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(rprec), dimension(*), intent(out) :: r, z
      real(rprec), dimension(*), intent(in) :: 
     1  rmnc, zmns, rmns, zmnc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1._dp/(ntheta - 1)
      piz = 1._dp/(nfp*nplots)
      nrth = ns1*ntheta
      r(:nrth) = 0.
      z(:nrth) = 0.
      do j = 1, ns1
         jes = ntheta*(j - 1)
         mes = mnmax*(j - 1)
         do mn = 1, mnmax
            xm0 = xm(mn)*pit
            xn0 = xn(mn)*piz
            do kt = 1, ntheta
               l = kt + jes
               m = mn + mes
               arg = xm0*(kt - 1) - xn0*(kz - 1)
               iarg = arg
               arg = twopi*(arg - iarg)
               r(l) =r(l) + xm(mn)*(rmns(m)*cos(arg) - rmnc(m)*sin(arg))
               z(l) =z(l) - xm(mn)*(zmnc(m)*sin(arg) - zmns(m)*cos(arg))
            end do
         end do
      end do

      end subroutine totzu


      subroutine tableau(ititlet, ititlex, ititley, fun1, fun2, nfuns, 
     1   npts, xpoint, grsize, lx, ly, nox, noy, igraph, itype, 
     2   runlabel)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vplotdata
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nfuns, npts, lx, ly, nox, noy, igraph, itype
      real(rprec) :: grsize
      real(rprec), dimension(npts) :: fun1, fun2, xpoint
      character*(*) :: ititlet, ititlex, ititley, runlabel
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: oneppg = 2
      integer :: i, nrunlabel
      real, dimension(npts,2) :: ydra
      real, dimension(npts) :: x_real
      real :: ymax, ymin
C-----------------------------------------------
c
c       to fix the axes ... call AGSETR('Y/MINIMUM.',YMIN)
c                           call AGSETR('Y/MAXIMUM.',YMAX)
c                           CALL AGSETR('Y/NICE.',0.)
c       cf pp 224 et al in Fundamentals Manual
c
      call findlen (runlabel, nrunlabel)
      call agseti ('FRAME.', 2)
      if (itype .eq. 1) call agseti ('Y/LOG.', 1)

      x_real = real(xpoint)
      ydra(:npts,1) = fun1(:npts)
      ydra(:npts,2) = fun2(:npts)
      ymax = maxval(fun1(:npts))
      ymin = minval(fun1(:npts))

      if (igraph .ne. oneppg) then
         call agsetp ('GRAPH WINDOW.', gwnd(1,nox+1), 4)
         call agsetc ('LABEL/NAME.', 'B')
         call agseti ('LINE/NUMBER.', -100)
         call agsetc ('LINE/TEXT.', ititlex)
         call agsetc ('LABEL/NAME.', 'L')
         call agseti ('LINE/NUMBER.', 100)
         call agsetc ('LINE/TEXT.', ititley)
         call agseti ('BACKGROUND TYPE.', 1)
         call agseti ('LABEL/CONTROL.', 2)
         if (nfuns .le. 1) then            
            if (nfuns .gt. 0) then
               call ezxy (x_real, ydra, npts, ititlet)
            else
               x_real(2) = x_real(npts)
               !!Forces off plot (so it will not be visible)
               ydra(1,1) = 1.E36      
               ydra(2,1) = 1.E36
               call ezxy (x_real, ydra, 2, ititlet)
            end if
            if (ndata .gt. 0) then
               call setclrpm ('green')
               call points (rdata, pldata, ndata, -4, 0)
            endif
            if (itype .eq. 1) call agseti ('Y/LOG.', 0)
         else
            call agseti ('ROW.', 1)
            call lineselect ('DASH', iselect)
            call ezmxy (x_real, ydra, npts, nfuns, npts, ititlet)
            if (itype .eq. 1) call agseti ('Y/LOG.', 0)
         endif
         nox = nox + 1
 
*           Define RUNLABEL Informational Line
*           Nox=4 signifies the 4th and last plot
*           Igraph=3 signifies forced last plot on page
*           Igraph=2 signifies one plot on page
         if (nox==4 .or. igraph==3 .or. igraph==2) then
            nox = 0
            call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
            call plchhq(0.5,0.01,runlabel(1:nrunlabel),0.008,0.,0.)
            call frame
         endif
      else 
         if (nox .ne. 0) call frame
         call agseti ('LINE/MAX.', 80)
         call agseti ('ROW.', 1)
         call agsetp ('GRAPH WINDOW.', gwnd(1,5), 4)
         call agsetc ('LABEL/NAME.', 'B')
         call agseti ('LINE/NUMBER.', -100)
         call agsetc ('LINE/TEXT.', ititlex)
         call agsetc ('LABEL/NAME.', 'L')
         call agseti ('LINE/NUMBER.', 100)
         call agsetc ('LINE/TEXT.', ititley)
         call agseti ('BACKGROUND TYPE.', 1)
         call agseti ('LABEL/CONTROL.', 2)
         if (nfuns .le. 1) then
            ymax = max(ymax,maxval(pldata(:ndata)))
            ymin = min(ymin,minval(pldata(:ndata)))
            if (index(ititley,'P-MID') > 0) ymin = 0
            if (index(ititley,'Q-MID') > 0) ymin = 0
            ymax = ymax + (ymax - ymin)/20
            call agsetr ('Y/MINIMUM.', ymin)
            call agsetr ('Y/MAXIMUM.', ymax)
            call agsetr ('Y/NICE.', 0.)
            if (nfuns .gt. 0) then
               call ezxy (x_real, ydra, npts, ititlet)
            else
               x_real(2) = x_real(npts)
               !!Forces off plot (so it will not be visible)
               ydra(1,1) = 1.E36      
               ydra(2,1) = 1.E36
               call ezxy (x_real, ydra, 2, ititlet)
            end if
            if (ndata .gt. 0) then
               call setclrpm ('green')
               call points (rdata, pldata, ndata, -4, 0)
            endif
            if (itype .eq. 1) call agseti ('Y/LOG.', 0)
         else
            call lineselect ('DASH', 2)
            call ezmxy (x_real, ydra, npts, nfuns, npts, ititlet)
         endif
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call plchhq (0.5, 0.01, runlabel(1:nrunlabel), 0.008, 0., 0.)
         call frame
         nox = 0
      endif

      end subroutine tableau
      

      subroutine prunlabel
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Vpname2
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrunlabel
C-----------------------------------------------
 
 
*
*     Define RUNLABEL Informational Line
*
      call findlen (runlabel, nrunlabel)
c                                                !at bottom
      call plchhq (0.5, 0.01, runlabel(1:nrunlabel), 0.008, 0., 0.)

      end subroutine prunlabel
      

      subroutine scalesiz(r1, r2, z1, z2, sizer, sizez)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real r1, r2, z1, z2, sizer, sizez
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real :: sr0, sz0, facr, sz1, fac1
C-----------------------------------------------
      sr0 = 7.9
      sz0 = 9.8
      sizer = sr0
      sizez = sz0
      facr = sr0/(r2 - r1)
      sz1 = facr*(z2 - z1)
      if (sz1 > sz0) then
         fac1 = sz1/sz0
         sizer = sizer/fac1
         sizez = sizez/fac1
      endif

      end subroutine scalesiz
      

      subroutine cklimpos(ios) 
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: rx1, rx2, zy1, zy2, condif, pfcspec, 
     1   nsets, nsetsn, nrgrid, nzgrid
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ios
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j
C-----------------------------------------------
c     Array PFCSPEC describes the PF Coil Set :: PFCSPEC(i1,i2,i3)
c     1st index : 1 = R; 2 = Z; 3 = nturns; 4 = cfac ! describes a loop
c     2nd index : labels the loop in a particular coil set
c     3rd index : labels the coil set
 
c     - - -  Auto Count Each Set - - -
 
      l20: do j = 1, nsets
         do i = 1, n2max
            if (pfcspec(1,i,j) .eq. 0.) then
               nsetsn(j) = i - 1
               cycle  l20
            endif
         end do
      end do l20
 
c     - - - Check for Overflows - - -
 
      if (nrgrid>ndim1 .or. nzgrid>ndim1) then
         print *, 
     1      ' Namelist values for nr,nz too large .. set to default'
         go to 990
      endif
 
      if (nsets > nsetsmax) then
         print *, ' Namelist values for nsets too large ... max is ', 
     1      nsetsmax
         call exit (2)
      endif
 
      do i = 1, nsets
         if (nsetsn(i) > n2max) then
            print *, 
     1         ' Namelist values for 1 nsetsn is too large ... max is '
     2         , n2max
            call exit (2)
         endif
      end do
 
  990 continue
 
c     - - - Default Grid Parameters - - -
c           Is LIMPOS file missing or
c       are the parameters unspecified ?
 
 
      if (nrgrid .eq. 0) nrgrid = nrgridsv
      if (nzgrid .eq. 0) nzgrid = nzgridsv
 
      if (rx1*rx2*zy1*zy2 .eq. 0.0) then          ! no LIMPOS File present
         rx1 = rx1sv
         rx2 = rx2sv
         zy1 = zy1sv
         zy2 = zy2sv
      endif
 
c     - - - Contour Plot Defaults - - -
 
      if (condif .eq. 0.) condif = condifsv
 
      ios = 0                                    ! default now

      end subroutine cklimpos
      
 
      subroutine gengrid
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: rx1, rx2, zy1, zy2, nrgrid, nzgrid
      use Vpname3
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: kr, kz, lenr, lenz
      real :: rminb, rmaxb, zminb, zmaxb, delrb, delzb
      character :: nrc*10, nzc*10, ier*10
C-----------------------------------------------
c     Define grid limits for poloidal field plot based on grid
c     in the LIMPOS file.
 
 
      delrb = (rx2 - rx1)/(nrgrid - 1)
      delzb = (zy2 - zy1)/(nzgrid - 1)
 
c      print *, ' Rmin,Rmax= ', rx1,rx2
c      print *, ' Zmin,Zmax= ', zy1,zy2
c      print *, ' nR,nZ= ', nrgrid,nzgrid
c      print *, ' delR,delZ= ',delrb,delzb
 
      do kr = 1, nrgrid
         rgrid(kr) = rx1 + (kr - 1)*delrb
      end do
 
      do kz = 1, nzgrid
         zgrid(kz) = zy1 + (kz - 1)*delzb
      end do
 
      end subroutine gengrid
      

      subroutine limplot
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: limitr, nlim, ns,
     1   zy1, zy2, rx1, rx2, nrgrid, nzgrid, condif      
      use Vpname2
      use Vpname3
      USE Vrzarray
      USE Vindat2
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ios, iostat, j, jes, l, icharg
      real, dimension(4) :: gwnd
      real :: orient, mr, mz, width = 10.0, sizex = 8.0, sizey = 9.8, 
     1   rx0 = 0.0, zy0 = -4.0, rxst = 0.1, zyst = 0.1, rxstp = 0.5, 
     2   zystp = 1.3, ticksz = 0.02, unit, file, status
C-----------------------------------------------
 
      namelist /namin17/rx1, rx2, zy1, zy2, nrgrid, nzgrid, condif
 
 
*     Read the LIMPOS.GRID namelist (/NAMIN17/) file to override
*     defacto defns made in Block Data; ignore values in 38 file
 
*     Use defacto values rather than 38 values
!      nrgrid = nrgrid_17
!      nzgrid = nzgrid_17
!      rx1 = rx1_17
!      rx2 = rx2_17
!      zy1 = zy1_17
!      zy2 = zy2_17
 
*     Save MGRID values of plotting specs
      nrgridsv = nrgrid
      nzgridsv = nzgrid
      rx1sv = rx1
      rx2sv = rx2
      zy1sv = zy1
      zy2sv = zy2
      condifsv = condif
 
      ios = 0
      if (initnml .eq. 0) open(unit=50, file='limpos.grid', 
     1    status='old', iostat=ios)
      if (ios .eq. 0) then
         write (*, *) 'Overriding r/z ', 
     1      'grid parameters with local LIMPOS.GRID file for plot 17'
         read (50, namin17)
         initnml = 1
      endif
 
      call cklimpos (ios)
      if (ios .ne. 0) return 
      call gengrid
 
*     Assume proper aspect ratio requires shortening of x axis.
*     Fix Mz (Z window-span) .eq. 1. Let Dz,Dr be units-span for R,Z, resp.
*     Then demand that the correct Mr (R window-span) is given by
*     Mr = Dr/Dz*Mz
*     If Mr > 1. then reverse everything.
 
      mr = (rx2 - rx1)/(zy2 - zy1)*1.0
      if (mr .le. 1.) then
         gwnd(1) = 0. + (1. - mr)/2.
         gwnd(2) = 1. - (1. - mr)/2.
         gwnd(3) = 0.0
         gwnd(4) = 1.0
      else
         mz = (zy2 - zy1)/(rx2 - rx1)*1.0
         gwnd(1) = 0.0
         gwnd(2) = 1.0
         gwnd(3) = 0. + (1. - mz)/2.
         gwnd(4) = 1. - (1. - mz)/2.
      endif
 
      call setclr ('yellow')
 
      call agseti ('LINE/MAX.', 80)
      call agseti ('FRAME.', 2)
      call agsetc ('LABEL/NAME.', 'B')
      call agseti ('LINE/NUMBER.', -100)
      call agsetc ('LINE/TEXT.', 'R')
      call agsetc ('LABEL/NAME.', 'T')
      call agseti ('LINE/NUMBER.', 100)
      call agsetc ('LINE/TEXT.', 'LIMITER & COIL POSITIONS')
      call agsetc ('LABEL/NAME.', 'L')
      call agseti ('LINE/NUMBER.', 100)
      call agsetc ('LINE/TEXT.', 'Z')
      call agseti ('BACKGROUND TYPE.', 1)
      call agsetp ('GRAPH WINDOW.', gwnd, 4)
      call agsetf ('X/MIN.', real(rx1))
      call agsetf ('X/MAX.', real(rx2))
      call agsetf ('X/NICE.', 0.)
      call agsetf ('Y/MIN.', real(zy1))
      call agsetf ('Y/MAX.', real(zy2))
      call agsetf ('Y/NICE.', 0.)
      call agseti ('LABEL/CONTROL.', 2)
c      call agseti('SET.',-1)       !suppress curve drawing
 
c                                                !dummy parameters
      call agstup (real(rx1), 1, 0, 0, 0, real(zy1), 1, 0, 0, 0)
      call agback
 
 
c     plot flux surface(s)
      call gspmci (61)                           ! white polymarkers
      do j = ns, 1, -3                           ! plot many
c      do j=ns,ns                            ! plot outer only
         jes = nthpts*(j - 1)
         l = 1 + jes
c                                                !plot as curve
         call agcurv (real(r(l)), 1, real(z(l)), 1, nthpts, 1)
c        call points(r(l),z(l),nthpts,0,0)    !plot as individual pts
      end do
 
      call plcoil
 
*
*       Allow for more than one limiter set
*
 
c      call agseti('DASH/PATTERNS/1.', 61680)
c      call agseti('DASH/SELECTOR.',1)
c      do l = 1,nlim
c       if(limitr(l).gt.1) then
c         call agcurv(xlim(1,l),1,ylim(1,l),1,limitr(l),1)
c        endif
c      enddo
 
!c     plot limiter sets as thick dashed lines
!      call gsln(2)
!      call gslwsc(width)
!      call gsplci(44)                    ! line color
!      do l = 1,nlim
!       if(limitr(l).gt.1) then
!c         plot as a dashed line [gsln(2) above] with width given by
!c         gslwsc above
!         call curved(xlim(1,l),ylim(1,l),limitr(l))
!       endif
!      enddo
!      call gslwsc(1.)
!      call gsln(1)
 
 
c     plot limiter sets as individual points
      call gspmci (61)                    ! white polymarkers for points
      do l = 1, nlim
         if (limitr(l) > 1) then
            if (l .eq. 1) j = -2
            if (l .eq. 2) j = -5
            call points(real(xlim(1,l)),real(ylim(1,l)),limitr(l),j,0)
         endif
      end do
 
c     mark points where boundary touches limiter (-4==circle)
      icharg = ichar('C')
      call points (real(xlimt), real(ylimt), ltouch, icharg, 0)
 
      call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
      call prunlabel
 
      call frame
 
      end subroutine limplot
      

      subroutine plcoil
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: pfcspec, nsets, nsetsn
      use Vpname3
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i1, i2, jj
      real(rprec) :: orient
C-----------------------------------------------
 
      do i1 = 1, nsets
         do i2 = 1, nsetsn(i1)
            orient = pfcspec(3,i2,i1)*pfcspec(4,i2,i1)
            if (orient .ne. 0) then
               if (orient < 0) then
                  jj = -5
               else 
                  jj = -4
               end if      
               call points(real(pfcspec(1,i2,i1)),
     1                     real(pfcspec(2,i2,i1)),1,jj,0)
            endif
         end do
      end do

      end subroutine plcoil


      subroutine pfluxeng(gsqrt, r12, z12, ru12, zu12)
c
c---------------------------------------------------------------------
c*  PFLUXENG.FOR
c   Dick Wieland 12-Mar-1993  PPL - TFTR only -
c   Package to contour plot total (plasma current + PF coil set)
c   poloidal flux on an (r,z) scale that encompasses the PF Coil set.
c   For speed, assume up-down symmetry in calculation contribution
c   due to plasma current.
c
c   Added following Contour Plot related vbls to the limpos namelist:
c       nrgrid, nzgrid : grid size
c       condif         : DISSPLA contour density
 
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: xn, pres, extcur, phip, bsubumn,
     1  bsubsmn, ns, mnmax, rprec,  
     1  iotas, xm, zy1, rx1, rx2, zy2, condif, nrgrid, nzgrid, nsets
      use Vpname1
      use Vpname2
      use Vpname3
      USE Vrzarray
      USE Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*) :: gsqrt, r12, z12, ru12, zu12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, js, kt, l, lj, lp, mn, ios, itype, iostat
      real :: t1, t2
      real(rprec), dimension(:), allocatable :: 
     1   dbsubsdu, dbsubuds, dummy1, dummy2, bsubsdu, bsubuds, taua
      character, dimension(4) :: pflhead*20
C----------------------------------------------- 

c***    Local Vbls   ***
  
c***    Candidates for New Include File  ***
 
      namelist /namin18/rx1, rx2, zy1, zy2, nrgrid, nzgrid, condif
 
****  ----- Read the LIMPOS.GRID namelist (/NAMIN18/) file to override -----
 
*     Use defacto values rather than 38 values
!      nrgrid = nrgrid_18
!      nzgrid = nzgrid_18
!      rx1 = rx1_18
!      rx2 = rx2_18
!      zy1 = zy1_18
!      zy2 = zy2_18
 
*     Save MGRID values of plotting specs
      nrgridsv = nrgrid
      nzgridsv = nzgrid
      rx1sv = rx1
      rx2sv = rx2
      zy1sv = zy1
      zy2sv = zy2
      condifsv = condif
 
      ios = 0
      if (initnml .eq. 0) open(unit=50, file='limpos.grid', 
     1    status='old', iostat=ios)
      if (ios .eq. 0) then
         write (*, *) 'Overriding r/z ', 
     1      'grid parameters with local LIMPOS.GRID file for plot 18'
         read (50, namin18)
         initnml = 1
      endif
 
      call cklimpos (ios)
      if (ios .ne. 0) return 
      call gengrid
 
c***        --------  Contour Plot Grid Definition  --------
 
c     Generate set of "Green's" Function arrays
c     for the PF Coil set, and
c     determine the size and resolution of the (R,Z) grid used to
c     generate the contour plot.
 
      allocate (pflcoilgr(nrgrid,nzgrid,nsets), pfltot(nrgrid,nzgrid), 
     1          stat=iostat)
      if (iostat .ne. 0) stop 'pflcoilgr allocation error in pfluxeng'
      pflcoilgr = 0.
      
      call gpfcoil

      allocate (dbsubsdu(mnmax*ns), dbsubuds(mnmax*ns), taua(nrt),
     1  dummy1(mnmax*ns), dummy2(nrt), bsubsdu(nrt), bsubuds(nrt))
 
 
c***        --------  Plasma Current Pre-Calculations --------
 
      do kt = 1, nthpts
         do j = 2, ns
            l = kt + nthpts*(j - 1)
            taua(l) = gsqrt(l)/r12(l)
         end do
      end do
 
c*     * ****************************
c*     * Compute Fourier coefficients of d Bu/ds and d Bs/du
c*     * on integer radial grid
c*     * ****************************
 
      hs = 1.0/(ns - 1)
      do js = 2, ns - 1
         lj = (js - 1)*mnmax
         do mn = 1, mnmax
            l = mn + lj
            if (mod(ixm(mn),2) .eq. 0) then
               dbsubsdu(l) = .5*xm(mn)*(bsubsmn(mn,js)+bsubsmn(mn,js+1))
               dbsubuds(l) = (bsubumn(mn,js+1)-bsubumn(mn,js))/hs
            else
               t1 = sqrt((js - 1.)/(js - 1.5))
               t2 = sqrt((js - 1.)/(js - 0.5))
               dbsubsdu(l) = .5*xm(mn)*(t1*bsubsmn(mn,js)
     1                                 +t2*bsubsmn(mn,js+1))
               dbsubuds(l) = (t2*bsubumn(mn,js+1)
     1                       -t1*bsubumn(mn,js))/hs + 0.25*(t1*
     1            bsubumn(mn,js)+t2*bsubumn(mn,js))/(hs*(js - 1.))
            endif
         end do
      end do
 
      call totz (nthpts, ns, 1, 1, bsubsdu, dummy2, 
     1   dbsubsdu, dummy1, dummy1, dummy1)
      call totz (nthpts, ns, 1, 1, bsubuds, dummy2, 
     1   dbsubuds, dummy1, dummy1, dummy1)
 
c***        -------- Calculate  & Plot Components of Poloidal Flux --------
 
c     Assume extcur are in the order specified below ...
c          {set1,set2,...,setn} - order determined by LIMPOS_FILE
c          {coh,cef,cvc,chf - for TFTR}
c     Offset EXTCUR by 1 to discard the TF coil current
c     Calculate components of poloidal flux
c     Itype=1 : use j-dot-gradphi ~= dbsubuds - dbsubsdu
c     Itype=2 : use j-dot-gradphi ~= ppsi(1-ar2) + 2pi*dI/dV*ar2
c     Itype=3 : use j-dot-gradphi ~= dbsubuds
c     Itype=4 : use j-dot-gradphi ~= ppsi + f*fprime/(u0*R**2)
 
      pflhead(1) = 'POLOIDAL FLUX - 1$'
      pflhead(2) = 'POLOIDAL FLUX - 2$'
      pflhead(3) = 'POLOIDAL FLUX - 3$'
      pflhead(4) = 'POLOIDAL FLUX - 4$'
 
      itype = 4
      call getpfltot (itype, ns, r, z, r12, z12, ru, zu, ru12, zu12, 
     1   pres, iotas, phip, gsqrt, taua, extcur(2), dbsubudstz,
     2   bsubsdu, bsubuds, bsubutz)
      call plot_pflux (real(pfltot), pflhead(itype))   !total
      call frame
 
 
      deallocate (dbsubsdu, dbsubuds, dummy1, dummy2, 
     1   bsubsdu, bsubuds, taua, pflcoilgr, pfltot)

      end subroutine pfluxeng

!-------------------------------------------------------------------------
!
!     "Generate the Green's Function for the PFC set"
!-------------------------------------------------------------------------
 
      subroutine gpfcoil
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Vpname3
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
 
      ipset(1) = 1.0
      ipset(2) = 0.0
      ipset(3) = 0.0
      ipset(4) = 0.0
 
      call lpfcoil (1)
 
      ipset(1) = 0.0
      ipset(2) = 1.0
      ipset(3) = 0.0
      ipset(4) = 0.0
 
      call lpfcoil (2)
 
      ipset(1) = 0.0
      ipset(2) = 0.0
      ipset(3) = 1.0
      ipset(4) = 0.0
 
      call lpfcoil (3)
 
      ipset(1) = 0.0
      ipset(2) = 0.0
      ipset(3) = 0.0
      ipset(4) = 1.0
 
      call lpfcoil (4)
 
      end subroutine gpfcoil
 
!-------------------------------------------------------------------------
!
!     Calc Poloidal Flux - looping over (R,Z) grid
!-------------------------------------------------------------------------
 
      subroutine lpfcoil(index)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nrgrid, nzgrid
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer index
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l1, l2
C-----------------------------------------------
 
      do l2 = 1, nzgrid
         do l1 = 1, nrgrid
            call xpfcoil (rgrid(l1), zgrid(l2), pflcoilgr(l1,l2,index))
         end do
      end do
 
      end subroutine lpfcoil


      subroutine xpfcoil(rr, zz, flux)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nsets
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: rr, zz, flux
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: f
C-----------------------------------------------
 
      flux = 0.
 
      do i = 1, nsets
         call loops (i, rr, zz, f)
         flux = flux + f
      end do
 
      end subroutine xpfcoil
      

      subroutine loops(ind, r, z, f)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: pfcspec, nsetsn
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ind
      real(rprec) :: r, z, f
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ini = 0, nw, iw
      real(rprec) :: pi, c, cw, df
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec), EXTERNAL :: flux
C-----------------------------------------------
c
      if (ini .eq. 0) then
         pi = acos(-1.d0)
         c = 4*pi*1E-7                           !unit current
         ini = 1
      endif
 
      nw = nsetsn(ind)
      cw = ipset(ind)
      f = 0.
 
      if (nw<1 .or. cw==0.) return 
      do iw = 1, nw
         if (pfcspec(1,iw,ind) <= 0.) exit 
         if (pfcspec(3,iw,ind) .ne. 0.) then
            df = (-2)*pi*flux(c,pfcspec(2,iw,ind),
     1                          pfcspec(1,iw,ind),z,r)
            f = f + pfcspec(3,iw,ind)*pfcspec(4,iw,ind)*df
c       ! times # of turns times cfac
         endif
      end do
      f = cw*f

      end subroutine loops
      

      subroutine loopscw(cw, nw, ind, r, z, f)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: pfcspec
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nw, ind
      real(rprec) :: cw, r, z, f
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ini = 0, iw
      real(rprec) :: pi, c, df

      save ini, pi, c
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: flux
C-----------------------------------------------
c
      if (ini .eq. 0) then
         pi = acos(-1.d0)
         c = 4*pi*1E-7                           !unit current
         ini = 1
      endif
 
      f = 0.
 
      if (nw<1 .or. cw==0.) return 
      do iw = 1, nw
         if (pfcspec(1,iw,ind) <= 0.) exit 
         if (pfcspec(3,iw,ind) .ne. 0.) then
            df = (-2.)*pi*flux(c,pfcspec(2,iw,ind),
     1           pfcspec(1,iw,ind),z,r)
            f = f + pfcspec(3,iw,ind)*pfcspec(4,iw,ind)*df
c       ! times # of turns times cfac
         endif
      end do
      f = cw*f

      end subroutine loopscw
    

      function flux (c, xc, yc, x, y)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) c, xc, yc, x, y
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: k, m1, m1s, d, e, t, sks, al, a0, a1, a2,
     1   b0, b1, b2, c1, c2, d1, d2, pi, flux
C-----------------------------------------------
 
      data a0, a1, a2/1.3862944, .1119723, .0725296/
      data b0, b1, b2/.5, .1213478, .0288729/
      data c1, c2/.4630151, .1077812/
      data d1, d2/.2452727, .0412496/
      data pi/3.1415926535898/
 
      d = (y + yc)**2 + (x - xc)**2
      t = c*yc/(pi*sqrt(d))
      sks = 4.*yc*y/d
      m1 = 1.d0 - sks
      m1s = m1**2
      al = log(1.d0/m1)
      k = a0 + a1*m1 + a2*m1s + (b0 + b1*m1 + b2*m1s)*al
      e = 1. + c1*m1 + c2*m1s + (d1*m1 + d2*m1s)*al
      flux = -t*((2. - sks)*k - 2.*e)/sks*y

      end function flux


      subroutine fcp(pcur, r, z, ntheta, ns, rext, zext, psiext, okay)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ntheta, ns
      real(rprec) :: rext, zext, psiext
      logical okay
      real(rprec), dimension(ntheta,*) :: r, z, pcur
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real, parameter :: u0 = 12.566371E-07
      real, parameter :: stata = 3.0E+09
      real, parameter :: hx = 0.01
      real, parameter :: hy = 0.01
      real, parameter :: xwh = 0.0
      real, parameter :: closex = .05
      real, parameter :: closey = .05
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ii, is, j, jj, ini = 0
      real :: area, aip, ajloc, cur, da, ds, dth, fac, xtau, xminus, 
     1   xplus, yplus, yminus, rloc, zloc, x, y, xc, yc, pi, flux, c, k
     2   , m1, m1s, d, e, t, sks, al, a0, a1, a2, b0, b1, b2, c1, c2, d1
     3   , d2
C-----------------------------------------------
 
      data c/1.0/
      data a0, a1, a2/1.3862944, .1119723, .0725296/
      data b0, b1, b2/.5, .1213478, .0288729/
      data c1, c2/.4630151, .1077812/
      data d1, d2/.2452727, .0412496/
c
c      if (ini.eq.0) then
      pi = acos(-1.0)
      dth = 2*pi/(ntheta - 1)
      ds = 1./(ns - 1)
      okay = .TRUE.
      ini = 1
c      endif
c
      psiext = 0.
      do j = 2, ns
         do i = 1, ntheta - 1
            c = 1.0
            xc = z(i,j)
            yc = r(i,j)
            x = zext
            y = rext
c           real function flux(c,xc,yc,x,y,ier) - in-line to save time
            d = (y + yc)**2 + (x - xc)**2
            t = c*yc/(pi*sqrt(d))
            sks = 4.*yc*y/d
            m1 = 1. - sks
            if (sks*d*m1 <= 0.) m1 = 4.0978193E-08
            m1s = m1**2
            al = alog(1./m1)
            k = a0 + a1*m1 + a2*m1s + (b0 + b1*m1 + b2*m1s)*al
            e = 1. + c1*m1 + c2*m1s + (d1*m1 + d2*m1s)*al
            flux = -t*((2. - sks)*k - 2.*e)/sks*y
            psiext = psiext + flux*pcur(i,j)
         end do                                  ! enddo i
      end do                                     ! enddo j
      psiext = psiext*2*pi*ds*dth*u0

      end subroutine fcp
      
 
      subroutine plot_pflux(pflcomp, pfhead)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: rx1, rx2, zy1, zy2, condif,
     1  nlim, limitr, nrgrid, nzgrid, ns  
      use Vpname2
      use Vpname3
      USE Vrzarray
      USE Vindat2
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character pfhead*(*)
      real, dimension(*) :: pflcomp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l2, l1, l, j, jes
      real :: mrmnc, mzmns, sizex = 8.0, sizey = 10.0, rx0 = 0.0, 
     1   zy0 = -4.0, rxst = 0.0, zyst = 0.1, rxstp = 0.5, 
     2   ticksz = 0.02, tk = 6.0, width = 10.
      real, dimension(4) :: grid
      real :: gmin, gmax, gwnd, gscale
C-----------------------------------------------

      call agseti ('FRAME.', 2)
      gmin = 1.E6
      gmax = -gmin
      do l2 = 1, nzgrid
         gmax = maxval(pflcomp((l2-1)*nrgrid+1:l2*nrgrid))
         gmin = minval(pflcomp((l2-1)*nrgrid+1:l2*nrgrid))
      end do
 
*
*       Disable the RUN Info Line
*
      call agsetc ('LABEL/NAME.', 'B')
      call agseti ('LINE/NUMBER.', -100)
      call agsetc ('LINE/TEXT.', 'R')
      call agsetc ('LABEL/NAME.', 'T')
      call agsetr ('LABEL/SUPPRESSION.', -1.)
      call agsetc ('LABEL/NAME.', 'L')
      call agseti ('LINE/NUMBER.', 100)
      call agsetc ('LINE/TEXT.', 'Z')
      call agseti ('BACKGROUND TYPE.', 1)
      call agsetp ('GRAPH WINDOW.', gwnd, 4)
      call agseti ('LABEL/CONTROL.', 2)
      call agsetf ('X/NICE.', 0.)
      call agsetf ('Y/NICE.', 0.)
 
      gscale = (gmax - gmin)/condif
 
 
      call contgraf2 (pflcomp, nrgrid, nzgrid, real(rx1), 
     1   real(rx2), real(zy1), real(zy2), 1, -2, gscale, 0, real(rgrid), 
     2   real(zgrid), 2, 0, pfhead, 'R$', 'Z$', runlabel, 6, 6, grid)
 
*
*       Plot Limiter Boundary: Allow for more than one limiter set
*
      call set(grid(1),grid(2),grid(3),grid(4),real(rx1),real(rx2),
     1     real(zy1), real(zy2),1)
      call gsln (2)                              !dashed
      call gslwsc (width)
      call gsplci (44)
 
      do l = 1, nlim
         if (limitr(l) > 1) call curved(real(xlim(1,l)), 
     1       real(ylim(1,l)),limitr(l))
      end do
 
c     Plot last closed flux surface
      call gslwsc (width/2.)
      call gsln (1)                              !solid
      do j = ns, ns
         jes = nthpts*(j - 1)
         l = 1 + jes
         call agcurv (real(r(l)), 1, real(z(l)), 1, nthpts, 1)
      end do
 
      call gslwsc (1.)
      call gsln (1)
 
*
*       Plot Coil Positions
*
      call plcoil
 
      end subroutine plot_pflux
      
 
      subroutine getpfltot(itype, ns, r, z, r12, z12, ru, zu, ru12, 
     1   zu12, pres, iotas, phip, gsqrt, taua, extcur, dbsubudstz, 
     2   bsubsdu, bsubuds, bsubutz)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nsets, nrgrid, nzgrid
      use Vpname1
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer itype, ns
      real(rprec), dimension(*) :: r, z, r12, z12, ru, zu, ru12, 
     1   zu12, pres, iotas, phip, gsqrt, taua, extcur, dbsubudstz, 
     2   bsubsdu, bsubuds, bsubutz
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real, parameter :: u0 = 12.566371E-07
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer::i,j,kr,kt,kz,l,lt,l0,l1,l2,l3,nrz,iii,ini,idbg
      real(rprec), dimension(:), allocatable :: bsubudszc, 
     1  gsqrtzb, jfb, 
     1  jfb1, jfb2, jphiffp1zc, jphiffp2zc, arsq, jfbzc, jfb1zc, 
     2  jfb2zc, ip9zb, jgradphi, jgradphizc, jphiffpzc
      real(rprec) :: sg, fac, emax, pfx, vol, pi, ds, dth, area
      real(rprec) :: twopi
      real(rprec) :: ipchk1, ipchk2
      real(rprec) :: vchk
      real(rprec) :: ipchk
      real(rprec) :: aipffpx, aipffp1x, aipffp2x, aip3x, aip31x, 
     1  aip32x
      real(rprec) :: gsqrte, psicur, psicoil, pflcache
      real(rprec), dimension(:,:), allocatable :: pflcur, pflcoil
      real(rprec), dimension(ns) :: aip, aip2, aip3, aip31, aip32
      real(rprec), dimension(ns) :: pprime, iprime
      real(rprec), dimension(ns) :: pprimezc, psiprimezc, fzc, fpzc
      real(rprec), dimension(ns) :: vprime
      real(rprec), dimension(ns) :: aguugsq, aguugsqzc, is, 
     1   iprime2, ispzc, iszc, ispzb
      real(rprec), dimension(ns) :: agsqrt, agsqrtzb, arm2gsqrt, 
     1   arm2gsqrtzb, absubutz, absubutzzb, iotazb
      real(rprec), dimension(ns) :: f, fp
      real(rprec), dimension(ns) :: aipffp, w, aipffp1
      real(rprec), dimension(ns) :: aipffp2
      real(rprec), dimension(ns) :: aarsq, szc, sp1, sp2, sp3
      real(rprec), dimension(ns) :: aip9zb, aip9zc
      logical :: okay
      character, dimension(4) :: jtype*35

!     save bsubudszc, jphiffpzc, jfbzc, jgradphizc, ini, idbg, jtype
C-----------------------------------------------
 
c                              {ctf,coh,cef,cvc,chf - for TFTR}
 
      data sg/ - 1./
      data fac/1000./
 
      data ini/0/
      data idbg/0/
 
      allocate (bsubudszc(nrt), gsqrtzb(nrt), jfb(nrt), jfb1(nrt), 
     1  jfb2(nrt), jphiffp1zc(nrt), jphiffp2zc(nrt), arsq(nrt), 
     2  jfbzc(nrt), jfb1zc(nrt), jfb2zc(nrt), ip9zb(nrt), jgradphi(nrt),
     3  jgradphizc(nrt), jphiffpzc(nrt))
 
 
*-  Imitiation:
 
      if (ini .eq. 0) then
         ini = 1
      else
         go to 11000
      endif

 
*-  Calculate poloidal field due to external coil set
*-  The minus sign was put in ad hoc to get agreement with pflset
*-  values used previously. Why not use pflset in makepflext?
 
c   Auto-scale the PF Coil Currents into Amps
 
      emax = abs(extcur(1))
      emax = max(emax,maxval(abs(extcur(2:nsets))))
      if (emax > 1000.) fac = 1.0                ! currents are in Amps
 
      do l2 = 1, nzgrid
         do l1 = 1, nrgrid
            pfx = 0.
            do l0 = 1, nsets
c                                                ! curr in AMPS
               pfx = pfx + sg*extcur(l0)*pflcoilgr(l1, l2, l0)*fac
            end do
            pflcoil(l1, l2) = pfx
         end do
      end do
 
 
 
*-  Calculate poloidal field due to plasma current
*-  Try many different approaches
 
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
c   calculate j-dot-gradphi using force balance eqtn
c   construct it on the full mesh; then interp back onto the half mesh
c   in order to send a zone centered current density to FCP
 
 
      pi = acos(-1.0)
      twopi = 2*pi
      dth = 2*pi/(nthpts - 1)
      ds = 1./(ns - 1)
 
c   calculate p-prime on half mesh
c                                                !zb
      pprime(2:ns-1) = (pres(3:ns)-pres(2:ns-1))/ds/u0
      pprime(1) = 2.*pprime(2) - pprime(3)
      pprime(ns) = 2.*pprime(ns-1) - pprime(ns-2)
      pprimezc(3:ns-1) = 0.5*(pprime(3:ns-1)+pprime(2:ns-2))
      pprimezc(2) = 1.5*pprime(2) - 0.5*pprime(3)
      pprimezc(ns) = 1.5*pprime(ns-1) - 0.5*pprime(ns-2)
 
c   calculate psi-prime as iota(s)*phi-prime on zc
c   do it directly on the zc grid
c   double check this
      ipchk2 = 0
      psiprimezc(2:ns) = iotas(2:ns)*phip(2:ns)*twopi
      ipchk2 = ipchk2 + sum(psiprimezc(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) '    Check the Psi : Psi = ', ipchk2
 
c   calculate <gsqrt*R**-2> on the half mesh
c     arm2gsqrt is zc
      do j = 2, ns
         arm2gsqrt(j) = 0.
         lt = (j - 1)*nthpts
         arm2gsqrt(j) = arm2gsqrt(j) - sum((gsqrt(lt+1:nthpts-1+lt)/
     1      r12(lt+1:nthpts-1+lt)**2)*dth/twopi)     !gsqrt wrong sign
      end do
 
c     construct F-Fprime on half grid
      fzc(1) = 0.
c                                                !zc
      fzc(2:ns) = phip(2:ns)*twopi/twopi/arm2gsqrt(2:ns)
      fpzc(1) = 0.
      fp(2:ns-1) = phip(2:ns-1)*twopi/twopi*(1./arm2gsqrt(3:ns)-1./
     1   arm2gsqrt(2:ns-1))/ds                   !zb
      fpzc(3:ns-1) = 0.5*(fp(3:ns-1)+fp(2:ns-2)) !zc
      fpzc(2) = 1.5*fp(2) - 0.5*fp(3)
      fpzc(ns) = 1.5*fp(ns-1) - 0.5*fp(ns-2)
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
c     Calculate j-dot-gradphi from p-psi + ff-prime

c     now calculate it directly on the half grid
      do j = 2, ns
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            jphiffp1zc(l) = sg*pprimezc(j)/psiprimezc(j)*twopi
            jphiffp2zc(l) = sg*fzc(j)*fpzc(j)/psiprimezc(j)*twopi/u0/r12
     1         (l)**2
            jphiffpzc(l) = jphiffp1zc(l) + jphiffp2zc(l)
         end do
      end do
 
      aipffpx = 0.
      aipffp1x = 0.
      aipffp2x = 0.
      aipffp(1) = 0.
      aipffp1(1) = 0.
      aipffp2(1) = 0.
      vol = 0.
      do j = 2, ns
         do kt = 1, nthpts - 1
            l = kt + nthpts*(j - 1)
            vol = vol - gsqrt(l)*ds*dth*twopi
            aipffp1x = aipffp1x - jphiffp1zc(l)*gsqrt(l)*ds*dth
            aipffp2x = aipffp2x - jphiffp2zc(l)*gsqrt(l)*ds*dth
            aipffpx = aipffpx - jphiffpzc(l)*gsqrt(l)*ds*dth
         end do
         aipffp1(j) = aipffp1x
         aipffp2(j) = aipffp2x
         aipffp(j) = aipffpx
      end do
      if (idbg .eq. 1) then
         write (*, *) 
     1      '    Sanity Check using p-psi + ffprime: Ip[tot]= ', aipffp(
     2      ns)
         write (*, *) '    Ip[ppsi]= ', aipffp1(ns), 'Ip[ffp]= ', 
     1      aipffp2(ns)
         write (*, *) '    Vol= ', vol
      endif
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
c   calculate I-prime by differentiating the phiprime-iota version  onto zc
      do j = 2, ns
         aguugsqzc(j) = 0.
         do kt = 1, nthpts - 1
            l = kt + (j - 1)*nthpts
            aguugsqzc(j) = aguugsqzc(j) - (ru12(l)**2+zu12(l)**2)/gsqrt(
     1         l)*dth/twopi                      !gsqrt wrong sign
         end do
      end do
 
c   Is - zc
      iszc(1) = 0.
      szc(1) = 0.
      do j = 2, ns
         iszc(j) = phip(j)*twopi/u0*iotas(j)*aguugsqzc(j)
         szc(j) = 1./real((ns - 1)*(j - 1.5))
      end do
      call splaan (ns, szc, iszc, sp1, sp2, sp3)
      ispzc(2:ns) = sp1(2:ns)
c      do j = 2,ns
c       write(*,*) szc(j),iszc(j),ispzc(j),iotas(j),aguugsqzc(j)
c      end do
c   double check this spline derivative
      if (idbg .eq. 1) 
     1   write (*, *) '    Check the Ip from is(ns): Ip = ',
     2         1.5*iszc(ns) - 0.5*iszc(ns-1)
      ipchk = 0
      ipchk = sum(ispzc(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the Ip from splined Ispzc: Ip = ', ipchk
  
c   Is - zc
!      iszc(1) = 0.
!      do j =2,ns
!       iszc(j) = phip(j)*twopi/u0*iotas(j)*0.5*(aguugsq(j)+aguugsq(j-1))
!      end do
 
c   Iprime1 - zc
      ispzb(2:ns-1) = (iszc(3:ns)-iszc(2:ns-1))/ds
cwie      ispzb(1) = 2.*ispzb(2) - ispzb(3)
      ispzb(1) = 0.
      ispzb(ns) = 2.*ispzb(ns-1) - ispzb(ns-2)
c                                                !zc
      ispzc(3:ns-1) = 0.5*(ispzb(3:ns-1)+ispzb(2:ns-2))
cwie      ispzc(2) = 1.5*ispzb(2) - 0.5*ispzb(3)
      ispzc(2) = .5*(ispzb(1)+ispzb(2))
      ispzc(ns) = 1.5*ispzb(ns-1) - 0.5*ispzb(ns-2)
 
c   double check this
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the Ip from is(ns): Ip = ',
     2   1.5*iszc(ns) - 0.5*iszc(ns-1)
      ipchk = 0
      ipchk = sum(ispzc(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1 '    Check the Ip from Ispzc: Ip = ', ipchk
 
c   try the trapezoidal rule on iprime
      w(2:ns-1) = 1.
      w(1) = 0.5
      w(ns) = 0.5
      ipchk = 0.
      ipchk = sum(w(:ns)*ispzb(:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the Ip from TraPEZOIDAL Ispzb: Ip = ', ipchk
 
c   try the simpsons rule on iprime
      w(2:ns-1:2) = 4.
      w(3:ns-2:2) = 2.
      w(1) = 1.
      w(ns) = 1.
      ipchk = 0.
      ipchk = sum(w(:ns)*ispzb(:ns)*ds/3.)
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the Ip from Simpsons Ispzb: Ip = ', ipchk
 
 
c   calculate <gsqrt> on the half mesh
      do j = 2, ns
         agsqrt(j) = 0.
         lt = (j - 1)*nthpts
c                                                !gsqrt wrong sign
         agsqrt(j) = agsqrt(j) - sum(gsqrt(lt+1:nthpts-1+lt)*dth/twopi)
      end do
 
c   calculate dV/ds on the half  mesh
 
c   double check this
      vchk = 0
      vprime(2:ns) = twopi**2*agsqrt(2:ns)
      vchk = vchk + sum(vprime(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the gsqrt part of Iprime: Vol = ', vchk
 
 
c    form arsq .eq. r**-2/<<r**-2>> on the half mesh
      do j = 2, ns
         lt = (j - 1)*nthpts
         arsq(lt+1:nthpts+lt) = 1./r12(lt+1:nthpts+lt)**2/(arm2gsqrt(j)/
     1      agsqrt(j))
      end do
 
c    Let us check it by multiplying by gsqrt  ...  should be .eq. 1
      do j = 2, ns
         aarsq(j) = 0.
         lt = (j - 1)*nthpts
         aarsq(j) = aarsq(j) - sum(arsq(lt+1:nthpts-1+lt)*gsqrt(lt+1:
     1      nthpts-1+lt)*twopi/vprime(j)*dth)
      end do
      if (idbg .eq. 1) write (*, *) 
     1   '    Sanity Check on 2pi*gsqrt*R2/(dV/ds) = ', aarsq(ns)
 
 
 
 
c     Let us calculate dI/ds by doing the ds 1st, and then the dtheta avg
 
      do j = 2, ns - 1                       !zc function diffed onto zb
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            ip9zb(l) = (iotas(j)*(ru12(l)**2+zu12(l)**2)/gsqrt(l)-iotas(
     1         j+1)*(ru12(l+nthpts)**2+zu12(l+nthpts)**2)/gsqrt(l+nthpts
     2         ))/ds*phip(j)/u0
         end do
      end do
      do kt = 1, nthpts
         l = kt
cwie      ip9zb(l) = 2.0*ip9zb(l+nthpts) - ip9zb(l+2*nthpts)
         ip9zb(l) = 0.5*ip9zb(l+nthpts)
      end do
      do kt = 1, nthpts
         l = kt + nthpts*(ns - 1)
         ip9zb(l) = 2.0*ip9zb(l-nthpts) - ip9zb(l-2*nthpts)
      end do
cwie      ip9zb(l) = 2.0*ip9zb(l+nthpts) - ip9zb(l+2*nthpts)
      do j = 1, ns
         aip9zb(j) = 0.
         do kt = 1, nthpts - 1
            l = kt + nthpts*(j - 1)
            aip9zb(j) = aip9zb(j) + ip9zb(l)*dth
         end do
      end do
 
c   double check this
      ipchk = 0.
      aip9zc(2:ns) = 0.5*(aip9zb(:ns-1)+aip9zb(2:ns))
      ipchk = ipchk + sum(aip9zc(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1 '    Check the Ip from aIp9zc: Ip = ', ipchk
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
c     Calculate j-dot-gradphi from p-psi + dI/dV
 
      do j = 2, ns                               ! directly on zc
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            jfb1zc(l) = twopi*pprimezc(j)/psiprimezc(j)*(1. - arsq(l))
            jfb2zc(l) = twopi*aip9zc(j)/vprime(j)*arsq(l)
c          jfb2zc(l) = twopi*ispzc(j)/vprime(j)*arsq(l)
            jfbzc(l) = jfb1zc(l) + jfb2zc(l)
         end do
      end do
 
c     Sanity Check #3
      aip3(1) = 0.
      aip31(1) = 0.
      aip32(1) = 0.
      aip3x = 0.
      aip31x = 0.
      aip32x = 0.
      do j = 2, ns
         do kt = 1, nthpts - 1
            l = kt + nthpts*(j - 1)
            aip3x = aip3x - jfbzc(l)*gsqrt(l)*ds*dth
            aip31x = aip31x - jfb1zc(l)*gsqrt(l)*ds*dth
            aip32x = aip32x - jfb2zc(l)*gsqrt(l)*ds*dth
         end do
         aip3(j) = aip3x
         aip31(j) = aip31x
         aip32(j) = aip32x
      end do
      if (idbg .eq. 1) then
         write (*, *) '    Sanity Check using ppsi + dI/dV: ', 
     1      ' Ip[tot]= ', aip3(ns)
         write (*, *) '    Ip[ppsi]= ', aip31(ns), ' Ip[dI/dV]= ', aip32
     1      (ns)
      endif
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
 
c   calculate I-prime by differentiating <Bu> onto the full mesh
c     absubutz - zc
      do j = 2, ns
         absubutz(j) = 0.
         lt = (j - 1)*nthpts
         absubutz(j) = absubutz(j) + sum(bsubutz(lt+1:nthpts-1+lt)*dth/
     1      twopi)
      end do
c     absubutzzb - zb
      absubutzzb(2:ns-1) = 0.5*(absubutz(2:ns-1)+absubutz(3:ns))
      absubutzzb(1) = 2.*absubutzzb(2) - absubutzzb(3)
      absubutzzb(ns) = 2.*absubutzzb(ns-1) - absubutzzb(ns-2)
c     iprime - zc
c   double check this
      ipchk = 0
      iprime(2:ns) = (absubutzzb(2:ns)-absubutzzb(:ns-1))/ds*twopi/u0
      ipchk = ipchk + sum(iprime(2:ns)*ds)
      if (idbg .eq. 1) write (*, *) 
     1   '    Check the Ip from iprime: Ip = ', ipchk
      ipchk1 = twopi/u0*(1.5*absubutz(ns)-0.5*absubutz(ns-1))
      ipchk2 = twopi/u0*absubutzzb(ns)
      if (idbg .eq. 1) then
         write (*, *) '    Check the Ip from <Bu>-prime: Ip = ', ipchk1
         write (*, *) '    Check the Ip from <Bu(ns)>: Ip = ', ipchk2
      endif
 
c     Alternatively, calculate J from Bu
c     b*,r on full mesh;avg gsqr onto full mesh;jgradphi on full mesh
 
c     >> interior values
      do j = 2, ns - 1
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            jgradphi(l) = .5*(1./gsqrt(l)+1./gsqrt(l+nthpts))*
     1         (bsubuds(l)-bsubsdu(l))
         end do
      end do
 
c     >> boundary values
      do kt = 1, nthpts
         l = kt + nthpts*(ns - 1)
         bsubuds(l) = 2.0*bsubuds(l-nthpts) - bsubuds(l-2*nthpts)
         bsubsdu(l) = 2.0*bsubsdu(l-nthpts) - bsubsdu(l-2*nthpts)
      end do
      do kt = 1, nthpts
         l = kt + nthpts*(ns - 1)
         gsqrte = 1.5*gsqrt(l) - 0.5*gsqrt(l-nthpts)
         jgradphi(l) = gsqrte*(bsubuds(l)-bsubsdu(l))
      end do
      do kt = 1, nthpts
         jgradphi(kt) = 2.*jgradphi(kt+nthpts) - jgradphi(kt+2*nthpts)
      end do
 
c     >> boundary values
      do j = 2, ns
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            jgradphizc(l) = .5/u0*(jgradphi(l)+jgradphi(l-nthpts))
         end do
      end do
 
c     Sanity Check #1
      pi = acos(-1.0)
      dth = 2*pi/(nthpts - 1)
      ds = 1./(ns - 1)
      aip(1) = 0.
      area = 0.
      do j = 2, ns
         do kt = 1, nthpts - 1
            l = kt + nthpts*(j - 1)
            area = area + abs(taua(l))*ds*dth
            aip(j) = aip(j-1) - jgradphizc(l)*gsqrt(l)*ds*dth
         end do
      end do
      if (idbg .eq. 1) write (*, *) 
     1   '    Sanity Check using jgradphizc: area= ', area, ' Ip= ', aip
     2   (ns)
 
c * - * - * * - * - * * - * - * * - * - * * - * - * * - * - * * - * - *
 
c     Sanity Check #2 (also another way to do the Green-s Fnc)
c     by calculating dbsubudstz on half mesh
      aip2(1) = 0.
      do j = 2, ns
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            bsubudszc(l)=sg*0.5/u0*(dbsubudstz(l)+dbsubudstz(l-nthpts))
         end do
      end do
      do j = 2, ns
         do kt = 1, nthpts - 1
            l = kt + nthpts*(j - 1)
            aip2(j) = aip2(j-1) + bsubudszc(l)*ds*dth
         end do
      end do
      if (idbg .eq. 1) write (*, *) 
     1   '    Sanity Check using just bsubu: Ip= ', aip2(ns)
 
c   Send the part in brackets to fcp ... Jda=[j*gsqrt]dsdu
      do j = 1, ns
         do kt = 1, nthpts
            l = kt + nthpts*(j - 1)
            jgradphizc(l) = jgradphizc(l)*abs(gsqrt(l))
            jfbzc(l) = jfbzc(l)*abs(gsqrt(l))
            jphiffpzc(l) = jphiffpzc(l)*abs(gsqrt(l))
            jfb1zc(l) = jfb1zc(l)*abs(gsqrt(l))
            jfb2zc(l) = jfb2zc(l)*abs(gsqrt(l))
         end do
      end do
 
c     Check Psi calc at flux loop locations
 
      jtype(1) = '(dBu/ds - dBs/du)'
      jtype(2) = '(dp/dpsi*(1-arsq) + dI/dV*arsq)'
      jtype(3) = '(dBu/ds)'
      jtype(4) = '(dp/dpsi + F*dF/dpsi/(u0*r**2))'
 
      if (itype==2 .and. idbg==1) call fluxcheck (jfbzc, ns, r12, z12, 
     1   extcur, sg, fac, jtype(itype))
      if (itype==4 .and. idbg==1) call fluxcheck (jphiffpzc, ns, r12, 
     1   z12, extcur, sg, fac, jtype(itype))
 
11000 continue
 
c     Calculate PF due to Plasma on R,Z Grid in UHP
 
      write (*, 2900) '  Generating PF plots on nr by ny (', nrgrid, 
     1   nzgrid, ') grid ... '
      write (*, *) '    using J type : ', jtype(itype)
 2900 format(a,i3,',',i3,a)
 
      allocate (pflcur(nrgrid,nzgrid), pflcoil(nrgrid,nzgrid), stat=i)
      if (i .ne. 0) stop 'pflcur allocation error in getpfltot'
      
      do l2 = nzgrid/2 + 1, nzgrid
         l3 = nzgrid - l2 + 1
         do l1 = 1, nrgrid
            if (itype .eq. 1) call fcp (jgradphizc, r12, z12, nthpts, 
     1         ns, rgrid(l1), zgrid(l2), psicur, okay)
            if (itype .eq. 2) call fcp (jfbzc, r12, z12, nthpts, ns, 
     1         rgrid(l1), zgrid(l2), psicur, okay)
            if (itype .eq. 3) call fcp (bsubudszc, r12, z12, nthpts, 
     1         ns, rgrid(l1), zgrid(l2), psicur, okay)
            if (itype .eq. 4) call fcp (jphiffpzc, r12, z12, nthpts, 
     1         ns, rgrid(l1), zgrid(l2), psicur, okay)
            pflcur(l1,l2) = psicur
            pflcur(l1,l3) = psicur
         end do
      end do
 
c     Total external poloidal flux = flux from plasma current +
c     flux from external PF coil set
      do l2 = 1, nzgrid
         pfltot(:nrgrid,l2) = pflcur(:nrgrid,l2) + pflcoil(:nrgrid,l2)
      end do
 
 
      if (allocated(bsubudszc))
     1  deallocate (bsubudszc, gsqrtzb, jfb, 
     2  jfb1, jfb2, jphiffp1zc, jphiffp2zc, arsq, jfbzc, jfb1zc, 
     3  jfb2zc, ip9zb, jgradphi, jgradphizc, jphiffpzc, pflcur, 
     4  pflcoil)

      end subroutine getpfltot
      

      subroutine fluxcheck(jphi, ns, r12, z12, extcur, sg, fac, jtype)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nsets, nsetsn
      use Vpname1
      use Vpname3
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns
      real(rprec) :: sg, fac
      character jtype*(*)
      real(rprec), dimension(nthpts,*) :: jphi, r12, z12
      real(rprec), dimension(*) :: extcur
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l0, nl, nw
      real(rprec), dimension(6) :: floopr, floopz, psitestpf, 
     1   psitestip      
      real(rprec) :: cw, pfx
      logical :: okay
C-----------------------------------------------
 
!R flux loop positions
      data floopr/1.709, 2.293, 3.693, 3.693, 2.293, 1.709/
!Z flux loop position
      data floopz/0.762, 1.187, 0.671, -0.671, -1.187,  - 0.762/
 
      write (*, *) 'Flux loop test for : ', jtype
      do nl = 1, 6
         psitestpf(nl) = 0.
         do l0 = 1, nsets
            cw = 1.0
            nw = nsetsn(l0)
            call loopscw (cw, nw, l0, floopr(nl), floopz(nl), pfx)
c                                                ! curr in AMPS
            psitestpf(nl) = psitestpf(nl) + sg*extcur(l0)*pfx*fac
         end do
         call fcp (jphi, r12, z12, nthpts, ns, floopr(nl), floopz(nl), 
     1      psitestip(nl), okay)
      end do
 
      write (*, 2908)
 2908 format(t22,'J-dot-gradphi',t42,'PFC',t56,'Total')
      do nl = 1, 5
         write (*, 2910) nl, nl + 1, psitestip(nl) - psitestip(nl+1), 
     1      psitestpf(nl) - psitestpf(nl+1), psitestip(nl) + psitestpf(
     2      nl) - psitestip(nl+1) - psitestpf(nl+1)
      end do
 2910 format(' psi(',i3,' - ',i3,')= ',3f15.6)
      nl = 6
      write (*, 2910) nl, nl - 5, psitestip(nl) - psitestip(nl-5), 
     1   psitestpf(nl) - psitestpf(nl-5), psitestip(nl) + psitestpf(nl)
     2    - psitestip(nl-5) - psitestpf(nl-5)
 
      end subroutine fluxcheck
      

      subroutine splint(xa, ya, y2a, n, x, y, yp, ndim)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ndim
      real(rprec), dimension(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, klo, khi, k
      real(rprec) :: c1o6, deriv, a, b, h, h2, a2, b2, y26lo, y26hi
C-----------------------------------------------
*
*       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
*       XA: ordered array of length N of ordinates at which function YA=F(XA)
*           is tabulated
*       YA: array of length N , = F(XA)
*       Y2A: array of second derivatives at XA points
*       computed from call to SPLINE
*       X : value at which Y = F(X) is to be computed from splines
*       YP = dY/dX at X
*       NDIM: dimension of X, Y, YP arrays
 
 
      c1o6 = 1.0/6.0
      deriv = yp(1)
      klo = 1
      khi = n
      do i = 1, ndim
 
         do while(khi - klo > 1)
            k = (khi + klo)/2
            if (xa(k) > x(i)) then
               khi = k
            else
               klo = k
            endif
         end do
 
         h = xa(khi) - xa(klo)
         a = xa(khi) - x(i)
         b = x(i) - xa(klo)
         h2 = h*h
         a2 = a*a
         b2 = b*b
         y26lo = c1o6*y2a(klo)
         y26hi = c1o6*y2a(khi)
         y(i) = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h
         if (deriv .ne. 0.0) yp(i) = (ya(khi)-ya(klo)+y26hi*(3.0*b2-h2)-
     1      y26lo*(3.0*a2-h2))/h
         if (i<ndim .and. x(i+1)>x(i)) then
            khi = n
         else
            klo = 1
         endif
      end do

      end subroutine splint
      
 
      subroutine icsel1(line, lastplot, explist, maxlist, len)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lastplot, maxlist, len
      character line*(*), explist*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: maxtokens = 20
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ntokens
      integer, dimension(maxtokens) :: ltokens
      integer :: i, ndt
      integer, dimension(2) :: ldtokens
      integer :: lim1, lim2, j, lint
      character, dimension(maxtokens) :: tokens*20
      character, dimension(2) :: dtokens*20
      character :: cint*20, addon*11
C-----------------------------------------------
 
c     Purpose: expand ICTRANS style plot list ( 1,3 5 ) to char list
c       e.g.,  1,3 into -1-2-3-
c       replace "$" with last plot number
 
c     input:  line [character*(*)]        line to be parsed
c             lastplot [integer]        substitute this for "$" in line
c             maxlist  [integer]        max len of list
c     output:
c             explist(1:len) [char]     generated list
c             len      [integer]        len of expanded list
 
 
 
 
c     Local Vbls
 
      explist = '-'
      len = 1
 
c     Break Up Token List By Whitespace First
      ntokens = maxtokens
      call wie_parse (line, ntokens, tokens, ltokens)
 
c     Check The Line For Dashes ... Cry Out And Die !
      if (index(line,'-') .ne. 0) then
         write (*, *) 
     1' Wrong format ... use COMMAs, not DASHES to indicate plot range !
     2'
         call exit (2)
      endif
 
c      do i = 1, ntokens
c       write (*,'(i3,2x,a)') i, tokens(i)
c      end do
 
c     Check Each Whitespace Delimited Token To See If It Is A
c     Triplet ("i,j") Or A Singlet ("i")
      do i = 1, ntokens
         ndt = 2
         call dlm_parse (tokens(i), ',', ndt, dtokens, ldtokens)
         if (ndt .eq. 2) then                      !  like 3,5 - Triplet -
            read (dtokens(1), 10) lim1
   10       format(i20)
            if (dtokens(2) .ne. '$') then
               read (dtokens(2), 10) lim2
            else
               lim2 = lastplot
            endif
         else if (ndt .eq. 1) then                 !  like 3   - Singlet -
            if (dtokens(1) .ne. '$') then
               read (dtokens(1), 10) lim1
            else
               lim1 = lastplot
            endif
            lim2 = lim1
         endif
         do j = lim1, lim2
            write (cint, 100) j
  100       format(i20)
            call str_strip (cint, lint)
            addon = cint(1:lint)//'-'
            if (len + lint + 1 > maxlist) then
               write (*, *) 
     1            ' ICSEL ??: # of elems expanded exceeds limit ', 
     2            maxlist
               call exit (2)
            endif
            explist = explist(1:len)//addon(1:lint+1)
            len = len + lint + 1
         end do
      end do
 
      end subroutine icsel1
       
 
      subroutine bgraph(phidiam, delphi, dsiobt, dsiext, plflux, nobd, 
     1   brbc, brcoil, plbrfld, bzbc, bzcoil, plbzfld, nbrbz, bbc, bcoil
     2   , plbfld, nmirnov, flmwgt, bcwgt)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nobd, nbrbz, nmirnov
      real(rprec) :: phidiam, delphi, flmwgt, bcwgt
      real(rprec), dimension(*) :: dsiobt, dsiext, plflux, brbc, 
     1   brcoil, plbrfld, bzbc, bzcoil, plbzfld, bbc, bcoil, plbfld
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nbcoils, i, ii, nunit, ilen, iflag, nclass, npts2, ndim
     1   , nwrk, nmeas
      real :: xleft, xtop, xwidth, xdrop
      real, dimension(4) :: gwnd
      real, dimension(2) :: barspace
      real, dimension(100,2) :: y
      real, dimension(100) :: class
      real, dimension(1000) :: wrk
      character :: xclass*100, title*80, xlabel*80, ylabel*80
C-----------------------------------------------
 
c    Subplot placement depends on whether nbcoils>0.
c    If >0, then this is the 2nd of 2 pages and the 1st
c    subplot goes in ULC.
c    If =0, then NADA
 
 
      data nwrk/1000/
      data ndim/100/
 
 
      nbcoils = nbrbz + nmirnov
 
*     Flux Loops
      if (nobd > 0) then
 
         i = 1
         y(i,1) = phidiam
         y(i,2) = delphi
         class(i) = 0
         do ii = 1, nobd
            i = i + 1
            y(i,1) = abs(dsiext(ii)+plflux(ii))*1.E3
            y(i,2) = abs(dsiobt(ii))*1.E3
            class(i) = ii
         end do
         nmeas = i
 
         nclass = nmeas
         npts2 = nmeas
         iflag = 3
 
         title = 'FLUX LOOPS: DATA (front) VS CODE (rear)'
         ylabel = 'ABS[FLUX] (mWb)'
         xlabel = 'LOOP INDEX (0=DMG)'
         xclass = ' 0  1  2  3  4  5  6 '
         nclass = nmeas
         nunit = 3
         barspace(1) = 2.0
         barspace(2) = -1.5
 
         if (nbcoils > 0) then
            gwnd(1) = 0.
            gwnd(2) = .5
            gwnd(3) = .5
            gwnd(4) = 1.
         else
            gwnd(1) = .5
            gwnd(2) = 1.
            gwnd(3) = .5
            gwnd(4) = 1.
         endif
 
         call bgraphn (title, ylabel, xlabel, class, xclass, nclass, 
     1      nunit, barspace, gwnd, y, ndim, npts2, iflag, wrk, nwrk)
 
      endif
!        call frame
      if (nbcoils .eq. 0) return 
 
      if (nbrbz > 0) then
 
*       BR Loops
         do ii = 1, nbrbz
            y(ii,1) = abs(brcoil(ii)+plbrfld(ii))*1.E3
            y(ii,2) = abs(brbc(ii))*1.E3
            class(ii) = ii
         end do
         nmeas = nbrbz
 
         nclass = nmeas
         npts2 = nmeas
         iflag = 3
 
         title = 'BR COILS: DATA (front) VS CODE (rear)'
         xlabel = 'COIL INDEX'
         ylabel = 'ABS[BR] (mT)'
         xclass = '1 2 3 4 5 6 7 8 9 10111213141516171819202122232425'
         nunit = 2
 
         barspace(1) = 2.0
         barspace(2) = -1.5
 
         gwnd(1) = 0.
         gwnd(2) = .5
         gwnd(3) = 0.
         gwnd(4) = .5
 
         call bgraphn (title, ylabel, xlabel, class, xclass, nclass, 
     1      nunit, barspace, gwnd, y, ndim, npts2, iflag, wrk, nwrk)
 
*       BZ Loops
 
         do ii = 1, nbcoils
            y(ii,1) = abs(bzcoil(ii)+plbzfld(ii))*1.E3
            y(ii,2) = abs(bzbc(ii))*1.E3
            class(ii) = ii
         end do
         nmeas = nbrbz
 
         nclass = nmeas
         npts2 = nmeas
         iflag = 3
 
         title = 'BZ COILS: DATA (front) VS CODE (rear)'
         xlabel = 'COIL INDEX'
         ylabel = 'ABS[BZ] (mT)'
         xclass = '1 2 3 4 5 6 7 8 9 10111213141516171819202122232425'
         nunit = 2
 
         barspace(1) = 2.0
         barspace(2) = -1.5
 
         gwnd(1) = .5
         gwnd(2) = 1.
         gwnd(3) = 0.
         gwnd(4) = .5
 
         call bgraphn (title, ylabel, xlabel, class, xclass, nclass, 
     1      nunit, barspace, gwnd, y, ndim, npts2, iflag, wrk, nwrk)
 
      endif
 
      if (nmirnov > 0) then
 
*       MIRNOV Loops
         do ii = 1, nmirnov
            y(ii,2) = abs(bcoil(ii)+plbfld(ii))*1.E3
            y(ii,1) = abs(bbc(ii))*1.E3
            class(ii) = ii
         end do
         nmeas = nmirnov
 
         nclass = nmeas
         npts2 = nmeas
         iflag = 3
 
         title = 'MIRNOV COILS: DATA (front) VS CODE (rear)'
         xlabel = 'COIL INDEX'
         ylabel = 'ABS[B-Theta] (mT)'
         xclass = '1 2 3 4 5 6 7 8 9 10111213141516171819202122232425'
         nunit = 2
 
         barspace(1) = 2.0
         barspace(2) = -1.5
 
         gwnd(1) = .5
         gwnd(2) = 1.
         gwnd(3) = .5
         gwnd(4) = 1.
 
         call bgraphn (title, ylabel, xlabel, class, xclass, nclass, 
     1      nunit, barspace, gwnd, y, ndim, npts2, iflag, wrk, nwrk)
 
      endif
 
      call frame
 
      return 
 
  300 format('Diamagnetic Flux = ',f7.3,' [Measured]')
      ilen = len_trim(title)
      call plchlq (xleft, xtop, title(1:ilen), xwidth, 0., -1.)
 
      write (title, 310) delphi*1.E3
  310 format('                     ',f7.3,' [Calculated]')
      ilen = len_trim(title)
      call plchlq (xleft, xtop - xdrop, title(1:ilen), xwidth, 0., -1.)
 
      write (title, 320) flmwgt*100
  320 format('Pct RMS Error in Saddle Loops =',f5.2)
      ilen = len_trim(title)
      call plchlq(xleft,xtop-3*xdrop,title(1:ilen),xwidth,0.,-1.)
 
      write (title, 330) bcwgt*100
  330 format('Pct RMS Error in B-Field =',f5.2)
      ilen = len_trim(title)
      call plchlq(xleft,xtop-4*xdrop,title(1:ilen),xwidth,0.,-1.)
 
 
      call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
      call prunlabel
 
      call frame
 
      end subroutine bgraph
      
 
      subroutine bgraphn(title, ylabel, xlabel, class, xclass, nclass, 
     1   nunit, barspace, gwnd, y, ndim, npts2, iflag, wrk, nwrk)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nclass, nunit, ndim, npts2, iflag, nwrk
      character title*(*), ylabel*(*), xlabel*(*), xclass*(*)
      real, dimension(*) :: class, barspace, gwnd
      real, dimension(ndim,2) :: y
      real, dimension(*) :: wrk
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(8) :: icolors
C-----------------------------------------------
 
 
      icolors(1) = 61
      icolors(2) = 61
      icolors(3) = 61
      icolors(4) = 61
      icolors(5) = 61
      icolors(6) = 61
      icolors(7) = 61
      icolors(8) = 61
 
      call hstopl ('DEF=ON')
      call hstopl ('SHA=OFF')
      call hstopi ('COL=ON', 3, 0, icolors, 8)
      call hstopc ('TIT=ON', title, 7, 3)
      call hstopl ('PRM=OFF')
      call hstopc ('FQN=ON', ylabel, 7, 3)
      call hstopc ('FOR=ON', '(F7.3)', 9, 3)
      call hstopc ('CHR=ON', xclass, nclass, nunit)
      call hstopc ('LAB=ON', xlabel, 7, 3)
      call hstopr ('SPA=ON', barspace, 2)
      call hstopr ('WIN=ON', gwnd, 4)
      call hstopl ('FRA=OFF')
      call histgr (y, ndim, npts2, iflag, class, nclass, wrk, nwrk)
 
      end subroutine bgraphn
      
 
      logical function makeplot (explist, lenlist, page, pagedesc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lenlist
      character*(*) :: explist, page, pagedesc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      logical :: flag, echo
      character :: id*10
C-----------------------------------------------
 
 
c     Local Vbls
      data echo/.TRUE./
 
 
c     Local Fncs
 
      id = '-' // trim(adjustl(page)) // '-'
 
      flag = index(explist(1:lenlist),id(1:len_trim(id))) > 0
 
      if (echo .and. flag) write (*, *) ' Making plot # ', 
     1   trim(adjustl(page)), ' : ', pagedesc
 
      makeplot = flag
 
      end function makeplot
       
 
      subroutine wrstab(ns, psi_in, pres, iota, runid, mnmax, rmnc, 
     1   zmns, rbt, itor, ffp, jdotb, betatot, betapol, betator, 
     2   betaxis)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Vtraneq
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns, mnmax
      real(rprec) :: itor
      real(rprec) :: rbt, betatot, betapol, betator, betaxis
      character runid*(*)
      real(rprec), dimension(*) :: psi_in, iota, ffp 
      real(rprec), dimension(*) :: rmnc, zmns, pres, jdotb
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: dmu0 = 1.256637d-06, zero = 0.d0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ishot, itime, mn, mjbal
      real(rprec), dimension(101) :: sp1, sp2, sp3
      real(rprec) :: err
C-----------------------------------------------
c======               START OF TRANEQ.BLK            =================
C
C   +++ See file EQOUT for definition of variables +++
C   +++                                            +++
C   +++ Note: There are two common blocks, EQCOM,  +++
C   +++  and EQFIL. The former passes numerical    +++
C   +++  data, while the latter passes character   +++
C   +++  data. The common block EQCOM has been     +++
C   +++  broken into parts for documentation       +++
C   +++  purposes. These blocks serve all of the   +++
C   +++  output codes which are interfaced.        +++
C   +++                                            +++
C   +++  Steve Sabbagh 9/11/90 Ver 1.0             +++
C   +++                1/02/91 Ver 1.1             +++
C   +++               12/30/93 Ver 1.2             +++
C   +++               08/30/94 Ver 1.3             +++
C   +++                 (80 character filenames)   +++
C   ++++++++++++++++++++++++++++++++++++++++++++++++++
C
c   +++ EQCOM common block +++
c   +++ 0) EQOUT input variables (number)
c
c   +++ 1) SHARED scalar variables +++
c
c   +++ 2) SHARED vector variables +++
c   +++      Note: r0b is a scalar +++
c
c   +++ 3) EQGRUM scalars +++
c
c   +++ 4) EQGRUM specific array variables +++
c
c   +++ 5) JSOLVER scalars +++
c
c   +++ 6) JSOLVER specific vector variables +++
c
c   +++ EQFIL common block +++
c
c   +++ J Profile & Beta Values
c
      read (runid(1:5), '(i5)', err=5) ishot
      shotnm = ishot
      go to 7
    5 continue
      shotnm = zero
      time = zero
      go to 20
    7 continue
      i = index(runid(9:),'t')
      if (i .eq. 0) time = zero
      if (i > 0) read (runid(9+i:9+i+2), '(i3)', err=10) itime
      time = itime/100
      go to 20
   10 continue
      time = zero
   20 continue
      intord = 0
      thdfil = 'wout.'//runid
      flnmo1 = runid(1:len_trim(runid))//'.eqi'
      pfilt1 = zero
      pfilt2 = zero
      pfilt3 = zero
      nequil = 1
      mjbal = ns + 2
      njav = mjbal
      mombnd = mnmax - 1
      mom = 0
      r0b = rmnc(1+mnmax*(ns-1))
      rmb(:mnmax-1) = rmnc(2+mnmax*(ns-1):mnmax*ns)
      ymb(:mnmax-1) = zmns(2+mnmax*(ns-1):mnmax*ns)
      rtor = r0b
      rtor = rtor + sum(rmb(2:mnmax:2))
      btor = abs(rbt)/rtor
      eqcamp = itor/dmu0
      ieqax = 2
      ieqedg = ns + 1
      psi(2:ns+1) = psi_in(:ns)
      q(2:ns+1) = 1./iota(:ns)
      press(2:ns+1) = pres(:ns)/dmu0
      jdotbc(2:ns+1) = jdotb(:ns)
      ggp(2:ns+1) = ffp(:ns)
      psi(1) = psi(3)
      q(1) = q(3)
      press(1) = press(3)
      jdotbc(1) = jdotbc(3)
      ggp(1) = ggp(3)
*
*     calculate derivatives using interpolating cubic splines
*     use Natural BC for q and pres
*
      call spline (ns, psi(ieqax), q(ieqax), sp1, sp2, sp3)
      dqdpsi(ieqax:ns-1+ieqax) = sp1(:ns)
      dqdpsi(ieqax-1) = dqdpsi(ieqax+1)
 
      call spline (ns, psi(ieqax), press(ieqax), sp1, sp2, sp3)
      dpdpsi(ieqax:ns-1+ieqax) = sp1(:ns)
      dpdpsi(ieqax-1) = dpdpsi(ieqax+1)
 
      batot = betatot
      bapol = betapol
      bator = betator
      baxis = betaxis
 
      call putstb
 
      end subroutine wrstab
      
 
      subroutine putstb
c
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     +++ PUTSTB:A subroutine that writes a transfer file FLNMO1    +++
c     +++        used by the EQGRUM equilibrium solver.             +++
c     +++                                                           +++
c     +++   Steve Sabbagh 1/ 3/94 Version 1.1                       +++
c     +++                                                           +++
c     +++        Small changes - status="new" changed to status=    +++
c     +++        "unknown" in open statement (for ATHENA) & now     +++
c     +++        supports 80 character filenames.                   +++
c     +++                                                           +++
c     +++   Steve Sabbagh 7/26/90 Version 1.0                       +++
c     +++                                                           +++
c     +++        Notes: Version 1.0: Writes "standard" stabil file  +++
c     +++           format in which profile are p, q, and psi.      +++
c     +++                                                           +++
c     +++        Adapted from the EQGRUM routine BALDIN by          +++
c     +++        M. Phillips. Serves as a BALDUR equilibrium        +++
c     +++        output file emulator.                              +++
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c     +++ Dimension and common blocks for EQGRUM variables +++
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Vtraneq
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ioutf, ishot, m, j, k
      real(rprec) :: pad
C-----------------------------------------------
c======               START OF TRANEQ.BLK            =================
C
C   +++ See file EQOUT for definition of variables +++
C   +++                                            +++
C   +++ Note: There are two common blocks, EQCOM,  +++
C   +++  and EQFIL. The former passes numerical    +++
C   +++  data, while the latter passes character   +++
C   +++  data. The common block EQCOM has been     +++
C   +++  broken into parts for documentation       +++
C   +++  purposes. These blocks serve all of the   +++
C   +++  output codes which are interfaced.        +++
C   +++                                            +++
C   +++  Steve Sabbagh 9/11/90 Ver 1.0             +++
C   +++                1/02/91 Ver 1.1             +++
C   +++               12/30/93 Ver 1.2             +++
C   +++               08/30/94 Ver 1.3             +++
C   +++                 (80 character filenames)   +++
C   ++++++++++++++++++++++++++++++++++++++++++++++++++
C
c   +++ EQCOM common block +++
c   +++ 0) EQOUT input variables (number)
c
c   +++ 1) SHARED scalar variables +++
c
c   +++ 2) SHARED vector variables +++
c   +++      Note: r0b is a scalar +++
c
c   +++ 3) EQGRUM scalars +++
c
c   +++ 4) EQGRUM specific array variables +++
c
c   +++ 5) JSOLVER scalars +++
c
c   +++ 6) JSOLVER specific vector variables +++
c
c   +++ EQFIL common block +++
c
c   +++ J Profile & Beta Values
c
c
c     +++ Open output file +++
      open(unit=20, file=flnmo1, status='unknown')
c
c     +++ I/O channels +++
      ioutf = 20                                 ! output file LUN
c
c     +++ Initialize +++
      pad = 0.0                                  !Pad variable
      ishot = int(shotnm)                  !integer value of shot number
c
c     +++ Write EQGRUM transfer file 'stabil'  +++
c     +++ The following is adapted from BALDIN +++
      write (ioutf, 150) nequil, time, nstep
  150 format(' equilibrium',i5,' at time',1pe15.7,' nstep',i5)
c
      write (ioutf, 9060) ishot, time        !These banner lines replace
      write (ioutf, 9070) intord           !the label arrays used in the
      write (ioutf, 9080) thdfil                 !previous routine.
      write (ioutf, 9090)
      write (ioutf, 9100) pfilt1, pfilt2, pfilt3
c
      write (ioutf, 154) njav, mombnd, mom
  154 format(' mjbal=',i5,' mombnd=',i5,' mom=',i5)
c
      write (ioutf, 155) r0b
  155 format(' r0b=',1pe15.7)
c
      write (ioutf, 156) (rmb(m),m=1,mombnd)
  156 format(' rmb=',1p5e15.7)
c
      write (ioutf, 157) (ymb(m),m=1,mombnd)
  157 format(' ymb=',1p5e15.7)
c
      write (ioutf, 158) btor, rtor, eqcamp
  158 format(' btor=',1pe15.7,' rtor=',1pe15.7,' eqcamp=',1pe15.7)
c
      write (ioutf, 160)
  160 format(t3,'j',t10,'psi',t25,'q',t40,'d q / d psi',t55,'p',t70,
     1   'd p / d psi',t83,'<jdotb>/<bdotgradv>',t109,'FF''')
c
c      --- Write profile arrays --- (jdotbc is NEW)
c         Note: First array element is a guard point and will be
c               ignored when read in from STABIL file
      do j = ieqax - 1, ieqedg
         write (ioutf, 162) j, psi(j), q(j), dqdpsi(j), press(j), dpdpsi
     1      (j), jdotbc(j), ggp(j)
      end do
c
c      --- Write dummy term on the end to satisfy BALDIN routine ---
      k = ieqedg + 1
      write (ioutf, 162) k, pad, pad, pad, pad, pad, pad
  162 format(1x,i4,1p7e15.7)
 
c      --- Write Beta values
      write (ioutf, 168) batot, bapol, bator, baxis
  168 format(' betatot=',1pe15.7,' betapol=',1pe15.7,' betator=',1pe15.7
     1   ,' betaxis=',1pe15.7)
 
c
c      +++ Close Files +++
      close(unit=20)
c
c      +++ Format Block +++
 9060 format(1x,'+++ TFTR Shot ',i6,' at ',1pe10.4,' seconds.  +++')
 9070 format(1x,'+++ TRANSP file:     Interpolation Order: ',i2,' +++')
 9080 format(1x,a80)
 9090 format(1x,'+++ Smoothing (%): Pressure     q      psi   +++')
 9100 format(1x,'+++                 ',f6.2,2x,f6.2,2x,f6.2,'   +++')

      end subroutine putstb

!       the codes (spline & seval) are taken from:
!       forsythe,malcolm and moler,
!       "computer methods for mathematical computations",
!       prentice-hall, 1977.
!
!       the codes (spleen,splaan & speval) are adaptations
!       by r.m. wieland for special cases ... see comments
 
 
      subroutine spline(n, x, y, b, c, d)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(rprec) :: t
C-----------------------------------------------
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n - 1
      if (n < 2) return 
      if (n >= 3) then
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
         d(1) = x(2) - x(1)
         c(2) = (y(2)-y(1))/d(1)
         d(2:nm1) = x(3:nm1+1) - x(2:nm1)
         b(2:nm1) = 2*(d(:nm1-1)+d(2:nm1))
         c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
         c(2:nm1) = c(3:nm1+1) - c(2:nm1)
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
         b(1) = -d(1)
         b(n) = -d(n-1)
         c(1) = 0
         c(n) = 0
         if (n .ne. 3) then
            c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(1) = c(1)*d(1)**2/(x(4)-x(1))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
         endif
         do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
         end do
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         end do
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
         b(n) = (y(n)-y(nm1))/d(nm1) + d(nm1)*(c(nm1)+2.*c(n))
         b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1+1) + 
     1      2*c(:nm1))
         d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
         c(:nm1) = 3*c(:nm1)
         c(n) = 3*c(n)
         d(n) = d(n-1)
         return 
c
      endif
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0
      d(1) = 0
      b(2) = b(1)
      c(2) = 0
      d(2) = 0

      end subroutine spline
      

      subroutine splaan(n, x, y, b, c, d)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(rprec) :: t
C-----------------------------------------------
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline for which s-prime(x1)=0.
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
c
      nm1 = n - 1
      if (n < 2) return 
      if (n >= 3) then
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
         d(1) = x(2) - x(1)
         c(2) = (y(2)-y(1))/d(1)
         d(2:nm1) = x(3:nm1+1) - x(2:nm1)
         b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
         c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
         c(2:nm1) = c(3:nm1+1) - c(2:nm1)
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
         b(1) = 2.*d(1)
         b(n) = -d(n-1)
         c(1) = 0.
         c(n) = 0.
         if (n .ne. 3) then
            c(1) = (y(2)-y(1))/d(1)
            c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
            c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
         endif
         do i = 2, n
            t = d(i-1)/b(i-1)
            b(i) = b(i) - t*d(i-1)
            c(i) = c(i) - t*c(i-1)
         end do
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         end do
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
         b(n) = (y(n)-y(nm1))/d(nm1) + d(nm1)*(c(nm1)+2.*c(n))
         b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1+1)+2.
     1      *c(:nm1))
         d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
         c(:nm1) = 3.*c(:nm1)
         c(n) = 3.*c(n)
         d(n) = d(n-1)
         return 
c
      endif
      b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.

      end subroutine splaan
