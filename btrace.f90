#!/bin/sh
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
      module btrace_mod
        use kind_spec
!**********************************************************************
!
!  namelist for controlling btrace
!
!
!  nfp  number of field periods.
!       0 assume no toroidal periodicity
!
!  lsym T  assume stellarator symmetry
!       F  assume no symmetry
!
!  field_source
!       0 Biot-Savart from coils-file
!       1 3D cubic splined from grid-file
!       2 dynamically generated splined-grid using Biot-Savart from coils
!
!  bspline
!       F field is integrated directly
!       T field is integrated on a grid, fit to bspline as a map.  
!         bspline coeffs are interpolated to subsequent tracing
!
!  rmin, rmax, zmin, zmax
!       bounds of computational domain
!
!  ntransits
!       number of period transits (or toroidal transits for nfp=0) to trace
!
!  phi_pl
!       toroidal angle (degrees) of poloidal planes to make plots at
!
!  step
!       horizontal step size in initial positions
!**********************************************************************
        integer, parameter :: max_planes = 100
        integer, parameter :: max_currents = 100

        integer :: nfp, ntransits, nplanes, field_source, 
     1             n_starts, max_starts, n_evals
        real(rprec) :: rmax, rmin, zmax, zmin, step, tol
        real(rprec), dimension(max_planes) :: phi_pl
        real(rprec), dimension(max_planes) :: phi_poin
        real(rprec), dimension(max_currents) :: coil_current
        character*(100) :: coils_file, mgrid_file 
        logical :: bspline, lsym

        real(rprec), dimension(:,:,:), allocatable :: r_poin,z_poin
        
        namelist /btrnml/ nfp, ntransits, field_source, 
     1    coils_file, mgrid_file, lsym, coil_current, tol,
     1    rmax, rmin, zmax, zmin, step, phi_pl, bspline, max_starts

      end module btrace_mod

!**********************************************************************
      module mgrid_mod
      use kind_spec 
!
!     grid dimensions:
!              ir  = no. radial (r) points in box
!              jz  = no. z points in box
!
!              suggest choosing hr == (rmax-rmin)/(ir-1) equal to
!                               hz == (zmax-zmin)/(jz-1)
!
      integer, parameter :: maxgroups = 100
      integer, parameter :: ir = 101, jz = 101
!                          ,kp = 16
!     real(rprec), parameter :: rmin = 0.09, rmax = 2.51,
!    1                          zmin =-1.21, zmax = 1.21
      integer :: nspul(maxgroups), ngroup(maxgroups)
      integer :: unit_parsed = 200
      integer :: nextcur
      real(rprec) :: extcur(maxgroups)
      character*8 :: dsilabel, bloopnames
      end module mgrid_mod

      module bcoils_mod
      use kind_spec
      integer :: nall, isteu
      real(rprec), dimension(:), allocatable :: xw, yw, zw, curre
      real(rprec) :: curcoils
      real(rprec), dimension(:), allocatable :: vx, vy, vz, dx, dy, dz
      end module bcoils_mod
EOF
cat > vbtrace.f << EOF
      program btrace
!
!     BTrace
!        This code traces the vacuum flux plot, given field input of
!        of various forms.  Forms input to date includ:
!        a) coils file
!
!        Future forms to be explored
!        b) Grid file
!        c) dynamically generated grid
!        d) splined field shape
!
!     Input is controlled by a namelist (above)
!
!     Using an improved version of Makegrid for the field calculation
!     M. Zarnstorff   July 2005
!
!     BOX DIMENSIONS: RMIN <= R <= RMAX,  ZMIN <= Z <= ZMAX
!                     KP = NO. TOROIDAL PLANES/FIELD PERIOD (MUST EQUAL VMEC VALUE)
      use mgrid_mod
      use bcoils_mod
      use btrace_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
      real(rprec) :: twopi
      integer, parameter :: igrid0=10, icontr=11
      character*(*), parameter :: NoCoil_String = 'QNotInAnyGroup'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: kp, istat, ig, k, j, i, numargs
      integer :: ifile, iextc, nfp2
      integer :: navg, nstarts, mxtransits
      integer, dimension(1) :: lowest
      real(rprec), allocatable, dimension(:, :, :) :: br, bz, bp
      real(rprec) :: phi_period, ravg, zavg, rstrt, ztop
      real(rprec) :: fperiod, delr, delz, delp, phi, zee, rad
      logical :: done
      character*30 :: curlabel(maxgroups)
      character*120 :: ext_file
      character*120 control_ext
      character*120 :: arg1

      external deriv_bs
C-----------------------------------------------
!
!    Input is via a namelist 'btrace' listed at top of this file
!
!
!
!     SET UP
!
      twopi = 8 * atan(one)

      call getcarg(1, arg1, numargs)
10    if (numargs .ge. 1) then
         arg1 = adjustl(arg1)
         control_ext = trim(arg1)

      else
         write (*,*) 'Enter extension of "btrace" file: '
         read (*, *) control_ext
      endif

!     initialize control variables
!
      phi_pl = 0
      nplanes = 1
      nfp = 0
      field_source = 0
      coils_file = ""
      mgrid_file = ""
      bspline = .false.
      rmin = 1.
      rmax = 2.
      zmin = -1.
      zmax = 1.
      ntransits = 1000
      coil_current = 0
      max_starts = 0

      ext_file = "btrace."//trim(control_ext)
      ifile = icontr
      call safe_open(ifile, istat, ext_file, 'old', 'formatted')
      if( istat .ne. 0) then
         write(*,*) '***Problems opening control file ***'
         numargs = 0   ! force question/answer
         goto 10
      endif

      read(nml=btrnml, unit=icontr, iostat=istat)
      close(unit=icontr)
      if( istat .ne. 0) then
         write(*,*) '***Problems reading control file ***'
         numargs = 0   ! force question/answer
         goto 10
      endif

      if (rmin .lt. 0.) stop ' rmin must be > 0'
      if (rmax .le. rmin) stop ' rmax must be > rmin'
      if (zmax .le. zmin) stop ' zmax must be > zmin'
      if (ntransits .le. 0) stop 'ntransits must be > 0'
      if (nfp .lt. 0) nfp = 0


      if( max_starts == 0) max_starts = (rmax-rmin)/step + 1.

      do while (nplanes.lt.max_planes .and. phi_pl(nplanes+1).ne. 0)
         nplanes = nplanes + 1
      enddo

      phi_period = 360.
      if (nfp .ne. 0) phi_period = 360._rprec/nfp

      do i=1, nplanes
         phi_pl(i) = modulo(phi_pl(i), phi_period)

         if (lsym .and. phi_pl(i) > phi_period/2 ) 
     1       phi_pl(i) = phi_period - phi_pl(i)
      enddo

!    make sure phi=0 is always first on the list
      i = 1
      phi_poin(1) = 0

      done = .false.
      do while (.not. done)
         lowest = minloc(phi_pl(:nplanes))
         if( phi_pl(lowest(1)) > phi_period) then
            done = .true.

         else 
            if( phi_poin(i) .ne. phi_pl(lowest(1))) then
               i = i+1
               phi_poin(i) = phi_pl(lowest(1))
            endif

            phi_pl(lowest(1)) = 2 * phi_period
         endif
      end do

      nplanes = i
      mxtransits = ntransits
      if( lsym) mxtransits = 2*mxtransits

      phi_poin(:nplanes) = phi_poin(:nplanes)*twopi/360

      allocate (r_poin(mxtransits,max_starts,nplanes),
     1          z_poin(mxtransits,max_starts,nplanes), stat=istat)
      if (istat .ne. 0) stop 'allocation error in btrace'
      
      if( field_source==0 .or. field_source==2) then
!    field is calculated from coils:  read in coils file         


!     PARSE FILAMENT FILE 
         call parse_coils_file(coils_file, curlabel, nfp2)

         if( nfp .ne. nfp2 .and. nfp.ne.0) then
            write(*,*) 'Mismatch in number of periods!!'
            write(*,*) 'Control namelist nfp=',nfp
            write(*,*) 'Coils file       nfp=',nfp2
            stop
         endif
            
         write(*,*)  'Coils and currents'
         do i=1, nextcur
            write(*,*) curlabel(i), coil_current(i)
         enddo
         write(*,*) '  '

         call readcoils(coil_current)

c      else if (field_source==1) then
!    field is calculated from pre-calculated grid-file: read it in

c         ext_file = 'mgrid.' // trim(mgrid_ext)
c         call get_grid(ext_file, coil_current, curlabels,)

      else
         stop 'Unimplemented source of magnetic field information'
      endif


! trace the fieldlines
      zavg = (zmin+zmax)/2
      if (lsym) zavg = 0
      i = 0
      rstrt = rmax
      ravg = rmin

      do while(i<max_starts .and. rstrt > ravg)
         n_evals = 0
         i = i+1

         call field_trace(deriv_bs, i, rstrt, zavg, istat)

         write(*,*) 'R, z, istat, evals =', rstrt, zavg, istat, n_evals
         if( istat == 0) then
            ztop = maxval(z_poin(:ntransits,i,1))
            if(.not. lsym) 
     1          zavg = sum(z_poin(:ntransits,i,1))/ntransits

            ravg = 0
            navg = 0
            do j=1, ntransits
               if( z_poin(j,i,1) < ztop/10) then
                  ravg = ravg + r_poin(j,i,1)
                  navg = navg + 1
               endif
            enddo

            ravg = ravg / navg
            write(*,*) 'ravg, zavg = ', ravg, zavg

         else if( istat < 0) then
            exit
         endif

         rstrt = rstrt - step
      enddo

      nstarts = i

      
      if (allocated(vx)) deallocate (vx, vy, vz, dx, dy, dz)
       
      end program btrace


      subroutine parse_coils_file(coil_file, curlabel, nfp)      
      use bcoils_mod
      use mgrid_mod, only: maxgroups, nspul, ngroup, unit_parsed,
     1                     nextcur
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nfp
      character*(*) :: coil_file, curlabel(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: max_entries=10000, icoil0=22
      real(rprec), parameter :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jgroup(0:max_entries), istat, igroup, jcount
      integer :: index1(1), j, icoil, line_parsed, parsed_group
      real(rprec), dimension(:), allocatable :: xwin, ywin, zwin, currin
      logical :: lcoil, lparse
      character*200 :: line, group
C-----------------------------------------------
!     This routine opens a scratch file (unit=unit_parsed) in which is
!     is written a list of the x,y,z,and current values for each coil
!     group. It obtains these by parsing the coil_file and stores the total
!     number of records written in jgroup(i).
!
!     Variable Definitions
!
!     nextcur             Number of external current groups (determined by parsing coils file)
!     ngroup(i)           Group ID of current group i (i = 1,nextcur)
!     jgroup(i)           Total number of xyz-cur entries written from k=0,i (jgroup(0)=0)
!                         = sum(k=1,i) nspul(k)
!
!     lcoil               set false when there are no new coil groups to parse
!                         used to terminate parsing of coil_file
!
C-----------------------------------------------
      icoil = icoil0
      call safe_open(icoil, istat, trim(coil_file), 'old', 'formatted')
      if (istat .ne. 0) stop 'Error opening input coil file'        

      call safe_open(unit_parsed, istat, 'null', 'scratch', 
     1     'unformatted')
      if (istat .ne. 0) stop 'Error opening parsed coil file'        

      allocate (xwin(max_entries), ywin(max_entries), 
     1          zwin(max_entries), currin(max_entries), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in parse_coils_file'

      nspul(:) = 0
      ngroup(:) = -1000
      nextcur = 0
      jgroup(:) = 0
      lcoil = .true.
      line_parsed = 2

      read (icoil, '(a)' , iostat=istat) line
      istat = index(line, 'periods')
      if (istat .eq. 0) 
     1   stop 'First line of coils file must contain periods'
      read (line, *, iostat=istat) group, nfp
      print *, 'nfp = ', nfp
      
      FILE_LOOP: do while (lcoil)

!        Eat lines already parsed...
         do j = 1, line_parsed
           read (icoil, *, err = 200)
         end do

!
!     Parse coil_file for current group labels and number of different coil types
!        
         jcount = 0
         lparse = .true.
         parsed_group = -max_entries

         PARSE_LOOP: do 
            read (icoil, '(a)', end = 100) line
            if (line(1:3) .eq. 'end') exit
            jcount = jcount + 1
            read (line, *, iostat=istat) xwin(jcount), ywin(jcount), 
     1           zwin(jcount), currin(jcount), igroup, group
         
            if (istat .eq. 0) then                                          !!Using new label format
!
!           Check if this group (with ID=igroup) already exists (has been previously added
!           to ngroup). If not, start a new coil group (add igroup to ngroup)
!           
               if (parsed_group .eq. -max_entries) parsed_group = igroup
               j = minval(abs(igroup - ngroup(:)))
               index1(1) = -1
               if (j .ne. 0) then
                  if (igroup .eq. parsed_group) then
                     nextcur = nextcur + 1
                     ngroup(nextcur) = parsed_group
                     index1(1) = nextcur
                     curlabel(nextcur) = trim(adjustl(group))
                  end if
               else 
                  index1 = minloc(abs(igroup - ngroup(:)))
               end if   

               if (line_parsed == 2) line_parsed = 3
               if (nextcur > maxgroups) 
     1            stop ' Number coil groups > maxgroups'

               if (igroup .eq. parsed_group) then
                  do j = 1, jcount
                     write (unit_parsed, IOSTAT=ISTAT) 
     1                 XWIN(J), YWIN(J), ZWIN(J), CURRIN(J)
                  end do
                  nspul(index1(1)) = nspul(index1(1)) + jcount
                  print '(a,i3,2a)', ' COIL GROUP = ', igroup,
     1             ' ID: ', trim(adjustl(group))
               else if (index1(1) .eq. -1) then
                  lparse = .false.                                       !!new coil group starts,will be parsed next time
               end if   

               if (lparse) line_parsed = line_parsed + jcount            !!line where parsing starts for next group

               if (jcount .gt. max_entries) 
     1            stop 'More than max_entries coil entries!'
      
               jcount = 0

            end if

         end do PARSE_LOOP

 100     continue        
                 
         if (nextcur .eq. 0) then           !!One coil in file...(old style)                
            nextcur = 1
            ngroup(1) = 1
         end if

!     After parsing coil parse_group=igroup, rewind coils file and eat first line_parsed lines of file
!     to start parsing next coil group

         rewind (icoil)
         if (parsed_group == -max_entries) lcoil = .false.
       
      end do FILE_LOOP

 200  continue
      
      close (unit=icoil)

      rewind (unit_parsed, iostat=istat)
      
      if (istat .ne. 0) then
         print *,' Error rewinding parsed coil file = ', istat
         stop
      end if
      
      print *
      print *,'Total no. current groups = ', nextcur
      print *

      end subroutine parse_coils_file


!
!     PUT BFIELD AND PFLUX SUBROUTINES HERE
!----------------------------------------------------------------------
!
!
! Lieber Stefan !
! These are the routines computing  the  vacuum field of the standard
! W7AS configuration  produced by modular coils and TF coils.
! Subroutine readcoil reads the coil data in W7ASCOILS.DAT file and has
! to be called once before calling the subroutine becoil.
! Becoil  computes the vacuum field by Biot-Savart.
!                      Viele Gruesse Peter
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      SUBROUTINE readcoils(currents)
! ----------------------------------------------------------------------
C
C     COMMON BPOL modified by jag:  curcoils
C     CURRENTS IN KILOAMPERES; B IN TESLA;  R,Z IN METERS
C
C ----------------------------------------------------------------------
      use bcoils_mod
      use mgrid_mod, only: ngroup, nspul, nextcur, unit_parsed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*), intent(in) :: currents
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
      integer :: nw, ierr, ic, ntot, nstart, nl
      real(rprec) :: curnorm
C ----------------------------------------------------------------------

!(SPH: ALL COILS READ IN, NOT JUST IN ONE FP)

      ntot = sum(nspul(1:nextcur))
      IF (ntot .LE. 0) STOP 'ERROR READING COILS IN READCOIL'

      IF (ALLOCATED(XW)) DEALLOCATE (XW, YW, ZW, CURRE)
      ALLOCATE (XW(ntot), YW(ntot), ZW(ntot), CURRE(ntot), STAT=IERR)
      IF (IERR .NE. 0) STOP 'ALLOCATION ERROR IN READCOIL'

      nl = 0
      do ic=1, nextcur

         NALL = NSPUL(IC)
         IF (NALL .LE. 0) STOP 'ERROR READING COILS IN READCOIL'
 
         NW = 0
         IERR = 0
!
!     UNIT_PARSED WAS PREVIOUSLY OPENED AND REWOUND IN PARSE_COIL ROUTINE
!     IT IS ASSUMED - AND REQUIRED - THAT THE groups are read 
!     IN THE ORDER THAT THE COIL GROUPS ARE WRITTEN OUT (OTHERWISE WE
!     WOULD HAVE TO DO RECORD SEARCHES IN A DIRECT-ACCESS FILE)
!
         DO WHILE (IERR.EQ.0 .AND. NW.LT.NALL)
            NW = NW+1
            READ (UNIT_PARSED, IOSTAT=IERR, END=1000) 
     1           XW(nl+NW),YW(nl+NW),ZW(nl+NW),CURRE(nl+NW)
         END DO

 1000 CONTINUE

C
C     NORMALIZE B-FIELD TO UNIT CURRENT (BASED ON FIRST CURRENT ELEMENT 
c     IN GROUP)
C     SAVE ACTUAL (FIRST) CURRENTS IN EXTCUR ARRAY
C
         curnorm = CURRE(nl+1)
         IF (CURRE(nl+1) .NE. ZERO) 
     1       CURRE(nl+1:nl+nw) = CURRE(nl+1:nl+nw)*currents(ic)/curnorm
      
         PRINT *,NW,' FILAMENTS READ IN FROM COIL GROUP ', IC,
     1       ' CURRENT = ', curnorm, currents(ic)
     

C     PRINT *,' XW(1) = ',xw(nl+1), ' XW(NW) = ',xw(nl+NW)
C     NALL   = NW*NP

         IF (IERR.NE.0 .OR. NW.NE.NALL) 
     1       STOP 'ERROR READING COILS FILE IN READCOIL'

         if( curnorm .ne. 0 .and. currents(ic) .ne. 0) nl = nl + nw

      enddo

      nall = nl
      
      END SUBROUTINE readcoils


! ----------------------------------------------------------------------
      subroutine field_trace (deriv, ii, rstart, zstart, istat)
! ----------------------------------------------------------------------
      use bcoils_mod, only: rprec, nall
      use btrace_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ii
      real(rprec), intent(in) :: rstart, zstart
      integer, intent(out) :: istat
      external deriv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, parameter :: neq = 2
      real(rprec), parameter :: one = 1
      integer :: iphi, iphi_st, n_poin, itask, istate, 
     1           lrw, liw, foo, mf
      real(rprec) :: rp, zp, phip, t, tout, tbase, twopi, rtol, atol
      real(rprec), dimension(2) :: y
      real(rprec), dimension(20+16*neq) :: rwork
      integer, dimension(20) :: iwork

      logical :: lforward
C-----------------------------------------------
      twopi = 8 * atan(one)

      istat = 0  ! all ok

      lrw = 20 + 16*neq
      liw = 20

      rtol = tol
      atol = tol
      mf = 10

      r_poin(1, ii, 1) = rstart
      z_poin(1, ii, 1) = zstart

!  initialize state & put in output
      y(1) = rstart
      y(2) = zstart
      itask = 1
      istate = 1

      n_poin = 0
      iphi_st = 2

      t = phi_poin(1)
      tbase = 0

      do while( n_poin < ntransits)
         n_poin = n_poin + 1

         if( iphi_st == 2) then   ! need to record symmetry plane
            r_poin(n_poin, ii, 1) = y(1)
            z_poin(n_poin, ii, 1) = y(2)
         endif


         do iphi = iphi_st, nplanes
            tout = tbase + phi_poin(iphi)

            call lsode(deriv, neq, y, t, tout, 1, rtol, atol, itask, 
     1              istate, 0, rwork, lrw, iwork, liw, foo, mf)

            if( istate .ne. 2 ) exit

c            write(*,*) 'Got to ',t

            r_poin(n_poin, ii, iphi) = y(1)
            z_poin(n_poin, ii, iphi) = y(2)
         enddo

         if( istate .ne. 2) exit

         if( nfp == 0) then
            iphi_st = 1
            tbase = tbase + twopi

         else if( .not. lsym) then
            iphi_st = 1
            tbase = tbase + twopi/nfp

         else 
            r_poin(n_poin+ntransits, ii, iphi) = y(1)
            z_poin(n_poin+ntransits, ii, iphi) = - y(2)

            tbase = tbase + twopi/nfp

            do iphi = nplanes-1, 1, -1
               tout = tbase - phi_poin(iphi)

               call lsode(deriv, neq, y, t, tout, 1, rtol, atol,itask, 
     1              istate, 0, rwork, lrw, iwork, liw, foo, mf)

               if( istate .ne. 2 ) exit

c               write(*,*) 'Backwards to ',t

               r_poin(n_poin+ntransits, ii, iphi) = y(1)
               z_poin(n_poin+ntransits, ii, iphi) = - y(2)
            enddo

            if( istate .ne. 2 ) exit

            iphi_st = 2
         endif

      enddo

      if( istate .ne. 2) then
         write(*,*) 'LSODE Error ',istate,' at phi=',t
         istat = istate
      
      else if( y(1)<rmin .or. y(1)>rmax .or. 
     1         y(2)<zmin .or. y(2)>zmax) then
         istat = 1

      else
         istat = 0

      endif
      end subroutine field_trace

! ----------------------------------------------------------------------
      subroutine deriv_bs (neq, phi, y, ydot)
! ----------------------------------------------------------------------
      use bcoils_mod, only: rprec, nall
      use btrace_mod, only: rmin, rmax, zmin, zmax, n_evals
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: neq
      real(rprec), intent(in) :: phi
      real(rprec), intent(in), dimension(neq) :: y
      real(rprec), intent(out), dimension(neq) :: ydot
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: rp, xp, yp, zp, bx, by, br, bz, bphi
C-----------------------------------------------
      rp = y(1)
      zp = y(2)

      if( rp<rmin .or. rp>rmax .or. zp<zmin .or. zp>zmax) then
         ydot = 0

      else
         XP = RP * COS(PHI)
         YP = RP * SIN(PHI)

         CALL BECOIL(XP,YP,ZP,BX,BY,BZ, nall)

         BR   = BX * COS(PHI) + BY * SIN(PHI)
         BPHI =-BX * SIN(PHI) + BY * COS(PHI)

         ydot(1) = br * rp / bphi
         ydot(2) = bz * rp / bphi

c         write (*,'(6es11.2)') phi, rp, zp, bx, by, bz

      endif

      n_evals = n_evals + 1
c      write (*,'(11x,2es11.2)') ydot(1), ydot(2)


      END SUBROUTINE deriv_bs


! ----------------------------------------------------------------------
      SUBROUTINE BECOIL (XP,YP,ZP,BX,BY,BZ,NALLi)
! ----------------------------------------------------------------------
      use bcoils_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(RPREC), INTENT(IN) :: XP, YP, ZP
      INTEGER, INTENT(IN) :: nalli
      REAL(RPREC), INTENT(OUT) :: BX, BY, BZ
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(RPREC), PARAMETER :: FAC = 1.E-07_DP, ZERO = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: I, NALL1
      REAL(RPREC) :: EPS, EPS1, AX, AY, AZ, R12, XP1, faf
c      REAL(RPREC), DIMENSION(:), ALLOCATABLE ::  X1, Y1, Z1, RW, FA
      REAL(RPREC), DIMENSION(nalli) ::  X1, Y1, Z1, RW
      REAL(RPREC), DIMENSION(nalli-1) ::  FA
      LOGICAL :: TOO_CLOSE
C-----------------------------------------------
C
C     COMMON /VBFIELD/ CURTOR,CURPOL,BZ0
C
C
      TOO_CLOSE = .FALSE.
      XP1 = XP
      NALL1 = NALL - 1
      IF(ISTEU .EQ. NALL) GO TO 20
         IF (ALLOCATED(VX)) DEALLOCATE (VX, VY, VZ, DX, DY, DZ, STAT=I)
         ALLOCATE (VX(NALL1), VY(NALL1), VZ(NALL1), DX(NALL1), 
     1             DY(NALL1), DZ(NALL1), STAT=I)
         IF (I .NE. 0) STOP 'ALLOCATION ERROR IN BECOIL SETUP'

      DO I=1,NALL1
         DX(I) = (XW(I+1)-XW(I))*CURRE(I)
         DY(I) = (YW(I+1)-YW(I))*CURRE(I)
         DZ(I) = (ZW(I+1)-ZW(I))*CURRE(I)
         VX(I) =  YW(I)*DZ(I)-ZW(I)*DY(I)
         VY(I) =  ZW(I)*DX(I)-XW(I)*DZ(I)
         VZ(I) =  XW(I)*DY(I)-YW(I)*DX(I)
      END DO

         ISTEU  = NALL
   20 CONTINUE
   
c      ALLOCATE (X1(NALL), Y1(NALL), Z1(NALL), RW(NALL),FA(NALL1),STAT=I)
c      IF (I .NE. 0) STOP 'ALLOCATION ERROR IN BECOIL'
   25 CONTINUE

      DO I=1,NALL
         X1(I)  = XP1- XW(I)
         Y1(I)  = YP - YW(I)
         Z1(I)  = ZP - ZW(I)
         RW(I)  = SQRT(X1(I)*X1(I)+Y1(I)*Y1(I)+Z1(I)*Z1(I))
      END DO

      DO I=1,NALL1
         R12    = RW(I+1)*RW(I)
         FA(I)  = R12*(R12 +X1(I+1)*X1(I)+Y1(I+1)*Y1(I)+Z1(I+1)*Z1(I))
      END DO

      WHERE (CURRE(1:NALL1) .EQ. ZERO) FA = 1
      EPS = EPSILON(BX)
      EPS1= EPS**2*(XP**2 + YP**2)
      IF (ANY(ABS(FA) .LE. EPS1)) THEN
         IF (.NOT.TOO_CLOSE) THEN
            TOO_CLOSE = .TRUE.
            XP1 = XP1 * (1 + SQRT(EPS))
            GO TO 25 
         ELSE
            write (6, 1000) XP, YP, ZP
            BX = 0;  BY = 0;   BZ = 0
            GO TO 50
         END IF
      END IF

      ax = 0; ay = 0; az = 0
      bx = 0; by = 0; bz = 0

      do i=1, nall1
         FA(i) = (RW(i) + RW(i+1)) / FA(i)
      enddo


      do i=1, nall1
         AX     = ax + FA(i)*DX(i)
         AY     = ay + FA(i)*DY(i)
         AZ     = az + FA(i)*DZ(i)
         BX     = bx + FA(i)*VX(i)
         BY     = by + FA(i)*VY(i)
         BZ     = bz + FA(i)*VZ(i)
      enddo

      BX     = FAC*(bx - YP*AZ+ZP*AY)
      BY     = FAC*(by - ZP*AX+XP*AZ)
      BZ     = FAC*(bz - XP*AY+YP*AX)

      IF (TOO_CLOSE) THEN
         write (6, 1000) XP, YP, ZP
         write (6, 1100) BX, BY, BZ
      END IF

 50   CONTINUE
       
c      DEALLOCATE (X1, Y1, Z1, RW, FA)

 1000 FORMAT(' OBSER. PT. HITS COIL AT: XP = ', 1PE12.4, ' YP = ', 
     1       1PE12.4, ' ZP = ', 1PE12.4)
 1100 FORMAT(' OBSER. PT. SHIFTED, BX = ', 1PE12.4, ' BY = ', 
     1        1PE12.4, ' BZ = ', 1PE12.4)     


      END SUBROUTINE BECOIL

EOF
EOC
