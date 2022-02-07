      program makegrid
!
!     THIS CODE (MAKEGRID) GENERATES B-FIELD COMPONENTS ON R,Z, PHI GRID
!     AND COMPUTES THE POLOIDAL FLUX AT NOBS OBSERVATION POINTS
!
!     NOTE TO USER: EXPERIENCE SHOWS THAT A GRID SPACING OF
!     DEL R = DEL Z <= 1/50 WORKS WELL.  HERE, DEL R = RMAX-RMIN.
!     GRID SPACING MUCH LARGER THAN THIS MAY ADVERSELY EFFECT CONVERGENCE
!     OF VACUUM CODE.
!
!     BOX DIMENSIONS: RMIN <= R <= RMAX,  ZMIN <= Z <= ZMAX
!                     KP = NO. TOROIDAL PLANES/FIELD PERIOD (MUST EQUAL VMEC VALUE)
      use mgrid_mod
      use becoil1_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
      integer, parameter :: igrid0=10
      character*(*), parameter :: NoCoil_String = 'QNotInAnyGroup'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: kp, nfp, istat, nextcur, ig, k, j, i, numargs
      integer :: igrid, iextc
      integer :: kp2, jz2, kp_odd, jz_odd
      real(rprec), allocatable, dimension(:, :, :) :: br, bz, bp
      real(rprec) :: rmin, rmax, zmin, zmax
      real(rprec) :: fperiod, delr, delz, delp, phi, zee, rad
      logical :: lstell_sym
      character*30 :: curlabel(maxgroups)
      character*120 :: mgrid_file, coil_file, extcur_file
      character*120 mgrid_ext
      character*120 :: arg1
C-----------------------------------------------
!
!     THIS CAN BE RUN EITHER FROM THE COMMAND LINE (INTERACTIVELY) OR FROM A COMMAND FILE:
!
!     xgrid < filename.cmd
!
!     where the command file (filename.cmd) has entries:
!     coils file extension
!     stell_sym (T/F)
!     rmin value
!     rmax value
!     zmin value
!     zmax value
!

!
!     SET UP GRID DIMENSIONS
!
      call getcarg(1, arg1, numargs)
      lstell_sym = .false.
      
      if (numargs .lt. 7) then
         write (*, 220) 'Enter extension of "coils" file: '
         read (*, *) mgrid_ext
         write (*, 220) 'Assume stellarator symmetry (Y/N)? '
         read (*, *) mgrid_file
         if (mgrid_file(1:1) == 'Y' .or. mgrid_file(1:1) == 'y') 
     1      lstell_sym = .true.
         write (*, 220) 'Enter Rmin (min radial grid dimension): '
         read (*, *) rmin
         write (*, 220) 'Enter Rmax (max radial grid dimension): '
         read (*, *) rmax
         write (*, 220) 'Enter Zmin (min vertical grid dimension): '
         read (*, *) zmin
         write (*, 220) 'Enter Zmax (max vertical grid dimension): '
         read (*, *) zmax
         write (*, 220) 
     1    'Enter number of toroidal planes/period for data: '
         read (*, *) kp
         write (*, 220) 
     1    'Enter number of radial grid points (e.g. 101): '
         read (*, *) ir
         write (*, 220) 
     1    'Enter number of vertical grid points (e.g. 101): '
         read (*, *) jz

      else                                                               !Command file
         arg1 = adjustl(arg1)
         mgrid_ext = trim(arg1)
         call getcarg(2, arg1, numargs)
         read (arg1, *) lstell_sym
         call getcarg(3, arg1, numargs)
         read (arg1, *) rmin
         call getcarg(4, arg1, numargs)
         read (arg1, *) rmax
         call getcarg(5, arg1, numargs)
         read (arg1, *) zmin
         call getcarg(6, arg1, numargs)
         read (arg1, *) zmax
         call getcarg(7, arg1, numargs)
         read (arg1, *) kp

c        stick with default values of ir and jz (101 each)

      end if   

      if (rmin .lt. 0.) stop ' rmin must be > 0'
      if (rmax .le. rmin) stop ' rmax must be > rmin'
      if (zmax .le. zmin) stop ' zmax must be > zmin'
      if (kp .le. 0) stop 'kp must be > 0'

      allocate (br(ir,jz,kp), bz(ir,jz,kp), bp(ir,jz,kp), stat=istat)
      if (istat .ne. 0) stop 'allocation error in xgrid'
      
      if (lstell_sym) then
         kp2 = kp/2;  jz2 = jz/2
         kp_odd = mod(kp,2)
         jz_odd = mod(jz,2)
!
!        Must be sure zmax = -zmin 
!
         if (abs(zmax) > abs(zmin)) then
            zmax = abs(zmax)
            zmin = -zmax
         else 
            zmin = -abs(zmin)
            zmax = -zmin
         end if
      else
         kp2 = kp;    jz2 = jz
         kp_odd = 0;  jz_odd = 0
      end if

      coil_file = 'coils.' // trim(mgrid_ext)
      mgrid_file = 'mgrid.' // trim(mgrid_ext)
      extcur_file = 'extcur.' // trim(mgrid_ext)
      if (lstell_sym) print *,' Stellarator symmetry IS assumed'
      if (.not.lstell_sym)print *,' Stellarator symmetry IS NOT assumed'
      print *,' RMIN = ', rmin,' RMAX = ', rmax      
      print *,' ZMIN = ', zmin,' ZMAX = ', zmax      
      print *,' kp = ',  kp
      print *
      print *,'Input  file: ',trim(coil_file)
      print *,'Mgrid  file: ',trim(mgrid_file)
      print *,'Extcur file: ',trim(extcur_file)

      igrid = igrid0
      call safe_open (igrid, istat, trim(mgrid_file), 
     1   'replace', 'unformatted')
      if (istat .ne. 0) then
         print *, ' XGRID could not create ', trim(mgrid_file)
         print *, ' IOSTAT = ', istat,' IUNIT = ', igrid
         stop 20
      end if

      iextc = igrid+1
      call safe_open(iextc, istat, trim(extcur_file), 
     1   'replace', 'formatted')
      if (istat .ne. 0) then
         print *, ' XGRID could not create ', trim(extcur_file)
         print *, ' IOSTAT = ', istat,' IUNIT = ', iextc
         stop 25
      end if

!
!     PARSE FILAMENT DUMP FILE FOR NUMBER OF FIELD PERIODS
!     SPLIT INTO COIL GROUPS. DETERMINE NEXTCUR
!     COMING OUT, IGROUP+100=UNIT NO IS OPENED AND READY TO READ
!     AFTER REWINDING
!
      call parse_coils_file(coil_file, curlabel, nextcur, nfp)      
        
      fperiod = (8*atan(one))/nfp
      delr = (rmax-rmin)/(ir-1)
      delz = (zmax-zmin)/(jz-1)
      delp = fperiod/kp


      write(igrid) ir, jz, kp, nfp, nextcur
      write(igrid) rmin, zmin, rmax, zmax
      write(igrid) (curlabel(i), i=1,nextcur)
!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "BFIELD" TO COMPUTE THESE VALUES
!
      do 100 ig = 1,nextcur

         k = 1                     ! this is always a symmetry plane
         phi = (k-1)*delp
         do j=1,jz2 + jz_odd
            zee = zmin + (j-1)*delz
            do i=1,ir
               rad = rmin + (i-1)*delr
               call bfield (rad, phi, zee, br(i,j,k),
     1                      bp(i,j,k), bz(i,j,k), ig)
               if (lstell_sym) then
                  br(i,jz+1-j,k) = -br(i,j,k)
                  bz(i,jz+1-j,k) =  bz(i,j,k)
                  bp(i,jz+1-j,k) =  bp(i,j,k)
               end if
            enddo
         enddo
         print *,' K = ',k,' (OUT OF ',KP,')'

         do k=2,kp2+kp_odd
            phi = (k-1)*delp
            do j=1,jz
               zee = zmin + (j-1)*delz
               do i=1,ir
                  rad = rmin + (i-1)*delr
                  call bfield (rad, phi, zee, br(i,j,k),
     1                         bp(i,j,k), bz(i,j,k), ig)

                  if (lstell_sym) then
                     br(i,jz+1-j,kp+2-k) = -br(i,j,k)
                     bz(i,jz+1-j,kp+2-k) =  bz(i,j,k)
                     bp(i,jz+1-j,kp+2-k) =  bp(i,j,k)
                  end if
               enddo
            enddo
            print *,' K = ',k
         enddo

         if ((kp_odd == 0) .and. lstell_sym) then       ! another symmetry plane
            k = kp2 + 1
            phi = (k-1)*delp
            do j=1,jz2 + jz_odd
               zee = zmin + (j-1)*delz
               do i=1,ir
                  rad = rmin + (i-1)*delr
                  call bfield (rad, phi, zee, br(i,j,k),
     1                         bp(i,j,k), bz(i,j,k), ig)

                  br(i,jz+1-j,k) = -br(i,j,k)
                  bz(i,jz+1-j,k) =  bz(i,j,k)
                  bp(i,jz+1-j,k) =  bp(i,j,k)
               enddo
            enddo
            print *,' K = ',k
         end if

         write(igrid)(((br(i,j,k),bz(i,j,k),bp(i,j,k),i=1,ir),
     1                   j=1,jz),k=1,kp)
         if (ig .lt. 10) write(iextc, 210) ig, extcur(ig)
         if (ig .ge. 10) write(iextc, 215) ig, extcur(ig)

 100  continue

      close (igrid)
      close (iextc)
      
      if (allocated(vx)) deallocate (vx, vy, vz, dx, dy, dz)
      deallocate (br, bp, bz)
      
      write (*, 200)
      DO IG = 1, NEXTCUR
         WRITE (*, '(a7,i2,a4,1pe22.14)') 'EXTCUR(',IG,') = ',EXTCUR(IG)
      END DO
 200  FORMAT(/, 
     1  ' USE THE FOLLOWING ARRAY FOR EXTERNAL CURRENTS IN INDATA',
     2  ' INPUT FILE:'/)

 210  format('EXTCUR(', i1,')  = ', 1pe22.14)
 215  format('EXTCUR(', i2,') = ', 1pe22.14)
 220  format(a)
 
      end program makegrid


      subroutine parse_coils_file(coil_file, curlabel, nextcur, nfp)      
      use bpol_mod
      use mgrid_mod, only: maxgroups, nspul, ngroup, unit_parsed
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nextcur, nfp
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

      SUBROUTINE BFIELD (RP, PHI, ZP, BR, BPHI, BZ, IG)
      use bpol_mod, only: isteu, rprec, nall
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ig
      real(rprec), intent(in) :: rp, phi, zp
      real(rprec), intent(out) :: br, bphi, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, save :: icount = 1
      real(rprec) :: xp, yp, bx, by
C-----------------------------------------------

      if (icount .eq. ig) then
        isteu = 0
        call readcoil(ig)
        icount = icount+1
      end if
        
      XP = RP * COS(PHI)
      YP = RP * SIN(PHI)

      CALL BECOIL(XP,YP,ZP,BX,BY,BZ, nall)

      BR   = BX * COS(PHI) + BY * SIN(PHI)
      BPHI =-BX * SIN(PHI) + BY * COS(PHI)

      END SUBROUTINE BFIELD


! ----------------------------------------------------------------------
      SUBROUTINE READCOIL(ICOUNT)
! ----------------------------------------------------------------------
C
C     COMMON BPOL modified by jag:  curcoils
C     CURRENTS IN KILOAMPERES; B IN TESLA;  R,Z IN METERS
C
C ----------------------------------------------------------------------
      use bpol_mod
      use mgrid_mod, only: ngroup, nspul, extcur, unit_parsed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: icount
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0
      integer :: nw, ierr
      real(rprec) :: curmod, curaux
C ----------------------------------------------------------------------

!(SPH: THESE ARE ONLY HERE FROM PKM - NOT USED)  
      CURMOD = 1
      CURAUX = 0
      curcoils = (120*CURAUX + 840*CURMOD)

!(SPH: ALL COILS READ IN, NOT JUST IN ONE FP)

      NALL = NSPUL(ICOUNT)
      IF (NALL .LE. 0) STOP 'ERROR READING COILS IN READCOIL'
      IF (ALLOCATED(XW)) DEALLOCATE (XW, YW, ZW, CURRE)
      ALLOCATE (XW(NALL), YW(NALL), ZW(NALL), CURRE(NALL), STAT=IERR)
      IF (IERR .NE. 0) STOP 'ALLOCATION ERROR IN READCOIL'
 
      NW = 0
      IERR = 0
!
!     UNIT_PARSED WAS PREVIOUSLY OPENED AND REWOUND IN PARSE_COIL ROUTINE
!     IT IS ASSUMED - AND REQUIRED - THAT THE CALL TO THIS SUBROUTINE IS 
!     MADE IN THE ORDER THAT THE COIL GROUPS ARE WRITTEN OUT (OTHERWISE WE
!     WOULD HAVE TO DO RECORD SEARCHES IN A DIRECT-ACCESS FILE)
!
      DO WHILE (IERR.EQ.0 .AND. NW.LT.NALL)
         NW = NW+1
         READ (UNIT_PARSED, IOSTAT=IERR, END=1000) 
     1       XW(NW),YW(NW),ZW(NW),CURRE(NW)
      END DO

 1000 CONTINUE

      PRINT *,NW,' FILAMENTS READ IN FROM COIL GROUP ', ICOUNT,
     1       ' CURRENT = ', CURRE(1)
     
C
C     NORMALIZE B-FIELD TO UNIT CURRENT (BASED ON FIRST CURRENT ELEMENT IN GROUP)
C     SAVE ACTUAL (FIRST) CURRENTS IN EXTCUR ARRAY
C
      EXTCUR(ICOUNT) = CURRE(1)
      IF (CURRE(1) .NE. ZERO) CURRE(:) = CURRE(:)/EXTCUR(ICOUNT)
      
C     PRINT *,' XW(1) = ',xw(1), ' XW(NW) = ',xw(NW)
C     NALL   = NW*NP

      IF (IERR.NE.0 .OR. NW.NE.NALL) 
     1   STOP 'ERROR READING COILS FILE IN READCOIL'
      
      END SUBROUTINE READCOIL

! ----------------------------------------------------------------------
      SUBROUTINE BECOIL (XP,YP,ZP,BX,BY,BZ,NALLi)
! ----------------------------------------------------------------------
      use bpol_mod
      use becoil1_mod
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

!        PI4   = 8*ASIN(1._DP)
!        FAC   =  1.0/PI4
!        CURPOL=  curcoils
!      WRITE(6,7800) curcoils,CURPOL,FAC
! 7800 FORMAT('   IN BECOIL curcoils,CURPOL,FAC=',1P3E12.4)
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

