!----
! based on NCARPLTS V1.1 \\ by Derek Fox
! Modified by Manish Vachharajani and Dick Wieland (Summer, 1993)
!
! Modified CONTGRAF to CONTGRAF1 & CONTGRAF2
!
! Warnings:
!     NCARG appears to conscript I/O unit 3. Beware of others!
! Changes:
!     1. added "gwnd"  arg to contgraf1 - window coordinates that allow
!        to position plot on page for 4 per page, e.g..
!        Use in place of "grph"
!     2. in CONTGRAF1, define a contracted "gwnd2" grid, based on "grid"
!        contraction, for use in contour plot calls. Use in place of "grid".
!     3. in CONTGRAF1, use "gwnd" in PLCHHQ calls to position labels.
!     4. in FINDLEN, modified the algorithm: before, look for the 1st
!        "$" as a signature of string length; now: use len_trim to
!        measure actual string_length. Last character must still be a "$"
!        (because it is discarded), but now embedded "$" can be in the
!        string, as may be desired for exotic fonts.
!     5. Changed the font-instruction character from "'" back to the
!        NCAR default "!". It is easier to use in text strings.
!     6. Added AGPWRT to enable the use of PLCHHQ with AUTOGRAPH routines.
!     7. Removed the Derek Fox title page
!     8. Added CONTGRAF2 for full page contour plots (no "gwnd" arg):
!        - no "gwnd" used in place of "grph" & "grid"
!        - no "gwnd2" used in place of "grid"
!     9. in CONTGRAF1 & CONTGRAF2, legend generation is commented out
!     10. replace external CPDRPL with modified local version; result
!         is to write through label boxes; remove contour labels;
!         all contours in solid now, no longer solid-dash alternating.
!
!----
!
!                      ********************
!                      *                  *
!                      *     ncarplts1    *
!                      *                  *
!                      ********************
!......................................................................
!
!
!  COMPILATION:  requires ncar and ncargks libraries, v3.0.  Compile
!     ncarplts.f and link to your program and to the libraries.
!
!        An easy way to do this is to use ncar-s special version of
!     the fortran compiler, ncargf90:
!
!             ncargf90 ncarplts.f program.f
!
!     where the name of your program replaces "program.f".  This
!     will produce the executable "a.out".
!
!        Duplication Cautions should be expected because of the
!     replacement of NCAR-s default CPMPXY routine.
!
!
!
!  Subsidiary routines (see the routines themselves for
!    further information):
!
!       colram:   slave routine for contgraf (coloring
!                 of areas between contour lines).
!       cpmpxy:   slave routine for contgraf (remapping
!                 of data for irregular grids).
!       findlen:  returns length of a string ending in
!                 '$'.
!       grafinit: initializes & finalizes graphics (see
!                 above, under USE).
!       intchr:   reads integer number into a string.
!       realchr:  reads a real number into a string.
!       setclr:   sets color for all ncar purposes, by
!                 name.
!       setpat:   sets pattern for all ncar line-
!                 drawing purposes, by name.
!
!
!  ARGUMENT LISTINGS:
!
!  Contgraf1:
! *********************************************************************
! *  zplot(i,j)     the m1dim by m2dim array to be plotted            *
! *  m1dim          first dimension of zplot                          *
! *  m2dim          second dimension of zplot                         *
! *  xmn            value of x at left boundary of frame              *
! *  xmx            value of x at right boundary of frame             *
! *  ymn            value of y at lower boundary of frame             *
! *  ymx            value of y at upper boundary of frame             *
! *  kscale         = 1: plot zplot                                   *
! *                 = 2: plot log10(zplot)                            *
! *  kchz           = 0: conpack chooses contour levels; chzcon ir-   *
! *                      relevant                                     *
! *                 =-1: chzcon contains number of contour levels     *
! *                      desired (as a real number)                   *
! *                 =-2: chzcon contains desired spacing of contour   *
! *                      levels (the interval as a real number)       *
! *                 > 0: chzcon is a real array with kchz elements    *
! *                      giving the values at which contour lines     *
! *                      should be drawn                              *
! *  chzcon         = real number or real array; meaning depends on   *
! *                   value of "kchz" (see above)                     *
! *  krect          = 0: data lies on evenly-spaced rectangular grid; *
! *                      xij and yij ignored                          *
! *                 = 1: NCAR ezmap mapping (x=longitude, y=latitude) *
! *                      xij(1) and xij(2) are min and max longitude; *
! *                      yij(1) and yij(2) are min and max latitude;  *
! *                      both coordinates measured in degrees.        *
! *                 = 2: Polar coordinate mapping (x=radius, y=theta) *
! *                      xij(1) and xij(2) are min and max radius;    *
! *                      yij(1) and yij(2) are min and max theta;     *
! *                      angle is measured in degrees.                *
! *                 = 3: x=x(i) and y=y(j); xij and yij contain one-  *
! *                      dimensional arrays for the mapping           *
! *                 = 4: x=x(i,j) and y=y(i,j); xij and yij contain   *
! *                      two-dimensional arrays of the mapping        *
! *  xij            =(real array):  meaning depends on value of krect *
! *  yij            =(real array):  meaning depends on value of krect *
! *  node!          number of decimal places to use for legend labels *
! *                  (contgraf decides whether exponential or not)    *
! *  kcolor         = 0: for output to mono printer (color lines)     *
! *                 <>0: for output to screen only or color printer   *
! *                      (color fill)                                 *
! *  labelplt(i)    title of plot                                     *
! *  labelx(i)      x-axis label                                      *
! *  labely(i)      y-axis label                                      *
! *  labrun(i)      run label                                         *
! *  ntty           logical unit number for writing error messages    *
! *  niore!         2nd logical unit number for writing error messages*
! *  gwnd           a 4 element real array specifying the section of  *
! *                 the screen on which to plot.
! *                                                                   *
! *  note:                                                            *
! *                                                                   *
! *  (1) labelplt, labelx, labely, labrun must end with the character *
! *      "$".                                                         *
! *********************************************************************
!
!end documentation
!*dk colram
!**********************************************************************
!beg                        colram                                    *
!**********************************************************************
      subroutine colram(xcra, ycra, ncra, iaia, igia, naia)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ncra, naia
      integer, dimension(*) :: iaia, igia
      real, dimension(*) :: xcra, ycra
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ifll, nocl, i
      real :: step
C-----------------------------------------------
c
c......................................................................
c        colram fills in the area between contour lines in color.
c......................................................................
c
      ifll = 1
cncar
      call cpgeti ('NCL - No. of Cnt Lines', nocl)
      step = 45./nocl
cncar
c
c  If any of the area identifiers is negative, then do not fill
c
      do i = 1, naia
         if (iaia(i) < 0) ifll = 0
      end do
c
c  Otherwise, fill the area with the color defined by its area
c    identifier relative to group 3, the contour line group
c
      if (ifll .ne. 0) then
         ifll = 0
         do i = 1, naia
            if (igia(i) .eq. 3) ifll = iaia(i)
         end do
         if (ifll>0 .and. ifll<=nocl+1) then
cncar
            call gsfaci (nint((ifll - 1)*step) + 1)
            call gfa (ncra - 1, xcra, ycra)
cncar
         endif
      endif

      end subroutine colram

!**********************************************************************
!beg                        contgraf1                                 *
!**********************************************************************
      subroutine contgraf1(zplot, m1dim, m2dim, xmn, xmx, ymn, ymx, 
     1   kscale, kchz, chzcon, krect, xij, yij, nodec, kcolor, 
     2   labelplt, labelx, labely, labrun, labhead, ntty, niorec, gwnd)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vcpmpcm1
      USE Vcpmpcm2
      USE Vcpmpinf
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: m1dim, m2dim, kscale, kchz, krect, nodec, kcolor, 
     1   ntty, niorec
      real xmn, xmx, ymn, ymx
      character labelplt*(*), labelx*(*), labely*(*), labrun*(*), 
     1   labhead*(*)
      real, dimension(m1dim*m2dim) :: zplot
      real, dimension(*) :: xij, yij
      real, dimension(*) :: chzcon
      real, dimension(4) :: gwnd
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer :: i1, i2
      integer, parameter :: ICPIWL=1000, ICPRWL=5000, IAREAL=20000
      integer, dimension(13), parameter :: iasf = (/(1, i1=1,13)/)
      real, dimension(10), parameter :: 
     1   tksa = (/1.E36,1.E36,6.,1.E36,0.010,0.,1.E36,1.E36,0.,0./),
     2   tksb = (/1.E36,1.E36,6.,1.E36,0.,0.015,1.E36,1.E36,0.,0.010/)
      real, dimension(4), parameter :: 
     1   grph = (/0.0, 1.0, 0.0, 1.0/),
     2   grid = (/0.15, 0.85, 0.15, 0.85/)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(ICPIWL) :: iwrk
      integer, dimension(:), allocatable :: iama
      integer, dimension(10) :: iaia, igia
      integer, dimension(:), allocatable :: lclr
      integer :: nlblplt, nlblx, nlbly, nlabrun, nlabhead, isolidpat, 
     1   idashpat, i, j, index, j12dim, nocl, nchclv, istat
      real, dimension(4) :: gwnd2, lgnd = (/0.905, 1., 0., 0.9/)
      real, dimension(2) :: xlplot, ylplot
      real, dimension(:), allocatable :: zzplot
      real, dimension(ICPRWL) :: rwrk, xcra, ycra
      real :: xc1, xcm, yc1, ycn, zmin, zmax, zploti, dx, dy, clrstp, 
     1   ciu, clv
      character :: kbackclr*8, chclv*50
      character*50, dimension(:), allocatable :: llbs
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real , EXTERNAL :: xmask, colram
C-----------------------------------------------
c
c       external cpdrpl
c                                                      Workspace arrays
 
c......................................................................
c       contgraf plots a contour plot of zplot using the conpack
c    utility of the ncar graphics library v3.0, frames the plot,
c    and writes a legend.  it was cannibalized from the contgraf
c    routine of degraf 60.0 written by derek fox under supervision
c    of daren stotler, 8/92.  further modified and improved by
c    derek fox 1/93.
c
c    modified by mvachharajan in 8/93 to pass in window coords (gwnd(1:4))
c    calling arguments
c    -----------------
c
c    zplot(i,j)     the m1dim by m2dim array to be plotted
c    m1dim          first dimension of zplot
c    m2dim          second dimension of zplot
c    xmn            value of x at left boundary of frame
c    xmx            value of x at right boundary of frame
c    ymn            value of y at lower boundary of frame
c    ymx            value of y at upper boundary of frame
c    kscale         = 1: plot zplot
c                   = 2: plot log10(zplot)
c    kchz           = 0: conpack chooses contour levels; chzcon ir-
c                        relevant
c                   =-1: chzcon contains number of contour levels
c                        desired (as a real number)
c                   =-2: chzcon contains desired spacing of contour
c                        levels (the interval as a real number)
c                   > 0: chzcon is a real array with kchz elements
c                        giving the values at which contour lines
c                        should be drawn
c    chzcon         = real number or real array; meaning depends on
c                     value of "kchz" (see above)
c    krect          = 0: data lies on evenly-spaced rectangular grid;
c                        xij and yij ignored
c                   = 1: NCAR ezmap mapping (x=longitude, y=latitude)
c                        xij(1) and xij(2) are min and max longitude;
c                        yij(1) and yij(2) are min and max latitude;
c                        both coordinates measured in degrees.
c                   = 2: Polar coordinate mapping (x=radius, y=theta)
c                        xij(1) and xij(2) are min and max radius;
c                        yij(1) and yij(2) are min and max theta;
c                        angle is measured in degrees.
c                   = 3: x=x(i) and y=y(j); xij and yij contain one-
c                        dimensional arrays for the mapping
c                   = 4: x=x(i,j) and y=y(i,j); xij and yij contain
c                        two-dimensional arrays of the mapping
c    xij            =(real array):  meaning depends on value of krect
c    yij            =(real array):  meaning depends on value of krect
c    nodec          number of decimal places to use for legend labels
c                     (contgraf decides whether exponential or not)
c    kcolor         = 0: for output to mono printer (color lines)
c                   <>0: for output to screen only or color printer
c                        (color fill)
c    labelplt(i)    title of plot
c    labelx(i)      x-axis label
c    labely(i)      y-axis label
c    labrun(i)      run label (bottom of page)
c    labhead(i)     run label (top of page)
c    ntty           logical unit number for writing error messages
c    niorec         2nd logical unit number for writing error messages
c    gwnd           window coords for this plot ( to allow 4 per page )
c
c    note:
c
c    (1) labelplt, labelx, labely, labrun, labhead must end with the character
c        "$".
c......................................................................
c
c                       Set all GKS aspect source flags to 'individual'
cncar                                              and force solid fill
      
      call gsasf (iasf)
      call gsfais (1)
cncar                                    Find the length of all strings
c                                        <nlblplt,nlblx,nlbly,nlabrun,nlabhead>
      call findlen (labelplt, nlblplt)
      call findlen (labelx, nlblx)
      call findlen (labely, nlbly)
      call findlen (labrun, nlabrun)
      call findlen (labhead, nlabhead)
c
c                                Set color for grid, background, labels
c                      (also see the legend plotting routines, however)
      kbackclr = 'cyan'
c                                 Set the dash patterns for solid, dash
c                                                  <isolidpat,idashpat>
      isolidpat = 65535
      idashpat = 52428
c                                        Read arrays into common blocks
c                                          for non-rectangular mappings
      if (krect.eq.1 .or. krect.eq.2) then
         xc1 = xij(1)
         xcm = xij(2)
         yc1 = yij(1)
         ycn = yij(2)
      endif
c
      if (krect.eq.3 .or. krect.eq.4) then
         mxdim = m1dim
         nydim = m2dim
         xc1 = 1.
         xcm = real(m1dim)
         yc1 = 1.
         ycn = real(m2dim)
      endif
c
      if (krect .eq. 3) then
         comxi(:m1dim) = xij(:m1dim)
         comyj(:m2dim) = yij(:m2dim)
      endif
c
      if (krect .eq. 4) then
         allocate (comxij(m1dim,m2dim), comyij(m1dim,m2dim), 
     1             stat=istat)
         if (istat .ne. 0) stop 'comxij allocation error in contgraf1'
         do j = 1, m2dim
            do i = 1, m1dim
               index = (j - 1)*m1dim + i
               comxij(i,j) = xij(index)
               comyij(i,j) = yij(index)
            end do
         end do
      endif
c
c
      j12dim = m1dim*m2dim
      if (j12dim .gt. 0) then
         if (xmn .ge. xmx) go to 1012
         if (ymn .ge. ymx) go to 1016
c                                                          clear zzplot
         allocate (zzplot(j12dim), stat=istat)
         if (istat .ne. 0) stop 'zzplot allocation error in contgraf1'
         
         zzplot(:j12dim) = 0.0
c
         select case (kscale) 
c                                                 kscale=1: linear plot
         case default
            zzplot(:j12dim) = zplot(:j12dim)
            zmin = minval(zzplot(:j12dim))
            zmax = maxval(zzplot(:j12dim))
            go to 50
c                                                    kscale=2: log plot
         case (2) 
            zmin = 1.E30
            zmax = -1.E30
            do i = 1, j12dim
               zploti = abs(zplot(i))
               if (zploti .eq. 0.0) then
                  zzplot(i) = 0.0
               else
                  zzplot(i) = alog10(zploti)
                  zmin = min(zmin,zzplot(i))
                  zmax = max(zmax,zzplot(i))
               endif
            end do
c
            if (zmin .ge. zmax) go to 1035
c
            where (zzplot(:j12dim) .eq. 0.0) zzplot(:j12dim) = zmin
         end select
c                                                           plot zzplot
cncar
   50    continue
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call setclr (kbackclr)
c                                                       Draw the labels
c
         dx = .20                             !make = proutf:(lx2-lx1)/2
         dy = .20                             !make = proutf:(ly2-ly1)/2
         call plchhq (dx + gwnd(1), gwnd(4) - .030, 
     1      labelplt(1:nlblplt), 0.015/2, 0., 0.)
         call plchhq (dx + gwnd(1), 0.01 + gwnd(3), labelx(1:nlblx), 
     1      0.01/1.5, 0., 0.)
         call plchhq (0.02 + gwnd(1), dy + gwnd(3), labely(1:nlbly), 
     1      0.01/1.5, 90., 0.)
c
c                                                   Initialize the plot
         call cprset
         call cpsetc ('HLT - Hi/Lo lable Text string', '''')
         call cpsetc ('ILT - Info Lable Text string', '''')
c**********
c  To turn off labels change the following digit '3' to '0'
c**********
         call cpseti ('LLP - Line Label Positioning', 3)
c**********
         call cpsetr ('T2D - 2D Smoothing', 2.5)
c                                                 Turn on special value
         call cpsetr ('SPV - SPecial Value', 1.E36)
         call cpseti ('PAI - Param Array Index', -2)
         call cpseti ('CLU - Contour Level Use flag', 1)
         call cpsetr ('CLL - Contour Level Line width', 2.)
         call cpseti ('MAP - Type of Mapping', krect)
c
c                         Adjust size of contour grid  wrt plot window

         gwnd2(1) = gwnd(1) + (gwnd(2)-gwnd(1))*grid(1)
         gwnd2(2) = gwnd(1) + (gwnd(2)-gwnd(1))*grid(2)
         gwnd2(3) = gwnd(3) + (gwnd(4)-gwnd(3))*grid(3)
         gwnd2(4) = gwnd(3) + (gwnd(4)-gwnd(3))*grid(4)
 
         if (krect .ne. 0) then
            call cpsetr ('XC1 - X Coord at i=1', xc1)
            call cpsetr ('XCM - X Coord at i=M', xcm)
            call cpsetr ('YC1 - Y Coord at j=1', yc1)
            call cpsetr ('YCN - Y Coord at j=N', ycn)
            call cpseti ('SET - do SET call?', 0)
            call set (gwnd2(1), gwnd2(2), gwnd2(3), gwnd2(4), xmn, xmx, 
     1         ymn, ymx, 1)
         else
            call cpsetr ('VPL - View Portal Left', gwnd2(1))
            call cpsetr ('VPR - View Portal Right', gwnd2(2))
            call cpsetr ('VPB - View Portal Bottom', gwnd2(3))
            call cpsetr ('VPT - View Portal Top', gwnd2(4))
            call cpsetr ('VPS - View Portal Scale', 0.0)
         endif
c
         call cprect(zzplot,m1dim,m1dim,m2dim,rwrk,ICPRWL,iwrk,ICPIWL)
c
         if (kchz > 0) then
            call cpseti ('CLS - Cont Lvl Selection', 0)
            nocl = kchz
            call cpseti ('NCL - Number of Cont Lvls', nocl)
            do i = 1, nocl
               call cpseti ('PAI - Param Array Index', i)
               call cpsetr ('CLV - Cont Level Value', chzcon(i))
            end do
         else
            if (kchz .eq. (-1)) call cpseti ('CLS - Cont Lvl Selection', 
     1         nint((-chzcon(1))))
            if (kchz .eq. (-2)) then
               call cpsetr ('CIS - Cont Interval Specifier', chzcon(1))
               call cpsetr ('CMN - Cont Minimum', zmin)
               call cpsetr ('CMX - Cont Maximum', zmax)
               call cpseti ('CLS - Cont Lvl Selection', 1)
            endif
            call cppkcl (zzplot, rwrk, iwrk)
            call cpgeti ('NCL - Number of Cont Lvls', nocl)
         endif
c                                               Adjust colors, patterns
         if (kcolor .eq. 0) then
            clrstp = 45./(nocl - 1)
            do i = 1, nocl
               call cpseti ('PAI - Param Array Index', i)
               call cpseti ('CLU - Cont Lvl Use', 3)
               call cpseti ('CLC - Cont Lvl Color', nint((i - 1)*clrstp)
     1             + 1)
c       Alternate solid and dashed ?
               call cpseti ('CLD - Cont Lvl Dash Ptn', isolidpat)
c            if(mod(i,2).eq.1) then
c              call cpseti('CLD - Cont Lvl Dash Ptn',idashpat)
c            else
c              call cpseti('CLD - Cont Lvl Dash Ptn',isolidpat)
c            endif
            end do
         endif
c
c
c                                       Prepare autograph to frame plot
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call agseti ('FRAME.', 2)
         call anotat ('$', '$', 0, 0, 0, 0)
         call agsetp ('GRAPH.', gwnd, 4)
         call agsetp ('GRID.', grid, 4)
         call agseti ('X/NICE.', 0)
         call agsetf ('X/MIN.', xmn)
         call agsetf ('X/MAX.', xmx)
         call agseti ('Y/NICE.', 0)
         call agsetf ('Y/MIN.', ymn)
         call agsetf ('Y/MAX.', ymx)
         if (kcolor .ne. 0) then
            call agsetp ('LEFT/TICKS.', tksa, 10)
            call agsetp ('BOTTOM/TICKS.', tksa, 10)
            call agseti ('RIGHT/CONTROL.', -1)
            call agseti ('TOP/CONTROL.', -1)
         endif
c                                                        Frame the plot
         call agstup (0., 0, 0, 0, 0, 0., 0, 0, 0, 0)
         call agback
c                                            Restore autograph defaults
         if (kcolor .ne. 0) then
            call agsetp ('LEFT/TICKS.', tksb, 10)
            call agsetp ('BOTTOM/TICKS.', tksb, 10)
            call agseti ('RIGHT/CONTROL.', 4)
            call agseti ('TOP/CONTROL.', 4)
         endif
         call agseti ('X/NICE.', -1)
         call agsetf ('X/MIN.', 1.E36)
         call agsetf ('X/MAX.', 1.E36)
         call agseti ('Y/NICE.', -1)
         call agsetf ('Y/MIN.', 1.E36)
         call agsetf ('Y/MAX.', 1.E36)
         call agseti ('FRAME.', 1)
c
c                                            Set up for line labels/color
         allocate (iama(IAREAL), stat=istat)
         if (istat .ne. 0) stop 'iama allocation error in contgraf1'
         call arinam (iama, IAREAL)
         if (krect .ne. 0) call set (gwnd(1), gwnd(2), gwnd(3), gwnd(4), 
     1      xmn, xmx, ymn, ymx, 1)
c
         if (kcolor .eq. 0) then
c                                                     Add labels to area map
            call cplbam (zzplot, rwrk, iwrk, iama)
c                                                     Draw the contours
            call cpcldm (zzplot, rwrk, iwrk, iama, xmask)
c                                                      Draw line labels
cnolabels          call cplbdr(zzplot,rwrk,iwrk)
c                                                        (kcolor.ne.0):
         else
c                                         Add contour lines to area map
            call cpclam (zzplot, rwrk, iwrk, iama)
c                                                         Color the map
            call arscam (iama, xcra, ycra, ICPRWL, iaia, igia, 
     1           10, colram)
c
         endif
c

         deallocate (iama)
         
c                                                         Set up legend
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         lgnd(3) = lgnd(4) - 0.02*nocl - 0.015
c
         if (kcolor .ne. 0) clrstp = 45./nocl
c
c                                   Find proper decimal form for levels
c
         if (kchz <= 0) then
            call cpgetr ('CIU - Cont Interval Used', ciu)
         else
            ciu = (zmax - zmin)/real(nocl)
         endif
c
         if (ciu>100 .or. ciu<0.01) nodec = -nodec
c
c                                       Set legend label & color arrays
         allocate (lclr(nocl+1), llbs(nocl))
         
         do i = 1, nocl
            call cpseti ('PAI - Param Array Index', i)
            call cpgetr ('CLV - Cont Lvl Value', clv)
            call realchr (clv, nodec, chclv, nchclv)
            llbs(i) = chclv(1:nchclv)
            lclr(i) = nint((i - 1)*clrstp) + 1
         end do
         lclr(nocl+1) = 46
c                                                     Legend with lines
         if (kcolor .eq. 0) then
            xlplot(1) = lgnd(1) + 0.005
            xlplot(2) = xlplot(1) + 0.02
c                                                Legend with color fill
         else
c                          (color 16 is cyan -- see setclr for details)
            call lbseti ('CLB - Color of LaBels', 16)
            call lbseti ('CBL - Color of Box Lines', 16)
            call lblbar (1, lgnd(1), lgnd(2), lgnd(3), lgnd(4), nocl + 1
     1         , 0.33, 1., lclr, 1, llbs, nocl, 1)
         endif
c                                                 Restore some defaults
         call cpseti ('CLS - Cont Lvl Selection', 16)
         call dashdb (ior(isolidpat,0))
c                                                        Plot run label
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call setclr (kbackclr)
         call plchhq (0.5, 0.05, labrun(1:nlabrun), 0.008, 0., 0.)
         call plchhq (0.5, 0.95, labhead(1:nlabhead), 0.016, 0., 0.)
cncar
c
         call cpgeti ('IWU', i1)
         call cpgeti ('RWU', i2)
c       write (ntty,*) ' IWU,RWU= ',i1,i2
 
c----------------------------------------------------------------------
c                                                                errors
      endif                !!End j12dim > 0 loop
      

      go to 1099
      
!     write (ntty, 10010) labelplt(1:40), j12dim
!     write (niorec, 10010) labelplt(1:40), j12dim
!     go to 1099
 1012 continue
      write (ntty, 10012)
      write (niorec, 10012)
      go to 1099
 1016 continue
      write (ntty, 10016)
      write (niorec, 10016)
      go to 1099
 1035 continue
      write (ntty, 10035)
      write (niorec, 10035)
      call agback
 1099 continue

      if (allocated(zzplot)) deallocate(zzplot)
      if (allocated(lclr))   deallocate (lclr, llbs)
      if (allocated(comxij)) deallocate (comxij, comyij) 
 
c----------------------------------------------------------------------
c                                                     format statements
10010 format(1x,a,1x,'error in contgraf: m1dim*m2dim=',1x,i7,1x,
     1   '.eq.0.or.gt.500000')
10012 format(' error in contgraf1: xmn.ge.xmx')
10014 format(' error in contgraf1: delx.le.0')
10016 format(' error in contgraf1: ymn.ge.ymx')
10018 format(' error in contgraf1: dely.le.0')
10035 format(' no contour plot made: zmin.ge.zmax')

      end subroutine contgraf1
      

      SUBROUTINE XMASK (XCS,YCS,NCS,IAI,IAG,NAI)
C
C Routine for masked drawing of contour and grid lines
C
C This version of DRAWCL draws the polyline defined by the points
C ((XCS(I),YCS(I)),I=1,NCS) if and only if none of the area identifiers
C for the area containing the polyline are negative.  The dash package
C routine CURVED is called to do the drawing.
C
      DIMENSION XCS(*),YCS(*),IAI(*),IAG(*)
C
C Turn on drawing.
C
      IDR=1
C
C If any area identifier is negative, turn off drawing.
C
!     DO 101 I=1,NAI
!        IF (IAI(I).LT.0) IDR=0
!101  CONTINUE
C
C If drawing is turned on, draw the polyline.
C
      IF (IDR.NE.0) CALL CURVED (XCS,YCS,NCS)
C
C Done.
C
      RETURN
C
      END SUBROUTINE XMASK

!**********************************************************************
!beg                        contgraf2                                 *
!**********************************************************************
      subroutine contgraf2(zplot, m1dim, m2dim, xmn, xmx, ymn, ymx, 
     1   kscale, kchz, chzcon, krect, xij, yij, nodec, kcolor, 
     2   labelplt, labelx, labely, labrun, ntty, niorec, grid)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vcpmpcm1
      USE Vcpmpcm2
      USE Vcpmpinf
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m1dim,m2dim,kscale,kchz,krect,nodec,kcolor,ntty,niorec
      integer :: istat
      real xmn, xmx, ymn, ymx
      character labelplt*(*), labelx*(*), labely*(*), labrun*(*)
      real, dimension(*) :: zplot, chzcon, xij, yij
      real, dimension(4) :: grid
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer :: i1, i2
      integer, parameter :: ICPIWL=1000, ICPRWL=5000, IAREAL=20000
      integer, dimension(13), parameter :: iasf = (/(1, i1=1,13)/)
      real, dimension(10), parameter :: 
     1   tksa = (/1.E36,1.E36,6.,1.E36,0.010,0.,1.E36,1.E36,0.,0./),
     2   tksb = (/1.E36,1.E36,6.,1.E36,0.,0.015,1.E36,1.E36,0.,0.010/)
      real, dimension(4), parameter ::  gr01 = (/0.00, 1.0, 0.0, 1.0/)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(ICPIWL) :: iwrk
      integer, dimension(:), allocatable :: iama
      integer, dimension(10) :: iaia, igia
      integer, dimension(50) :: lclr
      integer :: i, nlblplt, nlblx, nlbly, nlabrun, isolidpat, idashpat
     1   , j, index, j12dim, nocl, nchclv
      real, dimension(4) :: grph, lgnd = (/0.905, 1., 0., 0.9/)
      real, dimension(2) :: xlplot, ylplot
      real, dimension(:), allocatable :: zzplot
      real, dimension(ICPRWL) :: rwrk, xcra, ycra
      real :: mr, mz, x1, x2, y1, y2, amap, xc1, xcm, yc1, ycn, zmin, 
     1   zmax, zploti, clrstp, ciu, clv
      character :: kbackclr*8, chclv*50
      character, dimension(50) :: llbs*50
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real , EXTERNAL :: cpdrpl, colram
C-----------------------------------------------
c
c                                                      Workspace arrays

c......................................................................
c       contgraf2 plots a contour plot of zplot using the conpack
c    utility of the ncar graphics library v3.0, frames the plot,
c    and writes a legend.  it was cannibalized from the contgraf2
c    routine of degraf 60.0 written by derek fox under supervision
c    of daren stotler, 8/92.  further modified and improved by
c    derek fox 1/93.
c
c    calling arguments
c    -----------------
c
c    zplot(i,j)     the m1dim by m2dim array to be plotted
c    m1dim          first dimension of zplot
c    m2dim          second dimension of zplot
c    xmn            value of x at left boundary of frame
c    xmx            value of x at right boundary of frame
c    ymn            value of y at lower boundary of frame
c    ymx            value of y at upper boundary of frame
c    kscale         = 1: plot zplot
c                   = 2: plot log10(zplot)
c    kchz           = 0: conpack chooses contour levels; chzcon ir-
c                        relevant
c                   =-1: chzcon contains number of contour levels
c                        desired (as a real number)
c                   =-2: chzcon contains desired spacing of contour
c                        levels (the interval as a real number)
c                   > 0: chzcon is a real array with kchz elements
c                        giving the values at which contour lines
c                        should be drawn
c    chzcon         = real number or real array; meaning depends on
c                     value of "kchz" (see above)
c    krect          = 0: data lies on evenly-spaced rectangular grid;
c                        xij and yij ignored
c                   = 1: NCAR ezmap mapping (x=longitude, y=latitude)
c                        xij(1) and xij(2) are min and max longitude;
c                        yij(1) and yij(2) are min and max latitude;
c                        both coordinates measured in degrees.
c                   = 2: Polar coordinate mapping (x=radius, y=theta)
c                        xij(1) and xij(2) are min and max radius;
c                        yij(1) and yij(2) are min and max theta;
c                        angle is measured in degrees.
c                   = 3: x=x(i) and y=y(j); xij and yij contain one-
c                        dimensional arrays for the mapping
c                   = 4: x=x(i,j) and y=y(i,j); xij and yij contain
c                        two-dimensional arrays of the mapping
c    xij            =(real array):  meaning depends on value of krect
c    yij            =(real array):  meaning depends on value of krect
c    nodec          number of decimal places to use for legend labels
c                     (contgraf2 decides whether exponential or not)
c    kcolor         = 0: for output to mono printer (color lines)
c                   <>0: for output to screen only or color printer
c                        (color fill)
c    labelplt(i)    title of plot
c    labelx(i)      x-axis label
c    labely(i)      y-axis label
c    labrun(i)      run label
c    ntty           logical unit number for writing error messages
c    niorec         2nd logical unit number for writing error messages
c
c  Output:
c    grid           the viewport computed internally based on plot
c                   dimensions
c    note:
c
c    (1) labelplt, labelx, labely, labrun must end with the character
c        "$".
c......................................................................
c
*     Assume proper aspect ratio requires shortening of x axis.
*     Fix Mz (Z window-span) .eq. 1. Let Dz,Dr be units-span for R,Z, resp.
*     Then demand that the correct Mr (R window-span) is given by
*     Mr = Dr/Dz*Mz
*     If Mr > 1. then reverse everything.
*     Assume the aspect ratio is mapped x:y = 1:1 to the printer.
 
      x1 = .1
      x2 = .9
      y1 = .1
      y2 = .9
      amap = 1.
      mr = (xmx - xmn)/(ymx - ymn)/amap
      if (mr <= 1.) then
         grph(1) = (x1 + x2 - mr*(x2 - x1))/2
         grph(2) = (mr*(x2 - x1) + (x1 + x2))/2
         grph(3) = y1
         grph(4) = y2
      else
         mz = (ymx - ymn)/(xmx - xmn)*amap
         grph(1) = x1
         grph(2) = x2
         grph(3) = (y1 + y2 - mz*(y2 - y1))/2
         grph(4) = (mz*(y2 - y1) + (y1 + y2))/2
      endif
      grid = grph
 
c                       Set all GKS aspect source flags to 'individual'
cncar                                              and force solid fill
      call gsasf (iasf)
      call gsfais (1)
cncar                                    Find the length of all strings
c                                         <nlblplt,nlblx,nlbly,nlabrun>
      call findlen (labelplt, nlblplt)
      call findlen (labelx, nlblx)
      call findlen (labely, nlbly)
      call findlen (labrun, nlabrun)
c
c                                Set color for grid, background, labels
c                      (also see the legend plotting routines, however)
      kbackclr = 'cyan'
c                                 Set the dash patterns for solid, dash
c                                                  <isolidpat,idashpat>
      isolidpat = 65535
      idashpat = 52428
c                                        Read arrays into common blocks
c                                          for non-rectangular mappings
      if (krect==1 .or. krect==2) then
         xc1 = xij(1)
         xcm = xij(2)
         yc1 = yij(1)
         ycn = yij(2)
      endif
c
      if (krect==3 .or. krect==4) then
         mxdim = m1dim
         nydim = m2dim
         xc1 = 1.
         xcm = real(m1dim)
         yc1 = 1.
         ycn = real(m2dim)
      endif
c
      if (krect .eq. 3) then
         comxi(:m1dim) = xij(:m1dim)
c
         comyj(:m2dim) = yij(:m2dim)
      endif
c
      if (krect .eq. 4) then
         allocate (comxij(m1dim,m2dim), comyij(m1dim,m2dim), 
     1             stat=istat)
         if (istat .ne. 0) stop 'comxij allocation error in contgraf1'
         do j = 1, m2dim
            do i = 1, m1dim
               index = (j - 1)*m1dim + i
               comxij(i,j) = xij(index)
               comyij(i,j) = yij(index)
            end do
         end do
      endif
c
c
      j12dim = m1dim*m2dim
      if (j12dim .gt. 0) then
         if (xmn >= xmx) go to 1012
         if (ymn >= ymx) go to 1016
c                                                          clear zzplot

         allocate (zzplot(j12dim), stat=istat)
         if (istat .ne. 0) stop 'zzplot allocation error in contgraf2'
         zzplot(:j12dim) = 0.0
c
         select case (kscale) 
c                                                 kscale=1: linear plot
         case default
            zzplot(:j12dim) = zplot(:j12dim)
            zmin = minval(zzplot(:j12dim))
            zmax = maxval(zzplot(:j12dim))
            go to 50
c                                                    kscale=2: log plot
         case (2) 
            zmin = 1.E30
            zmax = -1.E30
            do i = 1, j12dim
               zploti = abs(zplot(i))
               if (zploti .eq. 0.0) then
                  zzplot(i) = 0.0
               else
                  zzplot(i) = alog10(zploti)
                  zmin = min(zmin,zzplot(i))
                  zmax = max(zmax,zzplot(i))
               endif
            end do
c
            if (zmin >= zmax) go to 1035
c
            where (zzplot(:j12dim) .eq. 0.0) zzplot(:j12dim) = zmin
         end select
c                                                           plot zzplot
cncar
   50    continue
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call setclr (kbackclr)
c                                                       Draw the labels
c
         call plchhq (0.5, .95, labelplt(1:nlblplt), 0.015, 0., 0.)
         call plchhq (0.5, grph(3) - .04, labelx(1:nlblx), 0.01, 0., 0.)
         call plchhq (grph(1) - .04, .5, labely(1:nlbly), 0.01, 90., 0.)
c
c                                                   Initialize the plot
         call cprset
 
         call cpsetc ('HLT - Hi/Lo Label Text string', '''')
         call cpsetc ('ILT - Info Label Text string', '''')
         call cpseti ('LLP - Line Label Positioning', 0)
         call cpsetr ('T2D - 2D Smoothing', 2.5)
c                                                 Turn on special value
         call cpsetr ('SPV - SPecial Value', 1.E36)
         call cpseti ('PAI - Param Array Index', -2)
         call cpseti ('CLU - Contour Level Use flag', 1)
         call cpsetr ('CLL - Contour Level Line width', 2.)
c
c                         For irregular mappings, outline the grid edge
c                         Wie: remove because of spurious midplane line
cwie        if(krect.eq.1 .or. krect.eq.2 .or. krect.eq.4)then
cwie          call cpseti('PAI - Param Array Index',-1)
cwie          call cpseti('CLU - Contour Level Use flag',1)
cwie          call cpsetr('CLL - Contour Level Line width',2.)
cwie        endif
c
         call cpseti ('MAP - Type of Mapping', krect)
c
         if (krect .ne. 0) then
            call cpsetr ('XC1 - X Coord at i=1', xc1)
            call cpsetr ('XCM - X Coord at i=M', xcm)
            call cpsetr ('YC1 - Y Coord at j=1', yc1)
            call cpsetr ('YCN - Y Coord at j=N', ycn)
            call cpseti ('SET - do SET call?', 0)
            call set (grid(1), grid(2), grid(3), grid(4), xmn, xmx, ymn
     1         , ymx, 1)
         else
            call cpsetr ('VPL - View Portal Left', grid(1))
            call cpsetr ('VPR - View Portal Right', grid(2))
            call cpsetr ('VPB - View Portal Bottom', grid(3))
            call cpsetr ('VPT - View Portal Top', grid(4))
            call cpsetr ('VPS - View Portal Scale', 0.0)
         endif
c
         call cprect(zzplot,m1dim,m1dim,m2dim,rwrk,ICPRWL,iwrk,ICPIWL)
c
         if (kchz > 0) then
            call cpseti ('CLS - Cont Lvl Selection', 0)
            nocl = kchz
            call cpseti ('NCL - Number of Cont Lvls', nocl)
            do i = 1, nocl
               call cpseti ('PAI - Param Array Index', i)
               call cpsetr ('CLV - Cont Level Value', chzcon(i))
            end do
         else
            if (kchz .eq. (-1)) call cpseti ('CLS - Cont Lvl Selection', 
     1         nint((-chzcon(1))))
            if (kchz .eq. (-2)) then
               call cpsetr ('CIS - Cont Interval Specifier', chzcon(1))
               call cpsetr ('CMN - Cont Minimum', zmin)
               call cpsetr ('CMX - Cont Maximum', zmax)
               call cpseti ('CLS - Cont Lvl Selection', 1)
            endif
            call cppkcl (zzplot, rwrk, iwrk)
            call cpgeti ('NCL - Number of Cont Lvls', nocl)
         endif
c                                               Adjust colors, patterns
         if (kcolor .eq. 0) then
            clrstp = 45./(nocl - 1)
            do i = 1, nocl
               call cpseti ('PAI - Param Array Index', i)
               call cpseti ('CLU - Cont Lvl Use', 3)
               call cpseti ('CLC - Cont Lvl Color', nint((i - 1)*clrstp)
     1             + 1)
               if (mod(i,2) .eq. 1) then
                  call cpseti ('CLD - Cont Lvl Dash Ptn', idashpat)
               else
                  call cpseti ('CLD - Cont Lvl Dash Ptn', isolidpat)
               endif
            end do
         endif
c
c
c                                       Prepare autograph to frame plot
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call agseti ('FRAME.', 2)
         call anotat ('$', '$', 0, 0, 0, 0)
         call agsetp ('GRAPH.', gr01, 4)
         call agsetp ('GRID.', grid, 4)
         call agseti ('X/NICE.', 0)
         call agsetf ('X/MIN.', xmn)
         call agsetf ('X/MAX.', xmx)
         call agseti ('Y/NICE.', 0)
         call agsetf ('Y/MIN.', ymn)
         call agsetf ('Y/MAX.', ymx)
         if (kcolor .ne. 0) then
            call agsetp ('LEFT/TICKS.', tksa, 10)
            call agsetp ('BOTTOM/TICKS.', tksa, 10)
            call agseti ('RIGHT/CONTROL.', -1)
            call agseti ('TOP/CONTROL.', -1)
         endif
c                                                        Frame the plot
         call agstup (0., 0, 0, 0, 0, 0., 0, 0, 0, 0)
         call agback
c                                            Restore autograph defaults
         if (kcolor .ne. 0) then
            call agsetp ('LEFT/TICKS.', tksb, 10)
            call agsetp ('BOTTOM/TICKS.', tksb, 10)
            call agseti ('RIGHT/CONTROL.', 4)
            call agseti ('TOP/CONTROL.', 4)
         endif
         call agseti ('X/NICE.', -1)
         call agsetf ('X/MIN.', 1.E36)
         call agsetf ('X/MAX.', 1.E36)
         call agseti ('Y/NICE.', -1)
         call agsetf ('Y/MIN.', 1.E36)
         call agsetf ('Y/MAX.', 1.E36)
         call agseti ('FRAME.', 1)
c
c                                          Set up for line labels/color
         allocate (iama(IAREAL), stat=istat)
         if (istat .ne. 0) stop 'iama allocation error in contgraf2'
         call arinam (iama, IAREAL)
         if (krect .ne. 0) call set (grid(1), grid(2), grid(3), grid(4), 
     1      xmn, xmx, ymn, ymx, 1)
c
         if (kcolor .eq. 0) then
c                                                Add labels to area map
            call cplbam (zzplot, rwrk, iwrk, iama)
c                                                     Draw the contours
            call cpcldm (zzplot, rwrk, iwrk, iama, cpdrpl)
c                                                      Draw line labels
cnolabels          call cplbdr(zzplot,rwrk,iwrk)
c                                                        (kcolor.ne.0):
         else
c                                         Add contour lines to area map
            call cpclam (zzplot, rwrk, iwrk, iama)
c                                                         Color the map
            call arscam (iama, xcra, ycra, ICPRWL, iaia, igia, 
     1          10, colram)
c
         endif

         deallocate (iama)
c
c                                                         Set up legend
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         lgnd(3) = lgnd(4) - 0.02*nocl - 0.015
c
         if (kcolor .ne. 0) clrstp = 45./nocl
c
c                                   Find proper decimal form for levels
c
         if (kchz <= 0) then
            call cpgetr ('CIU - Cont Interval Used', ciu)
         else
            ciu = (zmax - zmin)/real(nocl)
         endif
c
         if (ciu>100 .or. ciu<0.01) nodec = -nodec
c
c                                       Set legend label & color arrays
         do i = 1, nocl
            call cpseti ('PAI - Param Array Index', i)
            call cpgetr ('CLV - Cont Lvl Value', clv)
            call realchr (clv, nodec, chclv, nchclv)
            llbs(i) = chclv(1:nchclv)
            lclr(i) = nint((i - 1)*clrstp) + 1
         end do
         lclr(nocl+1) = 46
c                                                     Legend with lines
         if (kcolor .eq. 0) then
            xlplot(1) = lgnd(1) + 0.005
            xlplot(2) = xlplot(1) + 0.02
c                                                Legend with color fill
         else
c                          (color 16 is cyan -- see setclr for details)
            call lbseti ('CLB - Color of LaBels', 16)
            call lbseti ('CBL - Color of Box Lines', 16)
            call lblbar (1, lgnd(1), lgnd(2), lgnd(3), lgnd(4), nocl + 1
     1         , 0.33, 1., lclr, 1, llbs, nocl, 1)
         endif
c                                                 Restore some defaults
         call cpseti ('CLS - Cont Lvl Selection', 16)
         call dashdb (ior(isolidpat,0))
c                                                        Plot run label
         call set (0., 1., 0., 1., 0., 1., 0., 1., 1)
         call setclr (kbackclr)
         call plchhq (0.5, 0.01, labrun(1:nlabrun), 0.008, 0., 0.)
cncar
c
         call cpgeti ('IWU', i1)
         call cpgeti ('RWU', i2)
c        write (ntty,*) ' IWU,RWU= ',i1,i2
 
c----------------------------------------------------------------------
c                                                                errors
      endif          !!end j12dim .ne. 0
      
      goto 1099
      
!     write (ntty, 10010) labelplt(1:40), j12dim
!     write (niorec, 10010) labelplt(1:40), j12dim
!     go to 1099
 1012 continue
      write (ntty, 10012)
      write (niorec, 10012)
      go to 1099
 1016 continue
      write (ntty, 10016)
      write (niorec, 10016)
      go to 1099
 1035 continue
      write (ntty, 10035)
      write (niorec, 10035)
 1099 continue

      if (allocated(zzplot)) deallocate(zzplot)
      if (allocated(comxij)) deallocate (comxij, comyij) 

c----------------------------------------------------------------------
c                                                     format statements
10010 format(1x,a,1x,'error in contgraf2: m1dim*m2dim=',1x,i7,1x,
     1   '.eq.0.or.gt.500000')
10012 format(' error in contgraf2: xmn.ge.xmx')
10014 format(' error in contgraf2: delx.le.0')
10016 format(' error in contgraf2: ymn.ge.ymx')
10018 format(' error in contgraf2: dely.le.0')
10035 format(' no contour plot made: zmin.ge.zmax')

      end subroutine contgraf2
      
!**********************************************************************
!                           cpmpxy                                    *
!**********************************************************************
      subroutine cpmpxy(imap, xin, yin, xout, yout)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vcpmpcm1
      USE Vcpmpcm2
      USE Vcpmpinf
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer imap
      real xin, yin, xout, yout
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j
C-----------------------------------------------
c......................................................................
c  cpmpxy is the utility used by CONPACK to remap data points prior to
c      contouring.  This version allows for four types of special
c      plots:
c
c      imap=1:  NCAR EZMAP mapping.  x is longitude in degrees, and y
c               is latitude.
c      imap=2:  Polar coordinates.  x is radius and y is polar angle,
c               in degrees.
c      imap=3:  Orthogonal, unevenly-spaced.  x=comxi(i) and
c               y=comyj(j).  Calling routine must fill these arrays.
c      imap=4:  Generalized distortion.  x=comxij(i,j) and
c               y=comyij(i,j).  Calling routine must fill the arrays.
c......................................................................
c
c
      xout = xin
      yout = yin
      if (imap .eq. 1) call maptrn (yin, xin, yout, xout)
c
      if (imap .eq. 2) then
         xout = xin*cos(.017453292519943*yin)
         yout = xin*sin(.017453292519943*yin)
      endif
c
      if (imap .eq. 3) then
         i = max(1,min(mxdim - 1,int(xin)))
         j = max(1,min(nydim - 1,int(yin)))
         xout=(real(i+1)-xin)*comxi(i)+(xin-real(i))*comxi(i+1)
         yout=(real(j+1)-yin)*comyj(j)+(yin-real(j))*comyj(j+1)
      endif
c
      if (imap .eq. 4) then
         if (.not.allocated(comxij)) 
     1       stop 'comxij NOT allocated in cpmpxy'
         i = max(1,min(mxdim - 1,int(xin)))
         j = max(1,min(nydim - 1,int(yin)))
         xout = (real(j + 1) - yin)*((real(i + 1) - xin)*comxij(i,j)+(
     1      xin-real(i))*comxij(i+1,j)) + (yin - real(j))*((real(i + 1)
     2       - xin)*comxij(i,j+1)+(xin-real(i))*comxij(i+1,j+1))
         yout = (real(j + 1) - yin)*((real(i + 1) - xin)*comyij(i,j)+(
     1      xin-real(i))*comyij(i+1,j)) + (yin - real(j))*((real(i + 1)
     2       - xin)*comyij(i,j+1)+(xin-real(i))*comyij(i+1,j+1))
      endif
 
      end subroutine cpmpxy

!**********************************************************************
!beg                         findlen                                  *
!**********************************************************************
      subroutine findlen(string, lstring)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lstring
      character string*(*)
C-----------------------------------------------
c......................................................................
c Finds the "$" that signals termination of strings and returns its
c    position, so that strings can be centered appropriately.
c......................................................................
c
c
cwie      length = len(string)
c
cwie      do  10 i=1,length
cwie        if (string(i:i).eq.'$') go to 20
cwie   10 continue
c
cwie   20 lstring = i-1
C
      lstring = len_trim(string) - 1       ! drop mandatory trailing "$"

      end subroutine findlen

!**********************************************************************
!beg                       grafinit                                   *
!**********************************************************************
      subroutine grafinit(kinit, input_id)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kinit
      character input_id*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: IERRF=6, LUNIT=2, IWKID=1
      INTEGER, PARAMETER :: IWKID2=2, LUNITPS=4, IWKID3=3
      INTEGER, PARAMETER :: NCGM=1, X11=8, PSMONO=23
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, option = 4
      real :: hue, red, green, blue, sat, vat
      character :: fname*80, fdum*80
C-----------------------------------------------
c
c......................................................................
c     grafinit initializes (kinit=1) and finalizes (kinit=2) the
c   graphical output.
c......................................................................
c
      select case (kinit) 
c
c                                              kinit=1: initialize ncar
cncar
      case default
         fname = 'gmeta.'//trim(input_id)
         print *,' Select one of the following options: '// 
     1       '(in addition to a gmeta file)'
         print *,' (1) Mono Postscript (2) X-Window (3) 1 & 2 '//
     1       '(4) Only a gmeta file (default)'        
         READ (*, '(A80)') FDUM
         IF (FDUM .NE. ' ') THEN
           READ (FDUM,*) option
         END IF

         if (option > 4 .or. option < 1) option = 4
         
!        Open GKS
         CALL GOPKS (IERRF, 0)
         call gesc ((-1391), 1, fname, 1, 1, fdum)

!        Open and activate a CGM workstation with ID = IWKID
         CALL GOPWK (IWKID, LUNIT, NCGM)
         CALL GACWK(IWKID)
 
!        Open and activate a X11 workstation with ID = IWKID2
         if (option.eq.2 .or. option .eq.3) then
            CALL GOPWK(IWKID2,0,X11)
            CALL GACWK(IWKID2)
         end if   

!        Open and activate a Monochrome PS file with ID = IWKID3
#ifndef RISC
         if (option.eq.1 .or. option.eq.3) then
            CALL NGSETC('ME', trim(input_id) // ".ps")
            CALL GOPWK(IWKID3,LUNITPS,PSMONO)
            CALL GACWK(IWKID3)
         end if   
#else
         if (option .eq. 1) then
            print *,' This NCARG version does not support PS option'
            option = 4
         else if (option .eq. 3) then
            print *,' This NCARG version does not support PS option'
            option = 2
         end if   
#endif

         call gstxfp (12, 2)
c                                                    Set up color table
c
c  (0=blck,1=violet,6=blue,16=cyan,26=green,36=yellow,41=orng,46=red;
c   51=magenta,61=white,and other numbers 1-60 are in-between colors)
c
         do i = 1, 60
            if (i <= 46) then
               hue = 270. - (i - 1)*6.
            else
               hue = 360. - (i - 46)*6.
            endif
            if (hue .gt. 360.) hue = 360.
            if (hue. lt. 0.) hue = 0.
            sat = 0.99
            vat = 0.99
            call hsvrgb (hue, sat, vat, red, green, blue)
            call gscr (IWKID, i, red, green, blue)
            if (option.eq.2 .or. option .eq.3)
     1      call gscr (IWKID2,i, red, green, blue)
         end do
 
         call gscr (IWKID, 0, 0., 0., 0.)
         call gscr (IWKID, 61, 1., 1., 1.)
         if (option.eq.2 .or. option .eq.3) then
            call gscr (IWKID2, 0, 0., 0., 0.)
            call gscr (IWKID2, 61, 1., 1., 1.)
         end if   
 
c
c                                              kinit=2: finish graphics
cncar
      case (2) 
         CALL GDAWK(IWKID)
         CALL GCLWK(IWKID)
         if (option.eq.2 .or. option .eq.3) then
            CALL GDAWK(IWKID2)
            CALL GCLWK(IWKID2)
         end if   
         if (option.eq.1 .or. option.eq.3) then
            CALL GDAWK(IWKID3)
            CALL GCLWK(IWKID3)
         endif

         CALL GCLKS

         if (option.eq.2 .or. option.eq.4) then
            fname = 'gmeta.'//trim(input_id)
            fdum  = trim(input_id) // '.ps'
            print *, ' To create Post Script file, type:'
            print *, '     ctrans -d ps.mono ' // trim(fname) //
     1      ' > '//trim(fdum)
         end if
         
      end select

      end subroutine grafinit

!**********************************************************************
!beg                         intchr                                   *
!**********************************************************************
      subroutine intchr(intno, chint, lchint)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer intno, lchint
      character chint*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nintno, nascii0, i, j
C-----------------------------------------------
c......................................................................
c intchr takes an integer and returns a character string chint with
c     intno taking up the first lchint places.
c......................................................................
c
      nintno = intno
c
      if (nintno .eq. 0) then
         chint = '0'
         lchint = 1
         go to 20
      endif
c
      nascii0 = ichar('0')
c
      chint(1:1) = ' '
      if (nintno < 0) chint(1:1) = '-'
      nintno = iabs(nintno)
c
      lchint = int(log10(real(nintno))) + 1
c
      do i = lchint, 1, -1
         j = lchint - i + 2
         chint(j:j) = char(nascii0 + mod(nintno,10**i)/10**(i - 1))
      end do
c
      if (chint(1:1) .eq. '-') then
         lchint = lchint + 1
      else
         chint = chint(2:lchint+1)
      endif
c
   20 continue

      end subroutine intchr

!**********************************************************************
!beg                         realchr                                  *
!**********************************************************************
      subroutine realchr(realno, nopl, chreal, lchreal)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nopl, lchreal
      real realno
      character chreal*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: inopl, iexpflag, nascii0, iexp, lint, i, j, lexp
      real :: zrealno, zfuzz
      character :: chexp*3
C-----------------------------------------------
c......................................................................
c  realchr takes a real number realno and returns a string chreal with
c     the number in character form in the first lchreal positions.
c     nopl is the number of decimal places desired; if nopl.lt.0 the
c     number will be in exponential form.
c......................................................................
c
      inopl = nopl
      zrealno = realno
      iexpflag = 0
      nascii0 = ichar('0')
      zfuzz = 1.E-10
c
      if (inopl < 0) then
         inopl = iabs(inopl)
         if (zrealno .eq. 0.) go to 5
         if (abs(zrealno) < 1.) zfuzz = -zfuzz
         iexp = int(log10(abs(zrealno)) + zfuzz)
         if (abs(zrealno) < 1.) iexp = iexp - 1
         zrealno = zrealno/10.**iexp
         iexpflag = 1
      endif
c
      if (inopl .eq. 0) zrealno = nint(zrealno)
c                                              digits to left of dec pt
    5 continue
      call intchr (int(zrealno), chreal, lint)
      lchreal = lint
c
      if (inopl .ne. 0) then
c                                  (because intchr never returns '-0':)
c
         if (zrealno<0. .and. zrealno>(-1.)) then
            chreal(2:lint+1) = chreal(1:lint)
            chreal(1:1) = '-'
            lint = lint + 1
         endif
c                                                        decimal places
         zrealno = abs(zrealno)
         chreal(lint+1:lint+1) = '.'
         if (inopl .ne. 1) then
            do i = 1, inopl - 1
               j = lint + i + 1
               chreal(j:j) = char(nascii0 + mod(int(zrealno*10.**i),10))
            end do
         endif
         j = lint + inopl + 1
         chreal(j:j) = char(nascii0 + mod(nint(zrealno*10.**inopl),10))
         lchreal = j
c
         if (iexpflag .ne. 0) then
c                                               exponent for exp form
            call intchr (iexp, chexp, lexp)
c                                          (this assumes use of PLCHHQ)
cncar
            chreal(lchreal+1:lchreal+6) = '*10''S'''
cncar
            do i = 1, lexp
               j = lchreal + 6 + i
               chreal(j:j) = chexp(i:i)
            end do
            lchreal = j
c
         endif
      endif

      end subroutine realchr

!**********************************************************************
!beg                      setclr                                      *
!**********************************************************************
      subroutine setclr(color)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vthrint
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: color
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ishade, nascii0, i, index, ncolor
      character :: clr*3
C-----------------------------------------------
c......................................................................
c     Emulates DISSPLA routine for setting color for lines & text.
c......................................................................
c
      ishade = 1
      nascii0 = ichar('0')
      do i = 4, len_trim(color)
         index = ichar(color(i:i))
         if (index>nascii0 .and. index<=nascii0+5) then
            ishade = index - nascii0
            exit 
         endif
      end do
c
      clr = color(1:3)
      ncolor = 61
c
      if (clr=='bla' .or. clr=='blk') ncolor = 0
      if (clr=='whi' .or. clr=='wht') ncolor = 61
c
      if (clr=='pur' .or. clr=='prp') ncolor = 0 + ishade
      if (clr .eq. 'blu') ncolor = 5 + ishade
      if (clr=='tur' .or. clr=='trq') ncolor = 10 + ishade
      if (clr=='cya' .or. clr=='cyn') ncolor = 15 + ishade
      if (clr .eq. 'bgr') ncolor = 20 + ishade
      if (clr=='gre' .or. clr=='grn') ncolor = 25 + ishade
      if (clr .eq. 'ygr') ncolor = 30 + ishade
      if (clr=='yel' .or. clr=='ylw') ncolor = 35 + ishade
      if (clr=='ora' .or. clr=='orn') ncolor = 40 + ishade
      if (clr .eq. 'red') ncolor = 45 + ishade
      if (clr=='mag' .or. clr=='mgn') ncolor = 50 + ishade
      if (clr=='vio' .or. clr=='vlt') ncolor = 55 + ishade
cncar
      call sflush
      call gsplci (ncolor)
      call gsfaci (ncolor)
      call gstxci (ncolor)
      ithrmj = ncolor
      ithrmn = ncolor
      ithrtx = ncolor
cncar

      end subroutine setclr

!**********************************************************************
!beg                      setclrpm                                    *
!**********************************************************************
      subroutine setclrpm(color)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: color
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ishade, nascii0, i, index, ncolor
      character :: clr*3
C-----------------------------------------------
c......................................................................
c     Emulates DISSPLA routine for setting color for polymarkers
c......................................................................
c
      ishade = 1
      nascii0 = ichar('0')
      do i = 4, len_trim(color)
         index = ichar(color(i:i))
         if (index>nascii0 .and. index<=nascii0+5) then
            ishade = index - nascii0
            exit 
         endif
      end do
c
      clr = color(1:3)
      ncolor = 61
c
      if (clr=='bla' .or. clr=='blk') ncolor = 0
      if (clr=='whi' .or. clr=='wht') ncolor = 61
c
      if (clr=='pur' .or. clr=='prp') ncolor = 0 + ishade
      if (clr .eq. 'blu') ncolor = 5 + ishade
      if (clr=='tur' .or. clr=='trq') ncolor = 10 + ishade
      if (clr=='cya' .or. clr=='cyn') ncolor = 15 + ishade
      if (clr .eq. 'bgr') ncolor = 20 + ishade
      if (clr=='gre' .or. clr=='grn') ncolor = 25 + ishade
      if (clr .eq. 'ygr') ncolor = 30 + ishade
      if (clr=='yel' .or. clr=='ylw') ncolor = 35 + ishade
      if (clr=='ora' .or. clr=='orn') ncolor = 40 + ishade
      if (clr .eq. 'red') ncolor = 45 + ishade
      if (clr=='mag' .or. clr=='mgn') ncolor = 50 + ishade
      if (clr=='vio' .or. clr=='vlt') ncolor = 55 + ishade
cncar
      call sflush
      call gspmci (ncolor)
cncar

      end subroutine setclrpm

!**********************************************************************
!beg                      setpat                                      *
!**********************************************************************
      subroutine setpat(pattern)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character pattern*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ipat
C-----------------------------------------------
c
c......................................................................
c     Emulates DISSPLA routine for setting dashed-line patterns.
c......................................................................
c
      ipat = 65535
c
      if (pattern .eq. 'blank') ipat = 0
      if (pattern .eq. 'solid') ipat = 65535
      if (pattern .eq. 'longdash') ipat = 61680
      if (pattern .eq. 'dash') ipat = 52428
      if (pattern .eq. 'dot') ipat = 43690
      if (pattern .eq. 'chndot') ipat = 58596
cncar
      call dashdb (ior(ipat,0))
      call agseti ('DASH/PAT/1.', ipat)
cncar
      call agback

      end subroutine setpat


      subroutine agpwrt(xpos, ypos, chrs, nchs, isiz, iori, icen)
!
!-------------------------------------------------------------------------
!
!     cf /usr/local/doc/ncarg/ ncargks3v0.doc page 431
!     inserting this routine is a way of using the Version 3.00
!     PLOTCHAR routine in place of the version 2.00 PWRITX
!     routines used in AUTOGRAPH. Called by AUTOGRAPH (not by user).
!     This is how one writes Greek in EZXY, etc.
!     This also clears up the EZXY problem of zero-width blanks
!
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nchs, isiz, iori, icen
      real xpos, ypos
      character chrs*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real :: csfl
C-----------------------------------------------
      call pcgetr ('CS - CONSTANT SPACING FLAG', csfl)
      if (icen .ne. 0) then
         call pcsetr ('CS - CONSTANT SPACING FLAG', 1.25)
      else
         call pcsetr ('CS - CONSTANT SPACING FLAG', 0.)
      endif
      call plchhq (xpos, ypos, chrs(1:nchs), .8*real(isiz), real(iori), 
     1   real(icen))
      call pcsetr ('CS - CONSTANT SPACING FLAG', csfl)

      end subroutine agpwrt
      
 
      subroutine dlm_parse(string, delimiter, n_tokens, tokens, 
     1   len_tokens)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n_tokens
      character string*(*), delimiter
      integer, dimension(*) :: len_tokens
      character, dimension(*) :: tokens*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n_tokens_max, i, isl, iii, ipos, i_run
C-----------------------------------------------
 
c     Written by R. Wieland
c
c     break up the string "string" into its component tokens
c     delimiter given specifically
c
c     input:
c       string
c             delimiter delimiter character to be used
c       n_tokens: max no. of tokens allowed
c     output:
c       n_tokens: no. of tokens returned
c       tokens:   char array of individual tokens returned
c       len_tokens: integer array containing length of
c           each token
 
 
      n_tokens_max = n_tokens
 
      do i = 1, n_tokens
         tokens(i) = ' '
      end do
 
      isl = len_trim(string)
      if (isl .eq. 0) then
         n_tokens = 0
         iii = len(tokens(1))
cobsolete call str$copy_r (tokens(1), iii, %REF(string))
         len_tokens(1) = 0
         return 
      endif
 
c     Tag on a trailing delimiter; remove it at the end
      string(isl+1:isl+1) = delimiter
 
      i = 1
      ipos = 1
      do while(i .eq. 1)
         i = index(string(ipos:),delimiter)
         ipos = ipos + 1
      end do
      i_run = ipos - 1
 
      n_tokens = 0
      do while(i_run<=isl .and. n_tokens<n_tokens_max)
         i = index(string(i_run:),delimiter) + i_run - 1
         n_tokens = n_tokens + 1
         tokens(n_tokens) = string(i_run:i-1)
         len_tokens(n_tokens) = len_trim(tokens(n_tokens))
         i_run = i + 1
      end do
 
c     Remove trailing delimiter at the end
      string(isl+1:isl+1) = ' '
 
      end subroutine dlm_parse
      

      subroutine wie_parse(string, n_tokens, tokens, len_tokens)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n_tokens
      character string*(*)
      integer, dimension(*) :: len_tokens
      character, dimension(*) :: tokens*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nwsp = 3
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n_tokens_max, i, j, isl, iii, ipos, i_run, ix, sl
      logical :: whitespace
      character :: set*3
      character, dimension(nwsp) :: byteset
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      integer , EXTERNAL :: gen_find_first_in_set
      logical , EXTERNAL :: iswhitespace
C-----------------------------------------------
 
      equivalence (set, byteset)
 
      byteset(1) = char(32)                      !space
      byteset(2) = char(9)                       !tab
      byteset(3) = char(0)                       !null
 
      n_tokens_max = n_tokens
      do i = 1, n_tokens
         tokens(i) = ' '
         len_tokens(i) = 0
      end do
 
      isl = len_trim(string)
      if (isl .eq. 0) then
         n_tokens = 0
         iii = len(tokens(1))
         tokens(1) = string(1:iii)
         len_tokens(1) = 0
         return 
      endif
 
      i = 1
      ipos = 1
      do while(i .eq. 1)      !UNTIL SOMETHING OTHER THAN sp OR tab OR nul
         i = gen_find_first_in_set(string(ipos:),set)
         ipos = ipos + 1
      end do
      i_run = ipos - 1
 
      n_tokens = 0
      do while(i_run<=isl .and. n_tokens<n_tokens_max)
         i = gen_find_first_in_set(string(i_run:),set) + i_run - 1
         if (.not.iswhitespace(string(i_run:i_run),byteset,nwsp)) then
            n_tokens = n_tokens + 1
c ... added by alan {
            if (i - 1 >= i_run) then
               tokens(n_tokens) = string(i_run:i-1)
               i_run = i + 1
            else
               tokens(n_tokens) = string(i_run:)
               i_run = isl + 1
            endif
            len_tokens(n_tokens) = len_trim(tokens(n_tokens))
         else
            i_run = i + 1
         endif
c ... }
      end do
 
      end subroutine wie_parse
      
 
      integer function gen_find_first_in_set (string, set)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character string*(*), set*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, lstring, lset
C-----------------------------------------------
 
 
 
c     Mimic Str$Find_First_In_Set
 
      gen_find_first_in_set = 0
 
      lstring = len_trim(string)
      lset = len_trim(set)
 
      if (lstring==0 .or. lset==0) return 
 
      do i = 1, lstring
         do j = 1, lset
            if (string(i:i) .eq. set(j:j)) go to 100
         end do
      end do
 
      return 
 
  100 continue
 
c     Found a character in SET
 
      gen_find_first_in_set = i
 
      end function gen_find_first_in_set
       
 
      logical function iswhitespace (letter, byteset, nwsp)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nwsp
      character letter
      character, dimension(*) :: byteset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j
      logical :: whitespace
C-----------------------------------------------
      whitespace = .FALSE.
      j = 1
      do while( .not.whitespace .and. j<=nwsp)
         whitespace = letter .eq. byteset(j)
         j = j + 1
      end do
      iswhitespace = whitespace

      end function iswhitespace
 
 
      subroutine lineselect(pattern, lineid)
!---------------------
! This routine is extracted from the Derek Fox NCARPLTS routine.
! It provides a simple way to specify different line types.
! Input: pattern = pattern name (cf below)
!                  if pattern name not kosher, do nothing.
!        lineid  = integer identifier (e.g., 1 for principal line,
!                  2 for 2nd line in ezmxy call)
! Side Effects:    The DASH/SELECTOR. value gets set. So if you
!                  intend to specify multiple line types, do them
!                  in ascending order
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lineid
      character pattern*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ipat, len
      character :: clineid*7
C-----------------------------------------------
 
      write (clineid, 10, err=100) lineid
   10 format(i7)
      call str_strip (clineid, len)
 
      ipat = 65535
      if (pattern .eq. 'blank') ipat = 0
      if (pattern .eq. 'solid') ipat = 65535
      if (pattern .eq. 'longdash') ipat = 61680
      if (pattern .eq. 'dash') ipat = 52428
      if (pattern .eq. 'dot') ipat = 43690
      if (pattern .eq. 'chndot') ipat = 58596
      call dashdb (ior(ipat,0))
 
      call agseti ('DASH/PAT/'//clineid(1:len)//'.', ipat)
      call agseti ('DASH/SELECTOR.', lineid)
 
  100 continue

      end subroutine lineselect
 
!
!...       This routine compresses a character string, deleting all blanks
!...       the compressed string is returned in the pkg it came in [STRING]],
!...       together with a label [INBLEN] telling how long it is. (Dick
!
      subroutine str_strip(string, inblen)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer inblen
      character string*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ilen, ind, j
C-----------------------------------------------
c
 
      ilen = len(string)
      inblen = 0
      if (ilen .eq. 0) return 
c
      ind = 0
      do j = 1, ilen
         if (string(j:j) .ne. ' ') then
            ind = ind + 1
            string(ind:ind) = string(j:j)
         endif
      end do
      if (ind < ilen) string(ind+1:) = ' '
 
      inblen = ind
 
      end subroutine str_strip
