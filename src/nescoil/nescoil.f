! ----------------------------------------------------------------------
      subroutine nescoil (loopcount)
!................................................................
      USE kind_spec
      USE Vvacuum2, ONLY: iota_edge, phip_edge, curpol
      USE Vvacuum3
      USE Vprecal1, ONLY: np
      USE Vbiopre
      USE Vculine1
      USE Vprecal4
      use SvdCtrl, ONLY : noaccuracy
      use LoopCtrl
      use Vbnorm, ONLY : extension
      use NumParams, ONLY : inesc
      use safe_open_mod
      implicit none
!................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer :: loopcount
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: istat
      real(rprec) :: unit, file, t1, t2, tbeg, tend
      character*10 :: date, time

c................................................................
c
c      The NESCOIL code is based on a method described in
c         [1]  P.Merkel,Nucl.Fusion 27 (1987) 867.
c         [2]  P.Merkel, In: Theory of Fusion Plasmas,Eds. A.Bondeson
c              E.Sindoni,and F.Troyon, Varenna,Italy,EUR 11336 EN,25-46.
c
c      Introducing two nested toroidally closed surfaces a current
c      on the outer surface is determined such that the normal
c      component of the magnetic field on the inner surface is
c      minimized.
c
c      The current carrying surface is mapped onto the unit square:
c      0< u <1, 0< v <1.
c
c
c         ^
c       1.|
c         |                  cup - poloidal current per period
c         |                  cut - toroidal current
c         |
c   ^     |
c   |   u |
c   |     |
c         |
c  cup    |
c         |
c         |
c         |
c       0 ----------------------------->
c         0                v           1.
c                     <----- cut
c
c
c
c................................................................
c     PARAMETERS READ IN FROM INPUT FILE:
c
c     1. DIMENSIONS:
c        nu    - number of poloidal meshpoints for coil surface
c        nv    - number of toroidal meshpoints for coil surface
c
c        nu1   - number of poloidal meshpoints for plasma surface
c        nv1   - number of toroidal meshpoints for plasma surface
c
c        mf    - number of poloidal modes for current potential
c        nf    - number of toroidal modes for current potential
c
c        md    - number of poloidal modes for plasma and coil surface shape
c        nd    - number of toroidal modes for plasma and coil surface shape
c
c        npol  - number of segments of a modular or helical filament
c        ntor  - number of filaments per period
c
c        np    - number of field periods
c
c     2. NESCOIL CONTROL PARAMETERS:
c        cut   - Net toroidal current +-1 or 0
c        cup   - Net poloidal current +-1 or 0
c
c        ibex  - An  external field can be added by supplying
c                a subroutine bexter with IBEX = 1
c
c     3. INFORMATION FROM VMEC PLASMA:
c        iota_edge - iota at plasma edge
c        phip_edge - phi_prime at plasma edge
c        curpol - Total poloidal current Amps per field period from VMEC
c
c     4. SVD CONTROL PARAMETERS:
c        1. MSTRT = Method + svdscan start if >1
c                 >=0:Berr, <0:Xerr, unless MSTEP <=0: Least square
c        2. MSTEP = Method + svdscan stepsize:
c                 <=0: LeastSquare, =0: use old f04abe (now SOLVER), no svd
c        3. MKEEP = svd/scan control: =0: svdscan, else keep |nkeep| wgts
c                 <0: write all weights to output
c        4. MDSPW = 2+exponent of dsur multiplying bfn,ben:
c                 <0: post-calculate Xerr svdscan and write to output
c        5. CURWT = Weight for surface current minimization
c                   Works ONLY in LSQ branch
c        6. TRGWT : Not implemented yet (PMV)

c................................................................
c        Secondary parameters:
c
c        md,nd - dimension of arrays: md>=mf,nd>=nf
c        nmax  - dimension of array:  nmax> 3*(nu+2)*nv
c
c................................................................
c      Set loop counter iloop for internal use
      iloop = loopcount

c     Strategy based on iloop:
c      1. iloop = 0 or  1 : single run or first call, allocate
c      2. iloop = 0 or -1 : single run or last  call, deallocate
c................................................................
c     Read the input file only in the first loop
c     Read input data from file specified on command line:
c     xnesopt input 
c     NOTE: SINCE ALL DIMENSIONS ARE ALLOCATABLE NOW,
c           YOU MUST FIRST READ INPUT FILE (TO GET DIMENSIONS)
c           BEFORE DOING ANYTHING ELSE
      if (iloop == 0 .or. iloop == 1) call read_nescoil_input
c     This ends first-loop tasks, rest is done in each loop
c................................................................
c
      write (inesc, '(a,i3,a)') '-----  Begin Nescoil run num ', iloop,
     1    ' -----'      

c     Now that all inputs are here, begin nescoil run
      call date_and_time(date,time)
      write (inesc, 100) date(5:6),date(7:8),date(1:4),
     1  time(1:2),time(3:4),time(5:6)
 100  format('DATE = ',a2,'-',a2,'-',a4,' ',' TIME = ',2(a2,':'),a2)
      call second0(tbeg)   !This one for timing the whole run
c
c     setup constants, arrays, and allocations based on nwire
c     uncomment for NT      USE NUMERICAL_LIBRARIES
c
      nw = ntor*np*(npol+1)
      istat=0
      if (.not.allocated(vx))
     1 allocate (vx(nw), vy(nw), vz(nw), dx(nw), dy(nw), dz(nw),
     1           xw(nw), yw(nw), zw(nw), curre(nw),
     2           cok(0:np - 1), sik(0:np - 1), stat=istat)
      curre = 0; xw = 0; yw = 0; zw = 0                                !!MUST be initialized (SPH)
      if (istat .ne. 0) stop 'allocation error nw in NESCOIL'

c     Do all precalculations
      call precal
      call second0(t2)
      write (inesc,"('PRECAL took ',g12.3,' sec')") t2-tbeg

c     Compute quantities on plasma surface
      write (inesc, 10) '----- Calling Surface_Plasma -----'
      call second0(t1)
      call surface_plas
      call second0(t2)
      write (inesc,"('Time in Surface_Plasma: ',g12.3,' sec')") t2-t1

c     Compute quantities on coil surface
      write (inesc, 10) '----- Calling Surface_Coil -----'
      call second0(t1)
      call surface_coil
      call second0(t2)
      write (inesc,"('Time in Surface_Coil: ',g12.3,' sec')") t2-t1

c................................................................
c     Solve boundary value problem
      write (inesc, 10) '----- Calling Solver -----'
      call second0(t1)
      call solver_nescoil
      call second0(t2)
      write (inesc,"('Time in Solver: ',g12.3,' sec')") t2-t1

c................................................................
c     Post-process solution for various answers
      write (inesc, 10) '----- Calling Surfcur_Diag -----'
      call second0(t1)
      call surfcur_diag
      call second0(t2)
      write (inesc,"('Time in Surfcur_Diag: ',g12.3,' sec')") t2-t1

      if(noaccuracy .eq. 0) then
         write (inesc, 10) '----- Calling Accuracy -----'
         call second0(t1)
         call accuracy
         call second0(t2)
         write (inesc,"('ACCURACY took ',g12.3,' sec')") t2-t1
      endif

      call second0(tend)
      write (inesc,"('ONE NESCOIL RUN took ',g12.3,' sec')")
     1    tend-tbeg

 10   format(a)
c................................................................
c     Tasks to be done in last loop only
      if (iloop .le. 0) call nescoil_cleanup
c................................................................


      end subroutine nescoil

! ----------------------------------------------------------------------
      subroutine nescoil_help
      use NumParams, only: inesc
!................................................................
c     purpose:
c     Provide some help to nescoil user
c     Write sample input file for first-time user.
c     Note: Make sure to change the write statements here when
c           you change any "read (iunit,*)" statements
c................................................................
      write (inesc,*) '--------------------------------------------'
      write (inesc,*)
     1    'To run nescoil, type: xnesopt Infile '
      write (inesc,*)
     1   'Note: bnorm info is always read from file bnorm, so'
      write (inesc,*)
     1   'You should link bnorm to your bnorm file by following:'
      write (inesc,*)
     1   ' ln -s your-bnorm bnorm'
      write (inesc,*) 'See nesinp.demo file for proper infile syntax'

      write (inesc,*) '---------------------------------------------'

      write (inesc,*) 'Control settings:'
      write (inesc,*)
      write (inesc,*)'OUTPUT controls w_*** to set what nescoil writes:'
      write (inesc,*) ' > 0 : higher value for more detailed info'
      write (inesc,*) ' < 0 : for unique choice'
      write (inesc,*) ' w_psurf: Plasma surface info'
      write (inesc,*) ' w_psurf: Coil   surface info'
      write (inesc,*) ' w_bnuv : Bnorm field info'
      write (inesc,*) ' w_jsurf: J surface current info'
      write (inesc,*) ' w_xerr : X error (displacement) info'
      write (inesc,*) ' w_svd  : SVD solution info'
      write (inesc,*)
      write (inesc,*) 'SVD controls for svd calculations:'
      write (inesc,*) ' mstrt: start svd wght >=0:Berr, <0:Xerr'
      write (inesc,*) ' mstep: svd scan step <=0: LSQ, =0:no svd scan'
      write (inesc,*) ' mkeep: Max svd wght, 0=all'
      write (inesc,*) ' mdspw: 2+power of dsur'
      write (inesc,*) ' curwt: 0<wght of Jtarget < 1'
      write (inesc,*) ' trgwt: Not available yet:'
      write (inesc,*)

c     Write a sample input file
      open(unit=9,file='nesinp.demo')

      write (9,*) '------ Grid Spatial Dimensions ----'
      write (9,*) 'nu, nv, nu1, nv1, npol, ntor'
      write (9,*) '64, 64, 64,  64,  64,   7'
      write (9,*)

      write (9,*) '------ Fourier Dimensions ----'
      write (9,*) 'mf, nf, md, nd (max in surf and bnorm files)'
      write (9,*) '8,  8,  10, 10'
      write (9,*)

      write (9,*) '------ Plasma information from VMEC ----'
      write (9,*) 'np,     iota_edge,       phip_edge,       curpol'
      write (9,*) '3,    0.46888094303,    -9.961507888E-2 ,  1.2'
      write (9,*)

      write (9,*) '------ Current Controls ----'
      write (9,*) 'cut,  cup,  ibex(=1,use fixed background coils)'
      write (9,*) '0.0,  1.0,    0'
      write (9,*)

      write (9,*) '------ SVD controls -----'
      write (9,*) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (9,*) ' 0,     0,      0,     4,    0.0,   0.0'
      write (9,*)

      write (9,*) '------ Output controls -----'
      write (9,*) 'w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd'
      write (9,*) ' 0,       0,       0,      0,       0,      0'
      write (9,*)

      write (9,*) '------ Plasma Surface ---- '
      write (9,*) 'Number of fourier modes in following table'
      write (9,*) '2'
      write (9,*) 'Table of fourier coefficients:'
      write (9,*) 'm   n  R(m,n),  Z(m,n),    Lamda(m,n)'
      write (9,*) '0   0   1.42     0.00         0.00'
      write (9,*) '1   0   0.32     0.59         0.00'
      write (9,*)

      write (9,*) '------ Current Surface ---- '
      write (9,*) 'Number of fourier modes in following table'
      write (9,*) '2'
      write (9,*) 'Table of fourier coefficients'
      write (9,*) 'm  n  R(m,n), Z(m,n)'
      write (9,*) '0   0   1.42     0.00         0.00'
      write (9,*) '1   0   0.42     0.69         0.00'

      close(9)

      end subroutine nescoil_help

! ----------------------------------------------------------------------
      subroutine read_nescoil_input
!................................................................
c                                                             01/01/89
c     purpose:
c     Read input file with ALL nescoil settings
c     Write ALL inputs into output file so all settings will be with
c     the answers for later use, and also the output file can be used
c     as the input file to recreate an old run later.
c
c     Some of these dimensions used to be in meshes file in old version.
c     Now they are all dynamically allocated arrays
c................................................................
      USE kind_spec
      use safe_open_mod
      use SvdCtrl, ONLY: mstrt, mstep, mkeep, mdspw, curwt, trgwt
      use OutCtrl
      use LoopCtrl
      use Vmeshes
      USE Vprecal1, ONLY: np
      USE Vvacuum1, ONLY: cr, cz, ms, ns
      USE Vvacuum2, ONLY: ms1, ns1, cr1, cz1, cl1,
     1    iota_edge, phip_edge, curpol
      USE Vvacuum3
      USE Vprecal1
      USE Vvacuum4, ONLY: cr2, cz2, ms2, ns2
      USE Vvacuum5, ONLY: cr3, cz3, ms3, ns3
      USE Vvacuum6, ONLY: cf, sf, msf, nsf
      use Vbnorm, ONLY : extension
      use NumParams
      implicit none
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer, parameter :: nescoil0 = 7
      integer :: m, n, ntotal, k, mr, nr, istat, numargs, iunit
      real(rprec) :: crin, czin, clin
      character*200 :: arg1

c ----------------------------------------------------------------------
c read data from command line, xnescoil nescin.extension
C-----------------------------------------------
      call getcarg(1, arg1, numargs)
      if (numargs .lt.1 .or. arg1 .eq. '-h' .or. arg1 .eq. '/h') then
         print *,
     1   ' Syntax:  xnescoil nescin.extension'
         stop 'Invalid command line'
      endif
      
      numargs = index(arg1,'nescin.')
      if (numargs .le. 0) then
         extension = trim(arg1)
      else
         extension = arg1(numargs+7:len_trim(arg1))   
      end if   

      iunit = nescoil0
      call safe_open(iunit, istat, trim(arg1), 'old', 'formatted')
      if (istat .ne. 0) then
         print *,' Type xnescoil -h for proper syntax'
         stop 'Error opening input file in nescoil'
      endif   

      call safe_open(inesc, istat, 'nescout.'//extension,
     1   'unknown', 'formatted')
      if (istat .ne. 0) stop 'Error opening nescout file'

c................................................................
c read spatial dimensions
      read (iunit, *, err=99)
      read (iunit, *, err=99)
      read (iunit, *, err=99) nu, nv, nu1, nv1, npol, ntor
      read (iunit, *, err=99)
      write(inesc, 10) '----- Grid Spatial Dimensions -----'
      write(inesc, 10) 'nu, nv, nu1, nv1, npol, ntor'
      write(inesc,"(6i6)")  nu, nv, nu1, nv1, npol, ntor

c................................................................
c read fourier dimensions
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) mf, nf, md, nd
      read (iunit, *, err=99)
      write (inesc, 10) '----- Fourier Dimensions -----'
      write (inesc, 10) 'mf, nf, md, nd'
      write (inesc,"(4i6)")  mf, nf, md, nd

      nmax  = 3*nv*(nu+2)+1
      mnd   = (md + 1)*(2*nd + 1)
      nuv   = nu*nv
      nuv1  = nu1*nv1
      nuvh  = nuv/2 + nu
      nuvh1 = nuv1/2 + nu1

c................................................................
c read plasma parameters from vmec
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) np, iota_edge, phip_edge, curpol
      read (iunit, *, err=99)
      write (inesc, 10) '----- Plasma information from VMEC -----'
      write (inesc, 10) 'np, iota_edge, phip_edge, curpol'
      write (inesc,"(i6,3g25.16)")  np, iota_edge, phip_edge, curpol

c................................................................
c read currents and nescoil controls

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) cut, cup, ibex
      read (iunit, *, err=99)
      write (inesc, 10) '----- Current Controls -----'
      write (inesc, 10) 'cut, cup, ibex'
      write (inesc,"(2g25.16,i6)")  cut, cup, ibex

c................................................................
c read svd and resonance control switches
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) mstrt,mstep,mkeep,mdspw,curwt,trgwt
      read (iunit, *, err=99)
      write (inesc, 10) '----- SVD controls -----'
      write (inesc, 10) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (inesc,"(4i6,2g25.16)")  mstrt,mstep,mkeep,mdspw,curwt,trgwt

c................................................................
c read detailed output writing control switches
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99)
     1    w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd
      read (iunit, *, err=99)
      write (inesc, 10) '----- Output controls -----'
      write (inesc, 10) 
     1  'w_psurf, w_csurf, w_bnuv, w_jsurf, w_xerr, w_svd'
      write (inesc,"(6i6)") w_psurf,w_csurf,w_bnuv,w_jsurf,w_xerr,w_svd

c................................................................
c Allocate surface info arrays
      if (.not.allocated(cr)) 
     1 allocate (cr(0:md,-nd:nd), cz(0:md,-nd:nd), cr1(0:md,-nd:nd),
     1          cz1(0:md,-nd:nd), cr2(0:md,-nd:nd), cz2(0:md,-nd:nd),
     2          cr3(0:md,-nd:nd), cz3(0:md,-nd:nd), cf(0:md,-nd:nd),
     3          sf(0:md,-nd:nd), cl1(0:md,-nd:nd))
      cr1(0:md,(-nd):nd) = zero
      cz1(0:md,(-nd):nd) = zero
      cl1(0:md,(-nd):nd) = zero

c................................................................
c   read  plasma boundary data
      ms1 = 0
      ns1 = 0
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) ntotal
      write (inesc, 10) '----- Plasma Surface -----'
      write (inesc, 10) 'Number of fourier modes in table'
      write (inesc, *)  ntotal

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      write (inesc, 10) 
     1   '----- Plasma boundary fourier coefficients  -----'
      write (inesc, 10) 
     1   '   m    n        R(m,n)     Z(m,n)    Lamda(m,n)'

      do k = 1, ntotal
         read (iunit, *,err=99) mr, nr, crin, czin, clin
         if (mr.lt.0 .or. mr.gt.md .or. nr.lt.(-nd) .or. nr.gt.nd) then
            write (inesc, 10) 'error(s) in input plasma fourier coeffs:'
            write (inesc, '(2(a,i4))') 'Found m = ', mr, ' and n = ', nr
            if (mr.lt.0) then
               write (inesc, '(a,i4,a)') 'Cannot have m = ', mr, ' < 0'
            endif
            if (mr.gt.md) then
               write (inesc, 20) 'm > md = ', md
               write (inesc, 20) 'change md to ', mr
            else
               if (nr.lt.(-nd)) write (inesc, 20) 'n < -nd = ', -nd
               if (nr.gt.nd) write (inesc, 20) 'n > nd = ', nd
               write (inesc, 20) 'change nd to ', abs(nr)
            endif
            write (inesc, 10) 'Correct nescoil input file and rerun'
            stop 'error in nescoil input file plas rzmn'
         endif
         write (inesc,"(2i4,3g20.10)") mr, nr, crin, czin, clin
         cr1(mr,nr) = crin
         cz1(mr,nr) = czin
         cl1(mr,nr) = clin
         ms1 = max(ms1,abs(mr))
         ns1 = max(ns1,abs(nr))
      end do

      read (iunit, *, err=99)

c................................................................
c read coil current surface data
c
      cr(0:md,(-nd):nd) = zero
      cz(0:md,(-nd):nd) = zero
c
      ms = 0
      ns = 0
      read (iunit, *,err=99)
      read (iunit, *,err=99)
      read (iunit, *,err=99) ntotal
      write (inesc, 10) '----- Coil Surface -----'
      write (inesc, 10) 'Number of fourier modes in table'
      write (inesc,*)  ntotal

      read (iunit, *,err=99)
      read (iunit, *,err=99)
      write (inesc, 10) '----- Coil surface fourier coefficients -----'
      write (inesc, 10) '    m    n         R(m,n)         Z(m,n)'

      do k = 1, ntotal
         read (iunit, *,err=99) mr, nr, crin, czin
         if (mr.lt.0 .or. mr.gt.md .or. nr.lt.(-nd) .or. nr.gt.nd) then
            write (inesc, 10) 'error(s) in input plasma fourier coeffs:'
            write (inesc, '(2(a,i4))') 'Found m = ', mr, ' and n = ', nr
            if (mr.lt.0) then
               write (inesc,'(a,i4,a)') 'Cannot have m = ', mr, ' < 0'
            endif
            if (mr.gt.md) then
               write (inesc, 20) 'm > md = ', md
               write (inesc, 20) 'change md to ', mr
            else
               if (nr.lt.(-nd)) write (inesc, 20) 'n < -nd = ', -nd
               if (nr.gt.nd) write (inesc, 20) 'n > nd = ', nd
               write (inesc, 20) 'change nd to ', abs(nr)
            endif
            write (inesc, 10) 'Correct nescoil input file and rerun'
            stop 'error in nescoil input file coil surface rzmn'
         endif
         write (inesc,"(2i4,2g20.10)") mr, nr, crin, czin
         cr(mr,nr) = crin
         cz(mr,nr) = czin
         ms = max(ms,abs(mr))
         ns = max(ns,abs(nr))
      end do

c................................................................
c     Take care of n=0 coeffs
      cr1(0:ms1,0) = .5_dp*cr1(0:ms1,0)
      cz1(0:ms1,0) = .5_dp*cz1(0:ms1,0)
      cr(0:ms,0) = .5_dp*cr(0:ms,0)
      cz(0:ms,0) = .5_dp*cz(0:ms,0)

c................................................................
c     Write global info based on settings into output
      write (inesc, 10) 
     1  '----- end inputs, begin outputs. Nescoil Version 1.0 -----'
      if( is_zero(cup) .and. is_zero(cut) ) then
         write (inesc, 10) '----- Solving for Saddle coils -----'
      else
         write (inesc, 10) '----- Solving for Modular coils -----'
      endif

      if( ibex .eq. 0 ) then
         write (inesc, 10) '----- No background coils used -----'
      else
         write (inesc, 10) '----- Background coils used -----'
      endif

 10   format(a)
 20   format(a, i4)

      return

c................................................................
c     error branch, help user
 99   call nescoil_help
      stop 'error while reading input file in  read_nescoil_input'

      end subroutine  read_nescoil_input

!----------------------------------------------------------------------
      subroutine precal
!................................................................
c                                                             01/01/89
c     purpose:
c     Calculate everything that does not need pot.coeffs, e.g.
c     surface and fft stuff
c................................................................
      use Vmeshes
      use NumParams
      use LoopCtrl
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vprecal4
      USE Vprecal5
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10, ONLY: conv1, sinv1
      USE Vfourier2
      implicit none
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: k, m, ku, n, kv, i
c................................................................
c
      i=0
      if (.not.allocated(trigs)) 
     1 allocate (trigs(2*nu), trigsv(2*nv), conv1(nv1,0:nd),
     1   sinv1(nv1,0:nd), comu(nu,0:md), simu(nu,0:md), factor(nuvh1),
     2   coh(nv), sih(nv), conv(nv,0:nd), sinv(nv,0:nd), coh1(nv1),
     3   sih1(nv1), comu1(nu1,0:md), simu1(nu1,0:md), eps(0:md + nd),
     4   stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL precal'

      alp  = pi2/np
      alu  = pi2/nu
      alv  = pi2/nv
      alvp = pi2/(nv*np)
      fnuv = one/nuv
      fnv  = one/nv
      alu1 = pi2/nu1
      alv1 = pi2/nv1
      alvp1 = pi2/(nv1*np)
      fnuv1 = one/nuv1
      fnv1  = one/nv1
      call fftfax (nu, ifax, trigs)
      call cftfax (nv, ifaxv, trigsv)
      do 10 k=0,np-1
         cok(k) = cos(alp* k)
         sik(k) = sin(alp* k)
   10 continue
      do 20 m=0,md
      do 20 ku=1,nu
         comu(ku,m)   = cos(m*alu*(ku-1))
         simu(ku,m)   = sin(m*alu*(ku-1))
   20 continue
      do 30 n=0,nd
      do 30 kv=1,nv
         conv(kv,n)   = cos(n*alv*(kv-1))
         sinv(kv,n)   = sin(n*alv*(kv-1))
   30 continue
      do 40 kv=1,nv
         coh(kv)= cos(alvp*(kv-1))
         sih(kv)= sin(alvp*(kv-1))
   40 continue
      do 50 m=0,md
      do 50 ku=1,nu1
         comu1(ku,m)   = cos(m*alu1*(ku-1))
         simu1(ku,m)   = sin(m*alu1*(ku-1))
   50 continue
      do 60 n=0,nd
      do 60 kv=1,nv1
         conv1(kv,n)   = cos(n*alv1*(kv-1))
         sinv1(kv,n)   = sin(n*alv1*(kv-1))
   60 continue
      do 70 kv=1,nv1
         coh1(kv)= cos(alvp1*(kv-1))
         sih1(kv)= sin(alvp1*(kv-1))
   70 continue
      do 80 i=1,nu1
         factor(i) = one
         factor(i+nuv1/2) = one
   80 continue
      do 90 i=nu1+1,nuv1/2
         factor(i) = two
   90 continue
         eps(0) = one/two
      do 100 i=1,md+nd
         eps(i) = one
  100 continue

      end subroutine precal


! ----------------------------------------------------------------------
      subroutine nescoil_cleanup
!................................................................
c     Purpose:
c     deallocate everything at end of all nescoil runs
c................................................................
      USE Vsurface1, ONLY: x, y, z, r, nuv
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu, ruu
      USE Vsurface5, ONLY: xuv, yuv, zuv, ruv
      USE Vsurface6, ONLY: xvv, yvv, zvv, rvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9, ONLY: x1, y1, z1, r1
      USE Vsurfac10, ONLY: x1u, y1u, z1u
      USE Vsurfac11, ONLY: x1v, y1v, z1v
      USE Vsurfac12, ONLY: r1uu, z1uu, r1u
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1
      USE Vbiopre
      USE Vculine1
      USE Vprecal3, ONLY: factor, eps
      USE Vprecal4
      use Vprecal5
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10
      USE Vdiagno2, ONLY: pot
      USE Vdiagno3
      USE Vvacuum1, ONLY: cr, cz
      USE Vvacuum2, ONLY: cr1, cz1, cl1
      USE Vvacuum4, ONLY: cr2, cz2
      USE Vvacuum5, ONLY: cr3, cz3
      USE Vvacuum6, ONLY: cf, sf
      USE Vfourier2, ONLY: trigs, trigsv
      USE Vbnorm, ONLY: bn_ext
      use NumParams, ONLY: inesc
      implicit none
c................................................................

c      from nescoil:
      deallocate (vx, vy, vz, dx, dy, dz, xw, yw, zw, curre, cok, sik)

c      from read_nescoil_input:
      deallocate (x, y, z, r, xu, yu, xv, ru, rv, yv, zu, zv,
     1   xuu, yuu, zuu, ruu, xuv, yuv, zuv, ruv, xvv, yvv, zvv, rvv,
     2   snx, sny, snz, xcur, ycur, zcur)

c      from precal:
      deallocate (trigs, trigsv,
     1   conv1, sinv1, conv, sinv, coh, sih, comu, simu,
     2   factor, coh1, sih1, comu1, simu1, eps)

c      from surface_plas:
      deallocate (x1, y1, z1, r1, x1u, y1u, z1u, x1v, y1v, z1v,
     1   r1uu, z1uu, r1u, snx1, sny1, snz1, dsur1)

c      from surface_coil:
      deallocate (cr, cz, cr1, cz1,
     4   cl1, cr2, cz2, cr3, cz3, cf, sf)

c      from various places:
      deallocate (pote, pot, bn_ext)

      write (inesc, '(a)') '---- Deallocated all nescoil arrays ----'
      close (inesc)

      end subroutine nescoil_cleanup

! ----------------------------------------------------------------------
      subroutine surface_plas
!................................................................
c     purpose:                                                01/01/89
c     Set up plasma surface from fourier coeffs
c     Calculate dsur, normals etc
c
c................................................................
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_psurf
      use LoopCtrl
      USE Vvacuum1
      USE Vvacuum2, ONLY: cr1, cz1, ms1, ns1
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
c   following line removed: these are all from Vmeshes, some compilers can't cope
c     USE Vprecal3, ONLY: nuvh1, nuv1, nu1, nv1    
      USE Vprecal4
      USE Vprecal5
c     USE Vprecal6, ONLY: nu       !  same here
      USE Vprecal7
      USE Vprecal8, ONLY: coh1, sih1
      USE Vprecal9, ONLY: comu1, simu1
      USE Vprecal10, ONLY: conv1, sinv1
      USE Vsurface9, ONLY: x1, y1, z1, r1
      USE Vsurfac10, ONLY: x1u, y1u, z1u
      USE Vsurfac11, ONLY: x1v, y1v, z1v
      USE Vsurfac12, ONLY: r1uu, z1uu, r1u
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1, fla
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, kv, ku, k
      real(rprec) :: cm, cn, cofp, cofm, sifp, sifm, u, v
c ----------------------------------------------------------------------
c   plasma surface
c
      i=0
      if (.not.allocated(x1)) 
     1 allocate (x1(nuvh1), y1(nuvh1), z1(nuvh1), r1(nuvh1),
     1          x1u(nuvh1), y1u(nuvh1), z1u(nuvh1), x1v(nuvh1),
     2          y1v(nuvh1), z1v(nuvh1), r1uu(nuvh1), z1uu(nuvh1),
     3          r1u(nuvh1), snx1(nuvh1), sny1(nuvh1), snz1(nuvh1),
     4          dsur1(nuvh1), stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL surface_plas'

         r1(:)   = zero
         z1(:)   = zero
         r1u(:)  = zero
         y1v(:)  = zero
         z1u(:)  = zero
         z1v(:)  = zero
         r1uu(:)  = zero
         z1uu(:)  = zero

      do 20 m = 0,ms1
      do 20 n = 0,ns1
         cm     = m*pi2
         cn     = n*pi2
         i      = 0
      do 20 kv = 1 , nv1/2 + 1
      do 20 ku = 1 , nu1
         i      = i + 1
         cofp   = comu1(ku,m)*conv1(kv,n)-simu1(ku,m)*sinv1(kv,n)
         cofm   = comu1(ku,m)*conv1(kv,n)+simu1(ku,m)*sinv1(kv,n)
         sifp   = simu1(ku,m)*conv1(kv,n)+comu1(ku,m)*sinv1(kv,n)
         sifm   = simu1(ku,m)*conv1(kv,n)-comu1(ku,m)*sinv1(kv,n)
         r1(i)   = r1(i)    +       cr1(m,n)*cofp + cr1(m,-n)*cofm
         z1(i)   = z1(i)    +       cz1(m,n)*sifp + cz1(m,-n)*sifm
         r1u(i)  = r1u(i)   - cm *( cr1(m,n)*sifp + cr1(m,-n)*sifm )
         y1v(i)  = y1v(i)   - cn *( cr1(m,n)*sifp - cr1(m,-n)*sifm )
         z1u(i)  = z1u(i)   + cm *( cz1(m,n)*cofp + cz1(m,-n)*cofm )
         z1v(i)  = z1v(i)   + cn *( cz1(m,n)*cofp - cz1(m,-n)*cofm )
         r1uu(i) = r1uu(i) -cm*cm*( cr1(m,n)*cofp + cr1(m,-n)*cofm )
         z1uu(i) = z1uu(i) -cm*cm*( cz1(m,n)*sifp + cz1(m,-n)*sifm )
   20 continue
c
      do kv = 1, nv1/2 + 1
         do ku = 1, nu1
            i = nu1*(kv - 1) + ku
            x1(i) = coh1(kv)*r1(i)
            y1(i) = sih1(kv)*r1(i)
            x1u(i) = coh1(kv)*r1u(i)
            y1u(i) = sih1(kv)*r1u(i)
            x1v(i) = coh1(kv)*y1v(i) - alp*y1(i)
            y1v(i) = sih1(kv)*y1v(i) + alp*x1(i)
            snx1(i) = y1u(i)*z1v(i) - z1u(i)*y1v(i)
            sny1(i) = z1u(i)*x1v(i) - x1u(i)*z1v(i)
            snz1(i) = x1u(i)*y1v(i) - y1u(i)*x1v(i)
            dsur1(i) = sqrt(snx1(i)*snx1(i)+sny1(i)*sny1(i)+
     1         snz1(i)*snz1(i))
            snx1(i) = snx1(i)/dsur1(i)
            sny1(i) = sny1(i)/dsur1(i)
            snz1(i) = snz1(i)/dsur1(i)
         end do
      end do

c     Compute total plasma surface area
      fla = sum(dsur1(1:nu1)) + 2.*sum(dsur1(nu1+1:nuvh1-nu1)) +
     1   sum(dsur1(nuvh1-nu1+1:nuvh1))
      write (inesc, '(a,1pe10.3)') 'Total plasma surface area = ', fla

c.....Write plasma surface info to output if asked to do so
      if(w_psurf>0 .or. w_psurf==-1) then
         write (inesc, 320)
         write (inesc, 330) (i,i=0,ms1)
         write (inesc, 340)
         do k = -ns1, -1
            write (inesc, 360) k, (cr1(i,0),i=0,ms1)
         end do
         write (inesc, 360) k, (2*cr1(i,k),i=0,ms1)
         do k = 1, ns1
            write (inesc, 360) k, (cr1(i,k),i=0,ms1)
         end do
         write (inesc, 370)
         write (inesc, 330) (i,i=0,ms1)
         write (inesc, 340)
         do k = -ns1, -1
            write (inesc, 360) k, (cz1(i,0),i=0,ms1)
         end do
         write (inesc, 360) k, (2*cz1(i,k),i=0,ms1)
         do k = 1, ns1
            write (inesc, 360) k, (cz1(i,k),i=0,ms1)
         end do
 320     format(/1x,'     Plasma Surface cr1(m,n)'/)
 330     format(1x,'      m ',50(i5,:,5x))
 340     format(1x,'    n   ')
 360     format(4x,i2,4x,50(f9.6,:,1x))
 370     format(/1x,'     Plasma Surface cz1(m,n)'/)
      endif

      if(w_psurf>1 .or. w_psurf==-2) then
         write (inesc, '(a,i4,a)') '---- plasma x,y,z,r(i) for i = 1,',
     1       nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") x1(i),y1(i),z1(i),r1(i)
         enddo
         write (inesc, '(a)') '---- end plasma x,y,z,r(i) ----'
      endif

      if(w_psurf>2 .or. w_psurf==-3) then
         write (inesc, '(a,i4,a)')
     1       '---- plasma dsur,snx,sny,snz for i = 1,',nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") dsur1(i),snx1(i),sny1(i),snz1(i)
         enddo
         write (inesc,'(a)') '---- end plasma dsur,snx,sny,snz ----'
      endif

      if(w_psurf>3 .or. w_psurf==-4) then
         write (inesc, '(a,i4,a)') 
     1     '---- plasma xu,yu,xv,yv(i) for i = 1,',nuvh1,' ----'
         do i = 1, nuvh1
            write (inesc,"(4g16.6)") x1u(i),y1u(i),x1v(i),y1v(i)
         enddo
         write (inesc, '(a)') '---- end plasma xu,yu,xv,yv(i) ----'
      endif

c     Compute bn_ext (external normal field)
      call bnfld

      end subroutine surface_plas

! ----------------------------------------------------------------------
      subroutine surface_coil
!................................................................
c                                                             01/01/89
c     purpose:
c     Set up coil surface from fourier coeffs
c     Calculate dsur, normals etc
c
c................................................................
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_csurf
      use LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal5
      USE Vprecal6
      USE Vprecal7, ONLY: conv, sinv, nd
      USE Vvacuum1, ONLY: cr, cz, ms, ns
      USE Vsurface1, ONLY: x, y, z, r, nuv
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu, ruu
      USE Vsurface5, ONLY: xuv, yuv, zuv, ruv
      USE Vsurface6, ONLY: xvv, yvv, zvv, rvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, kv, ku, k
      real(rprec) :: cm, cn, cofp, cofm, sifp, sifm
c ----------------------------------------------------------------------
c   outer surface   1

      i=0
      if (.not.allocated(x))
     1 allocate (x(nuv), y(nuv), z(nuv), r(nuv), xu(nuv), yu(nuv),
     1   xv(nuv), ru(nuv), rv(nuv), yv(nuv), zu(nuv), zv(nuv),
     2   xuu(nuv), yuu(nuv), zuu(nuv), ruu(nuv), xuv(nuv), yuv(nuv),
     3   zuv(nuv), ruv(nuv), xvv(nuv), yvv(nuv), zvv(nuv), rvv(nuv),
     4   snx(nuv), sny(nuv), snz(nuv), xcur(nuv), ycur(nuv), zcur(nuv),
     5   stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL surface_coil' 
      x(:) = zero
      y(:) = zero
      z(:) = zero
      r(:) = zero

      xu(:) = zero
      yu(:) = zero
      zu(:) = zero
      ru(:) = zero

      xv(:) = zero
      yv(:) = zero
      zv(:) = zero
      rv(:) = zero

      xuu(:) = zero
      yuu(:) = zero
      zuu(:) = zero
      ruu(:) = zero

      xuv(:) = zero
      yuv(:) = zero
      zuv(:) = zero
      ruv(:) = zero

      xvv(:) = zero
      yvv(:) = zero
      zvv(:) = zero
      rvv(:) = zero

      snx(:) = zero
      sny(:) = zero
      snz(:) = zero

      xcur(:) = zero
      ycur(:) = zero
      zcur(:) = zero

      do 20 m = 0,ms
      do 20 n = 0,ns
         cm     = m*pi2
         cn     = n*pi2
         i      = 0
      do 20 kv = 1 , nv
      do 20 ku = 1 , nu
         i      = i + 1
         cofp   = comu(ku,m)*conv(kv,n)-simu(ku,m)*sinv(kv,n)
         cofm   = comu(ku,m)*conv(kv,n)+simu(ku,m)*sinv(kv,n)
         sifp   = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
         sifm   = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
         r(i)   = r(i)    +           cr(m,n)*cofp + cr(m,-n)*cofm
         z(i)   = z(i)    +           cz(m,n)*sifp + cz(m,-n)*sifm
         ru(i)  = ru(i)   - cm *    ( cr(m,n)*sifp + cr(m,-n)*sifm )
         rv(i)  = rv(i)   - cn *    ( cr(m,n)*sifp - cr(m,-n)*sifm )
         zu(i)  = zu(i)   + cm *    ( cz(m,n)*cofp + cz(m,-n)*cofm )
         zv(i)  = zv(i)   + cn *    ( cz(m,n)*cofp - cz(m,-n)*cofm )
         ruu(i) = ruu(i)  - cm * cm*( cr(m,n)*cofp + cr(m,-n)*cofm )
         ruv(i) = ruv(i)  - cm * cn*( cr(m,n)*cofp - cr(m,-n)*cofm )
         rvv(i) = rvv(i)  - cn * cn*( cr(m,n)*cofp + cr(m,-n)*cofm )
         zuu(i) = zuu(i)  - cm * cm*( cz(m,n)*sifp + cz(m,-n)*sifm )
         zuv(i) = zuv(i)  - cm * cn*( cz(m,n)*sifp - cz(m,-n)*sifm )
         zvv(i) = zvv(i)  - cn * cn*( cz(m,n)*sifp + cz(m,-n)*sifm )
   20 continue
      do 50 kv = 1 , nv
      do 50 ku = 1,nu
         i      = nu*(kv-1)+ku
         x(i)   =   coh(kv) * r(i)
         y(i)   =   sih(kv) * r(i)
         xu(i)  =   coh(kv) * ru(i)
         yu(i)  =   sih(kv) * ru(i)
         xv(i)  =   coh(kv) * rv(i) - alp*y(i)
         yv(i)  =   sih(kv) * rv(i) + alp*x(i)
         yuu(i) =   sih(kv) * ruu(i)
         yuv(i) =   sih(kv) * ruv(i) + alp*xu(i)
         yvv(i) =   sih(kv) * rvv(i) + alp*(2.*xv(i)+alp*y(i))
         xuu(i) =   coh(kv) * ruu(i)
         xuv(i) =   coh(kv) * ruv(i) - alp*yu(i)
         xvv(i) =   coh(kv) * rvv(i) - alp*(2.*yv(i)-alp*x(i))
   50 continue
      do 60  i = 1 , nuv
         snx(i) = yu(i)*zv(i)-zu(i)*yv(i)
         sny(i) = zu(i)*xv(i)-xu(i)*zv(i)
         snz(i) = xu(i)*yv(i)-yu(i)*xv(i)
         xcur(i) = cut*xv(i)-cup*xu(i)
         ycur(i) = cut*yv(i)-cup*yu(i)
         zcur(i) = cut*zv(i)-cup*zu(i)
   60 continue

c.....Write coilsurface info to output if asked to do so
      if(w_csurf>0 .or. w_csurf==-1) then
         write (inesc, 320)
         write (inesc, 330) (i,i=0,ms)
         write (inesc, 340)
         do k = -ns, -1
            write (inesc, 360) k, (cr(i,0),i=0,ms)
         end do
         write (inesc, 360) k, (2*cr(i,k),i=0,ms)
         do k = 1, ns
            write (inesc, 360) k, (cr(i,k),i=0,ms)
         end do
         write (inesc, 370)
         write (inesc, 330) (i,i=0,ms)
         write (inesc, 340)
         do k = -ns, -1
            write (inesc, 360) k, (cz(i,0),i=0,ms)
         end do
         write (inesc, 360) k, (2*cz(i,k),i=0,ms)
         do k = 1, ns
            write (inesc, 360) k, (cz(i,k),i=0,ms)
         end do
 320     format(/1x,'     Coil Surface cr(m,n)'/)
 330     format(1x,'      m ',50(i5,:,5x))
 340     format(1x,'    n   ')
 360     format(4x,i2,4x,50(f9.6,:,1x))
 370     format(/1x,'     Coil Surface cz(m,n)'/)
      endif

      if(w_csurf>1 .or. w_csurf==-2) then
         write (inesc, '(a,i4,a)') '---- coil x,y,z,r(i) for i = 1,',
     1       nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(4g16.6)") x(i),y(i),z(i),r(i)
         enddo
         write (inesc, '(a)') '---- end coil x,y,z,r ----'
      endif

      if(w_csurf>2 .or. w_csurf==-3) then
         write (inesc, '(a,i4,a)') 
     1      '---- coil snx,sny,snz(i) for i = 1,',nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(3g16.6)") snx(i),sny(i),snz(i)
         enddo
         write (inesc, '(a)') '---- end coil snx,sny,snz ----'
      endif

      if(w_csurf>3 .or. w_csurf==-4) then
         write (inesc, '(a,i4,a)') 
     1     '---- coil xcur,ycur,zcur(i) for i = 1,',nuvh,' ----'
         do i = 1, nuvh
            write (inesc,"(3g16.6)") xcur(i),ycur(i),zcur(i)
         enddo
         write (inesc, '(a)') '---- end coil xcur,ycur,zcur ----'
      endif

      end subroutine surface_coil

! ---------------------------------------------------------------------
      subroutine solver_nescoil
!................................................................
c                                                             01/01/89
c     purpose:
c     Set up g and h matrices, and invert g * phi = h to get phi
c     Uses matrix and fourier routines
c     Ref. P. Merkel, Nucl Fus 27, 5, 1987, pg 867-871
c................................................................
c     Modified by PMV and SPH:
c     SVD scan, Xerr target, weighted optimization
c
c ---------------------------------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      use Vmeshes
      use NumParams
      use SvdCtrl
      use OutCtrl, ONLY: w_svd, w_xerr
      use LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vsurface9, ONLY: x1, y1, z1
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1
      USE Vsurfac14, ONLY: bex, bey, bez, ben
      USE Vsolver1
      USE Vdiagno2, ONLY: pot
      use Vbnorm, ONLY: bn_ext
      use Vsurfcurdiag
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, mf1, ndim,
     1   nbdim, mfnv, l, m, nfl, n, kl, k1,
     2   mstrta, mstepa, mkeepa, mdspwa, mexp, ifail, nsvdw, icmplxmin,
     3   icurvamax, icdmaxmin, icdavemin, iemnabs, iemnrms, iemnmax,
     4   iwemnabs, iwemnrms, iwemnmax, j
      real(rprec) ::
     1   t1, t2, err, stat, epsilon, fla, eabsmin, ermsmin,
     2   emaxmin, weabsmin, wermsmin, wemaxmin, cmplxmin, curvamax,
     3   cdmaxmin, cdavemin, bmodav, bmodavinv, potpotm, potpotmm,
     4   potsq, complexity, errmax, errabs, errrms, werrmax, werrabs,
     5   werrrms, epsm
c................................................................
c   x1,y1,z1 hold co-ordinates of nuvh1 points on half surface
c   snx1,sny1,snz1 hold normal vector and dsur1 holds surface element
c   bex1,bey1,bez1 hold b field and ben holds its normal component
c   a(nmax) and work vectors are used in matrix and fourier
c................................................................

      ndim  = nf+mf*(2*nf+1)
      nbdim = 2*(mf+1)*(2*nf+1)

      i=0
      if (.not.allocated(a))
     1 allocate (a(nmax), work(nmax), ab(ndim,ndim+1), result(ndim),
     1          bfx(nuvh1,nbdim), bfy(nuvh1,nbdim), bfz(nuvh1,nbdim),
     2          bfn(nuvh1,nbdim), bex(nuvh1), bey(nuvh1), bez(nuvh1),
     3          bpx(nuvh1), bpy(nuvh1), bpz(nuvh1), dumv(nuvh1),
     4          ben(nuvh1), daccu(nuvh1), bmag(nuvh1),
     5        pot(0:mf,-nf:nf), svdv(ndim,ndim), svdw(ndim), stat=i)
      if (i .ne. 0) stop 'allocation error in nescoil solver'

c     Zero all allocated arrays
      a(:) = zero
      work(:) = zero
      result(:) = zero
      bex(:) = zero
      bey(:) = zero
      bez(:) = zero
      bpx(:) = zero
      bpy(:) = zero
      bpz(:) = zero
      dumv(:) = zero
      ben(:) = zero
      daccu(:) = zero
      bmag(:) = zero
      svdw(:) = zero

      bfx(:nuvh1,:nbdim) = zero
      bfy(:nuvh1,:nbdim) = zero
      bfz(:nuvh1,:nbdim) = zero
      bfn(:nuvh1,:nbdim) = zero
      ab(:ndim,:ndim+1) = zero

      pot(0:mf,-nf:nf) = zero
      svdv(:ndim,:ndim) = zero
c.........................
c     Either zero or get background b field at all points  depending on
      if (ibex .ne. 0) then
         do i = 1, nuvh1
            call bexter (x1(i), y1(i), z1(i), bex(i), bey(i), bez(i))
         end do
      endif

c.........................
      write (inesc, 10) '---- Calling MATRIX+FOURIER ----'
      call second0(t1)
c.........................
      do k = 1, nuvh1

c.....   Zero a(i) vector every time
         a(:nmax) = zero

c        Calculate g(u',v', u,v) and h(u,v) of Eqs. 7 and 8
         call matrix (a, k)

c        Integral over u',v' of {g(u v, u' v') * sin(m u'+ n v')} in eq. 6
         call fourier (3)

c        Set the info from a vector into bfx,bfy,bfz
         mf1 = 2*(mf + 1)
         mfnv = mf1*nv
         l = 0

         do m = 0, mf
            nfl = -nf
            if (m .eq. 0) nfl = 1
            do n = nfl, nf
               i = mf1*n + 2*(m + 1)
               if (n .lt. 0) i = i + mfnv
               l = l + 1
               bfx(k,l) = -a(i       )*fnv
               bfy(k,l) = -a(i+  mfnv)*fnv
               bfz(k,l) = -a(i+2*mfnv)*fnv
               bfn(k,l) = snx1(k)*bfx(k,l) + sny1(k)*bfy(k,l)
     1                   + snz1(k)*bfz(k,l)
            end do
         end do

      end do
c.....End of big loop over index k

      call second0(t2)
      write (inesc,"('Time in MATRIX+FOURIER: ',g12.3,' sec')") t2-t1
c.........................

c.....Calculate bnormal from bexternal dot normal
         ben(:nuvh1) = snx1(:nuvh1)*bex(:nuvh1) + sny1(:nuvh1)*
     1      bey(:nuvh1) + snz1(:nuvh1)*bez(:nuvh1) + bn_ext(:nuvh1)
c
c     At this point, g_m_n(u,v) and h(u,v) are known:
c     bfn(i,l) is g_m_n(u,v) (Green's function) and
c     ben(i)   is h(u,v) (External B + homogeneous part from Ip,It)
c
c     bfn(p,mn) * phi(mn) gives bnorm field at plasma point p
c         due to phi(mn) of coilcurrent, and
c     ben(p) is external+vmec normal field on plasma surface
c
c     Do not change bfn, ben if you want to use them later
c     If you do change them temporarily, make sure
c     you change them back.

c*****************************************
c
c     Begin Berror or Xerror targets modifications
c                     P. Valanju, (PMV) Mar 99
         mstrta = abs(mstrt)
         mstepa = abs(mstep)
         mkeepa = abs(mkeep)
c        Extract correct exponent of dsur from MDSPW
         if (mdspw == 0) mdspw = 3        !Default gives dsur to power 1
         mdspwa = abs(mdspw)
         mexp = mdspwa - 2                !exponent in dsur**mexp
         if (mdspwa >= 2) then
            write (inesc, '(a,i4)') 'Using fac*dsur**', mexp
         else
            write (inesc, 10) 'Using no fac or dsur**mexp'
         endif
c        Make sure LSQ branch is chosen if curr density is targeted
         if (curwt > 0) mstep = -iabs(mstep)
c.................................................................
c      NESCOIL control settings explained:
c
c 1. MSTRT = Method + svdscan start if >1
c            >=0:Berr, <0:Xerr, unless MSTEP <=0: Least square
c 2. MSTEP = Method + svdscan stepsize:
c            <=0: LeastSquare, =0: use old f04abe (now SOLVER), no svd
c 3. MKEEP = svd/scan control: =0: svdscan, else keep |nkeep| wgts
c            <0: write all weights to output
c 4. MDSPW = 2+exponent of dsur multiplying bfn,ben:
c            <0: post-calculate Xerr svdscan and write to output
c 5. CURWT = Weight for surface current minimization
c            Works ONLY in LSQ branch
c 6. TRGWT : Not implemented yet (PMV)
c           = 1 : use uniform weight,
c             else Weighted (over plasma surface) error minimization
c           = wtmn: file contains mn coeffs of weight
c           = wtuv: file contains uv values of weights
c......................................
c            MSTRT  MKEEP  MSTEP
c      LSQ     -      -      0    ->  LSQ  with f04abe (solver)
c              -      n     -s    ->  LSQ  with svd=n weights
c             m>1   +-n     -s    ->  LSQ  with svdscan m,n,s
c
c      Berr    1      n      s    ->  Berr with svd=n weights
c             m>1   +-n      s    ->  Berr with svdscan
c
c      Berr   -1      n      s    ->  Berr with svd=n weights
c            -m>1   +-n      s    ->  Berr with svdscan
c.................................................................
c      Note: LSQ cannot run with Xerr since that requires many extra
c            inverse fourier transforms inside min_xerr_phimn.
c            We can add that if needed later, but
c            Solution for phi(m,n) directly in fourier space is faster.
c......................................

c      Note about MDSPW: (controls factor(i)*dsur**mexp) :
c      We calculate all things over half+1 grid points nuvh1
c      This is ok for all nuv points due to stellrator symmetry.
c      What this means is when we want to sum anything over whole fs,
c      we MUST USE factor(i) which weighs the edge points correctly
c      in a trapezoidal sum. This is true EVEN IF WE DO NOT USE DSUR.
c
c      Due to this, the matrix in the linear (svd/Berr or Xerr) case:
c      bfn(:nuvh1,ndim)* phi(ndim) = ben(:nuvh1)
c      should be modified to
c       bfn(:nuvh1,:)*sqrt(factor(:nuvh1)) and
c       ben(:nuvh1)*sqrt(factor(:nuvh1))
c         if RMS ERROR (root of sum of square err) is minimized,
c      or should be modified to
c       bfn(:nuvh1,:)*factor(:nuvh1) and
c       ben(:nuvh1)*factor(:nuvh1)
c         if ABS ERROR (sum of |err|) is minimized.
c
c      We will always minimize RMS ERROR, as in original nescoil.
c      This is what Least-sq method does, so Linear will do same.
c      However, old accuracy reported ABS ERROR, which was not
c      eactly what was minimized.
c      We will report both.
c
c      This would make NO difference if bfn were a square matrix so
c      the error would be zero, but it DOES make a difference when
c      the there are more plasma points than fourier modes in phi
c      so that an exact fit is not posiible. This is similar to
c      a weighted least-square minimization.

c------------------------------------------------------
         if (mstep<=0 .or. mstrt==0) then
c        The least-square nescoil branch using f04abe (solver) or svd:
            call second0(t1)

c......     Calculate ab matrix:
c           use MDSPWA  to control power of dsur,
c             CURWT   to control amount of J minimization
            call ab_matrix (ab, bfn, dumv, ndim, mdspwa, curwt)

c......
            if (mstep == 0) then
c           Solve using u-v solver fo4abe as in old nescoil:
c           This is good choice since no svd or svdscan is done
               write (inesc, 10) 'Case: Standard Nescoil LSQ solver'

               result(:ndim) = ab(:ndim,ndim+1)
               call solver(ab, result, ndim)    !Do not use old f04abe solver

               call second0(t2)
               write (inesc,
     1           "('Berr Least SQ took ',g12.3,' sec')") t2-t1
c              There will be no svdscan, since svd is not used here
               go to 123            !Only this one avoids svd or svdscan

            else                                 !MSTEP not 0
c           Solve least-sq matrix eq using svd:
c           MSTRT, MKEEP, MSTEP will be used later for svdscan
               write (inesc, 10) 'Case: Standard Nescoil LSQ with svd'
               call svd_solve (ndim, ndim, ndim, ndim, ab,
     1             ab(1,ndim+1), svdv, svdw, nsvdw)
               call second0(t2)
               write (inesc,
     1             "('Berr LSQ SVD took ',g12.3,' sec')") t2-t1
            endif

         else
c        The Linear branch with Berr or Xerr:
c................................
            if (mstrt > 0) then
c           Target Berr, use svd to solve problem:
c           MSTRT, MKEEP, MSTEP will be used later for svdscan
               call second0(t1)
               write (inesc,10) 'Case: Targeting Berror Linear with svd'

c              Modify bfn, ben by proper factors of dsur and fac
c              so they agree with what is done in lsq method
c              Temporarily use dumv for this here and later
c              when un-modifying bfn, ben
               if (mdspwa >= 2) then
                  dumv(:nuvh1)=sqrt(factor(:nuvh1)*dsur1(:nuvh1)**mexp)
                  ben(:nuvh1) =ben(:nuvh1)*dumv(:nuvh1)
                  do i = 1, ndim
                     bfn(:nuvh1,i) = bfn(:nuvh1,i)*dumv(:nuvh1)
                  end do
               endif

c              Solve the linear svd problem
               call svd_solve (nuvh1, ndim, nuvh1, nbdim, bfn, -ben,
     1            svdv, svdw, nsvdw)

c              UnModify bfn, ben by proper factors of dsur and fac
c              so they can be used again
c              Use dumv calculated earlier
               if (mdspwa >= 2) then
                  ben(:nuvh1) = ben(:nuvh1)/dumv(:nuvh1)
                  do i = 1, ndim
                     bfn(:nuvh1,i) = bfn(:nuvh1,i)/dumv(:nuvh1)
                  end do
               endif
c              Bfn and Ben are again B fields at plasma points,
c              not b*ds. So they can be used to calc errors etc.

               call second0(t2)
               write (inesc,
     1             "('Berr Linear SVD took ',g12.3,' sec')") t2-t1
            endif

c................................
            if (mstrt < 0) then
c           Target Xerr, use svd to solve problem:
c           Note: we cannot use non-uniform weights in this case,
c                 since calculating xerr itself corresponds to using
c                 non-uniform weights on plasma surface
               call second0(t1)
               write (inesc,10) 'Case: Targeting Xerror Linear with svd'
               epsilon = 0

               call min_xerr_phimn (1, epsilon, nu1, nv1, nuvh1,
     1            nbdim, ndim, bfn, ben, dsur1, svdv, svdw, nsvdw)
               call second0(t2)
               write (inesc,
     1             "('Xerr Linear SVD took ',g12.3,' sec')") t2-t1
            endif

         endif

c------------------------------------------------------
c        At this point, the solution has been found (in some way)
c        All that remains is reporting the errors
c......................................................
c        Report svd errors, single if MSTRTA=1, scan if MSTRTA > 1:

c        Limit MKEEPA, use all weights if set to zero:
         if (mkeepa == 0) mkeepa = nsvdw
         mkeepa = min(nsvdw,mkeepa)

c        All branches arrive here with array of answers in svdv
c        At this point, we can do one of 3 three things:
c        1. Scan all svd weights and calculate err, jmax, complexity ...
c        2. Just use one of the weights (MKEEPA) specified by user
c        3. Find the optimum number of weights (isvdw) to keep
c           using ?? some ?? criterion and solutions svdv

         if (mstrta > 1) then
c        SVDSCAN: Find error, jmax, complexity, and anything else
c        for i=method,MKEEPA svd wgts kept

c....................................................................
c        Write all titles for scan (include identifiers for grep):
            write (inesc, '(4(a,i4),a)')
     1          '---- Svdscan from ',mstrta,' to ',mkeepa,' at ',
     2          mstepa,' steps out of ',ndim,' total wgts ----'

c           do not call accuracy later in nescoil (main) at each wght
            noaccuracy = 1

c           Calculate total area of plasma surface
            fla = sum(dsur1(:nu1)) + 2.*sum(dsur1(nu1+1:nuvh1-nu1))
     1          + sum(dsur1(nuvh1-nu1+1:nuvh1))

c           Set starting values for min/max error calculation in svdloop
            eabsmin = 1.0e14_dp                  !for min values of various errors
            ermsmin = 1.0e14_dp                  !min of rms err
            emaxmin = 1.0e14_dp                  !min of max err
            weabsmin = 1.0e14_dp                 !min of abs err
            wermsmin = 1.0e14_dp                 !weighted errors
            wemaxmin = 1.0e14_dp
            cmplxmin = 1.0e14_dp                 !min of complexity
            curvamax = zero                      !max of current curvature radius
            cdmaxmin = 1.0e14_dp                 !min of current density max
            cdavemin = 1.0e14_dp                 !min of current density ave

c...svdscan loop........................................................
            do i = mkeepa, mstrta, -mstepa     !scan in decreasing order
c           Go over all weights from best to worst, because
c           the calculation of bmod is done only first time, and
c           the best weight will give the best estimate for it.

c......        Calculate Bmag using green's functions rather than besucu
c              This is much faster since green's functions already known
               Bmag(:nuvh1) = sqrt(
     1         (bex(:nuvh1)+Matmul(bfx(:nuvh1,:ndim),svdv(:ndim,i)))**2
     2        +(bey(:nuvh1)+Matmul(bfy(:nuvh1,:ndim),svdv(:ndim,i)))**2
     3        +(bez(:nuvh1)+Matmul(bfz(:nuvh1,:ndim),svdv(:ndim,i)))**2)
c              Calculate fs average of Bmag
               dumv(:nuvh1) = dsur1(:nuvh1) * Bmag(:nuvh1)
               bmodav = ( sum(dumv(:nu1))
     1        + 2*sum(dumv(nu1+1:nuvh1-nu1))
     2        + sum(dumv(nuvh1-nu1+1:nuvh1)) )/ fla
               bmodavinv = one/bmodav

c.......       Calculate and write the maximum current density
c              Put svdv into phi(m,n) for sucudiag
c              Also calculate 1,2 moments for coil complexity
c              Save all pot(m,n) for later use by postprocessor.
               l = 0
               m = 0
               potpotm = zero       !for moments of pot(m,n) for complexity
               potpotmm = zero

c              Set and Write pot coeffs for this svd weight
               write (inesc, '(a,i4,a)') '---- Phi(m,n) for isvd = ', 
     1            i,'----'               
               do n = 1, nf
                  l = l + 1
                  pot(0,n)  = -svdv(l, i)
                  pot(0,-n) =  svdv(l, i)
               end do
               do n = -nf, nf
                  write (inesc,"(2i6,g25.16)") 0, n, pot(0,n)/two
               end do
               do m = 1, mf
                  do n = -nf, nf
                     l = l + 1
                     pot(m,n) = -svdv(l,i)
                     write (inesc,"(2i6,g25.16)") m, n, pot(m,n)
                     potsq = pot(m,n) * pot(m,n) * m
                     potpotm = potpotm + potsq
                     potpotmm = potpotmm + potsq*m
                  end do
               end do
               write (inesc, 10) '---- end Phi(m,n) ----'

               complexity = potpotmm/potpotm
               if (complexity < cmplxmin) then
                  cmplxmin = complexity
                  icmplxmin = i
               endif

               call surfcur_diag

               if (curvamax < curvra) then
                  curvamax = curvra
                  icurvamax = i
               endif
               if (cdmaxmin > cudema) then
                  cdmaxmin = cudema
                  icdmaxmin = i
               endif
               if (cdavemin > avcude) then
                  cdavemin = avcude
                  icdavemin = i
               endif

c......        Calculate errors for bnorm with <bmodavinv>, i.e.:
c              <|berr|>/<Bmod>, sqrt<|berr|^2>/<Bmod>, and
c              max<|berr|>/<Bmod>, where <> is fs average
c
c              1. Calculate unweighted residual |BFN*PHI+BEN|
               daccu(:nuvh1) = abs( ben(:nuvh1) +
     1            Matmul(bfn(:nuvh1,:ndim),svdv(:ndim,i)) )
c              2. Calc fs averages and max of residual
               dumv(:nuvh1) = daccu(:nuvh1)   !temporary storage
               call calc_percent_errs (dumv, dsur1, fla, nu1, nuvh1,
     1            errmax, errabs, errrms)
c              3. Calc globally weighted errors (with <Bmag>)
               errabs = errabs*bmodavinv
               errrms = errrms*bmodavinv
               errmax = errmax*bmodavinv

c.......       Calculate errors for locally wgted bnorm/bmod:
c              <|berr|/Bmod>, sqrt<|berr/Bmod|^2>>, and
c              max<|berr/Bmod|>, where <> is fs average
c
c              1. Calculate locally weighted residual |BFN*PHI-BEN|/Bmag
               dumv(:nuvh1) = daccu(:nuvh1)/bmag(:nuvh1)
c              2. Calc fs averages and max of locally weighted errors
               call calc_percent_errs (dumv, dsur1, fla, nu1, nuvh1,
     1            werrmax, werrabs, werrrms)

c.......       Write errors and j to output
               write (inesc, 10) 
     1            'Info for isvd, Berrors( abs, rms, max):'
               write (inesc, '(a,i4,1p3e10.3)') 
     1             '<|g*phi+h|> ',i, errabs, errrms, errmax
               write (inesc, '(a,i4,1p3e10.3)')
     1             '<|g*phi+h|/|B|> ', i, werrabs, werrrms, werrmax
               write (inesc, 10)
     1             'isvd, Jmax, Jave, J curv rad, Complexity'
               write (inesc, '(a,i4,1p4e10.3)')
     1             'J ',i, cudema, avcude, curvra, complexity

c.......       Find svd num and value of various min errors
               if (errabs < eabsmin) then
                  eabsmin = errabs
                  iemnabs = i
               endif
               if (errrms < ermsmin) then
                  ermsmin = errrms
                  iemnrms = i
               endif
               if (errmax < emaxmin) then
                  emaxmin = errmax
                  iemnmax = i
               endif

c              Find svd num and value of various weighted min errors
               if (werrabs < weabsmin) then
                  weabsmin = werrabs
                  iwemnabs = i
               endif
               if (werrrms < wermsmin) then
                  wermsmin = werrrms
                  iwemnrms = i
               endif
               if (werrmax < wemaxmin) then
                  wemaxmin = werrmax
                  iwemnmax = i
               endif

               write (inesc, '(a,i4,a)') 
     1            '---- end info for isvd = ', i,' ----'
            end do
c...end of svdscan loop.....................................

            write (inesc, 10)
     1          '---- %err_U Minimum (Bmod at local) Table ----'
            write (inesc, 20) 'Abs_U ', iemnabs, ' : ', eabsmin
            write (inesc, 20) 'Rms_U ', iemnrms, ' : ', ermsmin
            write (inesc, 20) 'Max_U ', iemnmax, ' : ', emaxmin
            write (inesc, 10) 
     1          '---- %err_W Minimum (Bmod at global) Table: ----'
            write (inesc, 20) 'Abs_W ', iwemnabs, ' : ', weabsmin
            write (inesc, 20) 'Rms_W ', iwemnrms, ' : ', wermsmin
            write (inesc, 20) 'Max_W ', iwemnmax, ' : ', wemaxmin
            write (inesc, 10) '---- Current density Minimum Table: ----'
            write (inesc, 20) 'Jave ', icdavemin, ' : ', cdavemin
            write (inesc, 20) 'max ', icdmaxmin, ' : ', cdmaxmin
            write (inesc, 20) 'JRad MAX at ', icurvamax, ' : ', curvamax
            write (inesc, 20) 'Complexity Min ',icmplxmin,' : ',cmplxmin

         endif                                   !End of svdscan if

c.....   Special cases to write various things:
         if( w_svd > 1 .or. w_svd == - 2 ) then
c        write all svd weights to output
            m = abs(mkeep)
            write (inesc, '(2(a,i4),a)') '---- svd weights: ', m,
     1        ' out of ', nsvdw,' ----'
            write (inesc, *) (svdw(i),i=1,m)
            write (inesc, 10) '---- end svd weights ----'
         endif

         write (inesc, 10) '---- end SVD Scan ----'

c.....   Calc Xerr with full svd scan from given solutions svdv
c        write all xerrors to output
         if (mdspw<0 .and. mstrt>=0 .and. mstep/=0) call
     1      min_xerr_phimn (0, epsilon, nu1, nv1, nuvh1, nbdim, ndim
     2      , bfn, ben, dsur1, svdv, svdw, nsvdw)

c....................................................
c        Always fill "result" array for postprocessing with MKEEPA wgts
         result(:ndim) = svdv(1,mkeepa:ndim-1+mkeepa)

c        Throw away svd arrays now
  123    continue
c..................................................................
c        All branches end up here with "result"
c        which is then used by the standard nescoil same way as always

c*********** End of Berr and Xerr modifications PMV ***********

c        Postprocessing:
c        All braches (original, Berr, Xerr) arrive here with "result"
c        By this time the phi(m,n) is known in array result(ndim),
c        all that needs to be done is unpack it into phi(m,n) etc

c        first zero all fields b
         bpx(:nuvh1) = zero
         bpy(:nuvh1) = zero
         bpz(:nuvh1) = zero
c...................................
c        Calculate b fields at all points from phi_m_n, i.e., pot(m,n)
         do l = 1, ndim
            bpx(:nuvh1) = bpx(:nuvh1) + bfx(:nuvh1,l)*result(l)
            bpy(:nuvh1) = bpy(:nuvh1) + bfy(:nuvh1,l)*result(l)
            bpz(:nuvh1) = bpz(:nuvh1) + bfz(:nuvh1,l)*result(l)
         end do
c..........................
c        First set pot(m,n) for m=0, n=1,nf. Array index l is just n here
         l = 0
         m = 0
         pot(0,1:nf) = -result(1:nf)
         pot(0,-1:(-nf):(-1)) = result(1:nf)
         l = nf + l
c.............................
c        Next set pot(m,n) for m=1,mf and n=1,nf. Array index l keeps incrementing
         do m = 1, mf
            pot(m,(-nf):nf) = -result(l+1:l+nf*2+1)
            l = l + nf*2 + 1
         end do
C.................................
C    Write pot.coeffs to output file if this was NOT an svdscan
         if (mstrta < 1) then
            potpotm = zero       !for moments of pot(m,n) for complexity
            potpotmm = zero
            write (inesc, 10) '---- Phi(m,n) for least squares ---'
            do m = 0, mf
               epsm = one
               if (m == 0) epsm = 0.5_dp
               do n = -nf, nf
                  write (inesc,'(i3,2x,i3,2x,g25.16)') m,n,epsm*pot(m,n)
                  potsq = pot(m,n) * pot(m,n) * m
                  potpotm = potpotm + potsq
                  potpotmm = potpotmm + potsq*m
               end do
            end do
            write (inesc, 10) '---- end Phi(m,n) for least squares ----'
            complexity = -1
            if (potpotm .ne. zero) complexity = potpotmm/potpotm
            write (inesc, '(a,1pe16.8)')'Complexity = ', complexity
         endif
         
      if (iloop .le. 0)
     1  deallocate (a, work, svdv, svdw, ab, result, bfx, bfy, bfz,
     1   dumv, bpx, bpy, bpz, daccu, bmag, bfn, bex, bey, bez, ben)

 10   format (a)
 20   format (a,i4,a,1pe16.8) 

      end subroutine solver_nescoil


!-------------------------------------------------------------
      subroutine calc_percent_errs(err, dsur1, fsarea, nu1, nuvh1,
     1   errmax, errabs, errrms)
!................................................................
c    Purpose: Calculate % errors err array on dsur1
c................................................................
C   M o d u l e s
c................................................................
      use Vmeshes, ONLY: rprec
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer nu1, nuvh1
      real(rprec) fsarea, errmax, errabs, errrms
      real(rprec), dimension(nuvh1) :: err, dsur1
c................................................................
c    1. Calculate the max err
      errmax = 100.*maxval(err)

c    2. Calculate the average of abs err: (Int = sum)
c       = Int[ds*|ERR|]/Int[ds]
      err(:nuvh1) = err(:nuvh1) * dsur1(:nuvh1)
c       and its surface average
      errabs = 100.*(sum(err(:nu1)) + 2.*sum(err(nu1+1:nuvh1-nu1))
     1   + sum(err(nuvh1-nu1+1:nuvh1)))/fsarea

c    3. Calculate the root mean square (RMS) err: (Int = sum)
c       = sqrt{ Int[ds*(BFN*PHI-BEN)^2]/Int[ds] }
      err(:nuvh1)=err(:nuvh1)*err(:nuvh1)/dsur1(:nuvh1)
      errrms = 100.*sqrt((sum(err(:nu1)) + 2.*sum(err(nu1+1:nuvh1-nu1))
     1       + sum(err(nuvh1-nu1+1:nuvh1)))/fsarea)

      end subroutine calc_percent_errs


!-------------------------------------------------------------
      subroutine ab_matrix(ab, bfn, dumv, ndim, mdspwa, curwt)
!................................................................
c     Purpose:
c     Calculate ab matrix:
c     use MDSPWA  to control power of dsur,
c         CURWT   to control amount of J minimization
c................................................................
C   M o d u l e s
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal3
      USE Vprecal6, ONLY: comu, simu
      USE Vprecal7, ONLY: conv, sinv
      USE Vsurfac14, ONLY: ben
      USE Vsurfac13, ONLY: dsur1
      USE Vsurface2, ONLY: xu, yu, xv, ru, rv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface7, ONLY: snx, sny, snz
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ndim, mdspwa
      real(rprec) :: ab(ndim, ndim+1), bfn(nuvh1, ndim),
     1   dumv(nuvh1), curwt
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer::mexp,l,lp,i,j,n,m,ndimp,jv,ju,m1,n1,m2,n2
      integer :: icm(ndim), icn(ndim)
      real(rprec) :: curwt2, cupnp, eux, euy, euz, evx, evy, evz,
     1   cmunv1, f1x, f1y, f1z, ex, ey, ez, cmunv2, f2x, f2y, f2z
C-----------------------------------------------

      ab(:ndim,:ndim+1) = zero

c     Calculate the dsur factor for integral over plasma surface
c     Also include the muliplicity factor(i) (1 or 2) in dumv(i)
      if (mdspwa >= 2) then
         mexp = mdspwa - 2                       !exponent in dsur**mexp
         dumv(:nuvh1) = factor(:nuvh1)*dsur1(:nuvh1)**mexp
      else
         dumv(:nuvh1) = one             !here just sum|berr|^2 minimized
      endif

c...1.First calculate the G and H parts for standard nescoil
c     Calculate ab matrix G part= Tr(G)*G for LSQ case
c     by integrating over nuvh1 on plasma surface
c     This makes the matrix square (ndim,ndim)
      do l = 1, ndim
         do lp = l, ndim
          ab(l,lp) = ab(l,lp) +
     1     sum( bfn(:nuvh1,l) *bfn(:nuvh1,lp) *dumv(:nuvh1) )
         end do
      end do

c     Calculate ab matrix h part (last column(ab)=Tr(G)*bnorm)
c     by integrating over nuvh1
      do l = 1, ndim
         ab(l,ndim+1) = ab(l,ndim+1) -
     1     sum( bfn(:nuvh1,l) *ben(:nuvh1) *dumv(:nuvh1) )
      end do

c...2.Add jminimization parts Fi and Ei to ab matrix
      if (curwt > 0) then
         write (inesc, '(a,1pe16.8)') 
     1     'Jmin targeted Using SVD Weight = ',curwt
         curwt2 = 2*curwt

c        Calculate the dsur factor for integral over coil surface
c        Also include the factor(i) (1 or 2) and curwt here into dumv(i)
c        Also fold in the 1/|eu x ev| = 1/dsur into dumv(i)
         if (mdspwa >= 2) then
            mexp = (-1) + .5_dp*(mdspwa - 2)    !divide by 2 to avoid sqrt

            do i = 1, nu
               dumv(i) = curwt*
     1            (snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i))**mexp
               j = i + nuv/2
               dumv(j) = curwt*
     1            (snx(j)*snx(j)+sny(j)*sny(j)+snz(j)*snz(j))**mexp
            end do

            do i = nu+1, nuv/2
               dumv(i)= CURWT2 *
     1         (snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i))**mexp
            enddo

         else
            dumv(:nuvh1) = one       !use just sum||^2 , not integral ds
         endif

c        Fill index arrays to get m,n from imn
c        ndim -> im,in map in fourier space
         i = 0
         do n = 1, nf                        !first m=0 cases (nf total)
            i = i + 1
            icm(i) = 0
            icn(i) = n
         end do
         do m = 1, mf             !next m>0 cases  ( mf*(2*nf+1) total )
            do n = -nf, nf
               i = i + 1
               icm(i) = m
               icn(i) = n
            end do
         end do

c        Add F and E to the G and H parts: Do not store F or E matrices
         cupnp = cup/np
         ndimp = ndim + 1
         i = 0
         do jv = 1, 1 + nv/2                     !half v + centerline
            do ju = 1, nu          !all of u. Total = nu*(1+nv/2) = nuvh
               i = i + 1
               eux = xu(i)                       !Tanget vector eu
               euy = yu(i)
               euz = zu(i)
               evx = xv(i)                       !Tanget vector ev
               evy = yv(i)
               evz = zv(i)

               do l = 1, ndim               !go over all ndim m,n values
                  m1 = icm(l)
                  n1 = iabs(icn(l))
                  cmunv1 = pi2*(comu(ju,m1)*conv(jv,n1)
     1                   -      simu(ju,m1)*sinv(jv,n1))
                  f1x = n1*eux - m1*evx
                  f1y = n1*euy - m1*evy
                  f1z = n1*euz - m1*evz

                  ex = cut*evx - cupnp*eux
                  ey = cut*evy - cupnp*euy
                  ez = cut*evz - cupnp*euz

c                 Add the F * E parts (x,y,z) to ab matrix:
                  ab(l,ndimp) = ab(l,ndimp) + dumv(i)*cmunv1*
     1               (f1x*ex + f1y*ey + f1z*ez)

                  do lp = l, ndim         !go over all ndim m',n' values
                     m2 = icm(lp)
                     n2 = iabs(icn(lp))
                     cmunv2 = pi2*
     1                 (comu(ju,m2)*conv(jv,n2)-simu(ju,m2)*sinv(jv,n2))
                     f2x = n2*eux - m2*evx
                     f2y = n2*euy - m2*evy
                     f2z = n2*euz - m2*evz

c                    Add the F * F parts (x,y,z) to ab matrix:
                     ab(l,lp) = ab(l,lp) + dumv(i)*cmunv1*cmunv2*
     1                  (f1x*f2x + f1y*f2y + f1z*f2z)
                  end do
               end do
            end do
         end do

      endif                       !End addition of jminimization section

c...3.Finally, Symmetrize ab matrix G part
c     no need to do center spine, so ndim-1
      do lp = 1, ndim - 1
         ab(1+lp:ndim,lp) = ab(lp,1+lp:ndim)
      end do

      end subroutine ab_matrix

! ----------------------------------------------------------------------
      subroutine matrix(a, ip)
! ----------------------------------------------------------------------
c                                                             01.01.89
c     purpose:
c
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vprecal1
      USE Vprecal4
      USE Vsurface1
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vsurface9
      USE Vsurfac14, ONLY: bex, bey, bez
      USE LoopCtrl
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ip
      real(rprec), dimension(3*nuv) :: a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: in, in1, in2, k, i
      real(rprec) :: sux, suy, suz
      real(rprec), allocatable, dimension(:) ::
     1     dax, day, daz, dbx, dby, dbz, dx, dy, dz, sq2, sq3, pn
c ----------------------------------------------------------------------

      i=0
      allocate (dax(nuv), day(nuv), daz(nuv), dbx(nuv),
     1          dby(nuv), dbz(nuv), dx(nuv), dy(nuv), dz(nuv),
     2          pn(nuv), sq2(nuv), sq3(nuv), stat = i)
      if (i .ne. 0) stop 'Allocation error in nescoil matrix'
      dax(:) = zero
      day(:) = zero
      daz(:) = zero
      dbx(:) = zero
      dby(:) = zero
      dbz(:) = zero
      dx(:) = zero
      dy(:) = zero
      dz(:) = zero
      pn(:) = zero
      sq2(:) = zero
      sq3(:) = zero

      in  = 0
      in1 = in + nuv
      in2 = in1 + nuv
      dz = z  - z1(ip)
      do k = 0, np - 1
         dx = x  - x1(ip)*cok(k) + y1(ip)*sik(k)
         dy = y  - y1(ip)*cok(k) - x1(ip)*sik(k)
         sq2 = one/(dx*dx + dy*dy + dz*dz)
         sq3 = sq2*sqrt(sq2)
         pn = three*sq2*(dx*snx + dy*sny + dz*snz )
         dax  = (snx  - pn*dx)*sq3
         day  = (sny  - pn*dy)*sq3
         daz  = (snz  - pn*dz)*sq3

         a(in +1:in1) = a(in +1:in1) + cok(k)*dax  + sik(k)*day
         a(in1+1:in2) = a(in1+1:in2) - sik(k)*dax  + cok(k)*day
         a(in2+1:in2+nuv) = a(in2+1:in2+nuv) + daz

         dbx  = (ycur*dz - zcur*dy)*sq3
         sux = sum(dbx(:nuv))*fnuv

         dby  = (zcur*dx - xcur*dz)*sq3
         suy = sum(dby(:nuv))*fnuv

         dbz  = (xcur*dy - ycur*dx)*sq3
         suz = sum(dbz(:nuv))*fnuv

         bex(ip) = bex(ip) + cok(k)*sux + sik(k)*suy
         bey(ip) = bey(ip) - sik(k)*sux + cok(k)*suy
         bez(ip) = bez(ip) + suz
      end do

c........Deallocate matrix arrays
      deallocate (dax, day, daz, dbx, dby, dbz, pn,
     1            sq2, sq3, dx, dy, dz)

      end subroutine matrix

! ----------------------------------------------------------------------
      subroutine fourier(ndim)
! ----------------------------------------------------------------------
c                                                             01.01.89
c     purpose:
c
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      USE Vfourier2
      USE Vsolver1, ONLY: a, work
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ndim
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nel, i, ni, nim, k, inc, jump, nuf, lot
C-----------------------------------------------
c
c   1. fourier transformation with respect to u
c
      nel = nv*ndim
      do i = nel, 1, -1
         ni = i*(nu + 2)
         nim = i*nu + 1
         a(ni)   = zero
         a(ni-1) = zero
         do k = 1, nu
            a(ni-1-k) = a(nim-k)
         end do
      end do
c
      inc = 1
      jump = nu + 2
      call fft991 (a, work, trigs, ifax, inc, jump, nu, nel, -1)
c
      nuf = 2*(mf + 1)
      do i = 1, nel
         ni = (i - 1)*(nu + 2)
         nim = (i - 1)*nuf
         do k = 1, nuf
            a(nim+k) = a(ni+k)
         end do
      end do
c
c   2.   fourier transformation with respect to v
c
      inc = mf + 1
      jump = inc*nv
      lot = ndim
      do i = 1, mf + 1
         call cfft99(a(2*i-1),work,trigsv,ifaxv,inc,jump,nv,lot,-1)
      end do

      end subroutine fourier

! ----------------------------------------------------------------------
      subroutine accuracy
! ----------------------------------------------------------------------
c                                                          05/01/89
c     purpose:
c
C-----------------------------------------------
      use Vmeshes
      use NumParams
      use OutCtrl, ONLY: w_bnuv
      USE LoopCtrl
      USE Vvacuum3
      USE Vsurface9
      USE Vsurfac13, ONLY: snx1, sny1, snz1, dsur1, fla
      USE Vbnorm
      USE Vcen
      USE Voptim1, ONLY: averro, ermax, ermin
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, istat
      integer, dimension(1) :: imax, imin
      real(rprec), dimension(:), allocatable :: derro, dvari, dme,
     1  spx, spy, spz
      real(rprec) :: bdf, berr, bb,
     1      bmod_ave, bmod_rms, varian
c ----------------------------------------------------------------------

      istat=0
      allocate (derro(nuvh1), dvari(nuvh1), dme(nuvh1),
     1          spx(nuv), spy(nuv), spz(nuv), stat=istat)
      if (istat .ne. 0) stop 'allocation error in accuracy'
      derro(:) = zero
      dvari(:) = zero
      dme(:) = zero
      spx(:) = zero
      spy(:) = zero
      spz(:) = zero

      call biowipre

      if(w_bnuv>0 .or. w_bnuv==-1)
     1   write(inesc,'(a)') '---- dB normal, Babs, dsur_plas table ----'
      bmod_ave = zero
      bmod_rms = zero
      do i=1,nuvh1
         xp   = x1(i)
         yp   = y1(i)
         zp   = z1(i)
         call besurfcur(spx, spy, spz)
         berr = snx1(i)*bx+sny1(i)*by+snz1(i)*bz + bn_ext(i)
         bdf  = abs( berr )
         dme  (i) = bdf*bmod1   !bmod1=1/bmod was calculated in besurfcur
         derro(i) = dme(i)*dsur1(i)
         dvari(i) = dme(i)*derro(i)
         bb = bmod * dsur1(i)
         bmod_ave = bmod_ave + bb
         bmod_rms = bmod_rms + bmod * bb
         if(w_bnuv>0 .or. w_bnuv==-1) write(inesc,*) berr,bmod,dsur1(i)
      enddo
      if(w_bnuv>0 .or. w_bnuv==-1) write(inesc, '(a)')
     1     '---- end dB normal, Babs, dsur_plas table ----'

      imax = maxloc(dme(1:nuvh1))
      imin = minloc(dme(1:nuvh1))
      ermax = dme(imax(1))
      ermin = dme(imin(1))
      averro = (sum(derro(1:nu1)) + 2*sum(derro(nu1+1:nuvh1-nu1))
     1       + sum(derro(nuvh1-nu1+1:nuvh1)))/fla
      varian = (sum(dvari(1:nu1)) + 2*sum(dvari(nu1+1:nuvh1-nu1))
     1       +  sum(dvari(nuvh1-nu1+1:nuvh1)))/fla
      varian = varian-averro*averro

      write(inesc, 20)
     1   'Berr ave, max, var = ', averro, ermax, sqrt(varian)
      write(inesc, 20)'Bmod ave, rms = ',bmod_ave/fla,sqrt(bmod_rms/fla)

 20   format (a,1p3e16.8)

      deallocate (derro, dvari, dme, spx, spy, spz)

      end subroutine accuracy

! ----------------------------------------------------------------------
      subroutine biowipre
! ----------------------------------------------------------------------
c                                                             05/01/89
c     purpose:
c
c
c ----------------------------------------------------------------------
      use Vmeshes
      USE Vsurfac13
      USE Vculine1
      USE Vbiopre, ONLY: vx, vy, vz, dx, dy, dz
      implicit none
C-----------------------------------------------
      dx(:nw-1) = (xw(2:nw)-xw(:nw-1))*curre(:nw-1)
      dy(:nw-1) = (yw(2:nw)-yw(:nw-1))*curre(:nw-1)
      dz(:nw-1) = (zw(2:nw)-zw(:nw-1))*curre(:nw-1)
      vx(:nw-1) = yw(:nw-1)*dz(:nw-1) - zw(:nw-1)*dy(:nw-1)
      vy(:nw-1) = zw(:nw-1)*dx(:nw-1) - xw(:nw-1)*dz(:nw-1)
      vz(:nw-1) = xw(:nw-1)*dy(:nw-1) - yw(:nw-1)*dx(:nw-1)
c ----------------------------------------------------------------------

      end subroutine biowipre

! ----------------------------------------------------------------------
      subroutine besurfcur(spx, spy, spz)
! ----------------------------------------------------------------------
c                                                             01/01/89
c     purpose:
c
c
c ----------------------------------------------------------------------
      use Vmeshes
      use NumParams
      USE Vprecal1
      USE Vprecal4
      USE Vsurface1, ONLY: x, y, z
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vsurface8, ONLY: xcur, ycur, zcur
      USE Vdiagno3, ONLY: pote
      USE Vvacuum3
      USE Vcen, ONLY: xp, yp, zp, bmod, bmod1, bx, by, bz
      USE LoopCtrl
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nuv) :: spx, spy, spz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ico = 0, i, k
      real(rprec), allocatable, dimension(:) ::
     1     dbx, dby, dbz, dx, dy, dz, sq2, sq3, dn
      real(rprec) :: bxk, byk, bxadd, byadd, bzadd
c ----------------------------------------------------------------------
c
      i=0
      if (.not.allocated(dbx))
     1 allocate (dbx(nuv), dby(nuv), dbz(nuv), dn(nuv), dx(nuv),
     1   dy(nuv), dz(nuv), sq2(nuv), sq3(nuv), stat=i)
      if (i .ne. 0) stop 'allocation error in NESCOIL besurfcur'
      dbx(:) = zero
      dby(:) = zero
      dbz(:) = zero
      dn(:) = zero
      dx(:) = zero
      dy(:) = zero
      dz(:) = zero
      sq2(:) = zero
      sq3(:) = zero

      if(ico .eq. 0) then
         ico    = 1
         spx = pote*snx
         spy = pote*sny
         spz = pote*snz
      endif

      bx     = zero
      by     = zero
      bz     = zero

      do k = 0,np-1
         dx  = xp*cok(k) + yp*sik(k) - x
         dy  = yp*cok(k) - xp*sik(k) - y
         dz  = zp  - z
         sq2 = one/(dx*dx + dy*dy + dz*dz)
         sq3 = sq2*sqrt(sq2)
         dn  = three*sq2*(spx*dx + spy*dy + spz*dz)
         dbx = (spx - dn*dx + ycur*dz - zcur*dy)*sq3
         dby = (spy - dn*dy + zcur*dx - xcur*dz)*sq3
         dbz = (spz - dn*dz + xcur*dy - ycur*dx)*sq3
c
         bxk    = sum(dbx(1:nuv))*fnuv
         byk    = sum(dby(1:nuv))*fnuv
         bx     = bx + (bxk*cok(k) - byk*sik(k))
         by     = by + (byk*cok(k) + bxk*sik(k))
         bz     = bz + sum(dbz(1:nuv))*fnuv
      end do
C
      bx = -bx
      by = -by
      bz = -bz
C
c     TURNED THIS OFF, ADD BEN(i) IN CALLING ROUTINE
      if (ibex .ne. 0) then
         call bexter (xp, yp, zp, bxadd, byadd, bzadd)
         bx = bx + bxadd
         by = by + byadd
         bz = bz + bzadd
      endif
c
      bmod   = sqrt(bx*bx+by*by+bz*bz)
      bmod1  = one/bmod

      if (iloop .le. 0)
     1 deallocate (dbx, dby, dbz, dx, dy, dz, sq2, sq3, dn)

      end subroutine besurfcur

! ----------------------------------------------------------------------
      subroutine bnfld
!.....................................................
      USE Vmeshes
      USE NumParams
      use OutCtrl, ONLY: w_bnuv
      USE Vbnorm
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal6
      USE Vprecal7
      USE Vsurfac13
      USE LoopCtrl
      use safe_open_mod
      implicit none
c.....................................................
C   L o c a l   V a r i a b l e s
c.....................................................
      integer, parameter :: ibnorm = 19
      integer :: ms, ns, iunit = ibnorm
      integer :: i, mm, nn, m, n, kv, ku
      real(rprec), allocatable, dimension(:,:) :: bnorm_mn
      real(rprec) :: bf, cm, cn, sifp, sifm, fact, rad2
      character*200 :: buffer
c.....................................................
      i=0
      if (.not.allocated(bn_ext))
     1  allocate (bn_ext(nuvh1), bnorm_mn(0:md,-nd:nd), stat = i)
      if (i .ne. 0) stop 'allocation error in NESCOIL bnfld'
      bn_ext(:nuvh1) = zero
      bnorm_mn(0:md,-nd:nd) = zero
c.....................................................
c     Open bnorm file and read bnorm fourier coeffs
      if (extension(1:1) .eq. ' ') then
         buffer = 'bnorm'
      else
         buffer = 'bnorm.' // extension
      end if      

      call safe_open(iunit, i, trim(buffer), 'old', 'formatted')
      if (i .ne. 0 ) then
         print *, 'No bnorm field file found. Assuming bnorm_ext = 0'
         return
      end if     
c
      ms = 0
      ns = 0
c
      do
        read(iunit,*,iostat = i) mm, nn, bf
        if (i .ne. 0) exit
        ms = max(mm, ms)
        ns = max(abs(nn), ns)
        if ( (ms > md) .or. (ns > nd) ) then
           write (inesc, '(a)') 
     1        'Not enough room to read bnorm_mn from file'
           write (inesc, '(2(a,i4))') 
     1        'Increase md, nd to >= ', ms, ',',ns
           write (inesc, '(a)') 'in nesinput file'
           stop 'Rerun nescoil'
        end if
        bnorm_mn(mm,nn) = bf
      end do
c
      write (inesc, '(2(a,i4))') 
     1   'Got bnorm_mn mmax = ', ms, ', nmax = ', ns

      close(iunit)
c.....................................................
c     Turn bnorm_mn into bn_ext(nuvh1) on plasma surface
      bnorm_mn(:ms,0) = .5*bnorm_mn(:ms,0)
      do m = 0, ms
         do n = 0, ns
c           cm = m*pi2
c           cn = n*pi2
            i = 0
            do  kv = 1 , nv1/2 + 1
               do  ku = 1 , nu1
                  i  = i + 1
                  sifp = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
                  sifm = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
                  bn_ext(i) = bn_ext(i) + bnorm_mn(m,n)*sifp
     1                                  + bnorm_mn(m,-n)*sifm
               enddo
            enddo
         end do
      end do
c
      if (iloop .le. 0) deallocate (bnorm_mn)

      fact = -2*pi2
      bn_ext(:nuvh1) = fact*bn_ext(:nuvh1)

      write (inesc, '(a)') 'Calculated bn_ext(u,v) on plasma surface'

c.....Write external B field to output if asked to do so
      if( w_bnuv > 1 .or. w_bnuv == -2 ) then
         write (inesc,'(a,i4,a)') 
     1     '---- bn_ext(i) for i = 1,',nuvh1,' ----'
         write (inesc,*) (bn_ext(i), i = 1, nuvh1)
         write (inesc, '(a)') '---- end bn_ext(i) ----'
      endif

      end subroutine bnfld

! -------------------------------------------------------------
      subroutine surfcur_diag
!................................................................
c     Purpose:
c................................................................

      USE Vmeshes
      USE NumParams
      use OutCtrl, ONLY: w_jsurf
      USE LoopCtrl
      USE Vvacuum3
      USE Vprecal1
      USE Vprecal2
      USE Vprecal3
      USE Vprecal6
      USE Vprecal7
      USE Vsurface1
      USE Vsurface2, ONLY: xu, yu, xv
      USE Vsurface3, ONLY: yv, zu, zv
      USE Vsurface4, ONLY: xuu, yuu, zuu
      USE Vsurface5, ONLY: xuv, yuv, zuv
      USE Vsurface6, ONLY: xvv, yvv, zvv
      USE Vsurface7, ONLY: snx, sny, snz
      USE Vdiagno2
      USE Vdiagno3
      use Vsurfcurdiag
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, n, ku, kv
      real(rprec), dimension(:), allocatable ::  wx, wy, wz,
     1     curvat, cude, potu, potv, potuu, potuv, potvv, vx, vy, vz
      real(rprec) :: cw,cwm,sw,swm,cm,cn,v2,w2,vw
c ----------------------------------------------------------------------
      i=0
      if (.not.allocated(wx))
     1 allocate (wx(nuv), wy(nuv), wz(nuv), curvat(nuv), cude(nuv),
     1   potuu(nuv), potuv(nuv), potvv(nuv), vx(nuv), vy(nuv),
     2   vz(nuv), potu(nuv), potv(nuv), pote(nuv), stat = i)
      if (i .ne. 0) return
c
      vx= zero
      vy= zero
      vz= zero
      curvat= zero
      cude= zero
      wx= zero
      wy= zero
      wz= zero

      pote = zero
      potu = cut
      potv = cup
      potuu = zero
      potuv = zero
      potvv = zero

      do 30 n= 0,nf
         cn      = pi2*n
      do 30 m= 0,mf
         cm      = pi2*m
         i       = 0
      do 30  kv = 1,nv
      do 30  ku = 1,nu
         i      = i+1
         cw     = comu(ku,m)*conv(kv,n)-simu(ku,m)*sinv(kv,n)
         cwm    = comu(ku,m)*conv(kv,n)+simu(ku,m)*sinv(kv,n)
         sw     = simu(ku,m)*conv(kv,n)+comu(ku,m)*sinv(kv,n)
         swm    = simu(ku,m)*conv(kv,n)-comu(ku,m)*sinv(kv,n)
         pote(i)  = pote(i)+ eps(m)*eps(n)*(pot(m,n)*sw + pot(m,-n)*swm)
         potu(i)  = potu(i) + cm   *eps(n)*(pot(m,n)*cw + pot(m,-n)*cwm)
         potv(i)  = potv(i) + cn   *eps(m)*(pot(m,n)*cw - pot(m,-n)*cwm)
         potuu(i) = potuu(i)- cm*cm*eps(n)*(pot(m,n)*sw + pot(m,-n)*swm)
         potuv(i) = potuv(i)- cm      *cn* (pot(m,n)*sw - pot(m,-n)*swm)
         potvv(i) = potvv(i)- cn*cn*eps(m)*(pot(m,n)*sw + pot(m,-n)*swm)
   30 continue

c     curvature  of the surface current lines
c
      do i=1,nuv
         vx(i)  = xu(i)*potv(i)-xv(i)*potu(i)
         vy(i)  = yu(i)*potv(i)-yv(i)*potu(i)
         vz(i)  = zu(i)*potv(i)-zv(i)*potu(i)
         wx(i)  = xu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + xv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  xuv(i)*potv(i)-  xvv(i)*potu(i))
     $          - potv(i)*(  xuv(i)*potu(i)-  xuu(i)*potv(i))
c
         wy(i)  = yu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + yv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  yuv(i)*potv(i)-  yvv(i)*potu(i))
     $          - potv(i)*(  yuv(i)*potu(i)-  yuu(i)*potv(i))
c
         wz(i)  = zu(i)  *(potuv(i)*potv(i)-potvv(i)*potu(i))
     $          + zv(i)  *(potuv(i)*potu(i)-potuu(i)*potv(i))
     $          - potu(i)*(  zuv(i)*potv(i)-  zvv(i)*potu(i))
     $          - potv(i)*(  zuv(i)*potu(i)-  zuu(i)*potv(i))
      end do

      do i=1,nuvh
         v2     = vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         w2     = wx(i)*wx(i)+wy(i)*wy(i)+wz(i)*wz(i)
         vw     = vx(i)*wx(i)+vy(i)*wy(i)+vz(i)*wz(i)
         curvat(i) = sqrt(w2-vw*vw/v2)/v2
         cude(i) = sqrt(v2/(snx(i)*snx(i)+sny(i)*sny(i)+snz(i)*snz(i)))
      end do

      curvra = maxval(curvat(1:nuvh))
      if (curvra .ne. zero) curvra = 1./curvra
      curvri = minval(curvat(1:nuvh))
      if (curvri .ne. zero ) curvri = 1./curvri
      cudema = maxval(cude(1:nuvh))
      cudemi = minval(cude(1:nuvh))
      avcude = (sum(cude(1:nu)) + 2.*sum(cude(1+nu:nuvh-nu))
     1       +  sum(cude(nuvh-nu+1:nuvh)))/float(nuv)

c.....Always Write this surface current info to output
      write (inesc, 20) 'J Surf max, min, ave = ',cudema, cudemi, avcude
      write (inesc, 20) 'Curv_R of J min, max = ',curvra, curvri


c.....Write extra surface current info to output if asked to do so
      if(w_jsurf>0 .or. w_jsurf==-1) then
         write (inesc, '(a,i4,a)') '---- potential phi(i) for i = 1,',
     1        nuv,' ----'
         write (inesc,*) (pote(i), i = 1, nuv)
         write (inesc, '(a)') '---- end  potential phi(i) ----'
      endif

      if(w_jsurf>1 .or. w_jsurf==-2) then
         write (inesc, '(a,i4,a)') 
     1      '---- current jsurf(i) for i = 1,',nuvh,' ----'
         write (inesc,*) (cude(i), i = 1, nuvh)
         write (inesc, '(a)') '---- end  current jsurf(i) ----'
      endif

      if (iloop .le. 0) deallocate (wx, wy, wz, curvat, cude,
     1   potuu, potuv, potvv, vx, vy, vz, potu, potv)

 20   format(a,1p3e16.8)

      end subroutine surfcur_diag

! ------------------------------------------------------------
      subroutine min_xerr_phimn(method, eps, nu, nv, nuvh,
     1   nbdim, ndim, bfn, ben, dsur, svdv,svdw,nsvdw )
!................................................................
c     Purpose:
c     Calculate phi(m,n) which minimize xerr (rather than berr)
c................................................................

      use kind_spec
      USE NumParams
      use OutCtrl, ONLY: w_xerr
      USE LoopCtrl
      use Vvacuum2, ONLY: cl1, iota_edge, phip_edge, ms1, ns1
      USE Vprecal1, ONLY: np
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      integer method, numsvd, nu, nv, nuvh, nbdim, ndim, nsvdw
      real(rprec) eps
      real(rprec), dimension(nuvh,nbdim) :: bfn
      real(rprec), dimension(nuvh) :: ben, dsur
      real(rprec), dimension(ndim) :: svdw
      real(rprec), dimension(ndim,ndim) :: svdv
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      integer :: mnl, mf, nf, nvh, mnf,
     1   ia, istat, i, m, n, k, kv, nwt, ku, j, nn, ns
      real(rprec) :: epsq, du, hval, u, res,
     1   rmsxerr, absxerr
      real(rprec), allocatable, dimension(:,:,:) :: smunv, cmunv
      real(rprec), allocatable, dimension(:,:) ::
     1    cu, cv, su, sv, xg
      real(rprec), allocatable, dimension(:) ::
     1    lmn, luv, resmn, xh, gval
c................................................................
c
c  Inside-nescoil processor to calculate phi(m,n) current potentials
c  by minimizing displacements by changing Nescoil Green functions to
c  those needed for calculating displacements from phimn.
c
c  Given ben, bfn in (theta,phi), and
c        fourier coeffs lmn of lambda and g from vmec,
c  Calculate the modified Green functions, and
c  solve in fourier space to get phi(m,n)
c...
c  t=theta, p=phi, u = t + lambda(t,p), v = p
c  Straight-line system is (u,v):
c  u   = theta + lambda( theta, phi ),  v = phi
c  t   = theta,  p = phi
c...
c  Inputs:
c  method  - 0:use given svdv to calc xerr, else solve for svdv
c  np  - number of field periods
c  eps  - Resonance broadening parameter
c
c  nu   - number of points used in u (poloidal) grid
c  nv   - number of points used in v (toroidal) grid in 1 period
c  nuvh - number of points in half period + edge = nu*(1+nv/2)
c  nbdim-
c  ben  -
c  bfn  -
c  dsur -
c...
c  Outputs:
c  phimn - Fourier coeffs of phi which minimize xerror
c  numsvd- Number of svd weights to use
c...
c  Local Arrays and variables:
c
c  Following come from VMEC:
c  ml   - number of m values for the lambda coefficients (0:ml)
c  nl   - number of n values for the lambda coefficients (0:nl)
c  mnl  - length of lmn array = 1+nl+ml*(2*nl+1))
c  lmn  - Fourier coefficients of the stream function lambda from VMEC
c  phip_edge - radial derivative of toroidal flux at this surface
c  iota_edge - iota on this plasma surface
c
c  Following are used for fourier transform calculations:
c  smunv, cmunv, su, cu, sv, cv
c
c................................................................
c
c  Allocate local storage for arrays on nuvh spatial halfgrid
c  and for cos/sin factors on the largest fourier grid mf,nf
      mf = nu/2
      nf = nv/2
      mnf = 1 + nf + mf*(2*nf + 1)
      nvh = 1 + nv/2                             !half+centerline
      nuvh = nu*nvh
      istat=0
      if (.not.allocated(smunv)) allocate(
     1   smunv(nuvh,0:mf,-nf:nf), cmunv(nuvh,0:mf,-nf:nf),
     2   cu(nu,0:mf),cv(nvh,0:nf), su(nu,0:mf),sv(nvh,0:nf),
     3   luv(nuvh), resmn(mnf), gval(ndim),
     4   xg(mnf,ndim), xh(mnf), stat=istat)
      if (istat .ne. 0) stop 'allocation error in min__xerr_phimn_nes'
      luv(:) = zero
      resmn(:) = zero
      gval(:) = zero
      xh(:) = zero
      smunv(:nuvh,0:mf,-nf:nf) = zero
      cmunv(:nuvh,0:mf,-nf:nf) = zero
      cu(:nu,0:mf) = zero
      cv(:nvh,0:nf) = zero
      su(:nu,0:mf) = zero
      sv(:nvh,0:nf) = zero
      xg(:mnf,:ndim) = zero

c  Calculate all the resonance factors 1/(m*iota+n*np) now, so
c  they can be used to turn dmn into xmn = dmn * resmn
       if (eps .lt. zero) then   !no resonance modification
          resmn(:mnf) = one
       else                    !standard resonance modification
          epsq = eps*eps
          resmn(1) = 1
          k = 0                !m=n=0 is never used here
          do n = 1, nf
             k = k + 1
            resmn(k) = real(n*np,rprec)/((n*np)**2 + epsq)
          enddo
          do m = 1, mf
             do n = -nf, nf
                k = k + 1
                res = (m*iota_edge + n*np)
                resmn(k) = res / (res**2 + epsq)
             enddo
          enddo
       endif

c............................................................
c  Calculate all the sin/cos factors needed on fourier grid mf,nf
c  for the theta, phi VMEC system, NOT the straight-line one
c  Note: For VMEC, m*u-n*np*v convention means np should be negative
c  These can be used on smaller lmn grids too, so saves time
      call trigfact(nu,nv,nvh,nuvh,mf,nf,cmunv,smunv,cu,su,cv,sv)

c  Following are all in VMEC co-ordinates, not the straight-line ones
c  so do them now before sin/cos(mu) are chaged later

c  Turn incoming lambda coeffs lmn into luv(theta,phi) with sine
c  in the theta, phi VMEC system, not the straight-line one
      mnl = 1+ns1 + ms1*(2*ns1 + 1)
      istat=0
      if (.not.allocated(lmn)) allocate (lmn(mnl), stat=istat)
      if (istat .ne. 0) stop 'allocation error in min_xerr_phimn'
      lmn(:) = zero
      k = 0
      do n = 0, ns1              !First the m = 0 coeffs
         k = k + 1
         lmn(k) = -cl1(0,-n)     !Because bnorm wrote only -n for m=0
      enddo                      !and we need only +n for m=0
      do m = 1, ms1              !Then the m > 0 ones
         do n = -ns1, ns1
            k = k + 1
            lmn(k) = cl1(m,n)    !here pack the cl1 array
         enddo
      enddo

      call fmn_to_uv(nu, nv, nuvh, luv, ms1, ns1, mnl, lmn,
     1     mf, nf, smunv)

c  Multiply ben, bfn by dsur/phip
      ben(:nuvh) = ben(:nuvh)*dsur(:nuvh)/phip_edge
      do i = 1, ndim
         bfn(:nuvh,i) = bfn(:nuvh,i)*(dsur(:nuvh)/phip_edge)
      end do
c............................................................
c  Note: Everything is calculated at (theta,phi) grid in VMEC
c
c  Integrate over u,v by doing xmn = xmn + f[ku,kv], i.e,
c  first calculate dmn from b dot ds
c  dmn = sum_ku_kv[sin(m[u+luv(kuv)]+n*np*v) * (bds[kuv])] *(du*dv/phip)
c  and then calculate xmn from dmn
c  xmn = dmn *(m*iota+n*np)/( (m*iota+n*np)^2 + eps^2 )
c
c  Note: You need to recaculate only sin/cos(m[u+luv(kuv)]), since
c        sin/cos(n*np*v) were already calculated above
c  Note: It is better to put the ku,kv loops outside the m,n loops
c        because at each u,v, u+lam(u,v) can be calculated once and
c        from it cos/sin(m*[u+lam(u,v)]) can be calculated by making
c        only one call to cos/sin(u+lam(u,v)) and then iterating
c  Note: sin( m[theta+lambda] + n*np*phi) and
c        cos( m[theta+lambda] + n*np*phi) will be calculated here.
c        They will replace the smunv, cmunv calculated earlier,
c         i.e., sin/cos( m*theta + n*np*phi)
c
      xg(:mnf,:ndim) = zero            !zero all fourier coeffs of modified
      xh(:mnf) = zero                  !Green functions
      du = pi2/nu                      !Step size in theta
      su(:nu,0) = zero                 !for m=0
      cu(:nu,0) = one
      i = 0                            !i is real space index 1,nuvh
c......
      do kv = 1, nvh                   !See surface_plas.f, this is the same
c
c       Check for centerline (v=1/2 or kv=1+nv/2)
c       use value of v because it tkaes care of both odd and even nv
         if ((kv.eq.1) .or. (2*(kv-1).eq.nv)) then
            nwt = 1                     !centerline gets only one weight
         else
            nwt = 2                   !all others carry twice the weight
         endif
c.......
         do ku = 1, nu
            i = i + 1                 !i goes from 1 to nuvh=nu*(1+nv/2)
            hval = nwt*ben(i)         !do this outside m,n loops to save time
            gval(:ndim) = nwt*bfn(i,:ndim)
c
c         Re-calculate sin/cos(m*u) for this u since it differs from
c          theta by u = theta + lambda: straight coordinate
c         Note: cos/sin(n*v) need not change, use precalculated value
c         Note: This calculation cannot be done outside the
c               v loop since lambda is a function of both u and v
c
c         Find the straight-line u at this theta,phi = theta+lambda:
            u = du*(ku - 1) + luv(i)             !u = theta + lambda
c
c         Make only one call per u,v to cos/sin and then iterate
c         calculate cos/sin(m*[theta+lambda]):
            su(ku,1) = sin(u)                    !sin/cos(theta+lambda)
            cu(ku,1) = cos(u)
            do m = 2, mf                     !Now calculate sin/cos(m*u)
               n = m - 1
               su(ku,m) = su(ku,1)*cu(ku,n) + cu(ku,1)*su(ku,n)
               cu(ku,m) = cu(ku,1)*cu(ku,n) - su(ku,1)*su(ku,n)
            end do
c
c         Now assemble everything at this (u,v) for all (m,n), i.e.,
c         1. Calculate sin(m*[u+luv(u,v)]+n*np*v)
c            at this u,v for all m,n in (mf,nf) fourier space
c         2. Then add the contribution to xmn from this u,v point
c
c         Note: Since luv(u,v) depends on both u and v, this cannot
c               be done outside the ku,kv loops, must be done here
c         Skip the m=0,n=0 case since 1/(m*iota+n) is singular, so
c         we set xmn(m=0,n=0) = 0, i.e., f.s.average of xerr = 0
c         Note: Only cu,cv, su,sv WILL be changed here, but
c               smunv and cmunv will NOT be changed. So from now on,
c               cu,cv, su,sv will be for the straight-line system, and
c               smunv and cmunv will be for theta, phi VMEC system
c         Note: - sign in summation since bmn's are derivs of cos
            j = 1                          !j = m,n space index 1 to mnf
            do n = 1, nf                         !first the m = 0 cases
               j = j + 1           !skips j=1, starts from 2 for m=0,n=1
c           for m=0, m*theta = m*u = 0 in both systems
c           Add to xg and xh integrals from this theta,phi point
               xh(j) = xh(j) + hval*sv(kv,n)
               xg(j,:ndim) = xg(j,:ndim) - gval(:ndim)*sv(kv,n)
            end do
c
            do m = 1, mf                         !Next the m>0 cases
               do n = -nf, nf
                  j = j + 1
c           Save cmunv and smunv in straight-line system (u,v)
c           for use in calculating xtp and dxdt, dxdp later.
                  nn = abs(n)     !Need this since su/cu are over n=0:nf
                  ns = sign(1,n)
                  cmunv(i,m,n)=cu(ku,m)*cv(kv,nn)-su(ku,m)*sv(kv,nn)*ns
                  smunv(i,m,n)=su(ku,m)*cv(kv,nn)+cu(ku,m)*sv(kv,nn)*ns
c            Add to dmn integral from this theta,phi point
                  xh(j) = xh(j) + hval*smunv(i,m,n)
                  xg(j,:ndim) = xg(j,:ndim) - gval(:ndim)*smunv(i,m,n)
               end do
            end do
c
         end do                                  !ku loop
      end do                                     !kv loop
c.......
c      Multiply xg, xh by resonance and overall normalization factor,
c      Note: the 1/(2*pi)^2 is obsent since u,v go over 0 to 1
c                                             !2 since only m>0 are used
      resmn(:mnf) = 2*resmn(:mnf)/(nu*nv)
      xh(:mnf) = xh(:mnf)*resmn(:mnf)
      do i = 1, ndim
         xg(:mnf,i) = resmn(:mnf)*xg(:mnf,i)
      end do

c............................................................
c....  Calculate phi(m,n) by SVD solving xg * phimn = -hmn
c      Note: negative sign in xh has already been absorbed
c      Note: this will be a square matrix unless lamdamn
c            had more coeffs than mf=nu/2 and nf=nv/2, where
c            nu,nv are points on plasma surface 1 period
c            Since lamda is a function defined olny on these
c            points, it should NOT have more fourier coeffs.
c            mnf should be = nuvh, but ndim can be smaller
       write(inesc, '(a)') 'mf, nf, mnf, nuvh, ndim ='
       write(inesc,*)  mf, nf, mnf, nuvh, ndim
       if(method.ne.0) call svd_solve(mnf, ndim, mnf, ndim,
     1                    xg, xh, svdv, svdw, nsvdw)

c............................................................
c      Calculate RMS Xerror scan over all weights kept.
c      svdv(:,i) = phi_i = ith solution (with i svd weights kept,
c      xg * phi_i - xh = error due to ith solution
c      Report its RMS value when summed over all mnf (m,n) modes
c      Use the no longer needed resmn(:mnf) array for temp storage
       write (inesc, '(a)') '---- Xerror table ----'
       write (inesc, '(a)')
     1'   isvd     weight(i),     rms Xerr(i),      abs Xerr(i)'
       do i = nsvdw, 1, -1
          resmn(:mnf)=Matmul(xg(:mnf,:ndim),svdv(:ndim,i))-xh(:mnf)
          rmsxerr=sqrt( sum(resmn(:mnf)*resmn(:mnf)) /mnf )
          absxerr=sum( abs(resmn(:mnf)) ) / mnf
          write(inesc,"(i5,3g25.15)") i, svdw(i), rmsxerr, absxerr
       enddo
       write (inesc, '(a)') '---- end Xerror table ----'

c............................................................
c      Restore ben, bfn by multiplying by phip/dsur
c      for later Berr calculation using phimn solutions found here
       ben(:nuvh) = phip_edge * ben(:nuvh) / dsur(:nuvh)
         do i = 1, ndim
            bfn(:nuvh,i) = phip_edge * bfn(:nuvh,i) / dsur(:nuvh)
         enddo

c............................................................
c      Deallocate everything
      if (iloop .le. 0) deallocate( smunv,cmunv, cu,cv, su,sv,
     1  luv, lmn, resmn, gval, xg, xh, stat=istat)

      end subroutine min_xerr_phimn

!-------------------------------------------------------------
      subroutine bexter(xp, yp, zp, bx, by, bz)

!................................................................
c     Purpose:
c     Calculate B field from background coils at xp,yp,zp
c................................................................

      USE Vmeshes
      USE NumParams
      USE Vprecal1, ONLY: np
      USE Vvacuum3, ONLY: cup
      implicit none
c................................................................
C   D u m m y   A r g u m e n t s
c................................................................
      real(rprec) ::  xp, yp, zp, bx, by, bz
c................................................................
C   L o c a l   V a r i a b l e s
c................................................................
      real(rprec) :: bxy
c................................................................
c
c     FOR NOW, JUST SIMULATE EXTERNAL 1/R FIELD

      bxy =  cup / (xp*xp + yp*yp)
      bx  =  yp * bxy
      by  = -xp * bxy
      bz  =  zero

      end subroutine bexter

!---------------------------------------------------------------
