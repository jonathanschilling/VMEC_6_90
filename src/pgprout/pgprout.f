      program pgpplotout

C-----------------------------------------------   
c  Ed Lazarus Jan. 2000

c  A version of prout which uses PGPLOT
c  http://astro.caltech.edu/~tjp/pgplot/

c  the 3D package PGXTAL is also required
c  http://www.isis.rl.ac.uk/dataanalysis/dsplot/

c  Aside from a new module in avmodules.f, subroutine and file
c  names are the same as prout with "pgp" prepended.

c  The plots having to do with reconstruction of an experiment are
c  likely to fail at this time (1/26/00)
C-----------------------------------------------   

C-----------------------------------------------   
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
c      use read_old_wout
      use Vpname1
      use Vpname2
      use Vpname3
      use system_mod, only: system
      use Vindat2
      use optim, only: rbc_vv,zbs_vv,mpol_vv,ntor_vv,vv
      use vmec_input, only: rbc,rbs,zbc,zbs, lfreeb
      use read_namelist_mod
      use safe_open_mod

#ifdef OSF1
      use decdebug
#endif
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: lastplot = 18
      integer, parameter :: maxlist = 100
      logical, parameter :: lwstbf = .TRUE.
      character*(*), parameter :: version = 
     1   ' VMEC Plotter Version 6.10pg  y2k '

      character*(*), parameter :: tempfile = 'QqZzXLftemp'

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: id1, id2, pgopen
      integer :: numargs, numchars, lenlist, i, 
     1   js, lj, mn, l, j, n, ii, k, m, iunit,  
     2   ntor0, lp, ierr
      real(kind=rprec), dimension(:), allocatable :: szc, szb, dummy
      real(kind=rprec) :: r0
      character :: input_id*60, threed2_file*60, ictrans_sel*60, 
     1   explist*100, device*20
      character ::device_type*11
      character, dimension(0:lastplot) :: pagedir*50
      character*120 :: stripped_input_file, in_file

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

#ifdef OSF1
      call decdebugset
      numargs = iargc()
#endif

!!!!!!!!!! the above gets ladebug to accept input 
      pagedir=' '
      pagedir(0) = 'Outer magnetic flux surface'
      pagedir(1) = 'Flux, Mod-B, Mesh, J Contours'
      pagedir(2) = 'Jacobian  Contours'
      pagedir(3) = 'Free Boundary - Boundary Error'
      pagedir(4) = 'Convergence, Magnetic Axis'
      pagedir(5) = 'Stability and Curvatures'
      pagedir(6) = 'Fluxes, Iota, Well'
      pagedir(7) =  'Mass, Pressure, Beta, Pol. Cur vs phi'
      pagedir(8) =   'Jtor , J|| & Itor vs phi'
      pagedir(9) = 'Force Balance'
      pagedir(10) = '|B| Spectrum'
      pagedir(16) = '{R,Z,L} set - the DEFAULT is to skip these'
c      pagedir(17) = 'Limiter & Coil Positions'
c      pagedir(18) = 'Poloidal Flux'

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
 
*       Get the Command Line Arguments (Case & Plot Selection)

      numargs = iargc()
      if (numargs >= 1) then
         call getarg (1, input_id)
         ictrans_sel = '1 4 5 6 7 8 $'     ! the default is ALL PLOTS PLEASE
                                 ! except for the {R,Y,L} set
      endif
!!!!!!!!  bypass for ladebug and take from terminal window
      if (input_id .eq. '-h' ) then
         write (*, 104)
  104    format(/' ** To generate all plots %xprout fileid'/
     1' ** To generate selected plots %xprout fileid "1,5 12 14,$"'/
     2' ** choose plotting device %xprout fileid "1,5 12 14,$" /device'/
     3' ** typically choose /xw (default) or /cps for/device'/
     4      ' ** To generate surface plot %xprout fileid 0')
         write (*, 105) (i,pagedir(i),i=0,lastplot)
  105    format(/,' *** Plot Directory ***',2/,(i5,3x,a))
         stop 
      endif
      if (numargs .ge. 2) call getarg (2, ictrans_sel)
      if (numargs .eq. 3) call getarg (3, device)
      if (numargs .lt. 3) device='/xw'
      if(numargs .eq. 0) then
         write (6, 104)
         write (6, 105) (i,pagedir(i),i=0,lastplot)
         write(6,*)
     .'enter fileid <cr> , selections (no quotes) <cr>, device <cr>'
         read(5,fmt='(a)')input_id
         read(5,fmt='(a)')ictrans_sel
         read(5,fmt='(a)')device
         if(device .eq. ' ') device='/xw'
         write(6,*)input_id
         write(6,*)ictrans_sel           
         write(6,*)device
       endif

 
*       Parse the Plot Selection Part, If Present
      call icsel1 (ictrans_sel, lastplot, explist, maxlist, lenlist)

*
*       Does the Input File Exist ?
*
      numchars = len_trim(input_id)
      if (numchars .eq. 0) then
         write(6,*) ' MUST ENTER FILE SUFFIX ON COMMAND LINE'
         stop 
      endif
      if (index(input_id,'wout.') .eq. 1) then
         input_id = input_id(6:)
         numchars = numchars - 6
      endif
      runlabel = input_id(1:numchars+1)//' '//cdate(1:len_trim(cdate)+1)
      threed2_file = 'threed2.'//input_id
      gmeta_file = 'gmeta.'//input_id(1:numchars)
      open(unit=39, file=threed2_file, status='unknown')


      call read_wout_file('wout.'//trim(input_id),ierr)
      if (ierr .ne. 0 ) stop "call read_wout_file"

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
          enddo
        endif
      endif

      nmirnovset = 0
      nmirnov = 0
      do n = 1, nbsets
         if (n .eq. nmirnov_id) then
            nmirnovset = n
            nmirnov = nbfld(nmirnovset)
            exit 
         endif
      enddo
*****************************
* Read input file for RBC, etc.
* Needed for LFREEB=T
*****************************
C  need to strip comments from input file producved by stellopt
*****************************
!------------------------------------------------------------------------------

!
!     read input (optimum) namelist
!     first, strip any comments from input file (F95-not necessary)
!
      in_file='input.'//trim(input_id)
      call strip_comments(trim(in_file)) 
!!Produces clean file input_file//'.stripped'      
      stripped_input_file = trim(in_file) // '.stripped'
      id1 = 99
      call safe_open (id1, id2, trim(stripped_input_file)  ,
     .  'old', 'formatted')
      if (id2 .ne. 0) then
         print *,  trim(stripped_input_file)  ,
     .   ': input file open error: iostat = ', 
     .      id2
      else
         call read_namelist (id1, id2, 'indata')
      endif
      close(id1) !  close and reopen to avoid error on hecate
      if (id2 .ne. 0) then
         write (*, *) ' indata namelist READ error: iostat = ', id2
         ierr = 3
      endif
      ierr=0; id1 = 99; id2=0
      call safe_open (id1, id2, trim(stripped_input_file)  ,
     .  'old', 'formatted')
      if (id2 .ne. 0) then
         print *,  trim(stripped_input_file)  ,
     .   ': input file open error: iostat = ', 
     .      id2
      else
         call read_namelist (id1, id2, 'optimum')
      endif
      if (id2 .ne. 0) then
         write (*, *) ' optimum namelist READ error: iostat = ', id2
         ierr = ierr+6
      endif
      vv=lfreeb .and. ierr.lt.6 .and. rbc_vv(0,0).gt.0 
       print *,'vv=',vv, lfreeb,ierr,rbc_vv(0,0)
      vv=lfreeb .and. rbc_vv(0,0).gt.0
       print *,'vv=',vv,rbc_vv(0,0)
      if(id2.eq.0) close(id1,status='delete')
      if(id2.ne.0) close(id1)

*****************************

      allocate (sqrt_phimod(mnmax*ns), phimod(mnmax*ns), dbzcods(ns), 
     1   ffp(ns), ub(ns), iotaf(ns), darea(ns), iotazb(ns), psi(ns),
     2   itors(ns), ipols(ns), szc(ns), szb(ns), dummy(ns),
     3   ixm(mnmax), ixn(mnmax), stat=mn)
      if (mn .ne. 0) stop 'Allocation error in plotout'        

      ixm(:mnmax) = nint(xm(:mnmax))
      ixn(:mnmax) = nint(xn(:mnmax))

      do mn = 1, mnmax
         if (ixm(mn).eq.0 .and. ixn(mn).eq.0) mn0 = mn
      enddo

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
      enddo

      do mn = 1, ntor
         gmn(mn,1) = 1.5*gmn(mn,2) - 0.5*gmn(mn,3)
      enddo

      currvmn(1+ntor0:mnmax,1) = 0.
      bmn(1+ntor0:mnmax,1) = 0.
      gmn(1+ntor0:mnmax,1) = 0.
      darea(2:ns) = vp(2:ns)*overr(2:ns)
      overr(2:ns-1) = .5*(overr(2:ns-1)+overr(3:ns))
      iotaf(2:ns-1) = .5*(iotas(2:ns-1)+iotas(3:ns))
      jcuru(1) = 2.*jcuru(2) - jcuru(3)
      overr(1) = 2.*overr(2) - overr(3)
      jcuru(ns) = 2.*jcuru(ns-1) - jcuru(ns-2)
      overr(ns) = 2.*overr(ns-1) - overr(ns-2)
      iotas(1) = 1.5*iotas(2) - 0.5*iotas(3)
      iotaf(1) = iotas(1)
      iotaf(ns) = 1.5*iotas(ns) - 0.5*iotas(ns-1)
 
*     Spline Knots are on sqrt(s) grid ...
      if (ireconstruct .gt. 0) then
         do i = 2, ns
            szc(i) = sqrt(hs*(i - 1.5))
            szb(i) = sqrt(hs*(i - 1))
         enddo
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
cj      where (overr .ne. 0.) jcurv = jcurv/overr  !  obsolete
      hs = 1.0/(ns - 1)
      ohs = 1.0/hs
      r0 = rmnc(mn0,ns)
      call wrfcn (input_id)
      id2 = pgopen (trim(device)) 
      if (id2.le.0)
     .  write(*,*)'Fail to Open ',trim(device),' with code ',id2
      if (id2.le.0) stop 'id2'
      call pgscf(2)
      call pgpplotter (r0, explist, lenlist, input_id, trim(device))
      call wrfcn2 (input_id)
*
*                   CALL EQI STABILITY FILE GENERATOR
*
*     Poloidal Flux / TwoPi
      if (lwstbf .and. ireconstruct.gt.0) then
         psi(1) = 0.
         do j = 2, ns
            psi(j) = psi(j-1) + hs*phip(j)*iotas(j)
         enddo
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
8033  call pgclos

      end program pgpplotout






      
 
      subroutine dlm_parse(string, delimiter, n_tokens, tokens, 
     1   len_tokens)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  12:27:59  11/23/98  
C...Switches: -p4 -yb
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
      enddo
 
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
      enddo
      i_run = ipos - 1
 
      n_tokens = 0
      do while(i_run<=isl .and. n_tokens<n_tokens_max)
         i = index(string(i_run:),delimiter) + i_run - 1
         n_tokens = n_tokens + 1
         tokens(n_tokens) = string(i_run:i-1)
         len_tokens(n_tokens) = len_trim(tokens(n_tokens))
         i_run = i + 1
      enddo
 
c     Remove trailing delimiter at the end
      string(isl+1:isl+1) = ' '
 
      end subroutine dlm_parse
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
            if (string(i:i) .eq. set(j:j)) goto 100
         enddo
      enddo
 
      return 
 
  100 continue
 
c     Found a character in SET
 
      gen_find_first_in_set = i
 
      end function !gen_find_first_in_set

      subroutine get_lun(iu)
      integer*4  iun(20),ist(20)
      save    iun,ist
      data iun/  119,118,117,116,115,114,113,112,111,110,
     &             109,108,107,106,105,104,103,102,101,100/
      data ist/  20*0/
      do j = 1,20
        if (ist(j) .eq. 0) then
          iu = iun(j)
          ist(j) = 1
          return
        endif
      enddo
      iu = -1
      return
      entry    free_lun(iu)
      ist(120-iu) = 0
      return
      end subroutine get_lun

        subroutine graf1 (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in) :: x, y
         real ymin, ymax, xmin, xmax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         xmin=minval(x);xmax=maxval(x)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgslw(3)
         call pgline(n,x,y)
         call pgebuf
         call pgunsa
        return
        end subroutine graf1
        subroutine graf1x (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in) :: x, y
         real ymin, ymax, xmin, xmax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         xmin=minval(x);xmax=maxval(x)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,1)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgslw(3)
         call pgline(n,x,y)
         call pgebuf
         call pgunsa
        return
        end subroutine graf1x

        subroutine graf1pt (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer n
         real , dimension(n), intent(in)  :: x, y
         real ymin, ymax, siz
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgsci(1)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgpage
         call pgvstd
         call pgswin(minval(x),maxval(x),ymin,ymax)
         call pgbox('BCNT',0.,0,'BCNTP1',0.,0)
c         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(n,x,y,2)
         call pgunsa
        return
        end subroutine graf1pt

        subroutine midplpt (x1,y1,nx,ny,lx,ly,lt,runlbl)
         USE Vmagaxis
         implicit none
         integer nx,ny
         real , dimension(nx,ny), intent(in)  :: x1, y1
         real ymin, ymax, siz
         real, dimension(:), allocatable :: x,y
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgsci(1)
         allocate(x(2*ny),y(2*ny))
         x=(/x1(1,1:ny),x1((nx+1)/2,1:ny)/)
         y=(/y1(1,1:ny),y1((nx+1)/2,1:ny)/)
         ymin=minval(y);ymax=maxval(y)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(2*ny,x,y,1)
         call pgsci(4)
         call pgpt1 ( real(rmagaxis), ymin+(ymax-ymin)/3.,-4)
         call pgunsa
         deallocate(x,y)
        return
        end subroutine midplpt

        subroutine graf2 (x,y,yfit,n,lx,ly,lt,runlbl)
c  two lines on a single scale
         implicit none
         integer i, n
         real ymin, ymax, siz
         real *4, dimension(n) :: x,y,yfit
         character*(*) lx,ly,lt,runlbl
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=min(minval(y),minval(yfit))
         ymax=max(maxval(y),maxval(yfit))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(3)
         call pgpt(n,x,y,-2)
         call pgslw(3)
         call pgsls(1)
         call pgsci(2)
         call pgline(n,x,yfit)
         call pgsci(3)
         call pgpt(n,x,y,2)
         call pgebuf
         call pgunsa
        return
        end subroutine graf2

        subroutine graf2pt (x1,x2,y1,y2,n,lx,ly1,ly2,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x1, x2, y1, y2
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,lt,runlbl
         call pgsave
         call pgslw(1)
         call pgsci(1)
         xmin=minval(x1);xmax=maxval(x1)
         xmin=min(xmin,minval(x2))
         xmax=max(xmax,maxval(x2))
         ymin=minval(y1);ymax=maxval(y1)
         ymin=min(ymin,minval(y2))
         ymax=max(ymax,maxval(y2))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgpage
         call pgvstd
         call pgswin(xmin,xmax,ymin,ymax)
         call pgbox('BCNT',0.,0,'BCNTP1',0.,0)
c         call pgenv(xmin,xmax,ymin,ymax,0,0)
         call pglab(trim(lx),'(m)',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(4)
         call pgslw(1)
         call pgpt(n,x1,y1,-1)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgslw(3)
         call pgpt(n,x2,y2,-2)
         call pgmtxt('L',2.65,1.0,1.0,trim(ly2))
         call pgunsa
        return
        end subroutine graf2pt

        subroutine graf3pt (x,y1,y2,y3,n,lx,ly1,ly2,ly3,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x, y1, y2, y3
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,ly3,lt,runlbl
         call pgsave
         call pgsci(1)
         ymin=minval(y1);ymax=maxval(y1)
         ymin=min(ymin,minval(y2))
         ymax=max(ymax,maxval(y2))
         ymin=min(ymin,minval(y3))
         ymax=max(ymax,maxval(y3))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,0)
         call pglab(trim(lx),' ',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(1)
         call pgpt(n,x,y1,2)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgpt(n,x,y2,3)
         call pgmtxt('L',2.65,0.5,0.5,trim(ly2))
         call pgsci(4)
         call pgpt(n,x,y3,4)
         call pgmtxt('L',2.65,1.,1.,trim(ly3))
         call pgunsa
        return
        end subroutine graf3pt

        subroutine graf2x (x,y,yfit,n,lx,lyl,lyr,lt,runlbl)
c two lines in a box with y axis left and right
         implicit none
         integer i, n
         real ymin, ymax, dx, siz
         real *4, dimension(n) :: x,y,yfit
         character*(*) lx,lyl,lyr,lt,runlbl
         character*8 xopt,yopt
         dx=x(3)-x(1)
         if(n.lt.30)dx=x(2)-x(1)
         call pgsave
         call pgbbuf
         call pgsci(1)
         call pgpage
         call pgvstd
         call pgswin(minval(x),maxval(x),minval(y),maxval(y))
         xopt='BCNST';yopt='B'
         call pgbox(trim(xopt),0.,0,trim(yopt),0.,0)
c         call pgmtxt('B',1.5,0.5,0.5,trim(lx))
         call pglab(trim(lx),' ',trim(lt))
c         call pgmtxt('T',1.5,0.5,0.5,trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         yopt='BNST'
         call pgbox(' ',0.,0,trim(yopt),0.,0)
         call pgmtxt('L',2.25,0.5,0.5,trim(lyl))
         call pgsls(1)
         call pgslw(3)
         call pgline(n,x,y)
c         call pgarro(x(n/3)-dx,y(n/3),x(n/3-n/5)-dx,y(n/3))
c         call pgunsa
         xopt=' ';yopt='CMST'
         call pgsci(4)
         call pgswin(minval(x),maxval(x),minval(yfit),maxval(yfit))
         call pgbox(' ',0.,0,trim(yopt),0.,0)
         call pgmtxt('R',2.25,0.5,0.5,trim(lyr))
c         call pgsave
         call pgsls(4)
         call pgslw(8)
         call pgline(n,x,yfit)
c         call pgarro(x(2*n/3)+dx,yfit(2*n/3),x(2*n/3+n/5)+dx,yfit(2*n/3))
         call pgebuf
         call pgunsa
        return
        end subroutine graf2x

      subroutine plot2d(x,nx,y,ny,zin,n1,n2,w,
     &  size,iwidth,xlbl,ylbl,titl,iwid)
c     ---------------------------------------------------------------
      use coltabs
c
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nx, ny, n1, n2, iwid, iwidth
      real          x(*),y(*),zin(n1,n2),w(*),size
      character*(*) xlbl, ylbl, titl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iplot, ictab, ncband, icmax, icmin, maxbuf, ncb
      real :: difus, shin, polish, dofset, dlow, dhigh
      real, dimension(:,:), allocatable :: z
      real rgb(3,3), eye(3), light(3), latice(3,3)
      real lutusr(3,256), rgbbkg(3)
      character string*32, type*16, chr*16
      logical   ovrlay, lshin
C-----------------------------------------------
C   I n i t i a l i z a t i o n
C-----------------------------------------------
      data eye,light /0.0,0.0,1000.0,-1.0,-1.0,-1.0/
      data rgb /0.0,0.0,1.0,0.35,0.35,0.35,1.0,1.0,1.0/
      data rgbbkg /1.,1.,0./
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c purpose
c      this subroutine plots "data" defined on a regularly-spaced 
c   rectangular grid of points z(i,j). with the default choice for the 
c   pgcell routine that is linked, the output is a linearly-interpolated 
c   map (rather than coarse rectangular boxes).
c
c parameters
c   argument  type  i/o  dimension  description
c    z        r*4    i    n1 x n2   the recatangular "data"-array.
c    z        r*4    o    n1 x n2   a scaled, and clipped, version of 
c                                   the input array(!).
c    n1       i*4    i       -      the first dimension of array z.
c    n2       i*4    i       -      the second dimension of array z.
c    x        r*4    i      nx      array of x-coordinates.
c    nx       i*4    i       -      number of x-pixels to be plotted 
c                                  (usually = n1, but must be <= n1).
c    y        r*4    i      ny      array of y-coordinates.
c    ny       i*4    i       -      number of y-pixels to be plotted 
c                                  (usually = n2, but must be <= n2).
c    w        r*4    i       -      dummy parameter (back-compatibility).
c    size     r*4    i       -      character-size for plot (try 1.5).
c    iwidth   i*4    i       -      line-width for plot (try 2).
c    xlbl     a*1    i     *(*)     label for x-axis.
c    ylbl     a*1    i     *(*)     label for y-axis.
c    titl     a*1    i     *(*)     title for plot.
c
c globals
c    coltabs.inc
c
c history
c   initial release.                                    dss:  3 jul 1992
c   minor changes to conform with new pgcell.           dss:  6 feb 1995
c   put in option to over-lay contours.                 dss: 21 feb 1995
c   now has proper 3-d surface surface rendering.       dss: 27 aug 1997
c   fortran made linux-friendly!                        dss: 15 sep 1997
c   choose white background colour for postscript.      dss:  5 feb 1999
c-----------------------------------------------------------------------
c
c
   1  iplot=1
        if (n1.ne.nx .or. n2.ne.ny) then
          write(*,*)' sorry folks, surface option needs n1=nx & n2=ny!'
          return
        endif
        allocate(z(n1,n2))
        z=zin(1:n1,1:n2)
        ictab=0
        shin=0.
        ncband=15
        difus=.7
        ncb=27
        polish=1
        lshin=.false.
        call eulerfix(latice)
        dhigh=maxval(z)
        dlow=minval(z)
        dofset=max(-dlow,0.)
        z=z+dofset
        dhigh=maxval(z)
        dlow=minval(z)
        call pgqcol(icmin,icmax)
        ncband=min(ncband,icmax-17+1)
        call pgsch(size)
        call pgslw(iwidth)
c        call pgpaper(0.0,1.0)
        call pgvport(0.0,1.0,0.0,1.0)
        call pgswin(-1.,1.,-1.,1.)
        call pgsci(0)
        call pgbox('BC',0.0,0,'BC',0.0,0)
        call pgsci(1)
          call sbfint(rgbbkg,15,1,iwid,maxbuf)
        if (ictab.le.2 .or. ictab.ge.7) then
          call colint(rgb,17,icmax,difus,shin,polish)
        elseif (ictab.eq.3) then
          call colsrf(heat,256,1.0,17,icmax,ncb,difus,shin,polish)
        elseif (ictab.eq.4) then
          call colsrf(spectrum(1,2),255,1.0,17,icmax,ncb,difus,shin,
     *                polish)
        elseif (ictab.eq.5) then
          call colsrf(bgyrw,256,1.0,17,icmax,ncb,difus,shin,polish)
        elseif (ictab.eq.6) then
          call colsrf(serp,256,1.0,17,icmax,ncb,difus,shin,polish)
        endif
        call sb2srf(eye,latice,z,n1-1,n2-1,dlow,dhigh,1.0,17,icmax,ncb,
     *              light,lshin)
        call pgsci(1)
        call axes3d(eye,latice,x(1),x(n1),y(1),y(n2),xlbl,ylbl,size,
     *         dlow,dhigh,dofset,z(1,1),z(n1,1),z(n1,n2),z(1,n2))
        call sbtext(eye,0.,0.5,0.5,titl)
        call sbfcls(iwid)
        call pgmtxt("T",-1.72,0.5,0.5,titl)
        call pgebuf
        call pgupdt
        deallocate(z)
        return
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC DRIVER CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c2d        call pgsave
c2d        call pgsubp(1,1)
c2d        call plot2d(real(hshere(1:)),nx,
c2d     1   real(lka(nx,1:)),ny,
c2d     1   real(residual12(1:nx,1:ny)),nx,ny,
c2d     1   real(br),1.5,1,'hs','mn','residual',id2)
c2d
c2d          call pgqinf('DEV/TYPE',device_type,kz)
c2d          if((device_type(2:4).eq.'/XW').or.
c2d     &     (device_type(2:4).eq.'/XS'))then 
c2d           write(6,*)trim(device_type),' Type <RETURN> for next page:'
c2d           read(5,fmt='(a)')pchar
c2d          endif
c2d        call pgunsa
c2d
c2d
c2d        call pgsave
c2d        call pgsubp(1,1)
c2d        call plot2d(real(hshere(1:)),nznt1,
c2d     1   real(pang(1:)),ntheta2/2+1,
c2d     1   real(r3d(1:nznt1,1:ntheta2/2+1,1)),nznt1,ntheta2/2+1,
c2d     1   real(br),1.5,1,'hs','theta','residual',id2)
c2d
c2d          call pgqinf('DEV/TYPE',device_type,kz)
c2d          if((device_type(2:4).eq.'/XW').or.
c2d     &     (device_type(2:4).eq.'/XS'))then 
c2d           write(6,*)trim(device_type),' Type <RETURN> for next page:'
c2d           read(5,fmt='(a)')pchar
c2d          endif
c2d        call pgunsa
c2d
c2d
c2d        call pgsave
c2d        call pgsubp(1,1)
c2d        call plot2d(real(hshere(1:)),nznt1,
c2d     1   real(tang(1:)),nzeta/2,
c2d     1   real(r3d(1:nznt1,1,1:nzeta/2)),nznt1,nzeta/2,
c2d     1   real(br),1.5,1,'hs','zeta','residual',id2)
c2d
c2d          call pgqinf('DEV/TYPE',device_type,kz)
c2d          if((device_type(2:4).eq.'/XW').or.
c2d     &     (device_type(2:4).eq.'/XS'))then 
c2d           write(6,*)trim(device_type),' Type <RETURN> for next page:'
c2d           read(5,fmt='(a)')pchar
c2d          endif
c2d        call pgunsa
c2d        call pgsubp(2,2)
c2dc        call pgpage
c2d        call resetcolor
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC END DRIVER CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      end subroutine plot2d

      subroutine eulerfix(latice)
c     ------------------------
c
      real         latice(3,*)
      data         pirad /0.01745329252/
      data ia,ib/45,30/
c
        sina=-sin(float(ia)*pirad)
        cosa=cos(float(ia)*pirad)
        sinb=-sin(float(ib)*pirad)
        cosb=cos(float(ib)*pirad)
      call roty(-0.5,-0.5,+0.5,u,v,w,sina,cosa)
      call rotx(u,v,w,latice(1,1),latice(2,1),z1,sinb,cosb)
      latice(3,1)=z1-1.0
      call roty(+0.5,-0.5,+0.5,u,v,w,sina,cosa)
      call rotx(u,v,w,latice(1,2),latice(2,2),z2,sinb,cosb)
      latice(3,2)=z2-1.0
      call roty(-0.5,-0.5,-0.5,u,v,w,sina,cosa)
      call rotx(u,v,w,latice(1,3),latice(2,3),z3,sinb,cosb)
      latice(3,3)=z3-1.0
      end subroutine eulerfix
 
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
c      enddo
 
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
         enddo
      enddo
 
      end subroutine icsel1

      logical function iswhitespace (letter, byteset, nwsp)
C...Translated by Pacific-Sierra Research VAST-90 2.06G2  12:27:59  11/23/98  
C...Switches: -p4 -yb
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
      enddo
      iswhitespace = whitespace

      end function !iswhitespace

        subroutine loggraf1 (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer i, n
         real ymin, ymax, siz
         real *4, dimension(n) :: x,y
         real *4, dimension(:), allocatable :: logx,logy
         character*(*) lx,ly,lt,runlbl
         character*100 loglabel
         allocate(logy(n))
       logy=alog10(y)
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=minval(logy);ymax=maxval(logy)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly)//'\(2224)'
         call pglab(
     &     trim(lx),trim(loglabel),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgslw(3)
         call pgline(n,x,logy)
         call pgebuf
         call pgunsa
         deallocate(logy)
        return
        end subroutine loggraf1


        subroutine loggraf3pt (x,y1,y2,y3,n,lx,ly1,ly2,ly3,lt,runlbl)
         integer n
         real , dimension(n), intent(in)  :: x, y1, y2, y3
         real *4, dimension(:), allocatable :: logx,logy1,logy2,logy3
         real ymin, ymax, siz
         character*(*) lx,ly1,ly2,ly3,lt,runlbl
         character*100 loglabel
         allocate(logy1(n),logy2(n),logy3(n)   )
       logy1=alog10(y1)
       logy2=alog10(y2)
       logy3=alog10(y3)
         call pgsave
         call pgsci(1)
         ymin=minval(logy1);ymax=maxval(logy1)
         ymin=min(ymin,minval(logy2))
         ymax=max(ymax,maxval(logy2))
         ymin=min(ymin,minval(logy3))
         ymax=max(ymax,maxval(logy3))
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly2)//'\(2224)'
         call pglab(trim(lx),' ',trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(1)
         call pgpt(n,x,logy1,2)
         call pgmtxt('L',2.65,0.,0.,trim(ly1))
         call pgsci(2)
         call pgpt(n,x,logy2,3)
         call pgmtxt('L',2.65,0.5,0.5,trim(loglabel))
         call pgsci(4)
         call pgpt(n,x,logy3,4)
         call pgmtxt('L',2.65,1.,1.,trim(ly3))
         call pgunsa
        return
        end subroutine loggraf3pt
        subroutine loggrafpt (x,y,n,lx,ly,lt,runlbl)
         implicit none
         integer i, n
         real ymin, ymax, dx, siz
         real *4, dimension(n) :: x,y
         real *4, dimension(:), allocatable :: logx,logy
         character*(*) lx,ly,lt,runlbl
         character*100 loglabel
         allocate(logy(n))
         logy=alog10(y)
         call pgsave
         call pgbbuf
         call pgsci(1)
         ymin=minval(logy);ymax=maxval(logy)
         if(ymax.eq.ymin)then
          if(ymax.gt.0.)ymin=.99*ymax
          if(ymax.lt.0.)ymax=.99*ymin
          if(ymax.eq.0.)ymax=1.e-3;ymin=-1.e-3
         endif
         call pgenv(minval(x),maxval(x),ymin,ymax,0,20)
         loglabel='Log\d10\u\(2223)'//trim(ly)//'\(2224)'
         call pglab(
     &     trim(lx),trim(loglabel),trim(lt))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         call pgsci(2)
         call pgpt(n,x,logy,2)
         call pgsci(15)
         call pgslw(1)
         call pgline(n,x,logy)
         call pgebuf
         call pgunsa
         deallocate(logy)
        return
        end subroutine loggrafpt

      subroutine logint(contr,zmin,zmax,n,irr)
        implicit none
c logarithmic interpolation
        real zmin,zmax,zlmin,zlmax,amin, amax
        logical misign,allpos,zerol,zerou
        complex c,d,bmin,bmax
        real cs, conv, pi, el, da
        integer i,n,n1, irr
        real*4 contr(*)
        pi=4.*atan(1.)
        if(zmax-zmin .lt. 10.)goto 1001  !  linear
        zlmax=zmax
        zlmin=zmin
        zerol=zmin.eq.0.
        zerou=zmax.eq.0.
        if(zerol)zlmin=-zmax/1.e4
        if(zerou)zlmax=abs(zmin/1.e4)
        conv=log(10.)
        allpos=(zlmax.gt.0.).and.(zlmin.gt.0.)
        misign=zlmin.lt.0.
        bmin=cmplx(zlmin,.0)
        bmin=(log(bmin))
        if(misign)cs=real(bmin)
        if(misign)bmin=cmplx(cs,pi)
        amin=real(bmin)
        if(zerol)amin=0
        if(misign)amin=-amin
        misign=zlmax.lt.0.
        bmax=cmplx(zlmax,.0)
        bmax=(log(bmax))
        if(misign)cs=real(bmax)
        if(misign)bmax=cmplx(cs,pi)
        amax=real(bmax)
        if(zerou)amax=0
        if(amin.gt.amax)then
          cs=amin
          amin=amax
          amax=cs
        endif
        da=real(amax-amin)/(n-1)
        if(allpos) then
          do i=1,n
          c=cmplx(amin+(i-1)*da)
          contr(i)=real(exp(c)) 
          enddo
        else
          el=abs(zlmax)+abs(zlmin)
          n1=-(zlmin/el)*n
          n1=max(1,n1)
          do i=1,n1
          c=cmplx(amin+(i-1)*da,pi)
          contr(i)=real(exp(-c)) 
          enddo
          i=n1+1
          contr(i)=0
          do i=n1+2,n
          c=cmplx(amin+(i-1)*da)
          contr(i)=real(exp(c)) 
          enddo
        endif
        irr=0
        return

 1001   continue     
        da=real(zmax-zmin)/(n-1)
        do i=1,n
           contr(i)=zmin +da*(i-1)
        enddo
        irr= 1001
        return
        end subroutine logint

      logical function makeplot (explist, lenlist, page, pagedesc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lenlist
      character explist*(*), page*(*), pagedesc*(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      logical :: flag, echo
      character :: id*10
C-----------------------------------------------
 
 
c     Local Vbls
      data echo/.TRUE./
 
 
c     Local Fncs
 
      id = '-'//page(1:len_trim(page))//'-'
 
      flag = index(explist(1:lenlist),id(1:len_trim(id))) > 0
 
      if (echo .and. flag) write (*, *) ' Making plot # ', page, ' : ', 
     1   pagedesc
 
      makeplot = flag
 
      end function !makeplot




      subroutine pgpcontour(xplot, yplot, func, nx,ny,
     1    ncon, ivar, mxnp, nmxnp, runlabel,xlab,ylab,zlab,mtitl)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use Toexternal
      USE READ_WOUT_MOD
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
      integer ncon, ivar, nx, ny, ny1
      real(kind=rprec) :: grsz
      character xlab*(*), ylab*(*), zlab*(*), mtitl*(*), runlabel*(*)
      real(kind=rprec), dimension(nx*ny) :: xplot, yplot, func
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nstep, nrav, nzav, n, k, l, kt, j, mxnp, nmxnp
      real, dimension(ny) :: rthet, zthet
      real :: step, exstep, fac, xmin, xmax, ymin, ymax, 
     1  delmax, delmin, delf, siz
      real, dimension(:), allocatable :: f_real, x_real, y_real
      real, dimension(:,:), allocatable :: contourfunc
      logical logs
      data logs/.false./
      character*35 logscale
      data logscale/'Using Logarithmically Scaled Levels'/
                                               
C-----------------------------------------------
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real, EXTERNAL :: RZtrans
C-----------------------------------------------

*        
*     LOCATE PHYSICAL ORIGIN
*
*     ROUND-OFF AXIS
*
      mxnt=nx
      step = max(rmax - rmin,zmax - zmin)
      exstep = 2. - nint(log10(step))
      fac = 10.**exstep
      nstep = nint(.25*(fac*step)*1.05)
      nrav = nint(.50*(rmax + rmin)*fac)
      nzav = nint(.50*(zmax + zmin)*fac)
      step = nstep/fac


      allocate (f_real(nx*ny+nx), x_real(nx*ny+nx), 
     1          y_real(nx*ny+nx),contourfunc(nx,ny+1))
      f_real = (/real(func),real(func(1:nx))/) ! complete the period
      x_real = (/real(xplot),real(xplot(1:nx))/)
      y_real = (/real(yplot),real(yplot(1:nx))/)

      xmin = minval(x_real)
      xmax = maxval(x_real)
      ymin = minval(y_real)
      ymax = maxval(y_real)
      call pgsci(1)
      call  pgenv (xmin, xmax, ymin, ymax,1,0)      ! define the subplot area
      call pglab(trim(xlab),trim(ylab),trim(zlab))
      call pgmtxt('T',1.0,0.5,0.5,trim(mtitl))
      call pgqch(siz)
      call pgsch(siz-.5)
      if(logs)call pgmtxt('T',0.20,0.5,0.5,trim(logscale))
      call pgmtxt('B',2.85,0.5,0.5,trim(runlabel))
      call pgsch(siz)
      
c         Plot the Boundary -- Last Closed Flux Surface
      l = ny*(ny - 1) + 1
      call  pgscf (2  )                             ! select simple typeface
      call pgpt 
     & (nx,x_real(nx*(ny-1)+1:nx*ny), y_real(nx*(ny-1)+1:nx*ny),-1)
c plot contours 
         phit = 2.0_rprec*pi/real(nfp*nmxnp,kind=rprec)* 
     &    real((mxnp-1),kind=rprec)
      ny1=ny+1
      contourfunc=reshape(f_real,(/nx,ny1/))
      delmax = maxval(contourfunc)
      delmin = minval(contourfunc)
      delf = (delmax - delmin)/ncon
      contrs=0.
      lcontr=0
      klabel=0
      kbold=0
      logs=.false.
      ncont=ncon
      do j=1,ncon
        contrs(j)=delmin +delf*(j-1)
c          lcontr(j)=1+modulo(j,5)
c          if(lcontr(j).eq.2)klabel(j)=1
      enddo
ceal      if(abs((delmax-delmin)/delf).gt.10.)then
ceal         call logint(contrs,delmin,delmax,ncon,kt)
ceal         if(kt.eq.0)logs=.true.
ceal       endif
      call pgsave
      call pgconx(Transpose(contourfunc),
     &ny1,nx,1,ny1,1,nx,contrs,ncon,RZtrans)
      call pgunsa
c         Plot the Magnetic Axis ?
      if (ivar .eq. 0) then
         call pgpt1 ( real(rmagaxis), real(zmagaxis),3)
      endif
 
c         Plot Constant THETA Contours?
      if (ivar .eq. 2) then
         do kt = 1, nx, 6
            do j=1,ny
              rthet(j) = r(kt+nx*(j-1))
              zthet(j) = z(kt+nx*(j-1))
            enddo
            call pgline ( ny, rthet, zthet )
         enddo
      endif
 
      deallocate (f_real, x_real, y_real, contourfunc)
      end subroutine pgpcontour

      subroutine pgptdplotter(rmnb, zmnb, rdum, zdum,ns,iwid)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: nfp, rprec
      use Vpltcn2
      use coltabs
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns
      real(kind=rprec), dimension(ns), intent(in) :: 
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
        integer, parameter :: isrs = 3
c
c declare variables to hold labels.
c
        character*64 :: xnlb,ynlb,znlb,xilb,yilb,zilb
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real red(3,256), grey(3,256), green(3,256), blue(3,256) 
      real map(3,256)
        real, dimension(3) :: eye,rgbbkg,rgbsrf,p1,p2,light
        real  zshift,xw,shin,polish,difus,alpha
        integer :: k, k1, k2, iwid, maxbuf, pgopen
        integer :: ic1, ic2, ncb, ncol
c        data       rgbbkg /0.25,0.25,0.25/
        data       rgbbkg /1.,1.,0./
        data       rgbsrf /1.00,0.00,0.00/
        data       light /-1.0,-1.0,-1.0/
        integer :: jdim, mtri, ntri, jplt, j, i, jvvmi, jvvma,
     1     iuvmi, iuvma, lnlg, ltri
        integer, dimension(:), allocatable :: itwk
        real, dimension(:,:), allocatable :: rval, hval
        real(kind=rprec), dimension(:), allocatable :: rt,ht
c
c declare a local array to hold the triangle list and a couple of
c temporary variables to be used in sorting the list.
c
        real, dimension(:,:), allocatable :: rtri, rtwk, vert
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
        zval(r,p,h)=h-zshift  ! pgplot convention all z<0
        zshift=1.1*zmax

        jplt = 51/nfp        !!v points per period
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
        allocate (rtri(10,mtri),rtwk(mtri,2),itwk(mtri), vert(3,3),
     1   rval(idim,2), hval(idim,2), rt(idim+1),ht(idim+1),stat=i)
        if (i .ne. 0) stop 'allocation error in 3dplot'

c
c turn clipping off.
c
c
c double the line width.
c
c
c define colors to use.
c

c
c select font number 25, turn on the outlining of filled fonts, set the
c line width to 1, and turn off the setting of the outline color.
c
c
c make tdpack characters a bit bigger.
c
c
c define tdpack rendering styles 1 through 7, using black-and-white
c shading or colored shading, whichever is selected.  the indices
c 1-7 can then be used as final arguments in calls to tditri, tdstri,
c and tdmtri.
c

c              
c initialize the count of triangles in the triangle list.
c
        ntri=0
c
c for each box on a rectangular grid in the uv plane, generate two
c triangles and add them to the triangle list.  each triangle is
c transformed from cylindrical coordinates to cartesian coordinates.
c
        j=1;rval=0;hval=0
        call totz(idim,1,jplt,j,rt,ht,
     1              rmnb,zmnb,rdum,zdum)
        hval(1:idim,1)=real(ht(1:idim))
        rval(1:idim,1)=real(rt(1:idim))

        do 102 j=1,jdim-1
          call totz(idim,1,jplt,j+1,rt,ht,
     1              rmnb,zmnb,rdum,zdum)
          hval(1:idim,2)=real(ht(1:idim))
          rval(1:idim,2)=real(rt(1:idim))
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
            endif
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
            endif

  101     continue
          rval(:,1) = rval(:,2)
          hval(:,1) = hval(:,2)
  102   continue
           if(max(
     &        maxval(rtri(3,1:ntri)),
     &        maxval(rtri(9,1:ntri)),
     &        maxval(rtri(6,1:ntri))).gt.0)stop 'z>0'

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
c
c

c solids!
      xw=1.35*xmax
      eye=(/xeye,yeye,2*zeye/)
      light=(/-xeye,yeye,-2*zeye/)
      call pgsave
      call pgpaper(0.0,1.0)
      call pgvport(0.0,1.0,0.0,1.0)
c      call pgswin(-0.95*xw,xw/1.23,-0.95*xw,xw/1.23)
      call pgswin(-1.23*xw,xw/1.23,-1.23*xw,xw/1.23)
      call pgbox('bc',0.0,0,'bc',0.0,0)
      call colint(rgbsrf,16,39,0.5,0.0,1.0)
      alpha=1.0
      call pgqcir(ic1,ic2)
      ncb=5!!(ic2-ic1+1)/6*0
      ncol=256
      shin=0.0
      polish=.5
      difus=.25
      red=0.0
      green=0.0
      blue=0.0
      do i = 1,256
        grey(1,i)=1.0/real(256-1)*(i-1)
        grey(2,i)=1.0/real(256-1)*(i-1)
        grey(3,i)=1.0/real(256-1)*(i-1)
        red(1,i)=1.0/real(256-1)*(i-1)
        green(2,i)=1.0/real(256-1)*(i-1)
        blue(3,i)=1.0/real(256-1)*(i-1)
      enddo
      map = 0.8*red+0.8*green+0.5*blue


      call sbfint(rgbbkg,15,1,iwid,maxbuf)
      call colsrf(map,ncol,alpha,ic1,ic2,ncb,difus,shin,polish)
      do i=1,ntri
         vert=reshape(rtri(1:9,i),(/3,3/))
         call sbplan(eye,3,vert,ic1,ic2,light)
      enddo

      call sbfcls(iwid)
      call pgebuf
      call pgupdt
      call pgunsa
      deallocate (rtri,rtwk,itwk, rval, hval, stat=i)
      end subroutine pgptdplotter


      subroutine pgpplotter(r0, explist, lenlist, input_id, device)
C-----------------------------------------------
CSave as working, has axes3d, offset corrected
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, itor_w=>itor
      use Vpname2
      USE Vpltcn2
      USE Vpltcn6
      USE Vrzarray
      USE Vplotdata
      USE Vmagaxis
      USE Vindat2
      use Toexternal
      use system_mod, only: system
      use vmec_input, only: rbc,rbs,zbc,zbs, lfreeb
      use safe_open_mod

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer lenlist
      real(kind=rprec) :: r0
      character explist*(*), input_id*(*), device*(*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec) :: zero = 0.0_dp , one = 1.0_dp
      integer, parameter :: ntdu = 31, ntdv = 21
      integer, parameter :: nplot1 = 25, ntheta1 = 2, 
     1  nox = 1, noy = 3, lx = -1, ly = 0, iend = 3,
     2  ilog = 1, icart = 0, oneppg = 2, fourppg = 1
      integer, dimension(4) :: ivar = (/0, 1, 2, 3/),
     1  ncon = (/10, 25, 15, 20/)      
      real(kind=rprec), parameter :: grsize1 = 3.5, grsize2 = 2.5, 
     1  grfull = 7.5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: id1, id2, np, id, pgopen, nzeta, ntheta, ntheta2
      integer :: i, ii, nplots, les, kz, k, kt, j, l, js, kz1, nplt1,
     1   ndeg, ns1, mn, mp, nt, ichar, istat, nznt, nznt1, ipag
      integer :: nsurf , ntht , nphi ,jxbcode
      real(kind=rprec) :: small 
      real(kind=rprec), dimension(ns) :: br, bz, phin
      real(kind=rprec), dimension(2*ns) :: oqmid
      real(kind=rprec), dimension(100) :: time
      real(kind=rprec), dimension(nplot1) :: raxis_v, zaxis_v
      real(kind=rprec), dimension(nfloops + 1) :: delflm
      real(kind=rprec), dimension(nbloops) :: delbc
      real(kind=rprec), dimension(nfloops + 1) :: indflm
      real(kind=rprec), dimension(nbloops) :: indbc
      real(kind=rprec), dimension(:,:), allocatable :: dummy1
      real(kind=rprec), dimension(:), allocatable :: rbndy,zbndy
      real(kind=rprec), dimension(:), allocatable :: modb, sqrt_phim, 
     1  phim, gsqrt, torcur, dummy2, r12, z12, ru12, zu12
      real(kind=rprec) :: dth, denom, bdotgradv0, phiangle, offset,
     1   bdotgradvn, phinorm, t1, a0, an, cosphi, sinphi
      character :: page*10, pagedesc*50, nchar*100, mchar*100, 
     1   pchar*100, ititlet*100, fname*80, fdum*80,
     1   xlabel*100, ylabel*100, zlabel*100, mtitle*40
      character ::device_type*11
      real xspread,xxspread,yspread,yyspread,deltx,siz,
     1 xmin,xmax,ymin,ymax,yzero,ymmax,ymmin,xmmin,xmmax
      integer maxc,ic,nmarks,nspread
      real(kind=rprec), dimension(:), allocatable :: hiota, presb, 
     1   beta_volb, phipb,  phib, bvcob,  bucob , xmb, xnb
      integer , dimension(:), allocatable :: jlistb
      real(kind=rprec), dimension(:,:), allocatable :: 
     1   bmnb, rmnb, zmnb, pmnb
      logical first
      data first/.true./
      real(kind=rprec) :: aspectb, rmax_surfb, rmin_surfb, betaxisb
      integer :: nsb, nfpb, nboz, mboz, mnboz, versionb, 
     1 iunit, ierr, iread, jrad, mb, nb

C-----------------------------------------------
c jxbout columns
c    LK      JSUPU      JSUPV      JSUPS      BSUPU      BSUPV      J X B        
c    p'   J X B - p'     J * B      BSUBU      BSUBV      BSUBS   

      integer :: nx, ny
      integer, dimension(:,:), allocatable :: lka
      real(kind=rprec), dimension(:,:), allocatable :: 
     1  bsupu, bsupv, bsubu, bsubv, bsubs, residual12, xfres12
      real(kind=rprec), dimension(:,:), allocatable :: 
     1  jsupu, jsupv, jcrossb, pprime, bdotj, jsups
      real(kind=rprec), dimension(:), allocatable :: pang, tang,
     1  avforce, hshere, amaxfor, aminfor
      real, dimension(:), allocatable :: xplot,yplot
      logical, dimension(:), allocatable :: iplot
      real x1,x2,y1,y2
      integer itheta,izeta
      logical addon, alldone
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      logical , EXTERNAL :: makeplot
C-----------------------------------------------
      ipag=0


      mn = mnmax*ns
      allocate (dummy1(mnmax,ns), dummy2(nrt), r(nrt), z(nrt), 
     1   ru(nrt), zu(nrt), modb(nrt), sqrt_phim(nrt), phim(nrt), 
     2   gsqrt(nrt),torcur(nrt), dbsubudstz(nrt), bsubutz(nrt), 
     3   bsubvtz(nrt),r12(nrt), z12(nrt), ru12(nrt), zu12(nrt), 
     4   stat = istat)
      if (istat .ne. 0) stop 'allocation error in plotter'
      allocate(rbndy(nthpts),zbndy(nthpts), stat = istat)
      if (istat .ne. 0) stop 'allocation error in plotter'


      dummy1 = 0.
      dummy2 = 0.
 
      noxc = 0
      noyc = noy - 1
      lxc = lx
      lyc = ly
 
*
*                 SETUP PGPLOT GRAPHICS
*
      call pgpap(7.5,1.)
      call pgask(.true.)
      call pgsch(1.85)
      call pgslw(3)
      call pgsubp(1,1)

 
*
*                 COMPUTE R,Z SCALES
*

      nplots = 1
      if (ntor .ne. 0) nplots = 4
      les = 1 + nthpts*(ns - 1)
      do kz = 1, nplots
         call totz(nthpts,ns,nplots,kz,r,z,rmnc,zmns,rmns,zmnc)
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
      enddo
!
!     3-D Surface
!      
C###    PAGE # 0
      page = '0'
      pagedesc = 'Outer magnetic flux surface'

      if (makeplot(explist,lenlist,page,pagedesc)) then
        call pgsave
        call pgptdplotter(rmnc(1,ns),zmns(1,ns),bmn(1,ns),
     &   rmns(1,ns),zmnc(1,ns),ns,id2)

        call pgqinf('DEV/TYPE',device_type,kz)
        if((device_type(2:4).eq.'/XW').or.!pgplot does not know about
     &   (device_type(2:4).eq.'/XS'))then  ! sbfcls -- fix it here
          write(6,*)trim(device_type),' Type <RETURN> for next page:'
          read(5,fmt='(a)')pchar
        endif
        call pgunsa
      endif
      if(trim(explist).eq.'-0-')goto 8033  !  terminate plotting
      call resetcolor
*       ---Transform from (m,n) Space to (u,v) Space---
      do kz = 1, nplots
         call totz (nthpts,ns,nplots,kz,r,z,rmnc,zmns,rmns,zmnc)
         call totzu(nthpts,ns,nplots,kz,ru,zu,rmnc,zmns,rmns,zmnc)
         call totz (nthpts, ns, nplots, kz, modb, dummy2, bmn, 
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, sqrt_phim, dummy2, 
     1      sqrt_phimod, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, phim, dummy2, phimod, 
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, gsqrt, dummy2, gmn,
     1      zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, torcur, dummy2, 
     1      currvmn, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, bsubutz, dummy2, 
     1      bsubumn, zmns, dummy1, dummy1)
         call totz (nthpts, ns, nplots, kz, bsubvtz, dummy2, 
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
            enddo
         enddo
         rmagaxis = r12(1)
         zmagaxis = z12(1)
c         call bndypts (nthpts,nplots,kz,rbndy,zbndy)
 
         do kt = 1, nthpts
            do j = 2, ns
               l = kt + nthpts*(j - 1)
               ru12(l) = 0.5*(ru(l)+ru(l-nthpts))
               zu12(l) = 0.5*(zu(l)+zu(l-nthpts))
            enddo
         enddo
 
**      Finally ... <J-dot-B> / <B-dot-gradPhi>
        where (bdotgradv .ne. 0.) jdotb = jdotb/bdotgradv
      
         ndeg = nint(360.*(kz - 1)/nplots)
   30    format(a,i3,a)
         noxc = 0
 
         call pgsubp(2,2)

C###    PAGE # 1
         page = '1'

         pagedesc = 'Flux, Mod-B, Mesh, J Contours'
         if (makeplot(explist,lenlist,page,pagedesc)) then
           kz1=kz
           nplt1=nplots
           xlabel='R (m)'
           ylabel='Z (m)'
           zlabel='\(2267)\gF Contours'
           phit = 2.0_rprec*pi/real(nfp*nplt1,kind=rprec)* 
     &       real((kz1-1),kind=rprec)
           write(nchar,'(f5.1)')phit*180/real(pi)*nfp
           mtitle='N\dfp\u\(0647)= '//adjustl(nchar)
           labelformat='(2h L ,1pe9.2)'
           call pgpcontour(r(nthpts+1),
     &         z(nthpts+1), sqrt_phim(nthpts+1),
     &         nthpts,ns-1,ncon(1), ivar(1),
     &         kz1, nplt1, runlabel,
     &         xlabel,ylabel,zlabel,mtitle)
           zlabel='MOD-B CONTOURS'
           call pgpcontour (r(nthpts+1),
     &         z(nthpts+1),  modb(nthpts+1),
     &         nthpts,ns-1,ncon(3), ivar(2), 
     &         kz1, nplt1, runlabel,
     &         xlabel,ylabel,zlabel,mtitle)
           zlabel='\gF and \gH CONTOURS'
           call pgpcontour (r(nthpts+1),
     &         z(nthpts+1),  phim(nthpts+1),
     &         nthpts,ns-1,ncon(1), ivar(3),
     &         kz1, nplt1, runlabel,
     &         xlabel,ylabel,zlabel,mtitle)
         zlabel= 'CONTOURS OF RJ\.\(2266)[\gQ]'
         call pgpcontour (r(nthpts+1),
     &         z(nthpts+1),  torcur(nthpts+1),
     &         nthpts,ns-1,ncon(3), ivar(2),
     &         kz1, nplt1, runlabel,
     &         xlabel,ylabel,zlabel,mtitle)
         endif !makeplot  PAGE # 1
C###    PAGE # 2
         if (page .ne. '1')goto 8020
         page = '2'

         pagedesc = 'Jacobian Contours & Midplane Slices'
         if (.not.makeplot(explist,lenlist,page,pagedesc)) goto 8020
           zlabel='JACOBIAN CONTOURS'
           call pgpcontour (r12(nthpts+1), 
     &         z12(nthpts+1),  gsqrt(nthpts+1),
     &         nthpts,ns-1,ncon(1), ivar(1), 
     &         kz1, nplt1, runlabel,
     &         xlabel,ylabel,zlabel,mtitle)
 
           xlabel='R (m)'
           ylabel='MidPlane Values'
           zlabel= 'CONTOURS OF RJ\.\(2266)[\gQ]'
         call midplpt
     &      (real(r(1:nrt)),real(torcur(1:nrt)),
     &       nthpts,ns,xlabel,ylabel,zlabel,runlabel)
           zlabel='MOD-B CONTOURS'
         call midplpt
     &      (real(r(1:nrt)),real(modb(1:nrt)),
     &       nthpts,ns,xlabel,ylabel,zlabel,runlabel)
           zlabel='JACOBIAN CONTOURS'
         call midplpt
     &      (real(r(1:nrt)),real(gsqrt(1:nrt)),
     &       nthpts,ns,xlabel,ylabel,zlabel,runlabel)
 8020    continue  !  PAGE # 2
C###    PAGE # 3
         page = '3'
         pagedesc = 'Boundary Errors'

         if(.not.lfreeb)goto 8021
         if (.not.makeplot(explist,lenlist,page,pagedesc)) goto 8021

 8021    continue  !PAGE # 3
      enddo !kz
      call newframe
*
*                 COMPUTE EQUILIBRIUM PROFILES
*
 
      phin(:ns) = abs(phi(:ns))
 
C###    PAGE # 4
      page = '4'

      pagedesc = 'Convergence, Magnetic Axis'
      if (makeplot(explist,lenlist,page,pagedesc)) then
      wdot(1) = wdot(2)
      fsqt(1) = fsqt(2)
      itfsq = min(itfsq, 100)
      do j = 1, itfsq
         time(j) = real(niter*(j - 1))/real(itfsq - 1)
      enddo

      noxc = 1
      ndata = 0
      xlabel= 'ITERATIONS'
      ylabel='\(2268)F\u2\d dV '
      zlabel='RESIDUAL FORCE'
      call loggraf1
     &(real(time(1:)),real(fsqt(1:)),
     &itfsq,xlabel,ylabel,zlabel,runlabel)

      xlabel= 'ITERATIONS'
      ylabel='-dW/dt'
      zlabel='ENERGY MINIMIZATION'
      call loggraf1
     &(real(time(1:)),real(wdot(1:)),
     &itfsq,xlabel,ylabel,zlabel,runlabel)


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
         enddo


         zlabel='MAGNETIC AXIS (R)'
         xlabel='N\dfp\u\gF/2\gp' 
         ylabel='R\dM\u(\gF) (m)'
         call graf1
     &      (real(time(1:nplot1)),real(raxis_v(1:nplot1)),
     &       nplot1,xlabel,ylabel,zlabel,runlabel)

         zlabel='MAGNETIC AXIS (Z)'
         xlabel='N\dfp\u\gF/2\gp' 
         ylabel='Z\dM\u(\gF) (m)'
         call graf1
     &    (real(time(1:nplot1)),real(zaxis_v(1:nplot1)),
     &    nplot1,xlabel,ylabel,zlabel,runlabel)

 8022 continue
      endif
         endif !makeplot  PAGE # 4

C###    PAGE # 5
         page = '5'
         pagedesc = 'Stability and Curvatures'
      if (makeplot(explist,lenlist,page,pagedesc)) then

      specw(1) = 2.*specw(2) - specw(3)
      phinorm = maxval(abs(phi(:ns)))
      phin(:ns) = abs(phi(:ns)/phinorm)


      zlabel='SPECTRAL WIDTH'
      ylabel= '<M>'
      xlabel='\gF'
      call graf1
     &  (real(phin(1:)),real(specw(1:)),
     &  ns,xlabel,ylabel,zlabel,runlabel)
 8023 continue


        zlabel='MERCIER'
        ylabel= 'ArcSinh D\dM\u'
        xlabel='\gF'
        dummy2=0
        dummy2(3:ns)=log( dmerc(3:ns) + sqrt( dmerc(3:ns)**2+1 ) )!arcsinh
        call graf1x
     &   (real(phin(3:)),real(dummy2(3:)),
     &   ns-4,xlabel,ylabel,zlabel,runlabel)

        zlabel='SHEAR and WELL'
        ylabel= 'D\dS\u'
        fdum='D\dW\u'
        xlabel='\gF'
        call graf2x
     &  (real(phin(3:)),real(dshear(3:)),real(dwell(3:)),
     &  ns-4,xlabel,ylabel,fdum,zlabel,runlabel)

        zlabel='CURRENT and GEODESIC'
        ylabel= 'D\dcur\u'
        fdum='D\dgeo\u'
        xlabel='\gF'
        call graf2x
     &   (real(phin(3:)),real(dcurr(3:)),real(dgeod(3:)),
     &   ns-4,xlabel,ylabel,fdum,zlabel,runlabel)

      endif  ! makeplot C###    PAGE # 5


C###    PAGE # 6
      page = '6'
      pagedesc = 'Fluxes, Iota, Well'
      if (makeplot(explist,lenlist,page,pagedesc)) then
 
         do j = 1, ns
            br(j) = hs*(j - 1)
         enddo
         ndata = 0
         zlabel='TOROIDAL FLUX'
         ylabel= '\gF'
         xlabel='s'
         call graf1
     &   (real(br(1:)),real(phin(1:)),
     &   ns,xlabel,ylabel,zlabel,runlabel)
         t1 = 0.5*twopi*hs
         do j = 2, ns
            bz(j) = iotas(j)
            br(j) = br(j-1) + t1*phip(j)*(iotas(j)+iotas(j-1))
         enddo
         zlabel='POLOIDAL FLUX'
         ylabel= '\gx'
         xlabel='\gF'
         call graf1
     &    (real(phin(1:)),real(br(1:)),
     &    ns,xlabel,ylabel,zlabel,runlabel)
         a0 = 1.5*bz(2) - .5*bz(3)
         an = 1.5*bz(ns) - .5*bz(ns-1)
CDIR$   IVDEP
         bz(2:ns-1) = .5*(bz(2:ns-1)+bz(3:ns))
         br(2:ns-1) = 1./bz(2:ns-1)
         bz(1) = a0
         if (a0 .ne. 0.) br(1) = 1./a0
         bz(ns) = an
         br(ns) = 1./an
         br(:ns) = abs(br(:ns))
         bz(:ns) = abs(bz(:ns))
         zlabel='IOTA'
         ylabel= '\gi'
         fdum='Safety Factor'
         xlabel='\gF'
         call graf2x
     &     (real(phin(1:)),real(bz(1:)),real(br(1:)),
     &     ns,xlabel,ylabel,fdum,zlabel,runlabel)
         bz(:ns) = (ub(:ns)-ub(1))/ub(1)
         zlabel='MAGNETIC WELL ?? AND INVERSE'
         write(nchar,fmt='(1h'')')
         ylabel= '[V'//nchar(1:1)//nchar(1:1)//'(\gF)- V'
         ylabel=trim(ylabel)//nchar(1:1)//nchar(1:1)//'(0)]'
         ylabel=trim(ylabel)//' \(2770) V'//nchar(1:1)//'(0)'
         fdum='what is this?'
         xlabel='\gF'
         call graf2x
     &     (real(phin(1:)),real(br(1:)),real(bz(1:)),
     &     ns,xlabel,ylabel,fdum,zlabel,runlabel)
 
         endif !makeplot  PAGE # 6


C###    PAGE # 7
      page = '7'
      pagedesc = 'Mass, Pressure, Beta, Pol. Cur vs phi'
      if (makeplot(explist,lenlist,page,pagedesc)) then
C      print *,'VolAvgB, RBtor0, RBtor, Aminor, Rmajor, Volume,itor'
C      print *,VolAvgB, RBtor0, RBtor, Aminor, Rmajor, Volume,itor

            ndata = 0

         zlabel='MASS'
         ylabel= 'M(\gF)'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(mass(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('R\(0210)B\dt\u/R=',f5.2,1hT,' A=',f5.2)"
         write(ylabel,fmt=xlabel)RBtor/ Rmajor,  Rmajor/ Aminor
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))
 

         zlabel='PRESSURE'
         ylabel= 'P(\gF)'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(pres(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel=
     &  "('R\(0210)B\d1\u=',f5.2,x,'R\(0210)B\d0\u=',f5.2,4h m-T)"
         write(ylabel,fmt=xlabel)RBtor, RBtor0
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))

         zlabel='BETA'
         ylabel= '<\gb>'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(beta_vol(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('<\gb>=',f5.2,1h%,x,'<B>=',f5.2,1hT)"
         write(ylabel,fmt=xlabel)wp/wb*100.,volavgb
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))

         zlabel='POLOIDAL CURRENT DENSITY'
         ylabel= '<J\.\(2266)\gh>'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(jcuru(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
         xlabel="('<R>=',f5.2,1hm,x,'V=',f5.2,'m\u3\d')"
         write(ylabel,fmt=xlabel)Rmajor,Volume
         call pgmtxt('T',0.9,0.5,0.5,trim(ylabel))
 
         endif !makeplot  PAGE # 7

 

C###    PAGE # 8
      page = '8'
      pagedesc = 'Jtor , J|| & Itor vs phi'
      if (makeplot(explist,lenlist,page,pagedesc)) then
 
         ndata = 0
         zlabel='FLUX-AV. JTOR (A/m\u2\d)'
         ylabel= '<J\.\(2266)\(0647)>\(2770)<R\u-1\d>'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(jcurv(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)


         zlabel='<J\.B>\(2770)<B\.\(2266)\(0647)>'
         ylabel= 'Amps/m'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(jdotb(1:)),
     &      ns,xlabel,ylabel,zlabel,runlabel)


         zlabel='TOROIDAL CURRENT'
         ylabel= 'I\dTOR\u'
         xlabel='\gF'
         call graf1
     &     (real(phin(1:)),real(itors(1:)),
     &     ns,xlabel,ylabel,zlabel,runlabel)
 
      endif !makeplot  PAGE # 8
      call newframe
 
C###    PAGE # 9
      page = '9'
      pagedesc = 'Local Force Balance'
      if (makeplot(explist,lenlist,page,pagedesc)) goto 3001 ! outdated format
      if (makeplot(explist,lenlist,page,pagedesc)) then
         nzeta=ntor*2+4           ! guess sizes and then check
         ntheta=mpol*2+6
         ntheta2=1+(2*(ntheta/2))/2
         nznt=(1+ntheta2/2)*(nzeta/2)
         nx=ns/2; ny=nzeta*ntheta2/4 +4
         ntht=ntheta2/2+1
         nsurf=ns/2
         nphi=nzeta/2  ! assume old defaults!
         jxbcode=2
         fname='jxbout.'//trim(input_id)
         call read_jxbout(lka, jsupu, jsupv, jsups, 
     1    bsupu, bsupv, jcrossb, pprime, residual12, 
     2    bdotj, bsubu, bsubv, bsubs, 
     3    hshere, avforce, aminfor, amaxfor,
     4    ns, nznt, mpol, ntor, trim(fname),
     5    nsurf , ntht , nphi , jxbcode)
         if (jxbcode.eq.1) then
          nznt1=ntht*nphi
          allocate(bsupu(nsurf,nznt1),bsupv(nsurf,nznt1), 
     1     bsubu(nsurf,nznt1), bsubv(nsurf,nznt1), bsubs(nsurf,nznt1), 
     2     residual12(nsurf,nznt1), bdotj(nsurf,nznt1), 
     3     jsupu(nsurf,nznt1), lka(nsurf,nznt1), jcrossb(nsurf,nznt1),
     4     jsupv(nsurf,nznt1), xfres12(nsurf,nznt1),
     5      pprime(nsurf,nznt1), jsups(nsurf,nznt1))
         allocate(avforce(nsurf), hshere(nsurf),iplot(nsurf), 
     1     amaxfor(nsurf), aminfor(nsurf),xplot(nsurf),yplot(nsurf))
         nx=nsurf; ny=nznt1
         fname='jxbout.'//trim(input_id)
         call read_jxbout(lka, jsupu, jsupv, jsups, 
     1    bsupu, bsupv, jcrossb, pprime, residual12, 
     2    bdotj, bsubu, bsubv, bsubs, 
     3    hshere, avforce, aminfor, amaxfor,
     4    ns, nznt, mpol, ntor, trim(fname),
     5    nsurf , ntht , nphi , jxbcode)
         aminfor=avforce*aminfor/100
         amaxfor=avforce*amaxfor/100

        elseif(jxbcode.eq.3)then
         nsurf=ns/2; nznt1=(1+ntheta2/2)*(nzeta/2)
          allocate(bsupu(nsurf,nznt1),bsupv(nsurf,nznt1), 
     1     bsubu(nsurf,nznt1), bsubv(nsurf,nznt1), bsubs(nsurf,nznt1), 
     2     residual12(nsurf,nznt1), bdotj(nsurf,nznt1), 
     3     jsupu(nsurf,nznt1), lka(nsurf,nznt1), jcrossb(nsurf,nznt1),
     4     jsupv(nsurf,nznt1), xfres12(nsurf,nznt1),
     5      pprime(nsurf,nznt1), jsups(nsurf,nznt1))
         allocate(avforce(nsurf), hshere(nsurf),iplot(nsurf), 
     1     amaxfor(nsurf), aminfor(nsurf),xplot(nsurf),yplot(nsurf))
         fname='jxbout.'//trim(input_id)
         call read_old_jxbout(lka, jsupu, jsupv, jsups, 
     1    bsupu, bsupv, jcrossb, pprime, residual12, 
     2    bdotj, bsubu, bsubv, bsubs, 
     3    hshere, avforce, aminfor, amaxfor,
     4    ns, nznt, ntheta2, nzeta, trim(fname),
     5    nsurf , ntht , nphi , jxbcode)
         nx=nsurf; ny=nznt1
        endif
        if (jxbcode.ne.0)goto 3001
         nzeta=ntor*2+4           ! retuned from read_jxbout
         ntheta=mpol*2+6
         ntheta2=1+(2*(ntheta/2))/2
         nznt=(1+ntheta2/2)*(nzeta/2) ! want this many and have nznt1
        kt=1+(nznt1-1)/((1+ntheta2/2)*(nzeta/2))
C        kt=1+(nznt1-1)/(ntht*nphi)
        kz=max(1+2*(kt/2),5)
        allocate(pang(ntheta2+1),tang(nzeta))
        do i=1,nznt1       
           itheta = 1 + lka(nx/2,i)/nzeta
           izeta = 1 + mod(lka(nx/2,i)-1,nzeta)
           pang(itheta)=pi*float((itheta-1))/ntheta2
           tang(izeta)=2.*pi*(izeta-1)/nzeta/nfp
        enddo
         zlabel='<j\xB Force>'
         ylabel= '\(2229)JxB-p''\(2229)'
         xlabel='s'
         call loggraf3pt
     &   (real(hshere(1:)),
     &    real(abs(aminfor(1:))),
     &    real(abs(avforce(1:))),
     &    real(abs(amaxfor(1:))),
     &    nx,xlabel,'\(2229)\fiMIN\fn\(2229)',ylabel,
     &    '\(2229)\fiMAX\fn\(2229)',zlabel,runlabel)



         zlabel='<j\xB Force>'
         ylabel= 'J\xB-p'''
         xlabel='s'
         call graf1pt
     &     (real(hshere(1:)),real(avforce(1:)),
     &     size(aminfor(1:)),xlabel,ylabel,zlabel,runlabel)
      print *,'hshere ',size(hshere), maxval(hshere), minval(hshere)
      print *,'avforce ',size(avforce), maxval(avforce), minval(avforce)

         zlabel='Min j\xB Force'
         ylabel= 'J\xB-p'''
         xlabel='s'
         call graf1pt
     &     (real(hshere(1:)),real(aminfor(1:)),
     &     size(aminfor(1:)),xlabel,ylabel,zlabel,runlabel)

         zlabel='Max j\xB Force'
         ylabel= 'J\xB-p'''
         xlabel='s'
         call graf1pt
     &     (real(hshere(1:)),real(amaxfor(1:)),
     &     size(aminfor(1:)),xlabel,ylabel,zlabel,runlabel)
CASCADE
        call pgsave
        call pgsubp(1,1)
        call pgpap(9.,8.5/11.)
        xfres12=abs(residual12)
        xmin=0
        xmax=1
        ymin=minval(xfres12)
        ymax=maxval(xfres12)
        xspread=1.
        yspread=maxval(xfres12)-minval(xfres12)
        xxspread=xspread/ntheta2/kt
        yyspread=yspread/ntheta2/kt
        xxspread=xspread/ntht/kt
        yyspread=yspread/ntht/kt
        deltx=real(hshere(2)-hshere(3))
        nspread=nint(xxspread/deltx)
        xxspread=1.*nspread*deltx
        maxc=nznt1
        ymmin=ymin
        ymmax=ymax+maxc*yyspread
        xmmin=0
        xmin=0
        xmax=1
        xmmax=1.6 + maxc*xxspread
        call pgsci(1)
        call pgpage
        call pgsvp(.1,.9,.1,.9) ! set original window size
        call pgswin(xmmin,xmmax,ymmin,ymmax) ! set original data limits
        call pgbox('bc',0.,0,'bc',0.,0) !display world box
        call xywtoxyan(xmin,ymin,x1,y1)  ! make coordinate labels
        call xywtoxyan(xmax,ymax,x2,y2)
        call pgsvp(x1,x2,y1,y2)  ! little box for world units
        call pgswin(xmin,xmax,ymin,ymax) ! limits for world units
        call pgsch(.79)
        call pgbox('BINT',.25,4,'BVINT',1.,2) !box for world units
        call pgsvp(.1,.9,.1,.9) ! reset to original size
        call pgswin(xmmin,xmmax,ymmin,ymmax) ! set world units for full box
        ylabel='Residual'
        zlabel= '2\.|JxB-p''|'
        zlabel=trim(zlabel)//' / [|JxB|+|p''|]'
        xlabel='s'
        call pgqch(siz)
        call pgsch(siz-.1)
        call pgmtxt('B',2.25,0.5,0.5,trim(runlabel))
        call pgsch(siz)
        call pglab(xlabel,ylabel,zlabel)
        write(pchar,fmt='(a,1pe9.2)')'Plane is at ',ymin
        call pgtext(xmmin+2*xxspread,ymmax-2*kt*yyspread,pchar)
        write(pchar,fmt='(a,1pe9.2)')'Maximum Residual is  ',ymax
        call pgtext(xmmin+2*xxspread,ymmax-4*kt*yyspread,pchar)
        write(pchar,fmt='(a)')'(\gh/\gp,\gz/\gp/N\dFP\u) pairs'
        call pgtext(xmmin+2*xxspread,ymmax-8*kt*yyspread,pchar)
        call pgsci(15)
        call pgsch(.68)
        xplot(1)=0.
        yplot(1)=ymin
        xplot(2)=xxspread*maxc
        yplot(2)=maxc*yyspread                         
        xplot(3)=1.+xplot(2)
        yplot(3)=yplot(2)
        xplot(4)=1
        yplot(4)=yplot(1)
        call pgpoly(4,xplot,yplot)  ! draw reference plan, full box
        call pgsci(2)
        call pgsci(1)
        ymmin=ymin
        ymmax=ymax
        do 2500 i=1,nznt1
         itheta = (1 + lka(nx,i)/nzeta)
         izeta = 1 + mod(lka(nx,i)-1,nzeta)
           write(xlabel,fmt='(a,f4.2,a1,f4.2,a1)')
     &    '(',
     &      1.*(itheta-1)/ntheta2,',',2.*(izeta-1)/nzeta/nfp,')'
             xlabel=' ..'//trim(xlabel)
         ic = i - ((i-1) / maxc) * maxc
         yzero = (maxc-i)*yyspread
         xplot(1:nx)=real(hshere(1:nx))+ (maxc-ic)*xxspread
         yplot(1:nx)=(real(xfres12(1:nx,i)+ (maxc-ic)*yyspread))
         call pgsci(1+modulo(i,kz))
         call pgline(nx,xplot,yplot)
         if(1+modulo(i,1+2*(kt/2)).eq.1)
     &     call pgtext(xplot(nx),yzero,trim(xlabel))
2500     continue
        xfres12=abs(residual12)
        ymin=minval(xfres12)
        ymax=maxval(xfres12)
        if(ymin.eq.0.)xfres12=xfres12+1.e-6
        xfres12=log10(xfres12)
        xspread=1.
        yspread=maxval(xfres12)-minval(xfres12)
        xxspread=xspread/ntheta2/kt
        yyspread=yspread/ntheta2/kt
        deltx=real(hshere(2)-hshere(3))
        nspread=nint(xxspread/deltx)
        xxspread=1.*nspread*deltx
        ymax=maxval(xfres12)
        ymin=minval(xfres12)
        ymin=ymax-3  ! set # decades
        maxc=nznt1

        ymmin=ymin
        ymmax=ymax+maxc*yyspread
        xmmin=0
        xmin=0
        xmax=1
        xmmax=1.6 + maxc*xxspread
        call pgsci(1)
        call pgpage
        call pgsvp(.1,.9,.1,.9) ! set original window size
        call pgswin(xmmin,xmmax,ymmin,ymmax) ! set original data limits
        call pgbox('bc',0.,0,'bc',0.,0) !display world box
        call  xywtoxyan(xmin,ymin,x1,y1)  ! make coordinate labels
        call  xywtoxyan(xmax,ymax,x2,y2)
        call pgsvp(x1,x2,y1,y2)  ! little box for world units
        call pgswin(xmin,xmax,ymin,ymax) ! limits for world units
        call pgsch(.79)
        call pgbox('BINT',.25,4,'BVSILNT',0.5,1) !box for world units
        call pgsvp(.1,.9,.1,.9) ! reset to original size
        call pgswin(xmmin,xmmax,ymmin,ymmax) ! set world units for full box
        ylabel='Log\d10\u Residual'
        zlabel= '2\.|JxB-p''|'
        zlabel=trim(zlabel)//' / [|JxB|+|p''|]'
        xlabel='s'
        call pgqch(siz)
        call pgsch(siz-.1)
        call pgmtxt('B',2.25,0.5,0.5,trim(runlabel))
        call pgsch(siz)
        call pglab(xlabel,ylabel,zlabel)
        write(pchar,fmt='(a,1pe9.2)')'Plane is at ',10**ymin
        call pgtext(xmmin+2*xxspread,ymmax-2*kt*yyspread,pchar)
        write(pchar,fmt='(a,1pe9.2)')'Maximum Residual is  ',10**ymax
        call pgtext(xmmin+2*xxspread,ymmax-4*kt*yyspread,pchar)
        write(pchar,fmt='(a)')'(\gh/\gp,\gz/\gp/N\dFP\u) pairs'
        call pgtext(xmmin+2*xxspread,ymmax-8*kt*yyspread,pchar)
        call pgsci(15)
        call pgsch(.68)
        xplot(1)=0.
        yplot(1)=ymin
        xplot(2)=xxspread*maxc
        yplot(2)=ymin+maxc*yyspread                         
        xplot(3)=1.+xplot(2)
        yplot(3)=yplot(2)
        xplot(4)=1
        yplot(4)=yplot(1)
        call pgpoly(4,xplot,yplot)  ! draw reference plan, full box
        call pgsci(2)
        call pgsci(1)
        ymmin=ymin
        ymmax=ymax
        do 3000 i=1,nznt1
         itheta = (1 + lka(nx,i)/nzeta)
         izeta = 1 + mod(lka(nx,i)-1,nzeta)
           write(xlabel,fmt='(a,f4.2,a1,f4.2,a1)')
     &    '(',
     &      1.*(itheta-1)/ntheta2,',',2.*(izeta-1)/nzeta/nfp,')'
             xlabel=' ..'//trim(xlabel)
         ic = i - ((i-1) / maxc) * maxc
         yzero = (maxc-i)*yyspread
         xplot(1:nx)=real(hshere(1:))+ (maxc-ic)*xxspread
         yplot(1:nx)=(real(xfres12(1:nx,i)+ (maxc-ic)*yyspread))
         call pgsci(1+modulo(i,kz))
         ymin=ymmin+(maxc-ic)*yyspread
         ymax=ymmax+(maxc-ic)*yyspread
         iplot=.false.
         where(yplot(1:nx).ge.ymin)
          iplot=.true.
         endwhere

         k=0
         l=k+1
         alldone=.false.
         segments:do
         alldone=.true.
         addon=k+1.lt.nx
         l=k+1
         do j=k+1,nx
            if(addon.and.iplot(j)) then
               addon=.false.
               l=j
            endif
         enddo
         addon=.true.
         do j=l+1,nx
            if(addon.and.iplot(j)) then
               k=j
            else
               addon=.false.
            endif
         enddo
         if(k.gt.l)call pgline(k-l+1,xplot(l:k),yplot(l:k))
         if(alldone.and.l.lt.nx)alldone=.false.
         if(alldone)exit segments
         if(iplot(l) .and. (.not.iplot(l+1)))then
           call pgpt1(xplot(l),yplot(l),-2)
           k=k+1
         endif
         l=k
         enddo segments
         if(1+modulo(i,1+2*(kt/2)).eq.1)
     &     call pgtext(xplot(nx),ymin,trim(xlabel))
3000   continue

CASCADE
         call pgunsa
         call pgsubp(2,2)
         call pgpap(8.5,1.)
         call pgsch(1.85)

      endif          !makeplot  PAGE # 9
      call newframe

 3001 continue

      if(allocated(  tang))deallocate(   tang)       
      if(allocated(  pang))deallocate(   pang)
      if(allocated(bsubu))deallocate(bsubu)
      if(allocated(bsubu))deallocate(bsubu)
      if(allocated(bsubv))deallocate(bsubv)
      if(allocated(bsubs))deallocate(bsubs)
      if(allocated(residual12))deallocate(residual12)
      if(allocated(bdotj))deallocate(bdotj)
      if(allocated(jsupu))deallocate(jsupu)
      if(allocated(lka))deallocate(lka)
      if(allocated(jcrossb))deallocate(jcrossb)
      if(allocated(jsupv))deallocate(jsupv)
      if(allocated(xfres12))deallocate(xfres12)
      if(allocated(pprime))deallocate(pprime)
      if(allocated(jsups))deallocate(jsups)
      if(allocated(avforce))deallocate(avforce)
      if(allocated(hshere))deallocate(hshere)
      if(allocated(iplot))deallocate(iplot)
      if(allocated(amaxfor))deallocate(amaxfor)
      if(allocated(aminfor))deallocate(aminfor)
      if(allocated(xplot))deallocate(xplot)
      if(allocated(yplot))deallocate(yplot)

C###    PAGE # 10
      page = '10'
      pagedesc = '|B| Spectrum'
      if (makeplot(explist,lenlist,page,pagedesc)) then
      iunit=110
      if(ns==0)ns=49
      call safe_open (iunit, istat, 'in_booz.' // trim(input_id), 
     1     'unknown','formatted')
      mboz=17
      nboz=9
      write(iunit,fmt='(2i4.3)')mboz,nboz
      write(iunit,fmt='(a)')trim(input_id)
      write(iunit,fmt='(4(i3.3,x))')2,ns/3,2*ns/3,ns-1
      close(iunit)
      fname='$HOME/bin/xbooz_xform in_booz.' // trim(input_id)//' F'
c      fname=trim(fname)//'; \rm in_booz.' // trim(input_id)
      call system(trim(fname),istat)
      if(istat/=0)stop 'xbooz_xform'
      iunit=120
      call safe_open (iunit, istat, 'boozmn.' // trim(input_id), 'old',
     1     'unformatted')
      if (istat .ne. 0) then
         write(6,*)' istat = ', istat
         stop 'Error opening boozmn file in XBOOZ_XFORM!'
      endif   
      read(iunit, iostat=ierr, err=1000) 
     1   nfpb, nsb, aspectb, rmax_surfb, rmin_surfb, betaxisb

      js=nsb;allocate(hiota(js), presb(js), jlistb(ns),
     1   beta_volb(js), phipb(js),  phib(js), bvcob(js),  bucob(js))

      do js = 2, nsb
         read(iunit, iostat=ierr, err=1000) hiota(js), presb(js), 
     1   beta_volb(js), phipb(js),  phib(js), bvcob(js),  bucob(js)
      enddo
      read(iunit, iostat=istat, err=1000) mboz, nboz, mnboz
      read(iunit, iostat=istat, err=1000) versionb
      allocate ( xmb(mnboz) , xnb(mnboz), bmnb(ns,mnboz), 
     &  rmnb(ns,mnboz), zmnb(ns,mnboz), pmnb(ns,mnboz))
      iread=0
999   read(iunit, iostat=istat, err=1000, end=1010) jrad
      iread=iread+1 ! # surfaces in data file at statement 1010
      jlistb(iread)=jrad
      if(first) then
         first = .false.
         do mn = 1,mnboz
            read(iunit, iostat=istat) nb, mb
            xmb(mn)=one * mb
            xnb(mn)=one * nb
            nb=nb/nfp
         enddo 
       endif
       do mn = 1, mnboz
         read(iunit, iostat=istat) bmnb(jrad,mn), rmnb(jrad,mn), 
     &       zmnb(jrad,mn), pmnb(jrad,mn)
         if(xnb(mn).eq.0..and.abs(xmb(mn)).lt.7.) write(6,
     &  fmt="('n=',i3,' m=',i3,' B=',1pe9.2,' s=',1pe9.2)")
     &  nint(xnb(mn)),nint(xmb(mn)),bmnb(jrad,mn), phi(jrad)/phi(ns)
       enddo
       goto 999
 1010 continue
      call pgpap(8.,1.)
      call pgsubp(2,2)
      do js=1,4
        write(pchar,fmt='(2(a,1pe9.2))')'\gi = ',hiota(jlistb(js)),
     &   ' s = ',phi(jlistb(js))/phi(ns)
        write(mchar,fmt="('B\d0,0\u=',1pe10.3)")bmnb(jlistb(js),1)
       call bdots(js,nfp,jlistb(js),mnboz,real(xnb(1:)),real(xmb(1:)),
     & real(bmnb(jlistb(js),1:)),'n','m',pchar,input_id,mchar)
      enddo
      goto 1111
 1000 continue
      write(6,*) 'reading error'
 1111 continue
      
      endif !          endif !makeplot  PAGE # 10

C###    PAGE # 16    ! differs considerably from the original
      page = '16'
      pagedesc = 'R,Z,Lamba Profiles'
      if (makeplot(explist,lenlist,page,pagedesc)) then
         dummy1(1:mnmax,1)=((/(i,i=1,mnmax)/))
         do i=1,mnmax
           dummy2(i)=abs(maxval(rmnc(i,1:ns))-minval(rmnc(i,1:ns)))
         enddo
         xlabel= 'mode #'
         ylabel='Max[R\dmn\uC\dmode\u]'
         zlabel='Maximum Amplitude'

         call loggrafpt
     &    (real(dummy1(1:mnmax,1)),real(dummy2(1:mnmax)),
     &    mnmax,xlabel,ylabel,zlabel,runlabel)

         dummy1(1:mnmax-1,1)=((/(i,i=2,mnmax)/))
         do i=2,mnmax
           j=i-1
           dummy2(j)=abs(maxval(zmns(i,1:ns))-minval(zmns(i,1:ns)))
         enddo
         if(minval(dummy2(1:mnmax-1)).eq.0.)then
          small=maxval(dummy2(1:mnmax-1))/1.e8
          dummy2(1:mnmax-1)=dummy2(1:mnmax-1)+small
         endif
         xlabel= 'mode #'
         ylabel='Max[Z\dmn\uS\dmode\u]'
         zlabel='Maximum Amplitude'
         call loggrafpt
     &    (real(dummy1(1:mnmax-1,1)),real(dummy2(1:mnmax-1)),
     &    mnmax-1,xlabel,ylabel,zlabel,runlabel)
         if (iasym .eq. 1) then
           dummy1(1:mnmax,1)=((/(i,i=1,mnmax)/))
           do i=1,mnmax
             dummy2(i)=abs(maxval(rmnc(i,1:ns))-minval(rmnc(i,1:ns)))
           enddo
           xlabel= 'mode #'
           ylabel='Max[R\dmn\uC]'
           zlabel='Maximum Amplitude'

           call loggrafpt
     &      (real(dummy1(1:mnmax,1)),real(dummy2(1:mnmax)),
     &       mnmax,xlabel,ylabel,zlabel,runlabel)

           dummy1(1:mnmax,1)=((/(i,i=1,mnmax)/))
           do i=1,mnmax
             dummy2(i)=abs(maxval(zmns(i,1:ns))-minval(zmns(i,1:ns)))
           enddo
           if(minval(dummy2(1:mnmax)).eq.0.)then
             small=maxval(dummy2(1:mnmax))/1.e8
             dummy2(1:mnmax)=dummy2(1:mnmax)+small
           endif
           xlabel= 'mode #'
           ylabel='Max[Z\dmn\uS]'
           zlabel='Maximum Amplitude'
           call loggrafpt
     &       (real(dummy1(1:mnmax,1)),real(dummy2(1:mnmax)),
     &       mnmax,xlabel,ylabel,zlabel,runlabel)
         endif
      endif  !makeplot  PAGE # 16

      if(allocated( dummy1))deallocate( dummy1)
      if(allocated( dummy2))deallocate( dummy2)
      if(allocated(  r))deallocate(  r)
      if(allocated( z))deallocate( z)
      if(allocated( ru))deallocate( ru)
      if(allocated( modb))deallocate( modb)
      if(allocated( sqrt_phim))deallocate( sqrt_phim)
      if(allocated( phim))deallocate( phim)
      if(allocated( gsqrt))deallocate( gsqrt)
      if(allocated( torcur))deallocate( torcur)
      if(allocated( dbsubudstz))deallocate( dbsubudstz)
      if(allocated( bsubutz))deallocate( bsubutz)
      if(allocated( bsubvtz))deallocate( bsubvtz)
      if(allocated(  r12))deallocate(  r12)
      if(allocated( z12))deallocate( z12)
      if(allocated( ru12))deallocate( ru12)
      if(allocated( zu12))deallocate( zu12)

 8033 return

  150 format(a)
      end subroutine pgpplotter
      subroutine bdots
     .  (ipos,nfp,js,mnboz,xnb,xmb,bmnb,lx,ly,lt,runlbl,l2)
      character*(*) lx,ly,lt,runlbl,l2
      character*8 xopt,yopt
      real xnb(*),xmb(*),bmnb(*),r1,r2
      data r1/6./
      integer js,mnboz,nfp,ipos
      if(ipos==1)then
         call pgpap(8.,1.)
      endif
      xnb(1:mnboz)=xnb(1:mnboz)/nfp
      xmin=minval(xnb(1:mnboz))
      xmax=maxval(xnb(1:mnboz))
      ymin=minval(xmb(1:mnboz))
      ymax=maxval(xmb(1:mnboz))
      bmax=maxval(bmnb(1:mnboz))
      bmin=minval(bmnb(1:mnboz))
      ymin=0
      ymax=(xmax-xmin)+1.
         call pgbbuf
         call pgsch(1.5)
         call pgsci(1)
         call pgpage
         call pgvstd
c         call pgswin(xmin,xmax,ymin,ymax)
         call pgwnad(xmin,xmax,ymin,ymax)
         xopt='BCNST';yopt='BCNST'
         call pgbox(trim(xopt),0.,0,trim(yopt),0.,0)
         call pglab(trim(lx),trim(ly),trim(lt))
         call pgmtxt('T',0.9,0.5,0.5,trim(l2))
         call pgqch(siz)
         call pgsch(siz-.5)
         call pgmtxt('B',2.85,0.5,0.5,trim(runlbl))
         call pgsch(siz)
         do mn=1,mnboz
         if(bmnb(mn) < 0.) then
           call pgsci(4)
           call pgsfs(1)
         else
           call pgsci(2)
           call pgsfs(1)
         endif
         r=(r1+log10(abs(bmnb(mn))/bmax))/r1/2
c         if(r > 0)call pgcirc(xnb(mn),xmb(mn),r)
         if(r > 0)call 
     .   pgrect(xnb(mn)-r/2,xnb(mn)+r/2,xmb(mn)-r/2,xmb(mn)+r/2)
         enddo
         call pgebuf
         if(ipos/=4)return
         call xywtoxyan(0.,ymin,x1,y1)
         call xywtoxyan(1.,ymax,x2,y2)
         call xywtoxyan(xmax,ymin,x3,y3)
         call pgsvp(.91,.91+x2-x1,y1,y2)
         call pgswin(0.,1.,ymin,ymax)
         do i=1,7
          r=(r1+log10(1./10.**(i-1)))/r1/2
c          call pgcirc(.5,ymax-(ymax-ymin)*i/8,r)
          call pgrect(
     .    .5-r/2,.5+r/2,
     .    ymax-(ymax-ymin)*i/8-r/2,ymax-(ymax-ymin)*i/8+r/2)
         enddo
         call pgsci(1)
         call pgsch(1.8)
         call pgptxt(1.,ymin+(ymax-ymin)/4,270.,1.,
     /    'decade scale (red is +)')
!         end subroutine bdots
         end

      subroutine putstb
c
c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!




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
      real(kind=rprec) :: pad
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
      enddo
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

      subroutine read_jxbout(lka, jsupu, jsupv, jsups,
     1 bsupu1, bsupv1, jcrossb, pprime, r12, 
     2 bdotj, bsubua, bsubva, bsubs, 
     3 hs, avforce, aminfor, amaxfor,
     4 ns, nznt, ntheta2, nzeta, jxbout_file,
     5    nsurf , ntht , nphi , jxbcode)
C-----------------------------------------------
C jxbcode:
C -1 -- file not found
C -2 -- nsurf > ns
C -2 -- ntht  > ntheta2
C -2 -- nphi  > nzeta
C +3 -- may be old style file; try read_old_jxbout
C +2 -- reading for table size set to 1 on output
C +1 -- reading table of values
C  0 -- successful reading of data

C-----------------------------------------------
      use kind_spec
      use Vpname1, only: dmu0
      use safe_open_mod

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer  ns, nznt , ntheta2, nzeta
      integer  nsurf , ntht , nphi , jxbcode
      real(kind=rprec), dimension(nsurf,*) :: 
     1 jsupu, jsupv, jsups, 
     2 bsupu1, bsupv1, jcrossb, pprime, r12, 
     3 bdotj, bsubua, bsubva, bsubs
      real(kind=rprec), dimension(*) :: 
     1 hs, avforce, aminfor, amaxfor
      integer, dimension(nsurf,*) :: lka
      character*(*) jxbout_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: p5 = 0.5d0, two = 2.0d0, c1p5=1.5d0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer lk, lz, lt, k, m, js, j, n, injxbout
      integer :: njxbout, n1, n2
      character*30 :: filename, dummytxt*72
      real(kind=rprec) dum1, dum2 
      real(kind=rprec), dimension(:,:), allocatable ::
     1 bsubuv, bsubvu
      logical newfiletype, exists
C-----------------------------------------------
      if(jxbcode.eq.1)goto 1001
      call get_lun(njxbout)
      call safe_open(njxbout, injxbout, trim(jxbout_file), 'old',
     1     'formatted')
      if (injxbout .ne. 0) then
         write(6,*)' Error opening JXBOUT file in read_jxbout'
         jxbcode=-1
         call free_lun(njxbout)
         return
      endif   
  7    format(/,19x,i3.3,24x,i3.3,24x,i3.3,/,19x,i3.3,19x,i3.3)
      read (njxbout,7,err=1002) nsurf, ntht, nphi, ntheta2, nzeta
      goto 1003
1002    jxbcode=3
        close(njxbout,iostat=injxbout)
        call free_lun(njxbout)
        return
1003    jxbcode=1
        close(njxbout,iostat=injxbout)
        call free_lun(njxbout)
        return  !  for allocation
1001  continue
      allocate(bsubuv(ns,nznt),bsubvu(ns,nznt))
      injxbout=0
      call get_lun(njxbout)
      call safe_open(njxbout, injxbout, trim(jxbout_file), 'old',
     1     'formatted')
      if(injxbout/=0)stop 'jxbreading error'
      do j=1,15
       read (njxbout,fmt='(a)')dummytxt
      enddo

!
!     NOTE: bsubua, bsubva are the ALIASED (truncated Fourier spectrum)
!           versions of bsubu, bsubv. bsubuv, bsubvu are used to compute
!           radial current (should be zero) 
  200 format(/17x,1pe12.3,18x,1p
     1   e12.3,15x,1pe12.3,21x,1pe12.3,/,
     2   44x,sp,0pf6.1,4x,f6.1)
  110 format(i5,12(x,en10.3e2))
      do js = 1, nsurf
            read (njxbout, 200) hs(js), avforce(js)
     1         ,dum1, dum2, amaxfor(js), aminfor(js)
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt

 
            lk=0
            do lt = 1, ntht
               do lz = 1, nphi
                  lk = lk + 1
                  read (njxbout, 110) lka(js,lk), jsupu(js,lk),
     1              jsupv(js,lk), dum1, bsupu1(js,lk),  
     2              bsupv1(js,lk),jcrossb(js,lk), pprime(js,lk),
     3              dum2, bdotj(js,lk), bsubua(js,lk),
     4              bsubva(js,lk), bsubs(js,lk)
                  jsups(js,lk)=dmu0*dum1
                  r12(js,lk)=2._rprec*dum2/  !  nomalized residual
     1             (abs(jcrossb(js,lk))+abs( pprime(js,lk)))
               enddo
            enddo
      enddo
 
      close (njxbout)
      call free_lun(njxbout)
      deallocate(bsubuv,bsubvu)
      jxbcode=0
      end subroutine read_jxbout

      subroutine read_old_jxbout(lka, jsupu, jsupv, jsups,
     1 bsupu1, bsupv1, jcrossb, pprime, r12, 
     2 bdotj, bsubua, bsubva, bsubs, 
     3 hs, avforce, aminfor, amaxfor,
     4 ns, nznt, ntheta2, nzeta, jxbout_file,
     5    nsurf , ntht , nphi , jxbcode)

      use kind_spec
      use Vpname1, only: dmu0
      use safe_open_mod

      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer  ns, nznt , ntheta2, nzeta
      integer nsurf , ntht , nphi , jxbcode
      real(kind=rprec), dimension(nsurf,*) :: 
     1 jsupu, jsupv, jsups, 
     2 bsupu1, bsupv1, jcrossb, pprime, r12, 
     3 bdotj, bsubua, bsubva, bsubs
      real(kind=rprec), dimension(*) :: 
     1 hs, avforce, aminfor, amaxfor
      integer, dimension(nsurf,*) :: lka
      character*(*) jxbout_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: p5 = 0.5d0, two = 2.0d0, c1p5=1.5d0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer lk, lz, lt, k, m, js, j, n, injxbout
      integer :: njxbout 
      character*30 :: dummytxt, filename
      real(kind=rprec) dum1, dum2 
      real(kind=rprec), dimension(:,:), allocatable ::
     1 bsubuv, bsubvu
      logical lxt,opnd
C-----------------------------------------------
      call get_lun(njxbout)
      call safe_open(njxbout, injxbout, trim(jxbout_file), 'old',
     1     'formatted')
      if (injxbout .ne. 0) then
         write(6,*)' Error opening JXBOUT file in read_jxbout'
         jxbcode=-1
         return
      endif   
      allocate(bsubuv(ns,nznt),bsubvu(ns,nznt))

       read (njxbout,fmt='(a10)')dummytxt
       read (njxbout,fmt='(a10)')dummytxt
       read (njxbout,fmt='(a10)')dummytxt
       read (njxbout,fmt='(a10)')dummytxt



!
!     NOTE: bsubua, bsubva are the ALIASED (truncated Fourier spectrum)
!           versions of bsubu, bsubv. bsubuv, bsubvu are used to compute
!           radial current (should be zero) 

      do js = 1, nsurf
            read (njxbout, 200) hs(js), avforce(js),
     1         dum1, dum2, amaxfor(js), aminfor(js)
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt
            read (njxbout,fmt='(a10)')dummytxt

*            write (njxbout, 90)
            lk=0
            do lt = 1, ntht
               do lz = 1, nphi
                  lk = lk + 1
                  read (njxbout, 110) lka(js,lk), jsupu(js,lk),
     1              jsupv(js,lk), dum1, bsupu1(js,lk),  
     2              bsupv1(js,lk),jcrossb(js,lk), pprime(js,lk),
     3              dum2, bdotj(js,lk), bsubua(js,lk),
     4              bsubva(js,lk), bsubs(js,lk)
                  jsups(js,lk)=dmu0*dum1
                  r12(js,lk)=2._rprec*dum2/  !  nomalized residual
     1             (abs(jcrossb(js,lk))+abs( pprime(js,lk)))
               enddo
            enddo
      enddo
      close(unit=njxbout)
      call free_lun(njxbout)
  200 format(/16x,1pe12.3,21x,1p
     1   e12.3,13x,1pe12.3,19x,1pe12.3,/,
     2   32x,1pe12.3,
     3   32x,1pe12.3)
c  110 format(i5,1p12e11.3)
  110 format(i5,12(x,en10.3e2))
      deallocate(bsubuv,bsubvu)
      jxbcode=0
      end subroutine read_old_jxbout

      subroutine resetcolor
       character  endchr*1,ename*50,evalue*50
       data endchr/';'/
       logical symbol
       integer ret, overwrite
       call pgscir(0,15)
       ename = 'PGPLOT_FOREGROUND'
       call getenv(ename,evalue)
       symbol = (trim(evalue) .ne. 'black').or.
     & (trim(evalue) .ne. 'BLACK')
       if(symbol) then
#ifndef IRIX64
         ret=putenv(trim(ename),'BLACK')
         if(ret.lt.0)write(6,*) 'ERROR setting ',trim(ename)
         ret=putenv(trim(ename),'WHITE')
         if(ret.lt.0)write(6,*) 'ERROR setting ',trim(ename)
#else
         call pxfsetenv(trim(ename),0,'BLACK',0,overwrite,ret)
         if(ret.ne.0)write(6,*) 'ERROR setting ',trim(ename)
         call pxfsetenv(trim(ename),0,'WHITE',0,overwrite,ret)
         if(ret.ne.0)write(6,*) 'ERROR setting ',trim(ename)
#endif
       endif
       call pgscr( 1   , 0.00, 0.00, 0.00)  !  black & white
       call pgscr( 0   , 1.00, 1.00, 1.00)  !  interchanged
       call pgscr( 2   , 1.00, 0.00, 0.00)
       call pgscr( 3   , 0.00, 1.00, 0.00)
       call pgscr( 4   , 0.00, 0.00, 1.00)
       call pgscr( 5   , 0.00, 1.00, 1.00)
       call pgscr( 6   , 1.00, 0.00, 1.00)
       call pgscr( 7   , 1.00, 1.00, 0.00)
       call pgscr( 8   , 1.00, 0.50, 0.00)
       call pgscr( 9   , 0.50, 1.00, 0.00)
       call pgscr(10   , 0.00, 1.00, 0.50)
       call pgscr(11   , 0.00, 0.50, 1.00)
       call pgscr(12   , 0.50, 0.00, 1.00)
       call pgscr(13   , 1.00, 0.00, 0.50)
       call pgscr(14   , 0.33, 0.33, 0.33)
       call pgscr(15   , 0.66, 0.66, 0.66)

      end subroutine resetcolor

      subroutine RZtrans(visble,x,y,z)
C-----------------------------------------------
C   WARNING do NOT use       call pgsave in this routine !
C-----------------------------------------------
      use kind_spec
      use read_wout_mod, only: ns,mnmax,xm,xn,rmnc,rmns,zmns,zmnc
      use Toexternal
      implicit none
      integer i, jc
      real(kind=rprec)   :: arg
      real x, y, z, rc, zc, rs, zs, xwrld, ywrld
      INTEGER VISBLE, ins, incntr
      character*20 flabel
      data incntr/0/
      do jc=1,ncont
        if(contrs(jc).eq.z)incntr=jc
      enddo
      call pgsci(1+modulo(incntr,5))
      if(kbold(incntr).eq.1) then
          call pgslw(3)
      else
          call pgslw(1)
      endif
      ins= int(x)
      xwrld = 0.0
      ywrld = 0.0
      do i = 1, mnmax
      arg = 2.0_rprec*3.141592654_rprec*((y-1.0)/real(mxnt-1))*xm(i)
     &    -phit*xn(i)
        if (ins .ne. ns) then
          rc=real((rmnc(i,ins+1)-rmnc(i,ins))*(x-real(ins))+rmnc(i,ins))
          rs=real((rmns(i,ins+1)-rmns(i,ins))*(x-real(ins))+rmns(i,ins))
          zc=real((zmnc(i,ins+1)-zmnc(i,ins))*(x-real(ins))+zmnc(i,ins))
          zs=real((zmns(i,ins+1)-zmns(i,ins))*(x-real(ins))+zmns(i,ins))
        else
          rc=rmnc(i,ins)
          rs=rmns(i,ins)
          zc=zmnc(i,ins)
          zs=zmns(i,ins)
        endif  
        xwrld=xwrld+rc*cos(arg)+rs*sin(arg)
        ywrld=ywrld+zc*cos(arg)+zs*sin(arg)
      enddo
      if (visble.eq.0) then
        if(klabel(incntr).ne.0)then
        klabel(incntr)=0.
        write(flabel,fmt=labelformat)contrs(incntr)
        call pgsch(.8)
        call pgsci(1.)
        call pgslw(1)
        call pgscf(2)
        call pgpt1(xwrld,ywrld,-4)!   make a diamond at the labelled contour
        call pgtext(xwrld,ywrld,trim(flabel))
        endif
        call pgmove (xwrld, ywrld)
      else
        call pgdraw (xwrld, ywrld)
      endif
      end

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
      real(kind=rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(kind=rprec) :: t
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
         enddo
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         enddo
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
 
      subroutine spline(n, x, y, b, c, d)
c       the codes (spline & seval) are taken from:
c       forsythe,malcolm and moler,
c       "computer methods for mathematical computations",
c       prentice-hall, 1977.
c
c       the codes (spleen,splaan & speval) are adaptations
c       by r.m. wieland for special cases ... see comments
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n
      real(kind=rprec), dimension(n) :: x, y, b, c, d
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nm1, ib, i
      real(kind=rprec) :: t
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
         b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
         c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
         c(2:nm1) = c(3:nm1+1) - c(2:nm1)
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
         b(1) = -d(1)
         b(n) = -d(n-1)
         c(1) = 0.
         c(n) = 0.
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
         enddo
c
c  back substitution
c
         c(n) = c(n)/b(n)
         do ib = 1, nm1
            i = n - ib
            c(i) = (c(i)-d(i)*c(i+1))/b(i)
         enddo
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
         b(n) = (y(n)-y(nm1))/d(nm1) + d(nm1)*(c(nm1)+2.*c(n))
         b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1) - d(:nm1)*(c(2:nm1+1) + 
     1      2.d0*c(:nm1))
         d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
         c(:nm1) = 3.d0*c(:nm1)
         c(n) = 3.d0*c(n)
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

      end subroutine spline

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
      real(kind=rprec), dimension(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, klo, khi, k
      real(kind=rprec) :: c1o6, deriv, a, b, h, h2, a2, b2, y26lo, y26hi
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
         enddo
 
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
      enddo
      end subroutine splint

      subroutine str_strip(string, inblen)
C
C...       This routine compresses a character string, deleting all blanks
c...       the compressed string is returned in the pkg it came in [STRING]],
c...       together with a label [INBLEN] telling how long it is. (Dick
c
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
      enddo
      if (ind < ilen) string(ind+1:) = ' '
 
      inblen = ind
 
      end subroutine str_strip

      subroutine bndypts (ntheta,nplots,kz,rb,zb)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, 
     .only: nfp, mnmax, xm, xn, rmnc, zmns, rmns, zmnc, rprec, ns
c      use read_wout_mod, only: nfp, rprec
      use vmec_input, only: rbc,rbs,zbc,zbs, lfreeb, ntor, mpol
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, nplots, kz, mpol_mn, ntor_mn
      real(kind=rprec), dimension(*), intent(out) :: rb, zb
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(kind=rprec) :: pit, xxm, xxn, arg, u, v, r, z, pi
      real(kind=rprec) :: cosa, sina, r1, z1, r0, z0
      integer :: j, m, n
C-----------------------------------------------
      mpol_mn=nint(maxval(xm))
      ntor_mn=nint(maxval(xn))/nfp
      pi=atan(1.d0)*4.d0
      pit = 2.d0*pi/(ntheta - 1)
      v=pi*(nfp-1)*(kz-1)/nplots
        r1=0;z1=0;m=0;r0=0;z0=0
         do n=-ntor,ntor
          xxm=m;xxn=n
          arg=xxm*u-xxn*v
          cosa=cos(arg);sina=sin(arg)
          r1=r1+rbc(n,m)*cosa+rbs(n,m)*sina
          z1=z1+zbc(n,m)*cosa+zbs(n,m)*sina
         enddo


         do n=1,ntor_mn
          arg=xn(n)/nfp*v
          cosa=cos(arg);sina=sin(arg)
          r0=r0+rmnc(n,ns)*cosa+rmns(n,ns)*sina
          z0=z0+zmnc(n,ns)*cosa+zmns(n,ns)*sina
         enddo
       write(6,fmt='(a,i3.3,a,/4(x,1pe11.4))')
     .  ' phi = ',nint(360*v/2/pi),
     .  ' compare [R,Z]cntrs {mn,bc} ',r0,r1,z0,z1
       do j=1,ntheta
        u=pit*(j-1)
        r=0
        z=0
        do m=0,mpol  
         do n=-ntor,ntor
          xxm=m;xxn=n
          arg=xxm*u-xxn*v
          cosa=cos(arg);sina=sin(arg)
          r=r+rbc(n,m)*cosa+rbs(n,m)*sina
          z=z+zbc(n,m)*cosa+zbs(n,m)*sina
         enddo
        enddo
        rb(j)=r
        zb(j)=z
       enddo

      return
!      end  subroutine bndypts
      end 

      subroutine totz(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, 
     .only: nfp, mnmax, xm, xn, rprec, mpol, ntor
      use Vpname1, only: nrt
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz = 1.d0/(nfp*nplots)
      nrth = ns1*ntheta
      if (nrth .gt. nrt) stop 'nrth > nrt in totz'

      r(:nrth) = 0.d0
      z(:nrth) = 0.d0
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
            enddo
         enddo
      enddo
      end subroutine totz

      subroutine totzu(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, only: rprec, nfp, mnmax, xm, xn
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1  rmnc, zmns, rmns, zmnc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz = 1.d0/(nfp*nplots)
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
            enddo
         enddo
      enddo

      end subroutine totzu
      subroutine totb(ntheta, ns1, nplots, kz, r, z, 
     1   rmnc, zmns, rmns, zmnc)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use read_wout_mod, 
     .only: nfp, mnmax, xm, xn, rprec, mpol, ntor
      use Vpname1, only: nrt
      use Vindat2, only: twopi
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, ns1, nplots, kz
      real(kind=rprec), dimension(*), intent(out) :: r, z
      real(kind=rprec), dimension(*), intent(in) :: 
     1   rmnc, zmns, rmns, zmnc

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nrth, l, j, jes, mes, mn, kt, m, iarg
      real(kind=rprec) :: pit, piz, xm0, xn0, arg
C-----------------------------------------------
      pit = 1.d0/(ntheta - 1)
      piz = 1.d0/(nfp*nplots)
      nrth = 1*ntheta
      if (nrth .gt. nrt) stop 'nrth > nrt in totz'

      r(:nrth) = 0.d0
      z(:nrth) = 0.d0
      j = ns1
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
               r(kt) = r(kt) + rmnc(m)*cos(arg) + rmns(m)*sin(arg)
               z(kt) = z(kt) + zmns(m)*sin(arg) + zmnc(m)*cos(arg)
            enddo
         enddo
      end subroutine totb

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
      enddo
 
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
      enddo
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
      enddo
 
      end subroutine wie_parse

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
      integer :: iloop, nmax1, mn, n1, j, jp, i
      real(kind=rprec) :: es, tb1, tm1, tv1, tp1, ti1, ub1
      character, dimension(3) :: ichar1*5, ichar2*5
C-----------------------------------------------
 
      data ichar1/'RmnC(', 'ZmnS(', 'LmnS('/
      data ichar2/'RmnS(', 'ZmnC(', 'LmnC('/
      twopi = 8.*atan(1.0)
*
*                 PRINT SYMMETRIC TERMS -- RmnC and ZmnS
*                 INTERPOLATE LAMBDA ONTO FULL MESH FIRST
*

      do mn = 1, mnmax
         if (ixm(mn) .ne. 0) then
            lmns(mn,1) = 0.
         else
            lmns(mn,1) = 1.5*lmns(mn,2) - 0.5*lmns(mn,3)
         endif
      enddo
               
      do j = 2, ns-1
         do mn = 1,mnmax
            lmns(mn,j) = 0.5*(lmns(mn,j) + lmns(mn,j+1))
         enddo         
      enddo   
      
      do mn = 1,mnmax
         lmns(mn,ns) = 2.0*lmns(mn,ns-1) - lmns(mn,ns-2)
      enddo      

      write (39, 200)
      do iloop = 1, 3
         nmax1 = 5
         do mn = 1, mnmax, 6
            if (mn > mnmax - 6) nmax1 = mnmax - mn
            write (39, 210) (ichar1(iloop),ixm(mn+n1),ixn(mn+n1),n1=0,
     1         nmax1)
            do j = 1, ns
               es = (j - 1)*hs
               select case (iloop) 
               case default
                  write (39, 220) es, (rmnc(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (2) 
                  write (39, 220) es, (zmns(mn+n1,j),n1=0,nmax1)
                  cycle 
               case (3) 
                  write (39, 220) es, (lmns(mn+n1,j),n1=0,nmax1)
               end select
            enddo
         enddo
      enddo
 
*
*                 PRINT ASYMMETRIC TERMS -- RmnS and ZmnC
*
      if (iasym .eq. 1) then
         write (39, 200)
         do iloop = 1, 3
            nmax1 = 5
            do mn = 1, mnmax, 6
               if (mn > mnmax - 6) nmax1 = mnmax - mn
               write (39, 210) (ichar2(iloop),ixm(mn+n1),ixn(mn+n1),
     1            n1=0,nmax1)
               do j = 1, ns
                  es = (j - 1)*hs
                  select case (iloop) 
                  case default
                     write (39, 220) es, (rmns(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (2) 
                     write (39, 220) es, (zmnc(mn+n1,j),n1=0,nmax1)
                     cycle 
                  case (3) 
                     write (39, 220) es, (lmnc(mn+n1,j),n1=0,nmax1)
                  end select
               enddo
            enddo
         enddo
      endif
*
*     DETERMINE RADIAL BETA PROFILE (SURFACE AVERAGED)
*
      phi(1) = 0.
      sqrt_phimod = 0.
      phimod = 0.
      do j = 1, ns
         sqrt_phimod(mn0+mnmax*(j-1)) = sqrt(abs(phi(j)))
         phimod(mn0+mnmax*(j-1)) = phi(j)
      enddo
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
      enddo
      ub(ns) = twopi*ub1
      vp(ns) = tv1
      mass(ns) = tm1
      pres(ns) = tp1
      beta_vol(ns) = tb1
      buco(ns) = 2.*buco(ns) - buco(ns-1)
      bvco(ns) = 2.*bvco(ns) - bvco(ns-1)
      write (39, 230)
*     Normalized FF-Prime : ffp(j)
      do j = 1, ns
         es = (j - 1)*hs
         itor = -twopi*buco(j)
         itors(j) = itor/dmu0
         ipol = twopi*(bvco(j)-bvco(ns))
         ipols(j) = ipol/dmu0
         ffp(j) = bvco(j)*dbzcods(j)/phip(j)/iotaf(j)/
     1      bvco(ns)**2
         write (39, 220) es, bvco(j), buco(j), itor, ipol, ffp(j)
      enddo
      write (39, 240)
      do i = 1, ns
         es = (i - 1)*hs
         write (39, 220) es, twopi**2*vp(i), ub(i), mass(i), pres(i), 
     1      iotaf(i), beta_vol(i)
      enddo
 
      return 
 
*  Do the above before the call to PLOTTER; do the below after
 
      entry wrfcn2 (input_id)
 
      write (39, 260)
      do i = 1, ns
         es = (i - 1)*hs
         write (39, 220) es, jcurv(i), jdotb(i)
      enddo
 
      write (39, 250) beta_vol(1)
 
  200 format(//,40x,'FOURIER COEFFICIENTS X(m,n)',/)
  210 format(//,9x,' S ',4x,6(4x,a5,i1,',',i3,')'),/)
  220 format(1p8e15.3)
  230 format(//,25x,'COVARIANT COMPONENTS OF B',
     1   ' AND INTEGRATED CURRENTS',2/,9x,' S ',10x,'<BZETA>',8x,
     2   '<BTHETA>',8x,'ITOR',11x,'IPOL',8x,'FF''',/)
  240 format(//,9x,' S ',11x,'VP',12x,'dV/dPHI',10x,'MASS',13x,'P',10x,
     1   'IOTA',12x,'BETA',/)
  250 format(//,'  BETA ON AXIS (SUM OVER MODES) = ',1pe10.3)
  260 format(//,9x,' S ',11x,'<JTOR>',2x,'<JdotB>/<bdotgradv>',/,9x,3x
     1   ,11x,'[A/M2]',10x,'[A/M]',/)
 
      end subroutine wrfcn

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
      real(kind=rprec) :: itor
      real(kind=rprec) :: rbt, betatot, betapol, betator, betaxis
      character runid*(*)
      real(kind=rprec), dimension(*) :: psi_in, iota, ffp 
      real(kind=rprec), dimension(*) :: rmnc, zmns, pres, jdotb
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(kind=rprec), parameter :: dmu0 = 1.256637d-06, zero = 0.d0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ishot, itime, mn, mjbal
      real(kind=rprec), dimension(101) :: sp1, sp2, sp3
      real(kind=rprec) :: err
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
      goto 7
    5 continue
      shotnm = zero
      time = zero
      goto 20
    7 continue
      i = index(runid(9:),'t')
      if (i .eq. 0) time = zero
      if (i > 0) read (runid(9+i:9+i+2), '(i3)', err=10) itime
      time = itime/100
      goto 20
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
*       calculate derivatives using interpolating cubic splines
*       use Natural BC for q and pres
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

      subroutine xywtoxyan(xw,yw,xso,yso)
       use pgplotinc, only: pgxscl,pgxorg,pgxsz,pgyscl,pgyorg,pgysz
c       include '/usr/local/pgplot/pgplot.inc' ! might not be there
       integer idx idc
       real xw,yw,xso,yso,xa,ya
         xa(xw,idc)=(xw*pgxscl(idc)+pgxorg(idc))
         ya(yw,idc)=(yw*pgyscl(idc)+pgyorg(idc))
         call pgqid(idx)
         xso=xa(xw,idx)/pgxsz(idx)
         yso=ya(yw,idx)/pgysz(idx)
         return
      end subroutine xywtoxyan

