      program vmec
      use vmec_input
      use vmec_seq
      use read_namelist_mod
      use safe_open_mod
      use vparams, only: nlog, nlog0
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nseq0 = 12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: numargs, ierr_vmec, index_end,
     1   iopen, isnml, iread, iseq, index_seq, 
     2   index_dat, ireset, iunit
      character*120 :: input_file, seq_ext, reset_file_name, arg
      character*120 :: log_file
      character*120, dimension(10) :: command_arg
      logical :: lfirst=.true., lreseta, lscreen
C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see subroutine BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       where phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For all but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' domain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!   
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s function table for the
!       vacuum magnetic field(s) and, if data analysis is to be done, flux and
!       field loops as well. The user provides a subroutine (BFIELD) which can be
!       called at an arbitrary spatial location and which should return the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required if
!       equilibrium reconstruction is to be done.)
!   
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, where 'EXT' is the command-line extension of the INPUT file.
!
!   
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to any one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!   
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these parameter statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the namelist INDATA.
!
!
!       Added features since last edition
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/sqrt(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!***

!    
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      call getcarg(1, command_arg(1), numargs)
      do iseq = 2, numargs
         call getcarg(iseq, command_arg(iseq), numargs)
      end do   

      lreseta = .true.            !!Default value: runvmec MUST be called this way the first time
      lscreen = .true.
      
      if (numargs .lt. 1) then
         stop 'Invalid command line'
      else if (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') then
         print *,
     1   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         print *
         print *,' For example: '
         print *,'    xvmec input.tftr OR xvmec tftr ',
     1           'OR xvmec ../input.tftr'
         print *
         print *,' Sequence files, containing a LIST of input files',
     1           ' are also allowed: '
         print *,'    xvmec input.tftr_runs'
         print *
         print *,' Here, input.tftr_runs contains a &VSEQ namelist',
     1           ' entry'
         print *
         print *,' Additional (optional) command arguments are',
     1           ' allowed:'
         print *
         print *,'    xvmec <filename> noscreen F'
         print *
         print *,' noscreen: supresses all output to screen ',
     1           ' (default, or "screen", displays output)'    
         print *,' F (or T): if "T", forces reset on',
     1           ' a coarse mesh (used for sequencing control)'   
         
         stop
      else if (numargs .gt. 1) then
         arg = command_arg(2)
         if (trim(arg).eq.'noscreen' .or. trim(arg).eq.'NOSCREEN')
     1      lscreen = .false.
      end if
      if (numargs .gt. 2) then
          arg = command_arg(3)
          if (arg(1:1).eq.'f' .or. arg(1:1).eq.'F') lreseta = .false.   
      end if    
      if (numargs .gt. 3) then
          reset_file_name = command_arg(4)
      end if    
 

!
!     Determine type of file opened (sequential or input-data)      
!     ARG1 (char var)
!          By default, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to open file ARG1 (full path + file name).
!                  Look for the VSEQ namelist to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i] 
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence namelist, then the data file 
!                  ARG1 or input.ARG1 is read directly, with NSEQ=1.
!
      arg = command_arg(1)
      index_dat = index(arg,'.')
      index_end = len_trim(arg)
      if (index_dat .gt. 0) then
         seq_ext  = arg(index_dat:index_end)
         input_file = trim(arg)
      else
         seq_ext = trim(arg)
         input_file = 'input.'//trim(seq_ext)
      end if   
 
      if (numargs .le. 3) reset_file_name = 'wout.' // seq_ext
      
      nseq = 1
      nseq_select(1) = 1
      extension(1) = input_file
!
!     READ IN NAMELIST VSEQ TO GET ARRAY
!     OF INPUT FILE EXTENSIONS AND INDEXING ARRAY, NSEQ_SELECT
!
      nlog = nlog0
      iunit = nseq0
      do iseq = 1, 2
         if (iseq .eq. 1) then
           arg = input_file
         else
           arg = seq_ext
         end if
         call safe_open(iunit, iopen, trim(arg), 'old', 'formatted')
         if (iopen .eq. 0) then
           call read_namelist (iunit, isnml, 'vseq')
           if (isnml.eq.0 .and. nseq .gt. nseqmax) stop 'NSEQ>NSEQMAX'
 
!
!       OPEN FILE FOR STORING SEQUENTIAL RUN HISTORY
!
           if (isnml .eq. 0) then
              log_file = 'log.'//seq_ext
 
              call safe_open(nlog, iread, log_file, 'replace', 
     1           'formatted')
              if (iread .ne. 0) then
                 print *, log_file, 
     1           ' LOG FILE IS INACCESSIBLE: IOSTAT= ',iread
                 stop 3
              else
                 exit        !!Break out of loop
              end if
           endif
        endif  

        close (iunit)

      end do  

!
!     CALL EQUILIBRIUM SOLVER 
!
!     nseq_select:      if sequence file (VSEQ namelist given with nseq >0)
!                       array giving indices into EXTENSION array prescribing
!                       the order in which the input files are run by VMEC
!     nseq:             number of sequential VMEC runs to make
!
!
!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
!
      do iseq = 1, nseq 
         index_seq = nseq_select(iseq)
         ireset = 0
         ierr_vmec = 0
         if (iseq .gt. 1) reset_file_name = 
     1       'wout.' // trim(extension(index_seq))
 100     continue
         call runvmec (extension(index_seq), iseq-1, lreseta, ierr_vmec, 
     1                 lfirst, lscreen, reset_file_name)
         lfirst = .false. 
         select case (ierr_vmec) 
!        case (1:2)    !BAD JACOBIAN AFTER 75 ITERATIONS...
!          ireset = ireset + 1
!          lreseta = .true.
!          if (ireset .le. 2) go to 100
         case (4)                                !Try a few more iterations
           ireset = ireset + 1
           lreseta = .false.
           if (ireset .le. 1) then
              if (lscreen) write (*, '(/,1x,a)') 
     1           'RUNNING A FEW MORE ITERATIONS THAN REQUESTED'
              go to 100
           else if (lscreen) then
              print *, 'DECREASE DELT OR INCREASE NITER'
           endif
         case (6)    !BAD JACOBIAN AFTER AXIS RESET: TRY DECREASING TO NS=3
           ireset = ireset + 1
           lreseta = .true.
           if (ireset .le. 1) go to 100             
         case default
           lreseta = .false.
         end select
      end do

!
!     FREE ANY LONG-TERM (PERSISTENT THROUGH ISEQ > 1, OR XC, SCALXC FOR
!     ITERATIVE OPTIMIZATION) POINTERS
!
      call free_persistent_mem
 
      close (nlog)
 
      end program vmec
      

      subroutine bcovar (lu, lv)
      use vmec_main
      use vmec_params, only: ns4, signgs
      use realspace, weight_f => sqrts, sqrts => sqrts
      use vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1   rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2   bsubu_e => clmn_e, bsubv_e => blmn_e, bsubu_o => clmn_o,
     3   bsubv_o => blmn_o, bsq => bzmn_o, phipog => brmn_o
      use vsvd, only: phifac, phifsave, imovephi
      use xstuff, only: xc
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt,0:1) :: lu, lv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
      real(rprec), parameter :: c1p5 = 1.5_dp, p5 = 0.5_dp,
     1    pdamp = 0.05_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, is, js, ndim
      real(rprec) :: r2, r3
      real(rprec), dimension(:), pointer :: r12sq, bsupu, bsupv,
     1    bsubuh, bsubvh
      real(rprec) :: luu12, luv12, lvv12
      real(rprec) :: arnorm, aznorm, volume, tcon_mul
      real(rprec), allocatable, dimension(:) :: luu, luv, lvv
      real(rprec), external :: dot_g
C-----------------------------------------------
      ndim = 1+nrzt
      r12sq => bsubu_o

      allocate (luu(ndim), luv(ndim), lvv(ndim), stat = is)

      if (is .ne. 0) stop 'allocation error in bcovar'

!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      guu = 0                           ! mcz
      guv = 0;  gvv = 0

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
!

CDIR$ IVDEP
      do l = 1, nrzt
         r2 = sqrts(l)*sqrts(l)
         guu(l)   =  ru(l,0)*ru(l,0) + r2*ru(l,1)*ru(l,1)
     1            +  zu(l,0)*zu(l,0) + r2*zu(l,1)*zu(l,1)
         luu(l)   = (ru(l,0)*ru(l,1) + zu(l,0)*zu(l,1))*2
         r12sq(l) = r1(l,0)*r1(l,0) + r2*r1(l,1)*r1(l,1)
         phipog(l)= 2*r1(l,0)*r1(l,1)
      end do

      if (lthreed) then
CDIR$ IVDEP
         do l = 1, nrzt
            r2 = sqrts(l)*sqrts(l)
            guv(l)   =  ru(l,0)*rv(l,0) + r2*ru(l,1)*rv(l,1)
     1               +  zu(l,0)*zv(l,0) + r2*zu(l,1)*zv(l,1)
            luv(l)   =  ru(l,0)*rv(l,1) + ru(l,1)*rv(l,0)
     1               +  zu(l,0)*zv(l,1) + zu(l,1)*zv(l,0)
            gvv(l)   =  rv(l,0)*rv(l,0) + r2*rv(l,1)*rv(l,1)
     1               +  zv(l,0)*zv(l,0) + r2*zv(l,1)*zv(l,1)
            lvv(l)   = (rv(l,0)*rv(l,1) + zv(l,0)*zv(l,1))*2
         end do
      end if
CDIR$ IVDEP
      do l = nrzt, 2, -1
         guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
         r12sq(l) = p5*(r12sq(l) + r12sq(l-1) + shalf(l)*
     1                (phipog(l) + phipog(l-1)))
      end do
      if (lthreed) then
CDIR$ IVDEP
         do l = nrzt, 2, -1
            guv(l) = p5*(guv(l) + guv(l-1) +
     1         shalf(l)*(luv(l) + luv(l-1)))
            gvv(l) = p5*(gvv(l) + gvv(l-1) +
     1         shalf(l)*(lvv(l) + lvv(l-1)))
         end do
      end if

      do js = 2, ns
         vp(js) = signgs*DOT_G(nznt,gsqrt(js),ns,wint(js),ns)
      end do
      if (iter2 .eq. 1) voli = twopi*twopi*hs*sum(vp(2:ns))

      gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)
      phipog(2:nrzt) = phip(2:nrzt)/gsqrt(2:nrzt)
      phipog(1)    = 0
      phipog(ndim) = 0
!
!     RECONSTRUCT IOTA, PRESSURE PROFILE FROM DATA
!
      if (lrecon) call newprofil (phipog)
      if (.not.lrecon .and. imovephi.gt.0) call newphi (phipog)

!
!     STORE CURRENT(0), MAGNETIC PITCH, PRESSURE DATA
!     STORE HALF-MESH VALUES OF LU, LV FOR ACCURATE CURRENT CALCULATION
!
      if (iequi .eq. 1) call storesvd (r1(1,0), r1(1,1), lu(1,0),
     1   lu(1,1), lv(1,0), phipog, zu0)

!
!     COMPUTE COVARIANT COMPONENTS OF B ON RADIAL FULL-MESH
!     NOTE: LU = 1+LAMU, LV = -LAMV COMING INTO THIS ROUTINE
!     WILL ADD IOTAF, PHIP/GSQRT FACTOR LATER...
!          LV == BSUPU = (PHIP/GSQRT)*(iota - LAMV),
!          LU == BSUPV = (PHIP/GSQRT)*(1    + LAMU)
!

      bsupu => bsubu_o                  !!Save these later in lv(l,0), lu(l,0)
      bsupv => bsubv_o

!
!     FIRST, PUT LAMBDA DERIVATIVES ON RADIAL HALF-MESH
!
CDIR$ IVDEP
      do l = 2, nrzt
         bsupv(l) = p5*phipog(l)*(lu(l,0) + lu(l-1,0) + shalf(l)*
     1                           (lu(l,1) + lu(l-1,1)))
         bsupu(l) = p5*phipog(l)*(lv(l,0) + lv(l-1,0) + shalf(l)*
     1                           (lv(l,1) + lv(l-1,1)))
      end do

      bsupv(1) = 0
      bsupu(1) = 0

!
!     COMPUTE UPDATED IOTA PROFILE
!
      call getiota(phipog, bsupu, bsupv)


!     ADD PRESENT VALUE OF IOTAF HERE.
      do js = 1, ns
         lv(js:nrzt:ns,0) = lv(js:nrzt:ns,0) + iotaf(js)
      end do

!
!     NEXT COMPUTE LAMBDA FORCES ON FULL MESH BY AVERAGING HALF-MESH METRICS
!
      luu = phipog(:ndim)*guu
      luv = phipog(:ndim)*guv
      lvv = phipog(:ndim)*gvv

CDIR$ IVDEP
      do l = 1,nrzt
         luu12 = p5*(luu(l) + luu(l+1))
         luv12 = p5*(luv(l) + luv(l+1))
         lvv12 = p5*(lvv(l) + lvv(l+1))
         bsubu_e(l) = luu12*lv(l,0) + luv12*lu(l,0)
         bsubv_e(l) = luv12*lv(l,0) + lvv12*lu(l,0)
      end do

      luu = luu*shalf
      luv = luv*shalf
      lvv = lvv*shalf

CDIR$ IVDEP
      do l = 1,nrzt
         luu12 = p5*(luu(l) + luu(l+1))
         luv12 = p5*(luv(l) + luv(l+1))
         lvv12 = p5*(lvv(l) + lvv(l+1))
         bsubu_e(l) = bsubu_e(l)+ luu12*lv(l,1) + luv12*lu(l,1)
         bsubv_e(l) = bsubv_e(l)+ luv12*lv(l,1) + lvv12*lu(l,1)
      end do

!
!     FINALLY, COMPUTE LAMBDA FORCES ON RADIAL HALF-MESH
!
      lu(:nrzt,0) = bsupv(:nrzt)
      lv(:nrzt,0) = bsupu(:nrzt)

      bsubuh => bsubu_o
      bsubvh => bsubv_o

      bsubuh(:nrzt) = guu(:nrzt)*lv(:nrzt,0) + guv(:nrzt)*lu(:nrzt,0)
      bsubvh(:nrzt) = guv(:nrzt)*lv(:nrzt,0) + gvv(:nrzt)*lu(:nrzt,0)

      bsubuh(ndim) = 0
      bsubvh(ndim) = 0

      do js = 1,ns
         fpsi(js) = DOT_G(nznt,bsubvh(js),ns,wint(js),ns)
      enddo

!
!     COMPUTE TOROIDAL AND POLOIDAL CURRENTS
!
      rbtor = c1p5*fpsi(ns) - cp5*fpsi(ns-1)
      rbtor0= c1p5*fpsi(2)  - cp5*fpsi(3)
      ctor = signgs*twopi*sum((c1p5*bsubuh(ns:nrzt:ns) -
     1   cp5*bsubuh(ns-1:nrzt:ns))*wint(ns:nrzt:ns))

!
!     COMPUTE KINETIC PRESSURE ON HALF-GRID
!     IF R*dp/ds NONVARIATIONAL FORM TO BE USED IN FORCES, UNCOMMENT ALL THE "RPRES" 
!     COMMENTS IN THIS SUBROUTINE, AS WELL AS FORCES, FUNCT3D SUBROUTINES
!
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*sum(vp(2:ns)*pres(2:ns))
      do js = 2,ns
         bsq(js:nrzt:ns) = pres(js)
      end do   
!RPRES     bsq = 0
!RPRES     do js = 2,ns
!RPRES        lu(js:nrzt:ns,1) = pres(js)
!RPRES     end do

!
!     COMPUTE MAGNETIC PRESSURE + KINETIC PRESSURE
!
      bsq(:nrzt) = bsq(:nrzt) + p5*(lv(:nrzt,0)*bsubuh(:nrzt) +
     1      lu(:nrzt,0)*bsubvh(:nrzt))
      wb = -wp + hs*abs(sum(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))
!RPRES      wb = hs*abs(sum(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))


!
!     AVERAGE LAMBDA FORCES ONTO FULL MESH
!
CDIR$ IVDEP
      do l = 1, nrzt
         r2 = (1 - weight_f(l)*weight_f(l))**2 * pdamp
         r3 = p5*(1 - r2)
         bsubu_e(l) = bsubu_e(l) * r2
     1              + r3*(bsubuh(l) + bsubuh(l+1))
         bsubv_e(l) = bsubv_e(l) * r2
     1              + r3*(bsubvh(l) + bsubvh(l+1))
      end do

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      if (iequi .eq. 1) then

!DBG     luu = bsubuh
!DBG     lvv = bsubvh

!RPRES         bsq(:nrzt) = bsq(:nrzt) + lu(:nrzt,1)

         do js = ns-1,2,-1
         do l = js, nrzt, ns
            bsubuh(l) = 2*bsubu_e(l) - bsubuh(l+1)
            bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
         end do
         end do

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)


!DBG     do js = 2, ns
!DBG        print *, 'JS = ', js
!DBG        print *,
!DBG 1      '    BSUBU(old)     BSUBU(new)    BSUBV(old)     BSUBV(new)'
!DBG        do l = js, nrzt, ns
!DBG           write(*,1223) luu(l), bsubu_e(l), lvv(l), bsubv_e(l)
!DBG        end do
!DBG     end do

 1223    format(1p4e14.5)

         go to 1000

      end if

      bsubu_o(:nrzt) = sqrts(:nrzt)*bsubu_e(:nrzt)
      bsubv_o(:nrzt) = sqrts(:nrzt)*bsubv_e(:nrzt)

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS EVERY NS4 STEPS.
!
      if(mod(iter2-iter1,ns4).eq.0)then
         phifsave = phifac
         phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt)
         call lamcal(phipog, guu, guv, gvv)
         call precondn(lu,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,crd)
         call precondn(lu,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                r1(1,1),azm,azd,bzm,bzd,crd)

         guu(:ndim) = guu(:ndim)*r12(:ndim)**2
         volume = hs*sum(vp(2:ns))
         r2 = max(wb,wp)/volume
         fnorm = one/(sum(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))
         fnorm1 = one/sum(xc(1+ns:2*irzloff)**2)

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!
!        OVERRIDE USER INPUT VALUE HERE
!
         r2 = ns
         tcon0 = min(abs(tcon0), one)                      !!ignore large tcon0 value from old-style file
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         do js = 2, ns-1
           arnorm = sum(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
           aznorm = sum(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
           if (arnorm .eq. zero .or. aznorm .eq. zero)
     1     stop 'arnorm or aznorm=0'

           tcon(js) = min(abs(ard(js,1)/arnorm),abs(azd(js,1)/
     1          aznorm))                 * tcon_mul * (32*hs)**2
         end do
         tcon(ns) = cp5*tcon(ns-1)

      endif

!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
      do l=2,nrzt
         guu(l) = lv(l,0)*lv(l,0)*gsqrt(l)
         guv(l) = lv(l,0)*lu(l,0)*gsqrt(l)
         gvv(l) = lu(l,0)*lu(l,0)*gsqrt(l)
         lv(l,0)  = bsq(l)*gsqrt(l)/r12(l)
         lu(l,0)  = bsq(l)*r12(l)
      enddo

 1000 continue

      deallocate (luu, luv, lvv, stat = l)

      end subroutine bcovar


      subroutine getiota(phipog, bsupu, bsupv)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) :: phipog, bsupv
      real(rprec), dimension(nrzt), intent(inout) :: bsupu
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l
      real(rprec) :: top, bot
C-----------------------------------------------

      if (ncurr .eq. 0) goto 100
      do js = 2, ns
         top = jcurv(js)
         bot = 0
         do l = js, nrzt, ns
            top = top - wint(l)*(guu(l)*bsupu(l) + guv(l)*bsupv(l))
            bot = bot + wint(l)*phipog(l)*guu(l)
         end do
         iotas(js) = top/bot
      end do

!     Do not compute iota too near origin
      iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)           !!zero gradient near axis
      iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
      do js = 2, ns-1
         iotaf(js) = p5*(iotas(js) + iotas(js+1))
      end do

 100  continue

      do js = 2, ns
         bsupu(js:nrzt:ns) = bsupu(js:nrzt:ns) + phipog(js:nrzt:ns)
     1                     * iotas(js)
      end do

      end subroutine getiota


      subroutine lamcal(phipog, guu, guv, gvv)
      use vmec_main
      use vmec_params, only: ntmax, jlam
      use realspace, only: sqrts
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt), intent(in) ::
     1   phipog, guu, guv, gvv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: damping_fac = -2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m,n,js
      real(rprec) :: tnn, tnm, tmm, power
C-----------------------------------------------

      blam(:ns) = sum(guu *phipog, dim=2)
      clam(:ns) = sum(gvv *phipog, dim=2)
      dlam(:ns) = sum(guv *phipog, dim=2)
      blam(1) =  blam(2)
      clam(1) =  clam(2)
      dlam(1) =  dlam(2)
      blam(ns+1) =  0
      clam(ns+1) =  0
      dlam(ns+1) =  0
      do js = 2, ns
        blam(js) = cp5*(blam(js) + blam(js+1))
        clam(js) = cp5*(clam(js) + clam(js+1))
        dlam(js) = cp5*(dlam(js) + dlam(js+1))
      end do

!
!       REDUCE FACLAM AT SOME FSQ THRESHOLD TO IMPROVE 3D CONVERGENCE
!
      do m = 0, mpol1
         tmm = m*m
         power = min(tmm/256, 8._dp)
         do n = 0, ntor
            if (m.eq.0 .and. n.eq.0) cycle
            tnn = (n*nfp)**2
            tnm = 2*m*n*nfp
            do js = jlam(m), ns
               faclam(js,n,m,1) = 2*damping_fac/
     1         ((blam(js)+blam(js+1))*tnn
     2         + sign((dlam(js)+dlam(js+1)),blam(js))*tnm
     3         + (clam(js) + clam(js+1))*tmm)
     4         * sqrts(js)**power                                   !Damps m > 16 modes
            end do
         end do
      end do

      do n = 2, ntmax
         faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
      end do

      end subroutine lamcal


      subroutine precondn(lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo,
     1   xodd, axm, axd, bxm, bxd, cx)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) ::
     1  lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
      real(rprec), dimension(ns+1,2), intent(out) ::
     1  axm, axd, bxm, bxd
      real(rprec), dimension(ns+1), intent(out) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, lk
      real(rprec), dimension(:,:), allocatable :: ax, bx
      real(rprec), dimension(:), allocatable :: ptau
      real(rprec) t1, t2, t3
C-----------------------------------------------
!
!     COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z
!     FORCE (ALL ARE MULTIPLIED BY 0.5)
!

      allocate (ax(ns+1,4), bx(ns+1,4), ptau(nznt))

      ax = 0
      bx = 0
      cx = 0

      do 20 js = 2,ns
!
!     COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING
!     MATRIX ELEMENTS
!
        lk = 0
        do l = js,nrzt,ns
          lk = lk + 1
          ptau(lk) = r12(l)*r12(l)*bsq(l)*wint(l)/gsqrt(l)
          t1 = xu12(l)*ohs
          t2 = cp25*(xue(l)/shalf(js) + xuo(l))/shalf(js)
          t3 = cp25*(xue(l-1)/shalf(js) + xuo(l-1))/shalf(js)
          ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
          ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
          ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
          ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
        end do
!
!       COMPUTE ORDER M**2 PRECONDITIONING MATRIX ELEMENTS
!
        lk = 0
        do l = js,nrzt,ns
          lk = lk+1
          t1 = cp5*(xs(l) + cp5*xodd(l)/shalf(js))
          t2 = cp5*(xs(l) + cp5*xodd(l-1)/shalf(js))
          bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
          bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
          bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
          cx(js) = cx(js) + cp25*lu1(l)**2*gsqrt(l)*wint(l)
        end do
 20   continue

      do js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)*sp(js)
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)*sp(js)
        cx(js)    = cx(js) + cx(js+1)
      end do

      axm(ns+1,:) = 0
      bxm(ns+1,:) = 0
      axd(ns+1,:) = 0
      bxd(ns+1,:) = 0
      cx(ns+1) = 0

      deallocate (ax, bx, ptau)

      end subroutine precondn

      subroutine forces
      use vmec_main
      use realspace, z1 => z1, gcon => z1
      use vforces, crmn_e => crmn_e, lv_e => crmn_e, 
     1   czmn_e => czmn_e, lu_e => czmn_e,
     2   czmn_o => czmn_o, lu_o => czmn_o
!RPRES , presg => czmn_o
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, lk, ndim
      real(rprec), dimension(:), allocatable :: 
     1    bsqr, gvvs, guvs, guus
      real(rprec) :: rcon1, zcon1, dphids
C-----------------------------------------------

      ndim = 1+nrzt 
      allocate (bsqr(ndim), gvvs(ndim), guvs(ndim), guus(ndim), 
     1           stat=l)
      if (l .ne. 0) stop 'Allocation error in VMEC FORCES'
      
!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p        
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
      dphids = p25/torflux_edge
      
      lu_e(1:ndim:ns) = 0; lv_e(1:ndim:ns) = 0
!RPRES      lu_o(1:ndim:ns) = 0
      guu(1:ndim:ns)  = 0; guv(1:ndim:ns)  = 0; gvv(1:ndim:ns) = 0
      guus = guu*shalf;    guvs = guv*shalf;    gvvs = gvv*shalf

CDIR$ IVDEP
      do l = 1, ndim
         armn_e(l)  = ohs*armn_e(l) * lu_e(l)
         azmn_e(l)  =-ohs*azmn_e(l) * lu_e(l)
         brmn_e(l)  = brmn_e(l) * lu_e(l)
         bzmn_e(l)  =-bzmn_e(l) * lu_e(l)
         bsqr(l)    = phip(l)*lu_e(l)/shalf(l)
      end do
      
CDIR$ IVDEP
      do l = 1, ndim
         armn_o(l)  = armn_e(l) *shalf(l)
         azmn_o(l)  = azmn_e(l) *shalf(l)
         brmn_o(l)  = brmn_e(l) *shalf(l)
         bzmn_o(l)  = bzmn_e(l) *shalf(l)
      end do   
! 
!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
CDIR$ IVDEP
      do l = 1, nrzt
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = dphids*(bsqr(l) + bsqr(l+1))
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
!RPRES         presg(l)= ohs*(presg(l+1) - presg(l))*
!RPRES     1             (r1(l,0) + sqrts(l)*r1(l,1))                         !R*dp/ds
      end do

CDIR$ IVDEP
      do l = 1, nrzt
         armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
     1             - gvv(l)*r1(l,0)
!RPRES     2             + presg(l)*zu0(l)
         azmn_e(l) = azmn_e(l+1) - azmn_e(l)
!RPRES                   - presg(l)*ru0(l)
         brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
         bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
      end do

CDIR$ IVDEP
      do l = 1, nrzt
        armn_e(l) = armn_e(l) - gvvs(l)*r1(l,1)
        brmn_e(l) = brmn_e(l) + bsqr(l)*z1(l,1) 
     1            - guus(l)*ru(l,1) - guu(l)*ru(l,0) 
        bzmn_e(l) = bzmn_e(l) - bsqr(l)*r1(l,1) 
     1            - guus(l)*zu(l,1) - guu(l)*zu(l,0)
      end do  
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT .5*(RU-odd*Zodd - ZU-odd*Rodd) TERM
!           IN JACOBIAN
      
CDIR$ IVDEP
      do l = 1, nrzt
         armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l) 
     1             + p5*(lv_e(l)*shalf(l) + lv_e(l+1)*shalf(l+1)) 
         azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l) 
         brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
         bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))         
      end do

!RPRES      presg(:nrzt)  = presg(:nrzt)*sqrts(:nrzt)
!RPRES      armn_o(:nrzt) = armn_o(:nrzt) + presg(:nrzt)*zu0(:nrzt)
!RPRES      azmn_o(:nrzt) = azmn_o(:nrzt) - presg(:nrzt)*ru0(:nrzt)

CDIR$ IVDEP
      do l = 1, nrzt
         lu_o(l)   = dphids*(lu_e(l)*phip(l) + lu_e(l+1)*phip(l+1))
         guu(l)    = guu(l)*sqrts(l)*sqrts(l)
         bsqr(l)   = gvv(l)*sqrts(l)*sqrts(l)
         armn_o(l) = armn_o(l) - zu(l,1)*lu_o(l)
     1             - bsqr(l)*r1(l,1) - gvvs(l)*r1(l,0)
         azmn_o(l) = azmn_o(l) + ru(l,1)*lu_o(l)
         brmn_o(l) = brmn_o(l) + z1(l,1)*lu_o(l)
     1             - guu(l)*ru(l,1) - guus(l)*ru(l,0)
         bzmn_o(l) = bzmn_o(l) - r1(l,1)*lu_o(l)
     1             - guu(l)*zu(l,1) - guus(l)*zu(l,0)
      end do      

      if (lthreed) then
CDIR$ IVDEP
         do l = 1, nrzt
            guv(l)  = p5*(guv(l) + guv(l+1))
            guvs(l) = p5*(guvs(l) + guvs(l+1))
            brmn_e(l) = brmn_e(l) - guv(l)*rv(l,0) - guvs(l)*rv(l,1)
            bzmn_e(l) = bzmn_e(l) - guv(l)*zv(l,0) - guvs(l)*zv(l,1)
            crmn_e(l) = guv(l) *ru(l,0) + gvv(l) *rv(l,0) 
     1                + gvvs(l)*rv(l,1) + guvs(l)*ru(l,1)
            czmn_e(l) = guv(l) *zu(l,0) + gvv(l) *zv(l,0)
     1                + gvvs(l)*zv(l,1) + guvs(l)*zu(l,1)
         end do

CDIR$ IVDEP
         do l = 1, nrzt
            guv(l) = guv(l) *sqrts(l)*sqrts(l)
            brmn_o(l) = brmn_o(l) - guvs(l)*rv(l,0) - guv(l)*rv(l,1)
            bzmn_o(l) = bzmn_o(l) - guvs(l)*zv(l,0) - guv(l)*zv(l,1)
            crmn_o(l) = guvs(l)*ru(l,0) + gvvs(l)*rv(l,0) 
     1                + bsqr(l)*rv(l,1) + guv(l) *ru(l,1) 
            czmn_o(l) = guvs(l)*zu(l,0) + gvvs(l)*zv(l,0) 
     1                + bsqr(l)*zv(l,1) + guv(l) *zu(l,1) 
         end do
      endif
!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      if (ivac .ge. 1) then
        lk = 0
CDIR$ IVDEP
        do l = ns,nrzt,ns
          lk = lk+1
          armn_e(l) = armn_e(l) + zu0(l)*rbsq(lk)
          azmn_e(l) = azmn_e(l) - ru0(l)*rbsq(lk)
          armn_o(l) = armn_o(l) + zu0(l)*rbsq(lk)
          azmn_o(l) = azmn_o(l) - ru0(l)*rbsq(lk)
        end do  
      endif

      deallocate (bsqr, gvvs, guvs, guus, stat=l)
!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
CDIR$ IVDEP
      do l = 1,nrzt
         rcon1   = (rcon(l,0) - rcon0(l)) * gcon(l,0)
         zcon1   = (zcon(l,0) - zcon0(l)) * gcon(l,0)
         brmn_e(l) = brmn_e(l) + rcon1
         bzmn_e(l) = bzmn_e(l) + zcon1
         brmn_o(l) = brmn_o(l)+ rcon1*sqrts(l)
         bzmn_o(l) = bzmn_o(l)+ zcon1*sqrts(l)
         rcon(l,0) =  ru0(l) * gcon(l,0)
         zcon(l,0) =  zu0(l) * gcon(l,0)
         rcon(l,1) = rcon(l,0) * sqrts(l)
         zcon(l,1) = zcon(l,0) * sqrts(l)
      end do   
 
      end subroutine forces

      subroutine tomnspa(frzl_array, armn, brmn, crmn, azmn, bzmn, 
     1   czmn, blmn, clmn, arcon, azcon)
      use vmec_main
      use vmec_params, only: jlam, jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: frzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: frcs = 3, frsc = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: fzcc, fzss,  flcc, flss
      integer :: jmax, m, mparity, i, jk, n, k, js
      real(rprec), dimension(ns,nzeta,12) :: work1
      real(rprec), dimension(ns*nzeta) :: temp1, temp3
C-----------------------------------------------
 
      fzcc = frcs + ntmax
      fzss = frsc + ntmax
      flcc = frcs + 2*ntmax
      flss = frsc + 2*ntmax
 
      jmax = ns
      if (ivac .lt. 1) jmax = ns1
!
!     BEGIN INVERSE FOURIER TRANSFORM
!     DO THETA (U) INTEGRATION FIRST
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = zero
         do i = 1, ntheta2
            do jk = 1, ns*nzeta
               temp1(jk) = armn(jk,i,mparity) + xmpq(m,1)*arcon(jk,i,
     1            mparity)
               temp3(jk) = azmn(jk,i,mparity) + xmpq(m,1)*azcon(jk,i,
     1            mparity)
               work1(jk,1,3) = work1(jk,1,3) + temp1(jk)*sinmui(i,m) + 
     1            brmn(jk,i,mparity)*cosmumi(i,m)
               work1(jk,1,5) = work1(jk,1,5) + temp3(jk)*cosmui(i,m) + 
     1            bzmn(jk,i,mparity)*sinmumi(i,m)
               work1(jk,1,9) = work1(jk,1,9) + blmn(jk,i,mparity)*
     1            sinmumi(i,m)
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  work1(jk,1,1) = work1(jk,1,1) + temp1(jk)*cosmui(i,m)
     1                + brmn(jk,i,mparity)*sinmumi(i,m)
                  work1(jk,1,7) = work1(jk,1,7) + temp3(jk)*sinmui(i,m)
     1                + bzmn(jk,i,mparity)*cosmumi(i,m)
                  work1(jk,1,11) = work1(jk,1,11) + blmn(jk,i,mparity)*
     1               cosmumi(i,m)
                  work1(jk,1,2) = work1(jk,1,2) - crmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,4) = work1(jk,1,4) - crmn(jk,i,mparity)*
     1               sinmui(i,m)
                  work1(jk,1,6) = work1(jk,1,6) - czmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,8) = work1(jk,1,8) - czmn(jk,i,mparity)*
     1               sinmui(i,m)
                  work1(jk,1,10) = work1(jk,1,10) - clmn(jk,i,mparity)*
     1               cosmui(i,m)
                  work1(jk,1,12) = work1(jk,1,12) - clmn(jk,i,mparity)*
     1               sinmui(i,m)
               end do
            endif
         end do
!
!        NEXT, DO ZETA (V) INTEGRATION
!
         do n = 0, ntor
            do k = 1, nzeta
 
               do js = jmin2(m), jmax
                  frzl_array(js,n,m,frsc) = frzl_array(js,n,m,frsc) + 
     1               work1(js,k,3)*cosnv(k,n)
                  frzl_array(js,n,m,fzcc) = frzl_array(js,n,m,fzcc) + 
     1               work1(js,k,5)*cosnv(k,n)
               end do
 
               do js = jlam(m), ns
                  frzl_array(js,n,m,flcc) = frzl_array(js,n,m,flcc) + 
     1               work1(js,k,9)*cosnv(k,n)
               end do
 
               if (lthreed) then
                  do js = jmin2(m), jmax
                     frzl_array(js,n,m,frsc) = frzl_array(js,n,m,frsc)
     1                   + work1(js,k,4)*sinnvn(k,n)
                     frzl_array(js,n,m,fzcc) = frzl_array(js,n,m,fzcc)
     1                   + work1(js,k,6)*sinnvn(k,n)
                     frzl_array(js,n,m,frcs) = frzl_array(js,n,m,frcs)
     1                   + work1(js,k,1)*sinnv(k,n) + work1(js,k,2)*
     2                  cosnvn(k,n)
                     frzl_array(js,n,m,fzss) = frzl_array(js,n,m,fzss)
     1                   + work1(js,k,7)*sinnv(k,n) + work1(js,k,8)*
     2                  cosnvn(k,n)
                  end do
                  do js = jlam(m), ns
                     frzl_array(js,n,m,flcc) = frzl_array(js,n,m,flcc)
     1                   + work1(js,k,10)*sinnvn(k,n)
                     frzl_array(js,n,m,flss) = frzl_array(js,n,m,flss)
     1                   + work1(js,k,11)*sinnv(k,n) + work1(js,k,12)*
     2                  cosnvn(k,n)
                  end do
               endif
 
            end do
         end do
      end do

      end subroutine tomnspa

      subroutine tomnsps(frzl_array, armn, brmn, crmn, azmn, bzmn, 
     1   czmn, blmn, clmn, arcon, azcon)
      use vmec_main
      use vmec_params, only: jlam, jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(out) :: frzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: frcc = 1, frss = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: fzcs, fzsc, flcs, flsc
      integer :: jmax, m, mparity, i, n, k, js, l
      real(rprec), dimension(ns*nzeta,12) :: work1
      real(rprec), dimension(ns*nzeta) :: temp1, temp3
C-----------------------------------------------

      fzcs = frcc+ntmax
      fzsc = frss+ntmax
      flcs = frcc+2*ntmax
      flsc = frss+2*ntmax

      frzl_array = zero
      
      jmax = ns
      if (ivac .lt. 1) jmax = ns1

!
!     BEGIN INVERSE FOURIER TRANSFORM
!     DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
!
!       FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
!       FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
!       FLmn =      - d(BLmn)/du + d(CLmn)/dv
!
!       NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = zero
         do i = 1, ntheta2
            temp1(:) = armn(:,i,mparity) + xmpq(m,1)*
     1            arcon(:,i,mparity)
            temp3(:) = azmn(:,i,mparity) + xmpq(m,1)*
     1            azcon(:,i,mparity)
            work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m) + 
     1            brmn(:,i,mparity)*sinmumi(i,m)
            work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m) + 
     1            bzmn(:,i,mparity)*cosmumi(i,m)
            work1(:,11) = work1(:,11) + blmn(:,i,mparity)*
     1            cosmumi(i,m)
            if (lthreed) then
               work1(:,2) = work1(:,2) - crmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,4) = work1(:,4) - crmn(:,i,mparity)*
     1               sinmui(i,m)
               work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m)
     1                + brmn(:,i,mparity)*cosmumi(i,m)
               work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m)
     1                + bzmn(:,i,mparity)*sinmumi(i,m)
               work1(:,6) = work1(:,6) - czmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,8) = work1(:,8) - czmn(:,i,mparity)*
     1               sinmui(i,m)
               work1(:,9) = work1(:,9) + blmn(:,i,mparity)*
     1               sinmumi(i,m)
               work1(:,10) = work1(:,10) - clmn(:,i,mparity)*
     1               cosmui(i,m)
               work1(:,12) = work1(:,12) - clmn(:,i,mparity)*
     1               sinmui(i,m)
            endif
         end do
!
!                NEXT, DO ZETA (V) INTEGRATION
!
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = jmin2(m), jmax
                  frzl_array(js,n,m,frcc) = frzl_array(js,n,m,frcc) + 
     1               work1(js+l,1)*cosnv(k,n)
                  frzl_array(js,n,m,fzsc) = frzl_array(js,n,m,fzsc) + 
     1               work1(js+l,7)*cosnv(k,n)
               end do
               do js = jlam(m), ns
                  frzl_array(js,n,m,flsc) = frzl_array(js,n,m,flsc) + 
     1               work1(js+l,11)*cosnv(k,n)
               end do
 
               if (lthreed) then
                  do js = jmin2(m), jmax
                     frzl_array(js,n,m,frcc) = frzl_array(js,n,m,frcc)
     1                   + work1(js+l,2)*sinnvn(k,n)
                     frzl_array(js,n,m,fzsc) = frzl_array(js,n,m,fzsc)
     1                   + work1(js+l,8)*sinnvn(k,n)
                     frzl_array(js,n,m,frss) = frzl_array(js,n,m,frss)
     1                   + work1(js+l,3)*sinnv(k,n) + work1(js+l,4)*
     2                  cosnvn(k,n)
                     frzl_array(js,n,m,fzcs) = frzl_array(js,n,m,fzcs)
     1                   + work1(js+l,5)*sinnv(k,n) + work1(js+l,6)*
     2                  cosnvn(k,n)
                  end do
                  do js = jlam(m), ns
                     frzl_array(js,n,m,flsc) = frzl_array(js,n,m,flsc)
     1                   + work1(js+l,12)*sinnvn(k,n)
                     frzl_array(js,n,m,flcs) = frzl_array(js,n,m,flcs)
     1                   + work1(js+l,9)*sinnv(k,n) + work1(js+l,10)*
     2                  cosnvn(k,n)
                  end do
               endif
 
            end do
         end do
      end do

      end subroutine tomnsps

      subroutine jacobian
      use vmec_main, only: ohs, nrzt, irst
      use vmec_params, only: meven, modd
      use realspace
      use vmec_dim, only: ns
      use vforces, r12 => armn_o, ru12 => azmn_e, z12 => blmn_e,
     1    zu12 => armn_e, rs => bzmn_e, zs => brmn_e, gsqrt => azmn_o
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec) :: taumax, taumin, dphids
C-----------------------------------------------
 
!
!     (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U)
!     AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
!     NOTE: LOOPS WERE SPLIT TO ALLOW EFFICIENT MEMORY BUS USAGE
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING MORE THAN ONE POINTER!
!
!
!     HERE, SQRT(G) = R( Ru * Zs - Rs * Zu ). THE DERIVATIVES OF SHALF = SQRT(PHI(s))
!     WERE COMPUTED EXPLICITLY AS: d Shalf/ds = .5/shalf * d(PHI)/ds, WHERE
!     d(PHI)/ds = phip(s)/torflux_edge
!
!
      dphids = p25/torflux_edge
      
CDIR$ IVDEP
      do l = 2,nrzt
        ru12(l) = p5*(ru(l,meven) + ru(l-1,meven) + 
     1       shalf(l)*(ru(l,modd)  + ru(l-1,modd)))
        zs(l)   = ohs*(z1(l,meven) - z1(l-1,meven) +
     1       shalf(l)*(z1(l,modd)  - z1(l-1,modd)))
        z12(l)  = p5*(z1(l,meven) + z1(l-1,meven) +
     1       shalf(l)*(z1(l,modd)  + z1(l-1,modd)))
        gsqrt(l) = ru12(l)*zs(l) + dphids*phip(l)*
     1  (ru(l,modd) *z1(l,modd) + ru(l-1,modd) *z1(l-1,modd) +
     2  (ru(l,meven)*z1(l,modd) + ru(l-1,meven)*z1(l-1,modd))/shalf(l))
      enddo


CDIR$ IVDEP
      do l = 2,nrzt
        zu12(l) = p5*(zu(l,meven) + zu(l-1,meven) +
     1       shalf(l)*(zu(l,modd)  + zu(l-1,modd)))
        rs(l)   = ohs*(r1(l,meven) - r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  - r1(l-1,modd)))
        r12(l)  = p5*(r1(l,meven) + r1(l-1,meven) +
     1       shalf(l)*(r1(l,modd)  + r1(l-1,modd)))
        gsqrt(l) = r12(l)*(gsqrt(l)- rs(l)*zu12(l) - dphids*phip(l)*
     1    (zu(l,modd) *r1(l,modd)+zu(l-1,modd) *r1(l-1,modd)
     2  + (zu(l,meven)*r1(l,modd)+zu(l-1,meven)*r1(l-1,modd))/shalf(l)))
      end do

!
!     TEST FOR SIGN CHANGE IN JACOBIAN
!
      gsqrt(1:nrzt:ns) = gsqrt(2:nrzt:ns)
      taumax = maxval(gsqrt(2:nrzt))
      taumin = minval(gsqrt(2:nrzt))
      if (taumax*taumin .lt. zero) irst = 2

      end subroutine jacobian
      

      subroutine allocate_funct3d
      use vmec_main
      use realspace
      use vforces
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1, ndim, ndim2
C-----------------------------------------------

      ndim  = 1+nrzt
      ndim2 = 2*ndim

      call free_mem_funct3d

      allocate( armn(ndim2), azmn(ndim2), brmn(ndim2), bzmn(ndim2), 
     1   crmn(ndim2), czmn(ndim2), blmn(ndim2), clmn(ndim2), 
     2   r1(nrzt,0:1), ru(nrzt,0:1), rv(nrzt,0:1),
     3   z1(nrzt,0:1), zu(nrzt,0:1), zv(nrzt,0:1),
     4   rcon(nrzt,0:1), zcon(nrzt,0:1), ru0(ndim), zu0(ndim),
     6   rcon0(ndim), zcon0(ndim), guu(ndim), guv(ndim), gvv(ndim),
     7   stat=istat1 )
      if (istat1.ne.0) stop 'allocation error #1 in funct3d'

      if (lfreeb) then
      allocate (brv(nznt), bphiv(nznt), bzv(nznt), 
     1   bpolvac(nznt), bsqvac(nznt), stat=istat1) 
      if (istat1.ne.0) stop 'allocation error #2 in funct3d'
      end if
!
!     Pointer alias assignments 
!     NOTE: In FORCES, X_e(nrzt+1) overlaps X_o(1), which should never be used...
!      
      armn_e => armn(:ndim)
      armn_o => armn(ndim:)
      armn(:ndim2) = zero
      brmn_e => brmn(:ndim)
      brmn_o => brmn(ndim:)
      brmn(:ndim2) = zero
      azmn_e => azmn(:ndim)
      azmn_o => azmn(ndim:)
      azmn(:ndim2) = zero
      bzmn_e => bzmn(:ndim)
      bzmn_o => bzmn(ndim:)
      bzmn(:ndim2) = zero
      crmn_e => crmn(:ndim)
      crmn_o => crmn(ndim:)
      crmn(:ndim2) = zero
      czmn_e => czmn(:ndim)
      czmn_o => czmn(ndim:)
      czmn(:ndim2) = zero
      blmn_e => blmn(:ndim)
      blmn_o => blmn(ndim:)
      blmn(:ndim2) = zero
      clmn_e => clmn(:ndim)
      clmn_o => clmn(ndim:)
      clmn(:ndim2) = zero
      rcon0(:ndim) = zero
      zcon0(:ndim) = zero
                  
      end subroutine allocate_funct3d
      

      subroutine allocate_ns (lreset, linterp, neqs2_old)
      use vmec_main
      use vmec_params, only: ntmax
      use realspace
      use vsvd
      use vspline
      use vforces
      use xstuff
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      integer neqs2_old
      logical :: lreset, linterp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ndim, nsp1, istat1
      real(rprec), dimension(:), allocatable :: xc_old, scalxc_old
      real(rprec) delr_mse
C-----------------------------------------------
!
!     FIRST STORE COARSE-MESH XC FOR INTERPOLATION
!

      ndim  = 1 + nrzt
      nsp1  = 1 + ns
      delr_mse = zero

!
!     Make consistency checks (just for extra protection)
!
      if (linterp .and. (neqs2_old .eq. neqs2)) then
         print *,' Setting linterp = F in allocate_ns. '
         linterp = .false.
      end if

!
!     Save old xc, scalxc for possible interpolation or if iterations restarted on same mesh...
!            
      if (neqs2_old .gt. 0 .and. allocated(scalxc) .and. linterp) then
         allocate(xc_old(neqs2_old), scalxc_old(neqs2_old), stat=istat1)
         if (istat1.ne.0) stop 'allocation error #1 in allocate_ns'
         xc_old(:neqs2_old) = xc(:neqs2_old)
         scalxc_old(:neqs2_old) = scalxc(:neqs2_old)
         if (lrecon) delr_mse = xc(neqs2_old)
      end if

!
!     ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
!     FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
!
      call free_mem_ns (lreset)

      allocate (phip(ndim), shalf(ndim), sqrts(ndim), wint(ndim), 
     1  stat=istat1)
      if (istat1.ne.0) stop 'allocation error #2 in allocate_ns'
      allocate( ireflect(ns*nzeta), indexr(2*ns) ,imid(2*ns),  
     1  stat=istat1)
      if (istat1.ne.0) stop 'allocation error #3 in allocate_ns'
      allocate( current(ns),rm2(ns),vrm2(ns),
     1  ovrm2(ns), ochip(ns), presph(ns), presint(ns),
     2  w_ia(ns), w1_ia(ns), u_ia(ns), u1_ia(ns),
     3  w_pa(ns), w1_pa(ns), u_pa(ns), u1_pa(ns),
     4  w_ib(ns), w1_ib(ns), u_ib(ns), u1_ib(ns),
     5  w_pb(ns), w1_pb(ns), u_pb(ns), u1_pb(ns),
     6  rmid(2*ns),datamse(2*ns),qmid(2*ns),
     7  shear(2*ns),presmid(2*ns),alfa(2*ns),curmid(2*ns),
     8  curint(2*ns),psimid(2*ns),ageo(2*ns),volpsi(2*ns),
     9  isplinef(ns),isplineh(ns),psplinef(ns),psplineh(ns),
     A  phimid(2*ns),pm(ns,0:nobser+nobd),
     B  im(ns,0:nobser+nobd),stat=istat1)
      if (istat1.ne.0) stop 'allocation error #4 in allocate_ns'

      if (nbsets.gt.0)
     1  allocate( pmb(ns,0:nbcoil_max,nbsets,2),
     2            imb(ns,0:nbcoil_max,nbsets,2),stat=istat1)
      if (istat1.ne.0) stop 'allocation error #5 in allocate_ns'

      allocate( ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),
     1          azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2), bzm(nsp1,2),
     2          sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns),
     3          stat=istat1)
      if (istat1.ne.0) stop 'allocation error #6 in allocate_ns'

      allocate( iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),
     1          jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),
     2          bdotgradv(ns), equif(ns), specw(ns), tcon(ns), fpsi(ns),
     3          psi(ns),yellip(ns),yinden(ns), ytrian(ns),yshift(ns),
     4          ygeo(ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),
     5          iotas(nsp1), phips(nsp1), pres(nsp1), vp(nsp1), 
     6          beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),
     7          blam(nsp1), clam(nsp1), dlam(nsp1), 
     7          stat=istat1)
      if (istat1.ne.0) stop 'allocation error #7 in allocate_ns'

      allocate( rmidx(2*ns), hmidx(2*ns), wmidx(2*ns), qmidx(2*ns),
     1          tenmidx(2*ns), ymidx(2*ns), y2midx(2*ns), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #8 in allocate_ns'
     
      allocate (gc(neqs2), xcdot(neqs2), xstore(neqs2), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #9 in allocate_ns'
      
      if (lreset) then
         allocate (xc(neqs2), scalxc(neqs2), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #10 in allocate_ns'
         xc(:neqs2) = zero
      end if         
   
      if (allocated(xc_old)) then
         xstore(1:neqs2_old) = xc_old(1:neqs2_old)
         scalxc(1:neqs2_old) = scalxc_old(1:neqs2_old)
         deallocate (xc_old, scalxc_old)         
      end if
      

      xc(neqs2) = delr_mse

!
!     Allocate nrzt-dependent arrays (persistent) for funct3d
!
      call allocate_funct3d

      end subroutine allocate_ns
      

      subroutine allocate_nunv
      use vmec_main
      use vmec_params, only: ntmax
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1
C-----------------------------------------------

      call free_mem_nunv

      allocate (bsubu0(nznt), rbsq(nznt), dbsq(nznt), stat=istat1) 
      if (istat1.ne.0) stop 'allocation error #1 in allocate_nunv'

      allocate (rmn_bdy(0:ntor,0:mpol1,ntmax),
     1          zmn_bdy(0:ntor,0:mpol1,ntmax), stat=istat1)
      if (istat1.ne.0) stop 'allocation error #2 in allocate_nunv'

!     PERSISTENT ARRAYS (DURATION OF PROGRAM)
      if (lfreeb) 
     1   allocate (amatsav(mnpd2*mnpd2),bvecsav(mnpd2), 
     2          bsqsav(nznt,3), potvac(2*mnpd), stat=istat1)
         if (istat1.ne.0) stop 'allocation error #3 in allocate_nunv'

      end subroutine allocate_nunv
        

      function aspectratio ()
      use vmec_main
      use realspace
      use vmec_io
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lk, l
      real(rprec) :: rb, zub, pi, t1, aspectratio
C-----------------------------------------------
!
!     routine for computing aspect-ratio
!     A = <R>/<a>
!
!     where pi <a>**2 = Area (toroidally averaged)
!           2*pi * <R> * Area = Volume
!     Use integration by parts to compute as surface integral (Stoke''s theorem)
!
 
      pi = 4*atan(one)
 
!
!     Compute Volume and Area
!
      volume_plasma = 0
      cross_area = 0
      do lk = 1, nznt
         l = ns*lk
         rb  = r1(l,0) + r1(l,1)
         zub = zu(l,0) + zu(l,1)
         t1  = rb*zub*wint(l)
         volume_plasma = volume_plasma + rb*t1
         cross_area = cross_area + t1
      end do
 
      volume_plasma = 2*pi*pi*abs(volume_plasma)
      cross_area = 2*pi*abs(cross_area)
 
      Rmajor_p = volume_plasma/(2*pi*cross_area)
      Aminor_p = sqrt(cross_area/pi)

      aspectratio = Rmajor_p/Aminor_p

      end function aspectratio 

                 
      subroutine open_output_files (extension, iseq, lmac, lscreen,
     1           lfirst)
      use safe_open_mod
      use vparams, only: nmac, nthreed, nmac0, nthreed0
      implicit none      
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq
      character*(*) :: extension
      logical :: lmac, lscreen, lfirst
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iread, inthreed=0, imac0=0
      character*120 :: mac_file, threed1_file
C-----------------------------------------------
 
!
!     OPEN FILES FOR READING, WRITING
!
      threed1_file = 'threed1.'//extension
      mac_file = 'mac.'//extension
 
      if (lscreen .and. lfirst) write (*, '(33('' -''))')
      nthreed = nthreed0
      call safe_open(nthreed, iread, threed1_file, 'new', 'formatted')
      if (iread .ne. 0) then
         if (iseq .eq. 0 .and. lscreen .and. lfirst) print *, 
     1   ' VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...'
         call safe_open(nthreed, inthreed, threed1_file, 'replace',
     1     'formatted')
      endif
   
 
      nmac = max(nmac0, nthreed)
      if (lmac) then
         call safe_open(nmac, imac0, mac_file, 'replace', 'formatted')
      end if
      if (inthreed.ne.0 .or. imac0.ne.0) then
         print *,' nthreed = ', nthreed, ' istat_threed = ', inthreed,
     1           ' nmac0   = ', nmac,' istat_mac0 = ', imac0 
         print *, 'Error opening output file in VMEC OPEN_OUTPUT_FILES'
         stop 10
      endif

      end subroutine open_output_files


      subroutine close_all_files
      use vparams, only: nmac, nthreed
      implicit none
C-----------------------------------------------

      if (nthreed .gt. 0) close (nthreed)
      if (nmac .gt. 0) close (nmac)

      end subroutine close_all_files
      

      subroutine storesvd(re, ro, lue, luo, lve, phipog, zu00)
      use vmec_main
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nzeta,*) :: re, ro, lue
      real(rprec), dimension(ns,*) :: luo
      real(rprec), dimension(*) :: lve, phipog
      real(rprec), dimension(ns,nzeta,*) :: zu00
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lk, l, noff, jmin, lt, js, j, ks, j1, j2
      real(rprec), dimension(ns,ntheta2) :: luef, luof
      real(rprec) :: signi, phipog0, fpsie, w1,
     1    fpsi2, fpsi1, diota
C-----------------------------------------------
!
!     COMPUTE FULL MESH VALUES OF LAMBDA AT THETA=0,PI
!
      call lamfull (luef, luof, lue, luo)
!
!     COMPUTE PHIP*GUU*(IOTA-LV)/GSQRT AS S->0 FOR CURRENT(0)
!
        do lk = 1,nznt
          l = 2+ns*(lk-1)
          phipog0 = 1.5_dp*phipog(l) - 0.5_dp*phipog(l+1)
          bsubu0(lk) = phipog0*( (iotas(2) + lve(l))*guu(l)
     1    + shalf(2)*luo(2,lk)*guv(l) )
        enddo
 
      if (.not.lrecon) return
      
      signi = sign(one,iotaf(2))
      noff = ntheta2
      jmin = 2
      do lt = 1, 2
         j2   = ns*nzeta*(noff-1)
         if (lt .eq. 1) then
            l = ns+1-jmin
            imid(l:1:(-1)) = (/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:1:(-1)) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:1:(-1)) = atan(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:1:(-1))*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:1:(-1)) = signi/iotaf(jmin:ns)
            presmid(l:1:(-1)) = presf(jmin:ns)/dmu0
         else
            l = ns+jmin-1
            imid(l:2*ns-1)=(/(j1,j1=jmin+j2,ns+j2)/)
            rmid(l:2*ns-1) = re(jmin:ns,1,noff) + sqrts(jmin:ns)
     1         *ro(jmin:ns,1,noff)
            datamse(l:2*ns-1) = atan(iotaf(jmin:ns)*zu00(jmin:ns
     1         ,1,noff)/(rmid(l:2*ns-1)*(luef(jmin:ns,noff)+
     2         sqrts(jmin:ns)*luof(jmin:ns,noff))))/dcon
            qmid(l:2*ns-1) = signi/iotaf(jmin:ns)
            presmid(l:2*ns-1) = presf(jmin:ns)/dmu0
         endif
         noff = 1
         jmin = 1
      end do
 
      call findphi (re, ro, rthom, delse2, delso2, rmid, indexs2, 
     1   indexu2, indexr, itse)
      pcalc(:itse) = (presf(indexs2(:itse))) + abs(delse2(:itse))*(presf
     1   (indexs2(:itse)+1)-(presf(indexs2(:itse))))
 
!
!       SORT ON RSORT(KS) ARRAY
!       INCLUDE EDGE PITCH MATCH TO TOTAL CURRENT
!
      fpsie = 1.5_dp*fpsi(ns) - 0.5_dp*fpsi(ns1)
 
      do j = 1, imse2
         ks = isorts(j)
         js = indexs1(ks)
         w1 = delse1(ks)
         if (js .eq. ns) then
            fpsi2 = fpsie
            fpsi1 = fpsie
         else if (js .eq. ns1) then
            fpsi2 = fpsie
         else
            fpsi2 = 0.5_dp*(fpsi(js+1)+fpsi(js+2))
         endif
         if (js.lt.ns .and. js.gt.1) fpsi1 = .5_dp*(fpsi(js)+fpsi(js+1))
         if (js .eq. 1) fpsi1 = fpsi(2)
         diota = (one - w1)*iotaf(js) + w1*iotaf(js+1)
         fpsical(ks) = (one - w1)*fpsi1 + w1*fpsi2
         if (stark_weight(ks) .ne. zero) qmeas(ks) = datastark(ks)/
     1      stark_weight(ks)
         qcalc(ks) = diota
         starkcal(ks) = diota*stark_weight(ks)
         rsort0(ks) = sqrt(hs*(js - 1 + w1))
      end do

      end subroutine storesvd

      
      subroutine bextrema(modb, bmin, bmax, nzeta, ntheta)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nzeta, ntheta
      real(rprec) :: modb(nzeta,ntheta), bmin(ntheta), bmax(ntheta)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer ku
C-----------------------------------------------
!
!     Computes max, min of |B| along v (zeta) between two angle lines (theta = 0,pi)
!
      do ku = 1,ntheta
         bmin(ku)  = minval(modb(:,ku))
         bmax(ku)  = maxval(modb(:,ku))
      enddo

      end subroutine bextrema


      subroutine convert(rmnc,zmns,lmns,rmns,zmnc,lmnc,rzl_array,js)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer js
      real(rprec), dimension(mnmax), intent(out) ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1    intent(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncc=1
      integer, parameter :: rmnss=2
      integer, parameter :: rmncs=3
      integer, parameter :: rmnsc=4
      real(rprec), parameter :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncs, zmnsc, zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
      integer :: mn, m, n, n1
      real(rprec) :: t1, sign0
C-----------------------------------------------
!
!     CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
!     FORM FOR OUTPUT (COEFFICIENTS OF cos(mu-nv), sin(mu-nv))
!
      zmncs = rmncc + ntmax
      lmncs = rmncc + 2*ntmax
      zmnsc = rmnss + ntmax
      lmnsc = rmnss + 2*ntmax
      zmncc = rmncs + ntmax
      lmncc = rmncs + 2*ntmax
      zmnss = rmnsc + ntmax
      lmnss = rmnsc + 2*ntmax

      mn = 0
!
!     DO M = 0 MODES SEPARATELY (ONLY KEEP N >= 0 HERE: COS(-NV), SIN(-NV))
!
      m = 0
      do n = 0, ntor
         t1 = mscale(m)*nscale(n)
         mn = mn + 1
         rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
         zmns(mn) =-t1*rzl_array(js,n,m,zmncs)
         lmns(mn) =-t1*rzl_array(js,n,m,lmncs)
         if (lasym) then
         rmns(mn) = t1*rzl_array(js,n,m,rmncs)
         zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
         lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
         end if
      end do

      if (js .eq. 1) then
         mn = 0
         do n = 0, ntor
            t1 = mscale(m)*nscale(n)
            mn = mn + 1
            lmns(mn) =-t1*(2._dp*rzl_array(2,n,m,lmncs)
     1               -           rzl_array(3,n,m,lmncs))            
         end do   
      end if
       
      do m = 1, mpol1
         do n = -ntor, ntor
            n1 = abs(n)
            t1 = mscale(m)*nscale(n1)
            mn = mn + 1
            if (n .eq. 0) then
               rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
               zmns(mn) = t1*rzl_array(js,n,m,zmnsc)
               lmns(mn) = t1*rzl_array(js,n,m,lmnsc)
               if (lasym) then
               rmns(mn) = t1*rzl_array(js,n,m,rmnsc)
               zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
               lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
               end if
            else if (js .gt. 1) then
               sign0 = n/n1
               rmnc(mn) = p5*t1*(rzl_array(js,n1,m,rmncc)+sign0*
     1            rzl_array(js,n1,m,rmnss))
               zmns(mn) = p5*t1*(rzl_array(js,n1,m,zmnsc)-sign0*
     1            rzl_array(js,n1,m,zmncs))
               lmns(mn) = p5*t1*(rzl_array(js,n1,m,lmnsc)-sign0*
     1            rzl_array(js,n1,m,lmncs))
               if (lasym) then
               rmns(mn) = p5*t1*(rzl_array(js,n1,m,rmncs)-sign0*
     1            rzl_array(js,n1,m,rmnsc))
               zmnc(mn) = p5*t1*(rzl_array(js,n1,m,zmncc)+sign0*
     1            rzl_array(js,n1,m,zmnss))
               lmnc(mn) = p5*t1*(rzl_array(js,n1,m,lmncc)+sign0*
     1            rzl_array(js,n1,m,lmnss))
               end if
            else if (js .eq. 1) then
               rmnc(mn) = zero
               zmns(mn) = zero
               lmns(mn) = zero
               if (lasym) then
               rmns(mn) = zero
               zmnc(mn) = zero
               lmnc(mn) = zero
               end if
            end if
         end do
      end do

      if( .not. lasym) then
         rmns = 0
         zmnc = 0
         lmnc = 0
      endif
 
      end subroutine convert
      

      subroutine eqsolve(nsval, interp_flag, ier_flag, 
     1    lreset, lfirst, lscreen, reset_file_name)
      use vmec_main
      use vmec_params, only: ntmax, ns4
      use realspace
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nsval, ier_flag
      character*(*) :: reset_file_name
      logical :: interp_flag, lreset, lfirst, lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: bad_jacobian_flag = 6
      real(rprec), parameter :: p98 = 0.98_dp, p96 = 0.96_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nsold=0, neqs2_old=0
      real(rprec) :: w1, r00s, w0, res0, wdota, r0dot
      real(rprec) :: delt0
      logical :: liter_flag, lreset_internal
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        hs      radial mesh size increment
!        iequi   counter used to call -EQFOR- at end of run
!        irzloff offset in xc array between R,Z,L components
!        ijacob  counter for number of times jacobian changes sign
!        irst    counter monitoring sign of jacobian; resets R, Z, and
!                Lambda when jacobian changes sign and decreases time step
!        signgs  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
!        iterj   stores position in main iteration loop (j=1,2)
!        itfsq   counter for storing FSQ into FSQT for plotting
!        ivac    counts number of free-boundary iterations
!        ndamp   number of iterations over which damping is averaged
!        meven   parity selection label for even poloidal modes of R and Z
!        modd    parity selection label for odd poloidal modes of R and
!        gc      stacked array of R, Z, Lambda Spectral force coefficients (see readin for stack order)
!        xc      stacked array of scaled R, Z, Lambda Fourier coefficients
 
!
!       INITIALIZE MESH-DEPENDENT SCALARS
!
      ns = nsval
      ns1 = ns - 1
      hs = one/ns1
      ohs = one/hs
      mns = ns*mnsize
      irzloff = ntmax*mns
      nrzt = nznt*ns
      neqs = 3*irzloff
      neqs1 = neqs + 1
      neqs2 = neqs1 + 1
      ijacob = 0
      delt0   = delt

      fsqt = 0                              ! mcz
      wdot = 0

      write (nthreed, 10) ns, mnmax, ftolv
      if (lscreen) print 10, ns, mnmax, ftolv
   10 format(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',
     1   1pe10.3)
      if (imovephi .gt. 0 .and. lscreen) print *, 
     1   'Turning on Ip Matching by varying Phi-Edge.'
 
!
!     ALLOCATE NS-DEPENDENT ARRAYS
!
      call allocate_ns(lreset, interp_flag, neqs2_old)

      lreset_internal = (ns .gt. nsold)

      if (lreset_internal .and. neqs2_old.gt.0) gc(1:neqs2_old) =
     1  scalxc(1:neqs2_old)*xstore(1:neqs2_old)

      lreset_internal = lreset_internal .or. lreset

 1000 continue
 
      itfsq = 0
      fsq     = one
      rsfac   = one
      w1      = zero
      r00s    = zero
      gphifac = zero
      grmse   = zero
 
!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
   20 continue
      iter2 = 1
      irst = 1
      call profil1d (xc, xcdot, lreset)
      if (lfirst .and. .not.lreset)               !!Command line with lreset=F
     1    call load_xc_from_wout(xc(1), xc(1+irzloff), xc(1+2*irzloff),
     2         ntor, mpol1, ns, reset_file_name, lreset_internal)     
      call profil3d (xc(1), xc(1+irzloff), lreset_internal)
 
!
!     INTERPOLATE FROM COARSE TO NEXT FINER RADIAL GRID
!
      if (interp_flag) call interp (xc, gc, scalxc, ns, nsold)
      nsold = ns
      neqs2_old = neqs2
 
!
!     Store XC,XCDOT for possible restart
!
      call restart(delt0)
      iter1 = iter2
      liter_flag = .true.
      ier_flag = 0

!
!     ENTER FORCE ITERATION LOOP
!
      iter_loop: do while (liter_flag)
!
!     ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
!
         call evolve (delt0, ier_flag, liter_flag, lscreen)
         if (ijacob.eq.0 .and. ier_flag.eq.1) then
            if (lscreen .and. irst.eq.2) print *,
     1      ' INITIAL JACOBIAN CHANGED SIGN: IMPROVING THE',
     2      ' GUESS FOR THE MAGNETIC AXIS'
            if (lscreen .and. irst.eq.4) print *,
     1      ' INITIAL FORCES LARGE: IMPROVING JACOBIAN ESTIMATE'
            call guess_axis (r1, z1, ru0, zu0, lscreen)
            lreset_internal = .true.
            ijacob = 1
            go to 20
         else if (ier_flag .ne. 0) then
            if (ier_flag .eq. 1) ier_flag = bad_jacobian_flag
            return
         endif

         w0 = wb + wp/(gamma - one)
 
!
!     ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)
!
         if (ijacob .eq. 25) then
            irst = 2
            call restart(delt0)
            delt0 = p98*delt
            if (lscreen) print 120, delt0
            go to 1000
         else if (ijacob .eq. 50) then
            irst = 2
            call restart(delt0)
            delt0 = p96*delt
            if (lscreen) print 120, delt0
            go to 1000
         else if (ijacob .ge. 75) then
            ier_flag = 2
            liter_flag = .false.
         else if (iter2.ge.niter .and. liter_flag) then
            ier_flag = 4
            liter_flag = .false.
         endif
 
!
!       TIME STEP CONTROL
!
         if (iter2 .eq. iter1) res0 = fsq
         res0 = min(res0,fsq)
!       Store current state (irst=1)
         if (fsq .le. res0 .and. iter2-iter1 .gt. 10) then
            call restart(delt0)
!       Residuals are growing in time, reduce time step
         else if (fsq .gt. 100.0_dp*res0 .and. iter2 .gt. iter1) then
            irst = 2
         else if (iter2 - iter1 .gt. ns4/2 .and. iter2 .gt. 2*ns4 
     1        .and. fsqr+fsqz .gt. c1pm2) then
            irst = 3
         endif
 
         if (irst .ne. 1) then
!       Retrieve previous good state
            call restart(delt0)
            iter1 = iter2
         else
!       Increment time step and Printout every nstep iterations
            if (mod(iter2,nstep).eq.0 .or. iter2.eq.1 .or.
     1         .not.liter_flag) call printout(iter2, delt0, w0, lscreen)
            iter2 = iter2 + 1
         endif
 
!       Store force residual, wdot for plotting
         wdota = abs(w0 - w1)/w0
         r0dot = abs(r00 - r00s)/r00
         r00s = r00
         w1 = w0
         if (ivac.eq.1 .and. lreset) then
            if (lscreen) print 110, iter2
            write (nthreed, 110) iter2
            ivac = ivac + 1
         endif
!
!       STORE FSQ FOR PLOTTING. EVENTUALLY, STORE FOR EACH RADIAL MESH
!
         if (mod(iter2,niter/nstore_seq + 1).eq.0 .and. ns.eq.
     1      ns_array(multi_ns_grid)) then
            if (itfsq .lt. nstore_seq) then
              itfsq = itfsq + 1
              fsqt(itfsq) = fsqr + fsqz
              wdot(itfsq) = max(wdota,c1pm13)
            end if  
         end if
 
      end do iter_loop
 
      write (nthreed, 60) wdota, r0dot
      write (nthreed, 70) 1.e-6_dp*ctor/dmu0, rbtor, rbtor0
      if (lrecon) write (nthreed, 65) r00*fsqsum0/wb
 
   60 format(/,' d(ln W)/dt = ',1pe10.3,' d(ln R0)/dt = ',1pe10.3,/)
   65 format(' Average radial force balance: Int[FR(m=0)]',
     1   '/Int(B**2/R) = ',1pe12.5,' (should tend to zero)'/)
   70 format(' TOROIDAL CURRENT = ',1pe10.2,' [MA] ',' R-BTOR(s=1) = ',
     1   1pe10.2,' R-BTOR(s=0) = ',1pe10.2)
  110 format(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)
  120 format(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,
     1  /,2x,'If this does NOT resolve problem, inform vendor ',
     2       'that lamcal scaling needs to be checked')

      end subroutine eqsolve


      subroutine guess_axis(r1, z1, ru0, zu0, lscreen)
      use vmec_main
      use vmec_params, only: nscale, signgs
      use realspace, only: sqrts
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nzeta,ntheta3,0:1),
     1     intent(in) :: r1, z1
      real(rprec), dimension(ns,nzeta,ntheta3), intent(in) :: ru0, zu0
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: limpts = 61
      real(rprec), parameter :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
      real(rprec), dimension(nzeta) :: rcom, zcom
      real(rprec), dimension(ntheta1) :: r1b, z1b, rub, zub
      real(rprec), dimension(ntheta1) :: r12, z12
      real(rprec), dimension(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
      real(rprec) :: rlim(limpts), zlim
      real(rprec) :: rmax, rmin, zmax, zmin, dzeta
      real(rprec) :: ds, mintau, mintemp
C-----------------------------------------------
!
!     COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
!     LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
!     SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
!     YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
!     CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED
!
      ns12 = (ns+1)/2

      planes: do iv = 1, nzeta
         if (.not.lasym .and. iv .gt. nzeta/2+1) then
            rcom(iv) = rcom(nzeta+2-iv)
            zcom(iv) =-zcom(nzeta+2-iv)
            cycle
         end if
         r1b(:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
         z1b(:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)
         r12(:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
         z12(:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)
         rub(:ntheta3) = ru0(ns,iv,:)
         zub(:ntheta3) = zu0(ns,iv,:)
         ru12(:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
         zu12(:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

         if (.not.lasym) then
!
!     USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
!     TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)
!
         ivminus = mod(nzeta + 1 - iv,nzeta) + 1           !!(twopi-v)
         do iu = 1+ntheta2, ntheta1
            iu_r = ntheta1 + 2 - iu
            r1b(iu) = r1(ns,ivminus,iu_r,0) + r1(ns,ivminus,iu_r,1)
            z1b(iu) =-(z1(ns,ivminus,iu_r,0) + z1(ns,ivminus,iu_r,1))
            r12(iu) = r1(ns12,ivminus,iu_r,0) +
     1                r1(ns12,ivminus,iu_r,1)*sqrts(ns12)
            z12(iu) =-(z1(ns12,ivminus,iu_r,0) +
     1                z1(ns12,ivminus,iu_r,1)*sqrts(ns12))
            rub(iu) =-ru0(ns,ivminus,iu_r)
            zub(iu) = zu0(ns,ivminus,iu_r)
            ru12(iu)=-p5*(ru0(ns,ivminus,iu_r) + ru0(ns12,ivminus,iu_r))
            zu12(iu)= p5*(zu0(ns,ivminus,iu_r) + zu0(ns12,ivminus,iu_r))
         end do

         end if
!
!        Scan over r-z grid for interior point
!
         rmin = minval(r1b);  rmax = maxval(r1b)
         zmin = minval(z1b);  zmax = maxval(z1b)
         rcom(iv) = (rmax + rmin)/2; zcom(iv) = (zmax + zmin)/2

!
!        Estimate jacobian based on boundary and 1/2 surface
!
         ds = (ns - ns12)*hs
         do iu = 1, ntheta1
            rs(iu) = (r1b(iu) - r12(iu))/ds + r1(1,iv,1,0)
            zs(iu) = (z1b(iu) - z12(iu))/ds + z1(1,iv,1,0)
            tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
         end do

         mintau = 0

         do nlim = 1, limpts
            zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
            if (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) then
               zlim = 0
               if (nlim .gt. 1) exit
            end if
            rlim(:) = ( / ((rmin + ((klim-1)*(rmax - rmin))/(limpts-1)),
     1              klim = 1, limpts) / )

!
!           Find value of magnetic axis that maximizes the minimum jacobian value
!
            do klim = 1, limpts
               tau = tau0 - ru12(:)*zlim + zu12(:)*rlim(klim)
               mintemp = minval(tau*signgs)
               if (mintemp .gt. mintau) then
                  mintau = mintemp
                  rcom(iv) = rlim(klim)
                  zcom(iv) = zlim
!           If up-down symmetric and lasym=T, need this to pick z = 0
               else if (mintemp .eq. mintau) then
                  if (abs(zcom(iv)).gt.abs(zlim)) zcom(iv) = zlim
               end if
            end do
         end do

      end do planes

!
!     FOURIER TRANSFORM RCOM, ZCOM
!
      dzeta = two/nzeta
      do n = 0, ntor
         raxis(n,1) = dzeta*sum(cosnv(:,n)*rcom(:))/nscale(n)
         zaxis(n,1) =-dzeta*sum(sinnv(:,n)*zcom(:))/nscale(n)
         raxis(n,2) =-dzeta*sum(sinnv(:,n)*rcom(:))/nscale(n)
         zaxis(n,2) = dzeta*sum(cosnv(:,n)*zcom(:))/nscale(n)
         if (n.eq.0 .or. n.eq.nzeta/2) then
            raxis(n,1) = p5*raxis(n,1)
            zaxis(n,2) = p5*zaxis(n,2)
         end if
!        if (lscreen) print 100, n, raxis(n), zaxis(n)
      end do

!  100 format(' n = ',i4,' raxis = ',1pe10.3,' zaxis = ',1pe10.3)

      end subroutine guess_axis
      
      subroutine heading(extension, time_slice, iseq_count, lmac, 
     1     lscreen, lfirst)
      use vmec_main, only: rprec
      use vparams, only: nthreed, nmac
      use vmec_params, only: version_
      use date_and_computer
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq_count
      real(rprec) :: time_slice
      character*(*) :: extension
      logical :: lmac, lscreen, lfirst
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(50), parameter ::
     1   banner = ' THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION '
      character*(*), parameter :: VersionID1 =
     1   ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.',
     2   VersionID2 = ' Forces Are Conservative'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: imon, nout
      character*(10) :: date0, time0, zone0
      character*(50) :: dateloc, Version
C----------------------------------------------- 
!
!     open output files
!
      call open_output_files (extension, iseq_count, lmac, lscreen, 
     1     lfirst)

c     FORTRAN-90 ROUTINE
      call date_and_time(date0,time0,zone0)
      read(date0(5:6),'(i2)')imon
      write(dateloc,100)months(imon),date0(7:8),date0(1:4),
     1  time0(1:2),time0(3:4),time0(5:6)
 100  format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

      if (lscreen .and. lfirst) write (*,'(a,i4,a,1pe12.4/2a)')
     1  '  SEQ = ', iseq_count+1,
     2  ' TIME SLICE',time_slice,'  PROCESSING INPUT.', trim(extension)

      Version = trim(adjustl(version_)) 
      write(nthreed,'(a,1x,a,/,2a,//,2a,2x,a)') trim(banner), 
     1     trim(Version), trim(VersionID1), trim(VersionID2),
     2     ' COMPUTER: ', trim(computer), trim(dateloc)
      if (lscreen .and. lfirst) 
     1   write (*,'(1x,a,1x,a,/,1x,2a,//,1x,2a,2x,a)') trim(banner),
     2   trim(Version), trim(VersionID1), trim(VersionID2), 
     3   ' COMPUTER: ', trim(computer), trim(dateloc)

      do nout = nthreed, nthreed+1
        imon = nout
        if (imon .eq. nthreed+1) imon = nmac
        if (imon.eq.nmac .and. .not.lmac) cycle
        write (imon,3) trim(extension),iseq_count,time_slice
      enddo

 3    format(' SHOT ID.: ',a,2x,'SEQ. NO.:',i4,/,
     1       ' TIME SLICE = ',f5.0,' ms')

      end subroutine heading      

      subroutine interp(xnew, xold, scalxc, nsnew, nsold)
      use vmec_main, only: dp, rprec, mnsize
      use vmec_params, only: ntmax
      use vmec_persistent, only: ixm
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nsnew, nsold
      real(rprec), dimension(nsnew,mnsize,3*ntmax),
     1  intent(out) :: xnew
      real(rprec), dimension(nsnew,mnsize,3*ntmax),
     1  intent(in) :: scalxc
      real(rprec), dimension(nsold,mnsize,3*ntmax) :: xold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, two = 2.0_dp
     1   ,zero = 0.0_dp, one = 1.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ntype, js, js1, js2
      real(rprec) :: hsold, sj, s1, xint
C-----------------------------------------------
 
      if (nsold .le. 0) return 
      hsold = one/(nsold - 1)
!
!       INTERPOLATE R,Z AND LAMBDA ON FULL GRID
!       (EXTRAPOLATE M=1 MODES,OVER SQRT(S), TO ORIGIN)
!       ON ENTRY, XOLD = X(COARSE MESH) * SCALXC(COARSE MESH)
!       ON EXIT,  XNEW = X(NEW MESH)   [ NOT SCALED BY 1/SQRTS ]
!
 
      do ntype = 1, 3*ntmax
 
         where (mod(ixm(:mnsize),2) .eq. 1) xold(1,:,ntype) =
     1       two*xold(2,:,ntype) - xold(3,:,ntype)

         do js = 1, nsnew
            sj = real(js - 1,rprec)/(nsnew - 1)
            js1 = 1 + ((js - 1)*(nsold - 1))/(nsnew - 1)
            js2 = min(js1 + 1,nsold)
            s1 = (js1 - 1)*hsold
            xint = (sj - s1)/hsold
            xint = min(one,xint)
            xint = max(zero,xint)
            xnew(js,:,ntype) = ((one - xint)*xold(js1,:,ntype)+xint*
     1         xold(js2,:,ntype))/scalxc(js,:,1)
         end do
 
 
!        Zero M=1 modes at origin
         where (mod(ixm(:mnsize),2) .eq. 1) xnew(1,:,ntype) = 0
 
      end do
 
      end subroutine interp
             

      subroutine load_xc_from_wout(rmn, zmn, lmn, ntor_in, mpol1_in,
     1    ns_in, reset_file, lreset)
      use read_wout_mod, only: rmnc, zmns, lmns, rmns, zmnc, lmnc, 
     1    xm, xn, ntor, ns, 
     2    nfp, mnmax, read_wout_file, read_wout_deallocate
      use vmec_params, only: rprec, mscale, nscale, ntmax
      use vmec_dim, only: mpol1
      use vparams, only: one, zero
      use vmec_input, only: lasym 
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns_in, mpol1_in, ntor_in
      real(rprec), dimension(ns_in,0:ntor_in,0:mpol1_in,ntmax), 
     1   intent(out) :: rmn, zmn, lmn
      character*(*) :: reset_file
      logical lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rbcc=1, rbss=2, rbcs=3, rbsc=4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ierr, mn, m, n, n1, zbcs, zbsc, zbcc, zbss
      real(rprec) :: t1, t2
C-----------------------------------------------


!
!     THIS ALLOWS SEQUENTIAL RUNNING OF VMEC FROM THE COMMAND LINE
!     i.e., WHEN VMEC INTERNAL ARRAYS ARE NOT KEPT IN MEMORY (compared to sequence file input)
!     THIS IS THE CASE WHEN VMEC IS CALLED FROM, SAY, THE OPTIMIZATION CODE
!
      call read_wout_file(reset_file, ierr)
      
      if (ierr .ne. 0) then
         print *,' Error opening/reading wout file in VMEC load_xc!'
         lreset = .true.
         return
      end if
      
      lreset = .false.   

      if (ns_in .ne. ns) then
         lreset = .true.
         print *, 'ns_in != ns in load_xc'
         return
      end if

      if (ntor_in  .ne. ntor ) stop 'ntor_in != ntor in load_xc'
      if (mpol1_in .ne. mpol1) stop 'mpol1_in != mpol1 in load_xc'
      if (nfp .eq. 0) stop 'nfp = 0 in load_xc'
      
      rmn = zero
      zmn = zero
      lmn = zero

      zbcs = rbcc
      zbsc = rbss
      zbcc = rbcs
      zbss = rbsc
      
      do mn = 1, mnmax
         m = nint(xm(mn))
         n = nint(xn(mn))/nfp
         n1 = abs(n)
         t1 = one/(mscale(m)*nscale(n1))
         t2 = t1
         if (n .lt. 0) t2 = -t2
         if (n .eq. 0) t2 = zero
         rmn(:ns, n1, m, rbcc) = rmn(:ns, n1, m, rbcc) + t1*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbsc) = zmn(:ns, n1, m, zbsc) + t1*zmns(mn,:ns)
         lmn(:ns, n1, m, zbsc) = lmn(:ns, n1, m, zbsc) + t1*lmns(mn,:ns)
         rmn(:ns, n1, m, rbss) = rmn(:ns, n1, m, rbss) + t2*rmnc(mn,:ns)
         zmn(:ns, n1, m, zbcs) = zmn(:ns, n1, m, zbcs) - t2*zmns(mn,:ns)
         lmn(:ns, n1, m, zbcs) = lmn(:ns, n1, m, zbcs) - t2*lmns(mn,:ns)
         if (lasym) then
         rmn(:ns, n1, m, rbcs) = rmn(:ns, n1, m, rbcs) - t2*rmns(mn,:ns)
         zmn(:ns, n1, m, zbcc) = zmn(:ns, n1, m, zbcc) + t1*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbcc) = lmn(:ns, n1, m, zbcc) + t1*lmnc(mn,:ns)
         rmn(:ns, n1, m, rbsc) = rmn(:ns, n1, m, rbsc) + t1*rmns(mn,:ns)
         zmn(:ns, n1, m, zbss) = zmn(:ns, n1, m, zbss) + t2*zmnc(mn,:ns)
         lmn(:ns, n1, m, zbss) = lmn(:ns, n1, m, zbss) + t2*lmnc(mn,:ns)
         end if
         if (m .eq. 0) then
            rmn(:ns, n1, m, rbss) = zero
            zmn(:ns, n1, m, zbsc) = zero
            lmn(:ns, n1, m, zbsc) = zero
            if (lasym) then
            rmn(:ns, n1, m, rbsc) = zero
            zmn(:ns, n1, m, zbss) = zero
            lmn(:ns, n1, m, zbss) = zero
            end if
         end if   
      end do   


      call read_wout_deallocate

      end subroutine load_xc_from_wout


      subroutine profil1d(xc, xcdot, lreset)
      use vmec_main
      use vmec_params, only: signgs
      use vsvd, torflux_edge => torflux
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(neqs2), intent(out) :: xc, xcdot
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: Itor, si, pedge, tflux, tfluxd
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec), external :: pcurr, pmass, piota, torflux, 
     1    torflux_deriv
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        ai       array of coefficients in phi-series for iota (ncurr=0)
!        ac       array of coefficients in phi-series for the quantity d(jcurv)/ds = toroidal
!                 current density * Vprime, so jcurv = Itor(s) (ncurr=1)
!        am       array of coefficients in phi-series for mass (NWT/m**2)
!        iotas    rotational transform , on half radial mesh
!        jcurv    (-)toroidal current inside flux surface (vanishes like s)
!        mass     mass profile on half-grid
!        phiedge  value of real toroidal flux at plasma edge (s=1)
!        phips    same as phip , one-dimensional array
!        presf    pressure profile on full-grid, mass/phip**gamma
!        spres_ped value of s beyond which pressure profile is flat (pedestal)
 
!
!       COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
!       COMPUTE MASS PROFILE ON HALF-GRID
!       BY READING INPUT COEFFICIENTS. PRESSURE CONVERTED TO
!       INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7
!
      torflux_edge = signgs * phifac * phiedge / twopi
      r00 = rmn_bdy(0,0,1)

      phips(1) = 0
      jcurv(1) = 0

      do i = 2,ns
         si = hs*(i - c1p5)
         tflux = torflux(si)
         phips(i) = torflux_edge * torflux_deriv(si)
         iotas(i) = piota(tflux)
         jcurv(i) = pcurr(tflux)
      end do

      do i = 1,ns
         si = hs*(i-1)
         tflux = torflux(si)
         iotaf(i) = piota(tflux)
      enddo
!
!     SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
!
      pedge = pcurr(one)
      Itor = 0
      if (abs(pedge) .gt. abs(epsilon(pedge)*curtor)) 
     1   Itor = signgs*currv/(twopi*pedge)
      jcurv(2:ns) = -signgs*Itor*jcurv(2:ns)

!
!     POSSIBLE PRESSURE PEDESTAL FOR S >= SPRES_PED
!
      spres_ped = abs(spres_ped)
      if (.not.lrecon) then
        do i = 2,ns
          si = hs*(i - c1p5)
          tflux = torflux(si)
          tfluxd = torflux_edge * torflux_deriv(si)
          if (si .gt. spres_ped) then 
             pedge = pmass(spres_ped)
          else
             pedge = pmass(tflux)
          end if   
          mass(i) = pedge*(abs(tfluxd)*r00)**gamma
        end do
 
      else
        iotas(:ns) = 0
        iotaf(:ns) = 0
        mass (:ns) = 0
        presf(:ns) = 0
      end if

      pres(:ns+1) = 0
      xcdot(:neqs2) = 0

      if (lreset) then
c        xc(:neqs1) = 0
        xc(:) = 0
        if (lrecon) iresidue = 0
      end if  
      if (lrecon) then
        if (iresidue .gt. 1) iresidue = 1        
!
!       COMPUTE INDEX ARRAY FOR FINDPHI ROUTINE
!
        do i = 1,ns
          indexr(i)    = ns + 1 - i                    !FINDPHI
          indexr(i+ns) = i + 1                         !FINDPHI
        enddo
        indexr(2*ns) = ns
      end if 

      end subroutine profil1d
      

      function pmass (xx)
      use kind_spec
      use vmec_input, only: am, bloat
      use vparams, only: dmu0
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: im
      real(rprec) :: xx, pmass, x
C-----------------------------------------------
!     NOTE: On entry, am is in pascals. pmass internal units are mu0*pascals (B**2 units)
      x = min (abs(xx * bloat), 1._dp)
      pmass = 0
      do im = size(am), 1, -1
         pmass = x*pmass + am(im-1)
      end do
      pmass = dmu0*pmass

      end function pmass

      
      function piota (x)
      use kind_spec
      use vmec_input, only: ai
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: iq
      real(rprec) :: x, piota
C-----------------------------------------------
      piota = 0
      do iq = size(ai), 1, -1
         piota = x*piota + ai(iq-1)
      end do

      end function piota

      
      function pcurr (xx)
      use kind_spec
      use vmec_input, only: ac, bloat, ac_form
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: ic
      real(rprec) :: xx, pcurr, x
C-----------------------------------------------
      x = min (abs(xx * bloat), 1._dp)
      pcurr = 0
      select case (ac_form)
      case(1)
         pcurr = ac(0)*x + (2.*ac(1)*x**1.5)/3
         do ic = 2, size(ac)
            pcurr = pcurr + (ac(ic)*x**ic)/ic
         end do

      case default
         do ic = size(ac), 1, -1
            pcurr = x*pcurr + ac(ic-1)/ic
         end do
         pcurr = x*pcurr
      end select
   
      end function pcurr
      

      function torflux (x)
      use kind_spec
!     use vmec_input, only: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, torflux
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     torflux = x*(af(1) + x*(af(2) + x*(af(3) + x*(af(4) + 
!    1     x*(af(5) + x*(af(6) + x*(af(7) + x*(af(8) + x*(af(9) + 
!    2     x*af(10))))))))))
      torflux = x

      end function torflux


      function torflux_deriv (x)
      use kind_spec
!     use vmec_input, only: af => aphi
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, torflux_deriv
C-----------------------------------------------
!     TEMPORARILY DISABLED....
!     torflux_deriv = af(1) + x*(2*af(2) + x*(3*af(3) + x*(4*af(4) + 
!    1           x*(5*af(5) + x*(6*af(6) + x*(7*af(7) + x*(8*af(8) + 
!    2           x*(9*af(9) + x*10*af(10)))))))))
      torflux_deriv = 1

      end function torflux_deriv
      

      subroutine profil3d(rmn, zmn, lreset)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      use vsvd, torflux_svd => torflux
      use vspline
      use realspace
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax),
     1    intent(inout) ::  rmn, zmn
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, lk, lt, lz, ntype, m, n, mn, ntype2
      real(rprec), dimension(0:ntor,ntmax) :: rold, zold
      real(rprec) :: sm0, t1, facj, si
      integer :: jcount, jk, k
      real(rprec), external :: torflux
C-----------------------------------------------

!
!                INDEX OF LOCAL VARIABLES
!
!        phip     radial derivative of phi/(2*pi) on half-grid
!        shalf    sqrt(s) ,two-dimensional array on half-grid
!        sqrts    sqrt(s), two-dimensional array on full-grid
!        wint     two-dimensional array for normalizing angle integrations
!        ireflect two-dimensional array for computing 2pi-v angle

      do js = 1, ns
         phip(js:nrzt:ns) = phips(js)
      end do

      phip(nrzt+1) = zero

      faclam(:ns,0:ntor,0:mpol1,:ntmax) = zero


      do js = 1, ns
         si = hs*abs(js - 1.5_dp)
         si = torflux(si)
         shalf(js:nrzt:ns) = sqrt(si)
         si = hs*(js - 1)
         si = torflux(si)
         sqrts(js:nrzt:ns) = sqrt(si)
      end do
      sqrts(ns:nrzt:ns) = one     !!Avoid round-off

      lk = 0
      do lt = 1, ntheta3
         do lz = 1, nzeta
            lk = lk + 1
            if (lasym) then
               do js = 2, ns
                  wint(js+ns*(lk-1)) = one/(nzeta*ntheta1)
               end do
            else
               do js = 2, ns
                  wint(js+ns*(lk-1)) = cosmui(lt,0)/mscale(0)
               end do
            end if
         end do
      end do

      shalf(nrzt+1) = 1
      wint(1:nrzt:ns) = 0

!
!       COMPUTE ARRAY FOR REFLECTING v = -v (only needed for lasym)
!
      jcount = 0
      do k = 1, nzeta
         jk = nzeta + 2 - k
         if (k .eq. 1) jk = 1
         do js = 1,ns
           jcount = jcount+1
           ireflect(jcount) = js+ns*(jk-1)           !Index for -zeta[k]
         enddo
      end do

      do js = 2,ns
         sm(js) = shalf(js)/sqrts(js)
         sp(js) = shalf(js+1)/sqrts(js)
      enddo
      sm(1) = 0
      sp(0) = 0
      sp(1) = sm(2)


!
!     COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
!     FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
!     (1/SQRTS FACTOR FOR ODD M VALUES)
!
!
!     ALLOW FOR lreset == F RESTART PERTURBATION OF FIXED-BDY (LFREEB=F)
!

         do js = 1, ns
            si = hs*(js - 1)
            sm0 = one - torflux(si)
            do ntype = 1, ntmax
               do m = 0, mpol1
                  do n = 0, ntor
                     t1 = one/(mscale(m)*nscale(n))
                     mn = n + ntor1*m
                     l = js + ns*mn + (ntype - 1)*mns
                     if (mod(m,2) .eq. 0) then
                        scalxc(l) = one
                        scalxc(l+2*irzloff) = one   !Lambda on full mesh
                     else
                        scalxc(l) = one/max(sqrts(js),sqrts(2))
                        scalxc(l+2*irzloff)=one/max(sqrts(js),sqrts(2))
                     endif
                     if (.not.lreset .and. lfreeb) cycle
                     if (m .eq. 0) then
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype) +
     1                     (one - sm0)*(rmn_bdy(n,m,ntype)*t1
     2                    - rmn(ns,n,m,ntype))
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype) +
     1                    (one - sm0)*(zmn_bdy(n,m,ntype)*t1
     2                    - zmn(ns,n,m,ntype))
                        if (mod(ntype,2).eq.1 .and. lreset) then
                           if (js .eq. 1) then
                              rold(n,ntype) = rmn(1,n,0,ntype)
                              zold(n,ntype) = zmn(1,n,0,ntype)
                           endif
                           ntype2 = (ntype+1)/2
                           rmn(js,n,m,ntype) = rmn(js,n,m,ntype) +
     1                       sm0*(raxis(n,ntype2)*t1 - rold(n,ntype))
                           zmn(js,n,m,ntype) = zmn(js,n,m,ntype) +
     1                       sm0*(-zaxis(n,ntype2)*t1 - zold(n,ntype))
                        endif

                     else
                        facj = sqrts(js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
!                       if (mod(m,2) .eq. 0) then
!                           facj = sqrts(js)*sqrts(js)
!                       else if (mod(m,2) .eq. 1) then
!                          facj = sqrts(js)**min(m,3)
!                       end if
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + (rmn_bdy
     1                     (n,m,ntype)*t1 - rmn(ns,n,m,ntype))*facj
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + (zmn_bdy
     1                     (n,m,ntype)*t1 - zmn(ns,n,m,ntype))*facj
                     endif

                  end do
               end do
            end do
         end do

         scalxc(1+irzloff:2*irzloff) = scalxc(:irzloff)              !Z-components

!
!     STORE PHIFAC IN XC(NEQS1) ARRAY ELEMENT
!     STORE DEL-RMSE IN XC(NEQS2) ARRAY ELEMENT
!
      xc(neqs1) = phifac
      scalxc(neqs1) = 1
      scalxc(neqs2) = 1

      if (lrecon) then
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF FUNCTIONS IN SPLININT
!       THESE ARE FIXED DURING THE ITERATION SEQUENCE FOR FIXED SK,PK-NOTS,
!       BUT MAY CHANGE IF SKNOTS, PKNOTS CHANGE FOR DIFFERENT DATA FILES
!
        call setup_int (sknots, shalf(2), hstark, w_ia, w1_ia, u_ia,
     1       u1_ia, nk_ia, isnodes, ns1)
c-08-96 call setup_int(pknots,shalf(2),hthom,w_pa,w1_pa,u_pa,u1_pa,
c-08-96 >  nk_pa,ipnodes,ns1)
!
!       COMPUTE SPLINE SEGMENTS FOR INTEGRALS OF DERIVATIVES IN SPLININT
!
        call setup_int (sknots, sqrts, hstark, w_ib, w1_ib, u_ib,
     1       u1_ib, nk_ib, isnodes, ns)
        call setup_int (pknots, sqrts, hthom, w_pb, w1_pb, u_pb,
     1       u1_pb, nk_pb, ipnodes, ns)
      endif

      end subroutine profil3d

      
      subroutine evolve(time_step, ier_flag, liter_flag, lscreen)
      use vmec_main
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: time_step
      integer ier_flag
      logical liter_flag, lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: fsq1, dtau, b1, fac
C-----------------------------------------------
 
!     COMPUTE MHD FORCES
 
      call funct3d (lscreen, ier_flag)
 
!     COMPUTE ABSOLUTE STOPPING CRITERION
 
      if (irst .ne. 1) then
         if (iter2 .eq. 1) ier_flag = 1
         return 
      else if (fsqr .le. ftolv .and. fsqz .le. ftolv .and. 
     1    fsql .le. ftolv) then
         liter_flag = .false.
         return 
      endif
 
!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
 
      if (iter2 .eq. iter1) otau(:ndamp) = cp15/time_step

      fsq1 = fsqr1 + fsqz1 + fsql1
      if (iter2 .gt. iter1) dtau = min(abs(log(fsq/fsq1)),cp15)
      fsq = fsq1
      if (iter2 .le. 1) return 
 
      otau(1:ndamp-1) = otau(2:ndamp)
 
      if (iter2 .gt. iter1) otau(ndamp) = dtau/time_step
      otav = sum(otau(:ndamp))/ndamp
      dtau = time_step*otav
      b1  = one - cp5*dtau
      fac = one/(one + cp5*dtau)
 
      do i = 1,neqs2
        xcdot(i) = fac*(xcdot(i)*b1 + time_step*gc(i))
        xc(i) = xc(i) + xcdot(i)*time_step
      end do  
 
      end subroutine evolve

      
      subroutine fileout(iseq, ier_flag, lscreen)
      use vmec_main
      use vac_persistent
      use realspace
      use vmec_params, only: mscale, nscale, signgs
      use vforces, czmn => czmn, lu => czmn, crmn => crmn, lv => crmn
      use vsvd
      use xstuff, ONLY: xc
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq, ier_flag
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      logical, parameter :: lreset_xc = .false.
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, istat1=0
      real(rprec) :: twouton, twoutoff, teqon, teqoff      
      character*(92), dimension(0:10), save :: werror
      character*(35), dimension(0:6), save :: form
      character*(*), parameter :: 
     1    Warning = " global memory deallocation error "
      logical log_open
C-----------------------------------------------
      data werror/'EXECUTION TERMINATED NORMALLY', 
     1   'INITIAL JACOBIAN CHANGED SIGN (IMPROVE BOUNDARY GUESS)', 
     2   'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)', 
     3   'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.',
     4   'FORCE RESIDUALS EXCEED FTOL: INCREASING NUMBER ITERATIONS', 
     5   'ERROR READING INPUT FILE OR NAMELIST',
     6   'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN',
     7   'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE',
     8   'NS ARRAY MUST NOT BE ALL ZEROES',
     9   'ERROR READING MGRID FILE',
     A   'VACUUM-VMEC MISMATCH IN TOROIDAL CURRENT: INITIAL BOUNDARY MAY
     A BE ENCLOSING EXTERNAL CURRENT'/
      data form/'  TOTAL COMPUTATIONAL TIME :       ',
     1   '  TIME IN VACUUM LOOP :            ',
     2   '  TIME TO READ IN DATA:            ',
     3   '  TIME TO WRITE DATA TO WOUT:      ',
     4   '  TIME IN EQFORCE                  ',
     5   '  TIME (REMAINDER) IN FUNCT3D:     ',
     6   '  TIME IN PROFILE RECONSTRUCTION:  '/


      INTERFACE

      subroutine eqfor(bsubu, bsubv, tau, rz_array)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax, signgs
      use realspace, br=>rcon, bz=>zcon
      use vforces, r12=>armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,
     2   lu => czmn, lv => crmn, bphi => czmn_o
      use vacmod
      use vsvd
      use vspline
      use csplinx
      use vmec_io
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt,0:1), intent(inout) :: 
     1  bsubu, bsubv
      real(rprec), dimension(nrzt), intent(out) :: tau
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1  target, intent(in) :: rz_array

      end subroutine eqfor

      END INTERFACE

C-----------------------------------------------
!
!     COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
!     CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
!     AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
!
      iequi = 1

      if (ier_flag.eq.0 .or. ier_flag.eq.4) then
!
!     The sign of the jacobian MUST multiply phi to get the PHYSICALLY 
!     correct toroidal flux
!      
         phi(1) = zero
         do js = 2, ns
            phi(js) = phi(js-1) + phip(js)
         end do
         phi = (signgs*twopi*hs)*phi
      
         call funct3d (lscreen, ier_flag)
      end if

      call second0 (teqon)
      if (ier_flag .eq. 0) call eqfor (clmn, blmn, rcon(1,1), xc)
      call second0 (teqoff)
      timer(4) = timer(4) + teqoff - teqon

!
!     MUST Call WROUT To Write Error, if nothing else
!
      call second0 (twouton)

      call wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o,
     1      czmn_e, crmn_e, azmn_e, ier_flag)
 
      call second0 (twoutoff)

      timer(3) = timer(3) + twoutoff - twouton
      timer(0) = timer(0) + timer(3) + timer(4)

      if (ier_flag .ne. 0) goto 1000

      timer(5) = timer(5) - timer(1) - timer(6)

      if (lscreen) print 10, ijacob
      write (nthreed, 10) ijacob
 10   format(/,'  NUMBER OF JACOBIAN RESETS = ',i4,/)

      if (lfreeb .and. lrecon) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1            form(3),timer(3), form(4),timer(4),
     1            form(1),timer(1), form(6),timer(6),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1            form(3),timer(3),form(4),timer(4), form(1),timer(1),
     2            form(4),timer(4), form(6),timer(6),form(5),timer(5)
      else if (lfreeb) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1               form(3),timer(3),form(4),timer(4), form(1),
     2               timer(1),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2), 
     1        form(3),timer(3),form(4),timer(4), form(1),timer(1),
     2        form(5),timer(5)
      else if (lrecon) then
         if (lscreen) print 20, form(0),timer(0),form(2),timer(2),
     1               form(3),timer(3),form(4),timer(4), 
     2               form(6),timer(6),form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1        form(3),timer(3),form(4),timer(4), form(6),timer(6),
     2        form(5),timer(5)
      else
         if (lscreen) print 20,form(0),timer(0),form(2),timer(2),
     1              form(3),timer(3),form(4),timer(4), form(5),timer(5)
         write (nthreed, 20) form(0),timer(0),form(2),timer(2),
     1              form(3),timer(3),form(4),timer(4), form(5),timer(5)
      end if 
   20 format(a35,f12.2,' SECONDS')

         inquire(unit=nlog,opened=log_open)
         if (lrecon .and. log_open) then
!
!     WRITE SEQUENCE HISTORY FILE
!
         if (iseq .eq. 0) write (nlog, 100)
         write (nlog, 110) iseq + 1, iter2, total_chi_square_n,
     1   1.e-6_dp*ctor/dmu0, 1.e-3_dp*ppeak/dmu0, torflux, r00,
     2   timer(0), input_extension
      
      endif

  100 format(' SEQ ITERS  CHISQ/N',
     1   '  TORCUR  PRESMAX  PHIEDGE     R00 CPU-TIME  EXTENSION')
  110 format(i4,i6,f8.2,3f9.2,f8.2,f9.2,2x,a20)

 1000 continue

      if (lscreen) print 120, trim(werror(ier_flag)), input_extension
      write (nthreed, 120) trim(werror(ier_flag)), input_extension
  120 format(2x,a,/,2x,'FILE : ',a/)

!
!     DEALLOCATE GLOBAL MEMORY
!
      if (allocated(cosmu))
     1  deallocate(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
     2  sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#1"

      if (allocated(xm)) deallocate (xm, xn, ixm, jmin3,
     1   mscale, nscale, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#2"
     
      if (allocated(tanu))
     1  deallocate(tanu, tanv, sin2v, sinper, cosper, sinuv, cosuv, 
     2  sinu, cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
     3  cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
      if (istat1 .ne. 0) print *, Warning // "#3"

      call free_mem_funct3d
      call free_mem_ns (lreset_xc)
      call free_mem_nunv

!
!     CLOSE OPENED FILES
!
      if (ier_flag .ne. 4) call close_all_files
 
      end subroutine fileout


      subroutine fixaray
      use vmec_main
      use vmec_params, only: jmin2, mscale, nscale, signgs
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: 
     1    two=2, three=3, five=5, pexp = 4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, m, j, n, mn, mn1, nmin0, istat1, istat2, mnyq, nnyq
      real(rprec):: argi, arg, argj, dnorm
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!       mscale   array for norming theta-trig functions (internal use only)
!                so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
!       nscale   array for norming zeta -trig functions (internal use only)
 
 
!
!    COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
!    NOTE: ARRAYS ALLOCATED HERE ARE GLOBAL AND ARE DEALLOCATED IN FILEOUT
!
      mnyq = max(0, ntheta1/2)           !!Nyquist (max) frequency: needed in jxbforce aliasing
      nnyq = max(0, nzeta/2)
      allocate(cosmu(ntheta2,0:mnyq), sinmu(ntheta2,0:mnyq),
     1           cosmum(ntheta2,0:mnyq), sinmum(ntheta2,0:mnyq), 
     2           cosmui(ntheta2,0:mnyq), cosmumi(ntheta2,0:mnyq),
     3           sinmui(ntheta2,0:mnyq), sinmumi(ntheta2,0:mnyq),
     4           cosnv(nzeta,0:nnyq),     sinnv(nzeta,0:nnyq), 
     5           cosnvn(nzeta,0:nnyq),    sinnvn(nzeta,0:nnyq),
     6           stat=istat1 )
      allocate( xm(mnmax),xn(mnmax), ixm(mnsize), jmin3(0:mnsize),
     1          mscale(0:mnyq), nscale(0:nnyq), stat=istat2)

      if (istat1.ne.0) stop 'allocation error in fixaray: istat1'
      if (istat2.ne.0) stop 'allocation error in fixaray: istat2'
      
      dnorm = one/(nzeta*(ntheta2 - 1))

      mscale(0) = osqrt2
      nscale(0) = osqrt2
      mscale(1:mnyq) = 1
      nscale(1:nnyq) = 1
      r0scale = mscale(0)*nscale(0)
 
      do i = 1, ntheta2
         argi = twopi*(i - 1)/ntheta1
         do m = 0, mnyq
            arg = argi*m
            cosmu(i,m) = cos(arg)*mscale(m)
            sinmu(i,m) = sin(arg)*mscale(m)
            cosmui(i,m) = dnorm*cosmu(i,m)
            sinmui(i,m) = dnorm*sinmu(i,m)
            if (i.eq.1 .or. i.eq.ntheta2) cosmui(i,m)=0.5_dp*cosmui(i,m)
            cosmum(i,m) = cosmu(i,m)*(m)
            sinmum(i,m) = -sinmu(i,m)*(m)
            cosmumi(i,m) = cosmui(i,m)*(m)
            sinmumi(i,m) = -sinmui(i,m)*(m)
         end do
      end do
 
      do j = 1, nzeta
         argj = twopi*(j - 1)/nzeta
         do n = 0, nnyq
            arg = argj*(n)
            cosnv(j,n) = cos(arg)*nscale(n)
            sinnv(j,n) = sin(arg)*nscale(n)
            cosnvn(j,n) = cosnv(j,n)*(n*nfp)
            sinnvn(j,n) = -sinnv(j,n)*(n*nfp)
         end do
      end do
 
!
!     R,Z,L / s**(m/2) ARE LINEAR NEAR ORIGIN
!
      mn = 0
      mn1 = 0
      jmin3 = 0                       ! mcz
      do m = 0, mpol1
         xmpq(m,1) = (m*(m - 1))
         xmpq(m,2) = (m**pexp)
         xmpq(m,3) = (m**(pexp+1))
         do n = 0, ntor
            jmin3(mn) = jmin2(m)
            mn = mn + 1
            ixm(mn) = m
         end do
         nmin0 = -ntor
         if (m .eq. 0) nmin0 = 0
         do n = nmin0, ntor
            mn1 = mn1 + 1
            xm(mn1) = (m)
            xn(mn1) = (n*nfp)
         end do
      end do

      if (mn1.ne.mnmax) stop 'mn1 != mnmax'
       
      faccon(0) = zero
      faccon(mpol1) = zero
c05-96    faccon(m) = -0.5_dp*signgs/((1+m)**4)
!!!         formfac = real(m**2,rprec)/(m+1)**2
c                                                !!!* formfac
      faccon(1:mpol1-1) = -0.25_dp*signgs/xmpq(2:mpol1,1)**2
 
      if (lrecon) call getgreen
      if (lfreeb) call precal                               !Fixed arrays for VACUUM

      end subroutine fixaray


      subroutine free_mem_funct3d
      use vmec_main
      use realspace
      use vforces
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0
C-----------------------------------------------
 
      if (allocated(armn))
     1   deallocate (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,
     1   r1, ru, rv, z1, zu, zv, rcon, zcon, ru0, zu0,
     2   rcon0, zcon0, guu, guv, gvv, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in funct3d'

      if (allocated(brv))
     1   deallocate (brv, bphiv, bzv, bpolvac, bsqvac, stat=istat1)
      if (istat1 .ne. 0) stop 'deallocation error in funct3d'

      end subroutine free_mem_funct3d

      
      subroutine free_mem_ns(lreset)
      use vmec_main
      use realspace
      use vforces
      use vsvd
      use xstuff
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      logical, intent(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0,
     >           istat5 = 0, istat6 = 0, istat7 = 0, istat8 = 0,
     >           istat9 = 0, istat10 = 0
C-----------------------------------------------
      if (allocated(phip)) 
     1  deallocate (phip, shalf, sqrts, wint, stat=istat3)
     
      if (allocated(ireflect))
     1  deallocate (ireflect,indexr,imid, stat=istat4)

      if (allocated(current))
     1  deallocate (current, rm2, vrm2, ovrm2, ochip, presph, presint,
     2      w_ia, w1_ia, u_ia, u1_ia, w_pa, w1_pa, u_pa, u1_pa,
     3      w_ib, w1_ib, u_ib, u1_ib, w_pb, w1_pb, u_pb, u1_pb,
     4      rmid, datamse, qmid, shear, presmid, alfa, curmid,
     5      curint, psimid, ageo, volpsi, isplinef, isplineh,
     6      psplinef, psplineh, phimid, pm, im, stat=istat5)


      if (allocated(pmb)) deallocate (pmb,imb,stat=istat6)

      if (allocated(ard))
     1  deallocate (ard,arm,brd,brm,crd,azd,azm,bzd,bzm, sm,sp,
     2        bmin, bmax,stat=istat7)

      if (allocated(iotaf))
     1  deallocate (iotaf,mass,phi,presf,jcuru,jcurv,jdotb,buco,bvco,
     2     bdotgradv,equif,specw,tcon,fpsi,psi,yellip,yinden,
     3     ytrian,yshift,ygeo,overr,faclam,iotas,phips,pres,vp,
     4     beta_vol, jperp2, jpar2, bdotb, clam, blam, dlam, 
     5     stat=istat8)

      if (allocated(rmidx))
     1  deallocate (rmidx,hmidx,wmidx,qmidx,tenmidx,ymidx,y2midx,
     2     stat=istat9)

      if (allocated(gc))
     1  deallocate (gc, xstore, xcdot, stat=istat10)
      if (allocated(xc) .and. lreset) deallocate (xc, scalxc)
     
      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or.
     1      istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0 .or.
     2      istat7.ne.0 .or. istat8.ne.0 .or. istat9.ne.0 .or.
     3      istat10.ne.0) then
          print *,' deallocation problem in free_mem_ns'
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3,' istat4 = ',istat4
          print *,' istat5 = ',istat5,' istat6 = ',istat6
          print *,' istat7 = ',istat7,' istat8 = ',istat8
          print *,' istat9 = ',istat9,' istat10= ',istat10
       endif  

      end subroutine free_mem_ns

        
      subroutine free_mem_nunv
      use vmec_main
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------

      if (allocated(bsubu0))
     1    deallocate (bsubu0, rbsq, dbsq, stat=istat1)
      if (allocated(rmn_bdy))
     1    deallocate (rmn_bdy, zmn_bdy, stat=istat2)

      if (allocated(amatsav))
     1    deallocate (amatsav, bvecsav, potvac, bsqsav, stat=istat3)
     

      if (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) then
          print *,' deallocation problem in free_mem_nunv'
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3
      endif  

      end subroutine free_mem_nunv

        
      subroutine free_mem_recon
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat1 = 0, istat2 = 0, istat3 = 0
C-----------------------------------------------
      if (allocated(nk_ia))
     1  deallocate (nk_ia,nk_ib, nk_pa, nk_pb,
     2      indexs2,indexu2, indexs1, indexu1,
     3      isortr, isorts, stat=istat1)
      if (allocated(hthom))
     1  deallocate( hthom, ythom, y2thom, pknots,
     2  hstark,y2stark,ystark, sknots, stat=istat2 )
      if (allocated(sthom))
     1  deallocate( sthom, delse2, delso2, pcalc,
     2      delse1, delso1, starkcal, qmeas, qcalc,
     3      fpsical, stark_weight, rsort,rsort0, stat=istat3)
      if ((istat1 .ne. 0) .or. (istat2 .ne. 0)
     1                    .or. (istat3 .ne. 0) )then
          print *,' in free_mem_recon, istat1 = ',istat1
          print *,' istat2 = ',istat2,' istat3 = ',istat3
        end if

      end subroutine free_mem_recon


      subroutine free_persistent_mem()
      use mgrid_mod
      use vmec_main
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat0 = 0
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0
c-----------------------------------------------
      if (allocated(xc)) deallocate (xc, scalxc, stat=istat0)
      if (allocated(btemp)) deallocate(btemp, stat=istat1)
      if (allocated(bvac)) deallocate( bvac,stat=istat2 )
      if (allocated(xobser))
     1  deallocate( xobser, xobsqr, zobser, unpsiext, dsiext,
     2      psiext,plflux, iconnect, needflx, needbfld, plbfld, 
     3      nbcoils, rbcoil, zbcoil, abcoil, bcoil, rbcoilsqr, dbcoil,
     4      pfcspec,dsilabel, bloopnames, curlabel, b_chi,stat=istat3 )
      if (allocated(rlim))
     1  deallocate( rlim,zlim, reslim,seplim,stat=istat4 )
      if (lrecon) call free_mem_recon

      if (istat1.ne.0 .or. istat2.ne.0 .or.
     1    istat3.ne.0 .or. istat4.ne.0) then
          print *,'problem in free_persistent_mem'
          print *,' istat0 = '
          print *,' istat1 = ',istat1,' istat2 = ',istat2
          print *,' istat3 = ',istat3,' istat4 = ',istat4
      endif
     
      end subroutine free_persistent_mem
        

      subroutine funct3d (lscreen, ier_flag)
      use vmec_main
      use vacmod
      use realspace, z1 => z1, gcon => z1
      use vforces, czmn => czmn, lu => czmn, crmn => crmn, lv => crmn
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(inout) :: ier_flag
      logical, intent(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l0pi, l, lk, ivacskip
      real(rprec), dimension(mnmax) :: 
     1   rmnc, zmns, lmns, rmns, zmnc, lmnc
      real(rprec), dimension(nznt) :: rax, zax
      real(rprec), allocatable, dimension(:,:) :: extra1
     1   , extra2, extra3, extra4
      real(rprec) :: tfunon, tvacon, presf_ns, tvacoff,
     1   tfunoff, delr_mse
C-----------------------------------------------
 
      call second0 (tfunon)

      allocate (extra1(nrzt,0:1), stat=l)
      if (lasym) allocate (
     1   extra2(nrzt,0:1), extra3(nrzt,0:1), extra4(nrzt,0:1),
     2   stat=l)
      if (l .ne. 0) stop 'Allocation error in funct3d'
!
!     CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
!
      call extrap(xc)
      gc(:neqs2) = xc(:neqs2)*scalxc(:neqs2)
 
!
!     RIGID BODY SHIFT OF RMNCC(JS.GT.1,0,0) BY DELR_MSE= R00-RAXMSE
!
      if (lrecon) then
         delr_mse = xc(neqs2)
         gc(1:ns) = gc(1:ns) + delr_mse
      endif

!
!     INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
!     FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
!     ON THE RANGE u = 0,pi  and v = 0,2*pi
!
      call totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

      if (lasym) then
!
!     NEXT, DO ANTI-SYMMETRIC PIECES
!
      call totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4, blmn, 
     1     clmn, extra1, extra2)

!     NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
!     TO GET "R's, Z's, L's" ON FULL RANGE OF u (0 to 2*pi)
 
      call symrzl (r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon, armn, 
     1   brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2)

      endif
 
      l0pi = ns*(1 + nzeta*(ntheta2 - 1))        !u = pi, v = 0, js = ns
      router = r1(ns,0) + r1(ns,1)
      rinner = r1(l0pi,0) + r1(l0pi,1)
      r00 = r1(1,0)
      z00 = z1(1,0)

!
!     COMPUTE CONSTRAINT RCON, ZCON
!
      do l = 1,nrzt
         rcon(l,0) = rcon(l,0) + rcon(l,1)*sqrts(l)
         zcon(l,0) = zcon(l,0) + zcon(l,1)*sqrts(l)
         ru0(l) = ru(l,0) + ru(l,1)*sqrts(l)
         zu0(l) = zu(l,0) + zu(l,1)*sqrts(l)
      end do 

! 
!     COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
!
      call jacobian
      if (irst.eq.2 .and. iequi.eq.0) goto 100 

!
!     COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
!     SCALE BY POWER OF SQRTS SO THAT RESTART FOR FIXED BOUNDARY DOES NOT
!     HAVE A DISCONTINUITY DUE TO NEW RCON0....
!
      if (iter2.eq.1 .and. .not.lfreeb) then
         do l = 1, ns
            rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
            zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
         end do
!        rcon0(:nrzt) = rcon(:nrzt,0)
!        zcon0(:nrzt) = zcon(:nrzt,0)
      endif

! 
!     COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
!     PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
!
      call bcovar (lu, lv)

!     COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
!     NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
!     AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
!     EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
!     EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
!     TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).
 
      if (lfreeb .and. iter2.gt.1) then
         if (fsqr + fsqz .le. 1.e2_dp) ivac = ivac + 1
         if (ivac .ge. 0) then
            call second0 (tvacon)
            ivacskip = mod(iter2 - iter1,nvacskip)
            if (ivac .le. 2) ivacskip = 0
            do lk = 1, nznt
               rax(lk) = r1(1 + ns*(lk - 1),0)
               zax(lk) = z1(1 + ns*(lk - 1),0)
            end do
            call convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)
            call vacuum (rmnc, rmns, zmns, zmnc, xm, xn, rax, zax, 
     1         ctor, rbtor, wint, ns, ivacskip, ivac, mnmax, 
     2         ier_flag, lscreen)
            if (ier_flag .ne. 0) return
            lk = 0
!
!           IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
!           UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
!           IF NON-VARIATIONAL FORCES ARE DESIRED 
!
            presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
!RPRES      if (iequi .ne. 1) presf_ns = 0

            do l = ns, nrzt, ns
               lk = lk + 1
               bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
               rbsq(lk) = (bsqvac(lk) + presf_ns)*ohs*(r1(l,0)+r1(l,1))
               dbsq(lk) = abs(bsqvac(lk) + presf_ns - bsqsav(lk,3))
            end do
            if (ivac .eq. 1) then
               bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
               bsqsav(:nznt,2) = bsqvac(:nznt)
            endif
            call second0 (tvacoff)
            timer(1) = timer(1) + (tvacoff - tvacon)
         endif
      endif

!
!     COMPUTE CONSTRAINT FORCE (GCON => Z1)
!
      extra1(:nrzt,1) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt)
     1    + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
      call alias (gcon, extra1(:,0), extra1(:,1), gc, gc(1+mns), 
     1   gc(1+2*mns), gc(1+3*mns))


      if (iequi .eq. 1) then
         if (lrecon) xc(:ns) = xc(:ns) + delr_mse
         goto 100
      end if

!
!     COMPUTE MHD FORCES ON INTEGER-MESH
!
      call forces

!
!     SYMMETRIZE FORCES (in u-v space) 
!
      if (lasym) call symforce (armn, brmn, crmn, azmn, bzmn,
     1     czmn, blmn, clmn, rcon, zcon, r1, ru, rv, z1, zu, zv, 
     2     extra3, extra4, extra1, extra2)
 
!
!     FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
!
      call tomnsps (gc, armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, 
     1     rcon, zcon)

      if (lasym) call tomnspa (gc, r1, ru, rv, z1, zu, zv,
     1     extra3, extra4, extra1, extra2)

!
!     COMPUTE FORCE RESIDUALS
!
      gc = gc * scalxc
      call residue (gc, gc(1+irzloff), gc(1+2*irzloff), faclam)
 
      gc(neqs1) = gphifac
      if (iopt_raxis .gt. 0) gc(neqs2) = grmse

 100  continue

      deallocate (extra1)
      if (lasym) deallocate (extra2, extra3, extra4)
      call second0 (tfunoff)
      timer(5) = timer(5) + (tfunoff - tfunon)
 
      end subroutine funct3d

              
      subroutine extrap(rzl_array)
      use vmec_main
      use vmec_params, only: jmin2, jlam, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: rzl_array
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: rmnmin, rmnmax
      integer :: zmnmax, lmnmin, lmnmax, m
      real(rprec) :: power
C-----------------------------------------------
      rmnmin = 1
      if (lasym) then
         rmnmax = 4
      else
         rmnmax = 2
      end if

      zmnmax = rmnmax + ntmax
      lmnmin = rmnmin + 2*ntmax
      lmnmax = rmnmax + 2*ntmax

!
!     EXTRAPOLATION AT JS=2 (PUT NONSYMMETRIC HERE, TOO)
!
      do m = 0, mpol1
         if (jmin2(m) .lt. 3) cycle
         power = osqrt2**m
         rzl_array(2,:,m,rmnmin:zmnmax) = 
     1   rzl_array(3,:,m,rmnmin:zmnmax)*power
      end do
      do m = 0, mpol1
         if (jlam(m) .lt. 3) cycle
         power = osqrt2**m
         rzl_array(2,:,m,lmnmin:lmnmax) = 
     1   rzl_array(3,:,m,lmnmin:lmnmax)*power
      end do

      end subroutine extrap




      subroutine alias(gcons, gcona, ztemp, gcs, gsc, gcc, gss)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
      real(rprec), dimension(ns*nzeta,ntheta3) ::
     1   gcons, gcona, ztemp
      real(rprec), dimension(ns,0:ntor,0:mpol1) ::
     1   gcs, gsc, gcc, gss
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, i, ir, jk, jka, n, k, js, l
      real(rprec), dimension(:,:), allocatable :: work
C-----------------------------------------------

      allocate (work(ns*nzeta,4))

      gcons = 0
      gcona = 0

      gcs(:,:ntor,:mpol1) = 0;  gsc(:,:ntor,:mpol1) = 0
      gcc(:,:ntor,:mpol1) = 0;  gss(:,:ntor,:mpol1) = 0

      do m = 1, mpol1 - 1
         work = 0
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i
            if (i .eq. 1) ir = 1
            do jk = 1, ns*nzeta
               work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui(i,m)
               work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
            end do
            if (lasym) then
            do jk = 1, ns*nzeta
               jka = ireflect(jk)
               work(jk,3) = work(jk,3) + ztemp(jka,ir)*cosmui(i,m)
               work(jk,4) = work(jk,4) + ztemp(jka,ir)*sinmui(i,m)
            end do
            end if
         end do

         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               if (.not.lasym) then
               do js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + tcon(js)*work(js+l,1)*
     1               sinnv(k,n)
                  gsc(js,n,m) = gsc(js,n,m) + tcon(js)*work(js+l,2)*
     1               cosnv(k,n)
               end do
               else
               do js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,1)-work(js+l,3))
                  gsc(js,n,m) = gsc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,2)-work(js+l,4))
                  gss(js,n,m) = gss(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,2)+work(js+l,4))
                  gcc(js,n,m) = gcc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,1)+work(js+l,3))
               end do
               end if
            end do
         end do
!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
         work = 0

         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = 2, ns
                  work(js+l,3) = work(js+l,3) + gcs(js,n,m)*sinnv(k,n)
                  work(js+l,4) = work(js+l,4) + gsc(js,n,m)*cosnv(k,n)
               end do
               if (lasym) then
               do js = 2, ns
                  work(js+l,1) = work(js+l,1) + gcc(js,n,m)*cosnv(k,n)
                  work(js+l,2) = work(js+l,2) + gss(js,n,m)*sinnv(k,n)
               end do
               end if
            end do
         end do
         do i = 1, ntheta2
            do jk = 1, ns*nzeta
               gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m)
     1                     + work(jk,4)*sinmu(i,m))*faccon(m)
            end do
            if (lasym) then
            do jk = 1, ns*nzeta
               gcona(jk,i) = gcona(jk,i) + (work(jk,1)*cosmu(i,m)
     1                     + work(jk,2)*sinmu(i,m))*faccon(m)
            end do
            end if
         end do
      end do

      if (lasym) then
!
!     EXTEND GCON INTO THETA = PI,2*PI DOMAIN
!
      do i = 1 + ntheta2, ntheta1
         ir = ntheta1 + 2 - i
         do jk = 1, ns*nzeta
            jka = ireflect(jk)
            gcons(jk,i) = -gcons(jka,ir) + gcona(jka,ir)
         end do
      end do
!
!     ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
!
      gcons(:,:ntheta2) = gcons(:,:ntheta2) + gcona(:,:ntheta2)

      end if

      deallocate (work)

      end subroutine alias

      
      subroutine bss(r12, rs, zs, ru12, zu12, bsubs, bsupu, bsupv,
     1               br, bphi, bz)
      use vmec_main
      use realspace
      use vsvd, only: torflux_edge => torflux
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt), intent(in) :: r12, rs, zs, 
     1     ru12, zu12, bsupu, bsupv 
      real(rprec), dimension(nrzt), intent(out) :: 
     1     br, bphi, bz, bsubs
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, p25 = p5*p5
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec) :: rv12, zv12, gsu, gsv, dphids, t1,
     1    rs12, zs12      
C-----------------------------------------------
      dphids = p25/torflux_edge
 
      do l = 2, nrzt
         t1 = dphids*phip(l)
         rv12 = p5*(rv(l,0)+rv(l-1,0) + shalf(l)*(rv(l,1) + rv(l-1,1)))
         zv12 = p5*(zv(l,0)+zv(l-1,0) + shalf(l)*(zv(l,1) + zv(l-1,1)))
         rs12 = rs(l) + t1*(r1(l,1) + r1(l-1,1))/shalf(l)
         zs12 = zs(l) + t1*(z1(l,1) + z1(l-1,1))/shalf(l)
         gsu  = rs12*ru12(l) + zs12*zu12(l)
         gsv  = rs12*rv12    + zs12*zv12
         br(l)    = bsupu(l)*ru12(l) + bsupv(l)*rv12
         bphi(l)  = bsupv(l)*r12(l)
         bz(l)    = bsupu(l)*zu12(l) + bsupv(l)*zv12
         bsubs(l) = bsupu(l)*gsu + bsupv(l)*gsv
      end do
      end subroutine bss
      

      subroutine eqfor(bsubu, bsubv, tau, rz_array)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax, signgs
      use realspace, br=>rcon, bz=>zcon
      use vforces, r12=>armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,
     2   lu => czmn, lv => crmn, bphi => czmn_o
      use vacmod
      use vsvd
      use vspline
      use csplinx
      use vmec_io
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt,0:1), intent(inout) :: 
     1  bsubu, bsubv
      real(rprec), dimension(nrzt), intent(out) :: tau
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1  target, intent(in) :: rz_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5=0.5_dp, c1p5=1.5_dp, two=2.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, icount, itheta, js, js1, l, loff,
     1   lpi, lt, n, n1, nchicur, nchiiota0, ncol, noff,
     2   nout, nsort, iv, iu, lk
      real(rprec), dimension(:,:), pointer ::
     1   rmags, zmags, rmaga, zmaga
      real(rprec), dimension(:,:,:), pointer :: rmncc,zmnsc
      real(rprec), dimension(ns) :: phipf, phi1, chi1, t3,
     1   bvcof, t11u, t21u, jPS2
      real(rprec) :: modb(nznt)
      real(rprec), dimension(:), allocatable ::
     1   bsup1u, dbts1u, dint1u, t12u, guu_1u, guus1u, r3v,
     2   redg1u, rbps1u
      real(rprec) :: aminr1, aminr2, aminr2in, anorm,
     1   area, aspectratio, betai, betpol, betstr,
     2   bminz2, bminz2in, bsq1, btor, iotamax, musubi,
     3   btorsq, bzcalc, bzin, chisq, chiwgt, circum, cur0,
     4   delphid_exact, delta1, delta2, delta3, denwgt, lambda,
     5   dlogidsi, dmusubi_meas, er, es, fac, facnorm, factor, fgeo,
     6   fmax, fmin, flao, fpsi0, pavg, pitchc, pitchm,
     7   pprime, qedge, qmin1, qmin2, qmin3, qzero,
     8   raxis0, rcalc, rcen, rcenin, rgeo, rjs,
     9   rjs1, rlao, rqmin1, rqmin2, rshaf, rshaf1, rshaf2, s11, s12,
     A   s13, s2, s3, sarea, sigr0, sigr1, sigz1, smaleli,
     B   splintx, splints, sqmin, sumbpo, sumbto, sumbtr, sump,
     C   sump2, sump20, t1, jpar_perp=zero, jparPs_perp=zero,
     D   tol, vnorm, volf, vprime, wght0, xmax,
     E   xmida, xmidb, xmin, rzmax, rzmin, zxmax, zxmin, zaxis0,
     F   zmax, zmin, yr1u, yz1u, waist(2), height(2)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
       external splintx,splints
C-----------------------------------------------

!
!     POINTER ASSOCIATIONS
!
      rmags => rz_array(:,:,0,1)
      zmags => rz_array(:,:,0,1+ntmax)
      rmaga => rz_array(:,:,0,3)
      zmaga => rz_array(:,:,0,3+ntmax)
      zmnsc => rz_array(:,:,:,2+ntmax)
      rmncc => rz_array(:,:,:,1)

!
!     NOTE: JXBFORCE ROUTINE MUST BE CALLED FIRST TO COMPUTE IZETA, JDOTB
!           ON OUTPUT, J, IZETA, JDOTB ARE IN MKS UNITS (HAVE 1/MU0)
!
      call bss (r12, bzmn, brmn, azmn, armn, lv(1+nrzt), bsupu, bsupv,
     1          br, bphi, bz)

!
!     STORE EDGE VALUES OF B-FIELD
!
      if (lfreeb .or. ledge_dump) then
         allocate (bredge(2*nznt), bpedge(2*nznt), bzedge(2*nznt))
         do iv = 1,nzeta
            do iu = 1,ntheta3
               lk = iv + nzeta*(iu-1)
               n1 = ns*lk
               bredge(lk) = 1.5_dp*br(n1,0)   - p5*br(n1-1,0)
               bpedge(lk) = 1.5_dp*bphi(n1)   - p5*bphi(n1-1)
               bzedge(lk) = 1.5_dp*bz(n1,0)   - p5*bz(n1-1,0)
            end do
         end do
      end if

      call jxbforce (bsupu, bsupv, bsubu, bsubv, lv(1+nrzt),
     1               br, bz, rcon0, zcon0, gsqrt, izeta, bsq)

!
!     HALF-MESH VOLUME-AVERAGED BETA
!
      tau(1) = zero
      tau(2:nrzt) = signgs*wint(2:nrzt)*gsqrt(2:nrzt)
      do i = 2, ns
         s2 = sum(bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i) - pres(i)
         beta_vol(i) = pres(i)/s2
      end do

      betaxis = c1p5*beta_vol(2) - p5*beta_vol(3)


      write (nthreed, 5)
    5 format(/,' NOTE: <RADIAL FORCE> = d(Ipol)/dPHI',
     1         ' - IOTA*d(Itor)/dPHI - dp/dPHI * dV/dPHI',/,
     1         ' (NORMED TO SUM OF INDIVIDUAL TERMS)',//,
     1         '     S      <RADIAL    TOROIDAL     IOTA     ',
     1         ' <JPOL>     <JTOR>     d(VOL)/',
     2         '   d(PRES)/    <M>     PRESF     <BSUBV> ',
     3         '     <J.B>     <B.B>',/,
     4         '             FORCE>      FLUX                ',
     5         ' (X mu0)    (X mu0)    d(PHI) ',
     6         '    d(PHI)                             ',/,137('-'),/)

      phipf(1) = twopi*signgs*(c1p5*phip(2) - p5*phip(3))
      presf(1) = c1p5*pres(2) - p5*pres(3)
      do i = 2,ns1
         presf(i) = p5*(pres(i) + pres(i+1))
         phipf(i) = p5*twopi*signgs*(phip(i) + phip(i+1))
      end do
      presf(ns) = c1p5*pres(ns)- p5*pres(ns-1)
      phipf(ns) = twopi*signgs*(c1p5*phip(ns) - p5*phip(ns1))

      phi1(1) = zero
      chi1(1) = zero
      do i = 2, ns
         buco(i) = sum(bsubu(i,:,0)*wint(i:nrzt:ns))
         bvco(i) = sum(bsubv(i,:,0)*wint(i:nrzt:ns))
         phi1(i) = phi1(i-1) + hs*phip(i)
         chi1(i) = chi1(i-1) + hs*(phip(i)*iotas(i))
      end do

      bvcof(1) = c1p5*bvco(2) - p5*bvco(3)

!
!     NOTE:  jcuru, jcurv on FULL radial mesh now.
!            They are local (surface-averaged) current densities (NOT integrated in s)
!
      do i = 2,ns1
         jcurv(i) = (signgs*ohs)*(buco(i+1) - buco(i))
         jcuru(i) =-(signgs*ohs)*(bvco(i+1) - bvco(i))
         t11u(i)  = jcurv(i)*iotaf(i)
         t21u(i)  = ohs*(pres(i+1) - pres(i))
         t3(i)    = p5*(vp(i+1)/phip(i+1) + vp(i)/phip(i))
         equif(i) = (jcuru(i) - t11u(i) - (t21u(i)*t3(i)))/(abs(t11u(i))
     1            +  abs(jcuru(i))+abs((t21u(i)*t3(i))))
         bvcof(i) = p5*(bvco(i) + bvco(i+1))
      end do

      bvcof(ns) = c1p5*bvco(ns) - p5*bvco(ns1)

      equif(1) = two*equif(2) - equif(3)
      jcuru(1) = two*jcuru(2) - jcuru(3)
      jcurv(1) = two*jcurv(2) - jcurv(3)
      t21u(1)  = two*t21u(2)  - t21u(3)
      t21u(ns) = two*t21u(ns1) - t21u(ns1-1)
      t3(1)  = two*t3(2) - t3(3)
      t3(ns) = two*t3(ns1) - t3(ns1-1)
      equif(ns) = two*equif(ns1) - equif(ns1-1)
      jcuru(ns) = two*jcuru(ns1) - jcuru(ns1-1)
      jcurv(ns) = two*jcurv(ns1) - jcurv(ns1-1)
      fac = twopi*signgs
      do js = 1, ns
         es = (js - 1)*hs
         cur0 = signgs*fac*t3(js)*phipf(js)
         write (nthreed, 30) es, equif(js), fac*phi1(js),
     1     iotaf(js), fac*jcuru(js)/cur0, fac*jcurv(js)/cur0,
     2     fac*t3(js), t21u(js)/phipf(js), specw(js), presf(js)/dmu0,
     3     bvcof(js), jdotb(js), bdotb(js)
      end do
 30   format(1p2e10.2,1p6e11.3,0pf7.3,1p4e11.3)

!
!     Calculate poloidal cross section area & toroidal flux
!
      anorm = twopi*hs
      vnorm = twopi*anorm
      area = sum(tau(2:nrzt)/r12(2:nrzt))
      do i = 2, ns
         overr(i) = sum(tau(i:(nznt-1)*ns+i:ns)
     1      /r12(i:(nznt-1)*ns+i:ns)) / vp(i)
      end do
      area = area*anorm
      volf = vnorm*sum(vp(2:ns))
      torflux = anorm * sum(bsupv(2:nrzt)*tau(2:nrzt))
      write (nthreed, 3) torflux, phifac
    3 format(/,' TOROIDAL FLUX    = ',f10.2,7x,'PHIFAC = ',f10.2)

!
!
!     OUTPUT BETAS, INDUCTANCES, SAFETY FACTORS, ETC.
!     (EXTRACTED FROM FQ-CODE, 9-10-92)
!
!     b poloidals (cylindrical estimates)
!
      rcen = p5*(router + rinner)               !geometric center
      n = 0
      n1 = n + 1
      rcenin = dot_product(rmncc(ns,n1,:mpol1+1:2),
     1                     mscale(:mpol1:2)*nscale(n))

      l = (mpol1+1)/2
      allocate (t12u(l))
      t12u(:l) = mscale(1:mpol1:2)*nscale(n)
      aminr2in = dot_product(rmncc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2in = dot_product(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2 = dot_product(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      deallocate (t12u)
      aminr1 = sqrt(two*volf/(twopi*twopi*r00))  !vol av minor radius
      aminr2 = p5*(router - rinner)             !geometric plasma radius
!
!       cylindrical estimates for beta poloidal
      sump = vnorm*sum(vp(2:ns)*pres(2:ns))
      pavg = sump/volf
      ppeak = presf(1)
      factor = two*pavg
!
!       delphid_exact = Integral[ (Bvac - B) * dSphi ]
!
      delphid_exact = zero                         !Eq. 20 in Shafranov
      musubi = zero

      rshaf1 = zero
      rshaf2 = zero
      allocate (bsup1u(nznt), dbts1u(nznt), dint1u(nznt))
      do js = 2, ns
         bsup1u(:nznt) = bsubvvac/(r12(js:nrzt:ns)*
     1      r12(js:nrzt:ns))
         delphid_exact = delphid_exact + sum((bsup1u(:nznt) -
     1      bsupv(js:nrzt:ns))*tau(js:nrzt:ns))
         dbts1u(:nznt) = bsup1u(:nznt)*bsubvvac -
     1      bsupv(js:nrzt:ns)*bsubv(js,:nznt,0)
         musubi = musubi + sum(dbts1u(:nznt)*tau(js:nrzt:ns))
         dint1u(:nznt) = (bsupu(js:nrzt:ns)*bsubu(js,:nznt,0)
     1      + dbts1u(:nznt) + two*pres(js))*tau(js:nrzt:ns)
         rshaf1 = rshaf1 + sum(dint1u(:nznt))
         rshaf2 = rshaf2 + sum(dint1u(:nznt)/r12(js:nrzt:ns))
      end do
      deallocate (bsup1u, dbts1u, dint1u)

      delphid_exact = anorm*delphid_exact
      rshaf = rshaf1/rshaf2

      if (lrecon) then
         if (apres .ne. aminr1) then
            write (*, 50) (apres/aminr1)**2
            write (nthreed, 50) (apres/aminr1)**2
         endif
   50 format(/'  Multiply final PHIEDGE by ',f7.3,
     1   ' to make apres = aminr1')
!
!       Output some of measured data
!

         raxis0 = sum(raxis(0:ntor,1))
         zaxis0 = sum(zaxis(0:ntor,1))
         fpsi0 = c1p5*bvco(2) - p5*bvco(3)
         b0 = fpsi0/r00
         cur0 = signgs*iotaf(1)*fpsi0/r00**2*(dkappa + one/dkappa)
         if (lpprof) write (nthreed, 60) presfac*pfac, (-100.*(r00 -
     1      rthompeak))
         write (nthreed, 65) b0, cur0/dmu0
   60    format(//,' Input pressure scaled by ',f6.3/,
     1      ' Pressure profile shifted ',f6.2,' cm.',
     2      ' relative to magnetic axis')
   65    format(//' B-PHI(R=Raxis,Z=0) = ',f8.3,' [Wb/M**2]',/,
     1   ' J-PHI(R=Raxis,Z=0) = iota(0)*Bt0/(mu0*R0)*(1+k**2)/k = ',1p
     2   e10.3,' [A/M**2]',2/,' Comparison of Input vs Output',8x,
     3   'Input',9x,'VMEC Output',/,1x,71('-'))
         write (nthreed, 70) 1.e-6_dp*currv/dmu0, 1.e-6_dp*ctor/dmu0,
     1      phiedge, torflux, 1.e3_dp*phidiam, 
     1      1.e3_dp*delphid, 1.e-3_dp*pthommax,
     2      1.e-3_dp*ppeak/dmu0, rcenin, rcen, 
     2      raxis0, r00, zaxis0, z00, rc0mse,
     3      apres, aminr1, aminr2in, aminr2, bminz2in, bminz2, rinner,
     4      router, 1.e3_dp*delphid_exact
   70    format(' Toroidal Current     =',2(10x,f10.3),'  [MA]',/,
     1   ' Edge Toroidal Flux   =',2(10x,f10.3),'  [Wb]',/,
     2   ' Diamagnetic Flux(*)  =',2(10x,f10.3),'  [mWb]',/,
     3   ' Peak Pressure        =',2(10x,f10.3),'  [KPa]',/,
     4   ' Geometric Center     =',2(10x,f10.3),'  [M]',/,
     5   ' Magnetic R Axis      =',2(10x,f10.3),'  [M]',/,
     6   ' Magnetic Z Axis      =',2(10x,f10.3),'  [M]',/,
     7   ' MSE R Axis           =',10x,f10.3,20x,'  [M]',/,
     8   ' Minor Radius (apres) =',2(10x,f10.3),'  [M]',/,
     9   ' Minor Radius (a)     =',2(10x,f10.3),'  [M]',/,
     .   ' Minor Radius (b)     =',2(10x,f10.3),'  [M]',/,
     1   ' Inboard  Midplane R  =',30x,f10.3,'  [M]',/,
     2   ' Outboard Midplane R  =',30x,f10.3,'  [M]',/,50('-')/,
     3   ' * Exact diamagnetic flux (not based on equilibrium) = ',f8.3,
     4   '  [mWb]')

!       Calculate components of total chi-squared
         total_chi_cur = (currv - ctor)**2/sigma_current**2
         nchicur = 1
         if (iphidiam .eq. 1) total_chi_dia = (phidiam - delphid)**2/
     1      sigma_delphid**2
         nchidia = 1

      else                               ! mcz
         b0 = 0
      end if   !!if(lrecon)

      rmax_surf = maxval(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      rmin_surf = minval(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      zmax_surf = maxval(z1(ns:nrzt:ns,0)+z1(ns:nrzt:ns,1))

!
!     Calculate poloidal circumference and normal surface area
!
      allocate (guu_1u(nznt), guus1u(nznt))
      guu_1u(:nznt) = ru0(ns:nrzt:ns)*ru0(ns:nrzt:ns) +
     1   zu0(ns:nrzt:ns)*zu0(ns:nrzt:ns)

c     *** old version ***
c     guus1u(:nznt) = twopi*wint(ns:nrzt:ns)*sqrt(guu_1u(:nznt))
c     circum = sum(guus1u(:nznt))
c     sarea = twopi*dot_product((r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1)),
c    1   guus1u(:nznt))

c    *** new version ***
      guus1u(:nznt) = wint(ns:nrzt:ns)*sqrt(guu_1u(:nznt))
      circum = twopi*sum(guus1u(:nznt))
      guus1u(:nznt) = wint(ns:nrzt:ns)*sqrt(
     1                (r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))**2
     2                    *guu_1u(:nznt) +
     3                 ((rv(ns:nrzt:ns,0)+rv(ns:nrzt:ns,1))*
     3                   zu0(ns:nrzt:ns) -
     4                  (zv(ns:nrzt:ns,0)+zv(ns:nrzt:ns,1))*
     5                   ru0(ns:nrzt:ns))**2 )
      sarea = twopi**2*sum(guus1u(:nznt))

      deallocate (guu_1u, guus1u)

      aspect = aspectratio()
      do js = 2, ns
         modb(:nznt) = sqrt(two*(bsq(js:nrzt:ns)-pres(js)))
         call bextrema (modb, bmin(1,js), bmax(1,js), nzeta, ntheta2)
      end do

!
!     output geometrical, |B| quantities
!
      call elongation (r1, z1, waist, height)

      write (nthreed, 75) bmin(1,ns), bmax(1,ns), bmin(ntheta2,ns), bmax
     1   (ntheta2,ns)
   75 format(/
     1   ' Magnetic field modulation (averaged over toroidal angle)',/,
     2   1x,71('-')/,' Bmin(u=0)             = ',f14.6/
     3   ' Bmax(u=0)             = ',f14.6/' Bmin(u=pi)            = ',
     4   f14.6/' Bmax(u=pi)            = ',f14.6/)

      sumbto = two*(vnorm*sum(bsq(:nrzt)*tau(:nrzt)) - sump)
      VolAvgB = sqrt(abs(sumbto/volf))
      IonLarmor = 0.0032_dp/VolAvgB
      jPS2(2:ns1) = (jpar2(2:ns1) - jdotb(2:ns1)**2/bdotb(2:ns1))
      s2 = sum(abs(jperp2(2:ns1))*(vp(2:ns1) + vp(3:ns)))
      jpar_perp = sum(jpar2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      jparPS_perp = sum(jPS2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      s2 = sum(jperp2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      if (s2 .ne. zero) then
         jpar_perp = jpar_perp/s2
         jparPS_perp = jparPS_perp/s2
      end if
      if (ntor .gt. 1) then
      write (nthreed, 80) aspect, volf, area, sarea, circum, Rmajor_p,
     1   Aminor_p, rmin_surf, rmax_surf, zmax_surf, waist(1), height(1),
     2   waist(2), height(2)
      else
      write (nthreed, 80) aspect, volf, area, sarea, circum, Rmajor_p,
     1   Aminor_p, rmin_surf, rmax_surf, zmax_surf, waist(1), height(1)
      end if
 80   format(/,' Geometric and Magnetic Quantities',/,1x,71('-')/,
     1   ' Aspect Ratio          = ',f14.6,/' Plasma Volume         = ',
     2   f14.6,' [M**3]',/' Cross Sectional Area  = ',f14.6,' [M**2]',/
     3   ' Normal Surface Area   = ',f14.6,' [M**2]',/
     4   ' Poloidal Circumference= ',f14.6,' [M]',/
     5   ' Major Radius          = ',f14.6,' [M]',
     6   ' (from Volume and Cross Section)',/
     7   ' Minor Radius          = ',f14.6,' [M]',
     8   ' (from Cross Section)',/
     9   ' Minimum (inboard)  R  = ',f14.6,' [M]',/
     A   ' Maximum (outboard) R  = ',f14.6,' [M]',/
     A   ' Maximum height     Z  = ',f14.6,' [M]',/
     B   ' Waist (v = 0)   in R  = ',f14.6,' [M]',/
     B   ' Full Height(v = 0)    = ',f14.6,' [M]',:,/
     B   ' Waist (v = pi)  in R  = ',f14.6,' [M]',:,/
     B   ' Full Height(v = pi)   = ',f14.6,' [M]')
      write (nthreed, 85) VolAvgB, IonLarmor, jpar_perp, jparPS_perp
 85   format(
     1   ' Volume Average B      = ',f14.6,' [T]',/
     2   ' Ion Larmor Radius     = ',f14.6,' [M] X Ti(keV)**0.5',/
     3   ' <J||**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/
     4   ' <JPS**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/)

      write (nthreed, 90)
   90 format(//' More Geometric and Physics Quantities',1x,71('-')/,5x,
     1   'j',3x,'psi-psiaxis',9x,'a [M]',3x,'ellipticity',3x,
     2   'indentation',7x,'d-shape',9x,'shift',6x,'<J||**2>/',4x,
     3   '<JPS**2>/',/,95x,
     4   '<J-perp**2>',3x,'<J-perp**2>'/,' -----',8(2x,12('-')))

      fac = twopi*hs*signgs
      psi(1) = zero
      allocate (r3v(ns-1))
      r3v(:ns-1) = fac*phip(2:ns)*iotas(2:ns)
      do i = 1, ns - 1
         psi(1+i) = psi(i) + r3v(i)
      end do
      deallocate (r3v)

      ygeo(1) = zero
      do js = 2, ns
         zmin =  HUGE(zmin)
         zmax = -HUGE(zmax)
         xmin =  HUGE(xmin)
         xmax = -HUGE(xmax)
         rzmax = zero

c                            !Theta = 0 to pi in upper half of X-Z plane
         noff = 1            !!nphi-plane, noff = 1,....,nzeta
         do icount = 1,2
            n1 = noff        !!nphi-plane, n1 = noff,...,nzeta
            if (icount .eq. 2)
     1      n1 = mod(nzeta + 1 - noff,nzeta) + 1           !!(twopi-v)
            loff = js + ns*(n1-1)
            t1 = one
            if (icount .eq. 2) t1 = -one
            do itheta = 1,ntheta2
               yr1u = r1(loff,0) + sqrts(js)*r1(loff,1)
               yz1u = z1(loff,0) + sqrts(js)*z1(loff,1)
               yz1u = t1*yz1u
               if (yz1u .ge. zmax) then
                  zmax = abs(yz1u)
                  rzmax = yr1u
               else if (yz1u .le. zmin) then
                  zmin = yz1u
                  rzmin = yr1u
               end if
               if (yr1u .ge. xmax) then
                  xmax = yr1u
                  zxmax = yz1u
               else if (yr1u .le. xmin) then
                  xmin = yr1u
                  zxmin = yz1u
               end if
               loff = loff + ns*nzeta
            end do
         end do


         lpi = ns*nzeta*(ntheta2 - 1)
         xmida = r1(js+lpi,0) + sqrts(js)*r1(js+lpi,1)
         xmidb = r1(js,0)     + sqrts(js)*r1(js,1)
!
         rgeo = p5*(xmidb + xmida)              !Geometric major radius
         ygeo(js) = p5*(xmidb - xmida)          !Geometric minor radius
c
         yinden(js) = (xmida - xmin)/(xmax - xmin) !Geometric indentation
         yellip(js) = (zmax - zmin)/(xmax - xmin)  !Geometric ellipticity
c
         ytrian(js) = (rgeo - rzmax)/(xmax - xmin) !Geometric triangularity
         yshift(js) = (r1(1,0)-rgeo)/(xmax - xmin) !Geometric shift
c
         if (jperp2(js) .eq. zero) jperp2(js) = epsilon(jperp2(js))
         jpar_perp = jpar2(js)/jperp2(js)
         if (js .lt. ns) then
            jparPS_perp = jPS2(js)/jperp2(js)
         else
            jparPS_perp = zero
         end if
         write (nthreed, 120) js, psi(js), ygeo(js), yellip(js),
     1      yinden(js), ytrian(js), yshift(js), jpar_perp, jparPS_perp

      end do
  120 format(1x,i5,6f14.5,1p3e14.2)

      write (nthreed, 130)
  130 format(//,' Magnetic Fields and Pressure',/,1x,71('-'))
      sumbpo = zero
      sumbtr = zero
      do i = 2, nrzt
         js = mod(i - 1,ns) + 1
         ncol = (i - 1)/ns + 1
         btorsq = (r12(i)*bsupv(i))**2
         bsq1 = bsupu(i)*bsubu(js,ncol,0) + bsupv(i)*bsubv(js,ncol,0)
         sumbpo = sumbpo + vnorm*tau(i)*(bsq1 - btorsq)
         sumbtr = sumbtr + vnorm*tau(i)*btorsq
      end do
      fac = p5/dmu0
      write (nthreed, 140) sump/dmu0, pavg/dmu0, fac*sumbpo, fac*sumbpo/
     1   volf, fac*sumbtr, fac*sumbtr/volf, fac*sumbto, fac*sumbto/volf,
     2   c1p5*sump/dmu0, c1p5*pavg/dmu0
  140 format(' Volume Integrals (Joules) and Volume ',
     1   'Averages (Pascals)',/,24x,'Integral',6x,'Average',/,
     2   ' pressure         = ',1p2e14.6,/,' bpol**2 /(2 mu0) = ',
     3   1p2e14.6,/,' btor**2/(2 mu0)  = ',1p2e14.6,/,
     4   ' b**2/(2 mu0)     = ',1p2e14.6,/,' EKIN (3/2p)      = ',
     5   1p2e14.6,/)

      write (nthreed, 800)
  800 format(/,' MAGNETIC AXIS COEFFICIENTS'/,
     1   '    n     rmag       zmag        rmag        zmag',/,
     2   '          (cc)       (cs)        (cs)        (cc)',/)
      n1 = 1
      zmags(1,n1) = zero                         !Used for pushing PHIEDGE
      do n = 0, ntor
         n1 = n + 1
         t1 = mscale(0)*nscale(n)
         if (lasym) then
            write (nthreed, 820) n, t1*rmags(1,n1), (-t1*zmags(1,n1)),
     1                             -t1*rmaga(1,n1),   t1*zmaga(1,n1)
         else
            write (nthreed, 820) n, t1*rmags(1,n1), (-t1*zmags(1,n1))
         end if
      end do
  820 format(i5,1p4e12.4)

      betpol = two*sump/sumbpo
      sump20 = two*sump
      sump2 = zero
      sump2 = sum(pres(2:ns)*pres(2:ns)*vp(2:ns)*vnorm)

      betstr = two*sqrt(sump2/volf)/(sumbto/volf)
      betatot = sump20/sumbto
      betapol = betpol
      betator = sump20/sumbtr

      write (nthreed, 150) betatot, betapol, betator
  150 format(' From volume averages over plasma, betas are',/,
     1   ' beta total    = ',f14.6,/,' beta poloidal = ',f14.6,/,
     2   ' beta toroidal = ',f14.6,/)

      write (nthreed, 160) 1.e-6_dp*ctor/dmu0, rbtor, betaxis, betstr
  160 format(' Toroidal Current     = ',f14.6,'  [MA]',/
     1   ' R * Btor-vac         = ',f14.6,' [Wb/M]',/,
     2   ' Peak Beta            = ',f14.6,/,' Beta-star            = ',
     3   f14.6,/)

      if (lrecon) then
!
!
!     Shafranov surface integrals s1,s2
!     Plasma Physics vol 13, pp 757-762 (1971)
!     Also, s3 = .5*S3, defined in Lao, Nucl. Fusion 25, p.1421 (1985)
!     Note: if ctor = 0, use Int(Bsupu*Bsubu dV) for ctor*ctor/R
!
      if (lfreeb) then
        factor = zero                             !Compute current-like norm
        do l = ns, nrzt, ns
           js = mod(l - 1,ns) + 1
           ncol = (l - 1)/ns + 1
           factor = factor + twopi*wint(l)*abs(bsubu(js,ncol,0))
        end do
        factor = one/factor**2
        facnorm = factor*twopi*twopi

        allocate (redg1u(nznt), rbps1u(nznt))
        redg1u(:nznt) = r1(ns:nznt*ns:ns,0) + r1(ns:nznt*ns:ns,1)
        rbps1u(:nznt) = two*facnorm*redg1u(:nznt)*(bpolvac(:nznt) +
     1    presf(ns))*wint(ns:nznt*ns:ns)
        sigr0 = dot_product(rbps1u(:nznt),zu0(ns:nznt*ns:ns))
        sigr1 = dot_product(rbps1u(:nznt)*zu0(ns:nznt*ns:ns),
     1                    redg1u(:nznt))
        sigz1 = -dot_product(rbps1u(:nznt)*ru0(ns:nznt*ns:ns),
     1           z1(ns:nznt*ns:ns,0) + z1(ns:nznt*ns:ns,1))
        deallocate (redg1u, rbps1u)

        er = sigr1 + sigz1
        rlao = volf/(twopi*area)               !LAO, NUCL.FUS.25(1985)1421
        flao = rshaf/rlao
        fgeo = rshaf/rcen
        factor = two*factor/rshaf

        smaleli = factor*sumbpo
        betai = two*factor*sump
        musubi = vnorm*factor*musubi
        dmusubi_meas = two*twopi*factor*delphid*rbtor
        lambda = p5*smaleli + betai
        s11 = (er - rshaf*sigr0)/rshaf         !Shafranov def. based on RT
        s12 = (er - rcen*sigr0)/rcen               !R = Rgeometric
        s13 = (er - rlao*sigr0)/rlao               !R = RLao
        s2 = sigr0
        s3 = sigz1/rshaf
        delta1 = zero
        delta2 = one - fgeo
        delta3 = one - flao
        write (nthreed, 170) rshaf, rcen, rlao, dmusubi_meas, delta1,
     1   delta2, delta3, s11, s12, s13, s2, s2, s2, s3, s3*fgeo, s3*flao
     2   , smaleli, smaleli*fgeo, smaleli*flao, musubi, musubi*fgeo,
     3   musubi*flao, betai, betai*fgeo, betai*flao, musubi + s11,
     4   musubi*fgeo + s12 + s2*(one - fgeo), musubi*flao + s13 + s2*(
     5   one - flao), lambda, fgeo*dlam, flao*dlam, 0.5*s11 + s2, 0.5*(
     6   s12 + s2*(one + fgeo)), 0.5*(s13 + s2*(one + flao)), (3*betai
     7    + smaleli - musubi)*0.5/(s11 + s2) - one, fgeo*(3*betai +
     8   smaleli - musubi)*0.5/(s12 + s2) - one, flao*(3*betai + smaleli
     9    - musubi)*0.5/(s13 + s2) - one, (betai + smaleli + musubi)*0.5
     .   /s2 - one, fgeo*(betai + smaleli + musubi)*0.5/s2 - one,
     .   flao*(betai + smaleli + musubi)*0.5/s2 - one
  170 format(' Integrals of Shafranov',/,1x,22('-'),/,
     1   ' RT (Flux-weighted)   = ',f14.6,' [M]',/,
     2   ' RG (Geometric)       = ',f14.6,' [M]',/,
     3   ' RL (Vol/2*pi*Area)   = ',f14.6,' [M]',/,
     4   ' Mui (diamagnetism)   = ',f14.6,2/,32x,'R = RT',12x,'R = RG',
     5   12x,'R = RL',/,20x,3(10x,8('-')),/,' delta = 1 - RT/R     = ',3
     6   (f14.6,4x),/,' s1                   = ',3(f14.6,4x),/,
     7   ' s2                   = ',3(f14.6,4x),/,
     8   ' s3                   = ',3(f14.6,4x),/,
     9   ' Li                   = ',3(f14.6,4x),/,
     .   ' Mui                  = ',3(f14.6,4x),/,
     1   ' Betai (Calculated)   = ',3(f14.6,4x),/,
     2   ' Betai (Mui + s1)     = ',3(f14.6,4x),/,
     3   ' Lambda (Calculated)  = ',3(f14.6,4x),/,
     4   ' Lambda (s1/2 + s2)   = ',3(f14.6,4x),/,
     5   ' 1st Shafr''v relation = ',3(f14.6,4x),/,
     6   ' (3*Betai + Li - Mui)/[2*(s1+s2)] - 1',/,
     7   ' Radial force balance = ',3(f14.6,4x),/,
     8   ' (Betai + Li + Mui)/(2*s2) - 1',/)

      end if
!
!     Safety Factors (q)
!
      qzero = one/abs(iotaf(1))
      qedge = one/abs(iotaf(ns))
      write (nthreed, 180) qzero, qedge, qzero/qedge
  180 format(' Safety Factors',/,1x,14('-'),/,' q (on axis) = ',f12.4,/,
     1   ' q (at edge) = ',f12.4,/,' q(0)/qedge  = ',f12.4,/)

!
!     PRINT OUT IOTA, PRESSURE SPLINE COEFFICIENTS
!     (PRESSURE IN MKS UNITS, NWT/M**2)
!
      write (nthreed, 190)
  190 format(/,' SUMMARY OF IOTA AND PRESSURE SPLINES'/,
     1   '  K   Spline Node       IOTA(K)      IOTA"(K)'/,
     2   '        sqrt(s)',/,3('-'),3(4x,10('-')))
      do i = 1, isnodes
         write (nthreed, 200) i, sknots(i), ystark(i), y2stark(i)
      end do
  200 format(i3,1p3e14.3)
      write (nthreed, 210)
  210 format(/,'  K   Spline Node       PRES(K)      PRES"(K)'/,
     1   '        sqrt(s)',/,3('-'),3(4x,10('-')))
      factor = pthommax
      do i = 1, ipnodes
         write (nthreed, 200) i, pknots(i), factor*ythom(i), factor*
     1      y2thom(i)
      end do

!
!     PRINT-OUT MAGNETIC PITCH DATA
!
      nchimse = imse
      total_mse_chi = zero

      if (imse .ne. 0) then

         write (nthreed, 220)
  220    format(//,4x,'N     R(data)      Spline',3x,
     1      'atan[BZ/BT] atan[BZ/BT] Chi-Sq-Err',4x,
     2      'Sigma       BZ(in)     BZ(calc)       BT        q(in)',/,
     3      12x,'[m]',5x,'sqrt(s)-node',2x,'(in-deg)',2x,'(calc-deg)',
     4      20x,3(9x,'[T]')/,2x,3('-'),1x,2(3x,9('-')),8(2x,10('-'))/)

         msewgt = zero
         denwgt = zero
         do n = 1, imse + 1
            nsort = isortr(n)
            pitchc = atan(starkcal(n))/dcon
            pitchm = atan(datastark(nsort))/dcon
            wght0 = atan(sigma_stark(nsort))/dcon
            js = indexs1(n)
            lt = indexu1(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            rcalc = (one - delso1(n))*rjs + delso1(n)*rjs1
            if (delse1(n) .lt. zero) rcalc = zero
            if (rcalc .ne. 0.) btor = fpsical(n)/rcalc
            bzin = tan(dcon*pitchm)*btor
            bzcalc = tan(dcon*pitchc)*btor
            chisq = zero
            if (abs(rcalc - router) .lt. epstan) then
               write (nthreed, 230) n, rcalc, rsort0(n), pitchm,
     1            pitchc, wght0, one/(qcalc(n) + epstan)
            else
               if (abs(rsort(n) - rcalc) .le. epstan) then
                  chisq = ((pitchc - pitchm)/wght0)**2
                  msewgt = msewgt + chisq
                  denwgt = denwgt + (pitchm/wght0)**2
               endif
               write (nthreed, 240) n, rcalc, rsort0(n), pitchm,
     1            pitchc, chisq, wght0, bzin, bzcalc, btor,
     2            one/(qmeas(n) + epstan)
            endif
         end do

         chiwgt = msewgt/(imse)
         msewgt = sqrt(msewgt/denwgt)
c                                             !total chi-squared for mse
         total_mse_chi = (imse)*chiwgt
         write (nthreed, 250) chiwgt, msewgt

  230    format(i5,' ',1p2e12.3,1p2e12.3,12x,1pe12.3,6x,
     1      '- Outer (phantom) Edge -',6x,1pe12.3)
  240    format(i5,' ',1p2e12.3,1p8e12.3)
  250    format(/' MEAN CHI-SQ ERROR IN STARK DATA MATCH : ',1pe10.3,/,
     1      ' RMS ERROR IN STARK DATA MATCH : ',1pe10.3,2/)

      endif

!
!     PRINT-OUT PRESSURE DATA
!
      nchipres = itse
      total_pres_chi = zero

      if (lpprof) then
         tswgt = zero
         denwgt = zero
         write (nthreed, 300) presfac*pfac
  300    format(4x,'N     R(in)      R(calc)',4x,f6.2,
     1      ' X Pres(in)      Pres(calc)','  Chi-Sq Err       Sigma'/,2x
     2      ,3('-'),2(2x,10('-')),2(6x,12('-')),2(3x,9('-')),/)

         do n = 1, itse
            js = indexs2(n)
            lt = indexu2(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            if (delso2(n) .eq. (-one)) delso2(n) = one
            rcalc = (one - delso2(n))*rjs + delso2(n)*rjs1
            wght0 = sigma_thom(n)
            chisq = ((datathom(n)*pfac-pcalc(n)/dmu0)/wght0)**2
            tswgt = tswgt + chisq
            denwgt = denwgt + (datathom(n)/wght0)**2
            write (nthreed, 310) n, rthom(n), rcalc, presfac*pfac*
     1         datathom(n), pcalc(n)/dmu0, chisq, wght0
         end do

         chiwgt = tswgt/(itse)
         tswgt = sqrt(tswgt/denwgt)
         total_pres_chi = (itse)*chiwgt !total
         write (nthreed, 320) chiwgt, tswgt
  310    format(i5,1p2e12.3,1p2e18.3,1p2e12.3)
  320    format(/' MEAN CHI-SQ ERROR IN PRESSURE DATA MATCH: ',1pe10.3,/
     1      ,' RMS ERROR IN PRESSURE DATA MATCH: ',1pe10.3/)
      endif

!
!     SUMMARIZE MAGNETICS DATA AND MATCH
!
      call magnetics_data

!
!     COMPUTE REAL TOROIDAL CURRENT ALONG MIDPLANE (MULTIPLY IZETA BY R/SQRT(G))
!     IN PHI = 0 PLANE
!
      call getcurmid (curmid, izeta, gsqrt, r12)

      do nout = nthreed, nmac, (nmac - nthreed)
         if (nout .eq. nmac) then
            if (.not.lmac) cycle
            write (nout, *)
     1      'FOLLOWING DATA EQUALLY SPACED IN TOROIDAL FLUX'
         end if
         write (nout, 700)
         if (nout .eq. nthreed) write (nout, 710)
         iotas(1) = two*iotas(2) - iotas(3)
         iotas(ns+1) = two*iotas(ns) - iotas(ns1)
         vp(1) = two*vp(2) - vp(3)
         vp(ns+1) = two*vp(ns) - vp(ns1)
         pres(1) = two*pres(2) - pres(3)
         pres(ns+1) = two*pres(ns) - pres(ns1)
         do icount = 1, 2*ns - 1
            js = mod(imid(icount) - 1,ns) + 1
            ageo(icount) = ygeo(js)
            phimid(icount) = torflux*(js - 1)/(ns1)
            psimid(icount) = psi(js)
            volpsi(icount) = vnorm*sum(vp(2:js))
            if (js .eq. 1) volpsi(icount) = zero
            dlogidsi = (iotas(js+1)-iotas(js))*ohs/iotaf(js)
            vprime = p5*vnorm*ohs*(vp(js)+vp(js+1))
            pprime = (pres(js+1)-pres(js))*ohs
            rgeo = rmid(icount) + ygeo(js)
            if (icount .gt. ns) rgeo = rgeo - two*ygeo(js)
            alfa(icount) = -two*pprime*vprime*sqrt(two*volpsi(icount)/
     1         rgeo/twopi)/(iotaf(js)*phipf(js))**2
            shear(icount) = -two*volpsi(icount)*dlogidsi/vprime
            write (nout, 720) rmid(icount), ygeo(js), psi(js),
     1         volpsi(icount), qmid(icount), shear(icount),
     2         presmid(icount),
     3         alfa(icount), curmid(icount), datamse(icount)
         end do
      end do
  700 format(/3x,'RMID[M]',6x,'AMID[M]',6x,'PSI[Wb]',4x,'VOL[M**3]',9x,
     1   'q(r)',2x,'SHEAR(Spsi)',6x,'P(PASC)',6x,'ALF(P'')',2x,
     2   'JTOR(A/M**2)',1x,'ATAN(Bz/Bt)[deg]')
  710 format(10('-'),9(3x,10('-')))
  720 format(1pe10.3,9(3x,1pe10.3))

!     Determine q-min on the s grid by min-splining on the sqrt-s grid
      nptsx = 2*ns - 1

      wmidx(:nptsx) = one
      tenmidx(:nptsx) = 0.1_dp
      rmidx(:nptsx) = rmid(:nptsx)
      qmidx(:nptsx) = qmid(:nptsx)

      tol = .005_dp
      sqmin = fmax(zero,one,splints,tol)
      iotamax = splints(sqmin)
      qmin3 = -99999.0_dp
      if (iotamax .ne. zero) qmin3 = one/iotamax
      sqmin = sqmin**2

!     Determine q-min on a fine r-midplane grid
      tol = .005_dp
!     outboard side
      rqmin1 = fmin(rmid(ns),rmid(nptsx),splintx,tol)
      qmin1 = splintx(rqmin1)
      rqmin2 = fmin(rmid(1),rmid(ns),splintx,tol)!outboard side only
      qmin2 = splintx(rqmin2)
      write (nthreed, 730) qmin1, rqmin1, qmin2, rqmin2, qmin3, sqmin
  730 format(//' MINIMUM Q :     OUTBOARD SIDE: QMIN= ',f6.3,'  AT R= ',
     1   f6.3,/,'                  INBOARD SIDE: QMIN= ',f6.3,'  AT R= '
     2   ,f6.3,/,'                    IN S SPACE: QMIN= ',f6.3,
     3   '  AT S= ',f6.3,/)

      nchiiota0 = 0
      total_chi_square = total_b_chi + total_saddle_chi + total_pres_chi
     1    + total_mse_chi + total_chi_cur + total_chi_dia + chisqerr(
     2   islope0)                        !Total CHISQ of equilibrium fit
      nchitot = nbfldn + nchiiota0 + nchisaddle + nchipres + nchimse +
     1   nchicur + nchidia
      total_chi_square_n = total_chi_square/max(1,nchitot)
      write (nthreed, 900)
      write (nthreed, 901) (bloopnames(n),nbfld(n),b_chi(n),b_chi(n)/
     1   max(1,nbfld(n)),n=1,nbsets)
      write (nthreed, 902) nchisaddle, total_saddle_chi,
     1   total_saddle_chi/max(1,nchisaddle), nchipres, total_pres_chi,
     2   total_pres_chi/max(1,nchipres), nchimse, total_mse_chi,
     3   total_mse_chi/max(1,nchimse), nchicur, total_chi_cur,
     4   total_chi_cur/max(1,nchicur), nchidia, total_chi_dia,
     5   total_chi_dia, nchitot, total_chi_square, total_chi_square_n

  900 format(/,' CHI-SQUARED SUMMARY:',t25,'Type',t50,'Nchi',t65,'ChiSq'
     1   ,t84,'ChiSq/N'/,t25,'----',t50,'----',t65,'-----',t84,'-------'
     2   )
  901 format(t25,'B-Loops-',a,t50,i3,t60,1pe12.4,t80,1pe12.4)
  902 format(t25,'Saddle',t50,i3,t60,1pe12.4,t80,1pe12.4,/,t25,
     1   'Pressure',t50,i3,t60,1pe12.4,t80,1pe12.4,/,t25,'MSE',t50,i3,
     2   t60,1pe12.4,t80,1pe12.4,/,t25,'Ip',t50,i3,t60,1pe12.4,t80,1p
     3   e12.4,/,t25,'Diamagnetic Flux',t50,i3,t60,1pe12.4,t80,1pe12.4,/
     4   ,t25,'TOTAL',t50,i3,t60,1pe12.4,t80,1pe12.4)

      end if          !!IF(LRECON)


      end subroutine eqfor


      subroutine elongation (r1, z1, waist, height)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), intent(out) :: waist(2), height(2)
      real(rprec), intent(in), dimension(ns,nzeta,ntheta3,0:1) :: r1, z1
      integer :: nv, n1
C-----------------------------------------------
!
!     Compute Waist thickness, Height in phi = 0, pi symmetry planes
!
      n1 = 0
      do nv = 1, nzeta/2+1
         if (nv.ne.1 .and. nv.ne.nzeta/2+1) cycle
         n1 = n1+1
         waist(n1) = (r1(ns,nv,1,0)       + r1(ns,nv,1,1)) -
     1               (r1(ns,nv,ntheta2,0) + r1(ns,nv,ntheta2,1))
         height(n1) = 2*maxval(z1(ns,nv,:,0) + z1(ns,nv,:,1))
      end do

      end subroutine elongation

      
      subroutine getcurmid (curmid, izeta, gsqrt, r12)
      use vmec_input, only: rprec, dp, nzeta
      use vmec_dim, only: ns, ns1, ntheta2
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: curmid(2*ns)
      real(rprec) :: izeta(ns,nzeta,*), gsqrt(ns,nzeta,*), 
     1    r12(ns,nzeta,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: midcur(ns)
C-----------------------------------------------
!     THETA = pi, PHI = 0      
      midcur(2:ns) = r12(2:ns,1,ntheta2)/gsqrt(2:ns,1,ntheta2)
      
      curmid(1) = izeta(ns,1,ntheta2)*midcur(ns)
      curmid(2:ns1) = 0.5_dp*izeta(ns1:2:-1,1,ntheta2)*
     1                   (midcur(ns1:2:-1) + midcur(ns:3:-1))
  
!     THETA = 0, PHI = 0      
      midcur(2:ns) = r12(2:ns,1,1)/gsqrt(2:ns,1,1)

      curmid(ns+1:2*ns-1) = 0.5_dp*izeta(2:ns1,1,1)*
     1                   (midcur(2:ns1) + midcur(3:ns))

      curmid(ns) = 0.5_dp*(curmid(ns-1) + curmid(ns+1))
      curmid(2*ns) = 2*curmid(2*ns-1) - curmid(2*ns-2)

      end subroutine getcurmid
            

      subroutine getfsq(gcr, gcz, gnormr, gnormz, gnorm, mprecon)
      use vmec_main
      use vmec_params, only: ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mprecon
      real(rprec) gnormr, gnormz, gnorm
      real(rprec), dimension(ns,mnsize,ntmax) :: gcr, gcz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jsmax
C-----------------------------------------------
      jsmax = ns1 + mprecon
      gnormr = gnorm * sum(gcr(:jsmax,:,:)**2)
      gnormz = gnorm * sum(gcz(:jsmax,:,:)**2)

      end subroutine getfsq
      

      subroutine jxbforce(bsupu, bsupv, bsubu, bsubv, bsubs, bsubsu,
     1   bsubsv, itheta, brho, gsqrt, izeta, bsq)
      use safe_open_mod
      use vmec_main
      use vmec_params, only: mscale, nscale, signgs
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,nznt), intent(in) ::
     1  bsupu, bsupv, bsq, gsqrt
      real(rprec), dimension(ns,nznt,0:1), intent(inout) ::
     1  bsubu, bsubv
      real(rprec), dimension(ns,nznt), intent(inout) :: bsubs
      real(rprec), dimension(ns,nznt), intent(out) ::
     1  itheta, brho, izeta
      real(rprec), dimension(ns,nznt,0:1), intent(out) ::
     1  bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      logical, parameter :: lbsubs = .true.       !!True  to use NEW bsubs calculation (from mag. diff. eq.)
                                                  !!False to use OLD bsubs calculation (from metrics)
      logical, parameter :: lprint = .false.      !!Prints out bsubs spectrum to fort.33
      real(rprec), parameter :: two=2, p5=0.5_dp, c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer lk, lz, lt, k, m, js, j, n, injxbout, mparity
      integer :: njxbout = jxbout0, mnyq, nnyq, nmin
      integer, parameter :: ns_skip = 1, nu_skip = 1, nv_skip = 1
      real(rprec), dimension(:,:), allocatable ::
     1    bdotj, bsubuv, bsubvu, brhomn
      real(rprec), dimension(:,:,:), allocatable :: bsubsmn
      real(rprec), dimension(:), allocatable     :: jperpu, jperpv, 
     2    sqgb2, jp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1
        real(rprec), dimension(:,:), allocatable   :: bsubua, bsubva
      real(rprec) ::
     1    bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,
     1    bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,
     2    dnorm1, tsini1, tsini2, tcosi1, tcosi2, tcosm1, tcosm2,
     3    tcosn1, tcosn2, tsinm1, tsinm2, tcos1, tcos2, tsin1, tsin2,
     4    tsinn1, tsinn2, avforce, pprime, amaxfor, aminfor, tjnorm,
     5    ovp, pnorm, brho00(ns)
      real(rprec), dimension(:,:), allocatable ::
     1    bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a
      character :: jxbout_file*100
      real(rprec), external :: dot_g
C-----------------------------------------------
      jxbout_file = 'jxbout.'//input_extension
      call safe_open(njxbout, injxbout, jxbout_file, 'replace',
     1    'formatted')
      if (injxbout .ne. 0) then
         print *,' Error opening JXBOUT file in jxbforce'
         return
      end if

      write (njxbout,6) (ns1-1)/ns_skip, ntheta2/nu_skip, nzeta/nv_skip,
     1    mpol, ntor
 6    format(/,' Radial surfaces = ',i3, ' Poloidal grid points = ',i3,
     1         ' Toroidal grid points = ',i3,/,
     2         ' Poloidal modes = ',i3,' Toroidal modes = ', i3)

!
!     PROGRAM FOR COMPUTING LOCAL JXB = grad-p FORCE BALANCE
!
!     Compute u (=theta), v (=zeta) derivatives of B sub s
!
      mnyq = max(0, ntheta1/2)
      nnyq = max(0, nzeta/2)
      write (njxbout, 5)
 5    format(' LEGEND:',/,
     1  " U = VMEC poloidal angle, V = VMEC (geometric) toroidal angle"/
     2  " SQRT(g') = SQRT(g-VMEC) / V': Jacobian based on VOLUME",/,
     3  " V' = dV/ds: differential volume element",/,
     4  " Es = SQRT(g') [grad(U) X grad(V)] : covariant radial",
     4  " unit vector based on volume",/,
     5  " BSUP{U,V} = B DOT GRAD{U,V}:  contravariant components of B"/,
     6  " JSUP{U,V} = SQRT(g') J DOT GRAD{U,V}",/,
     7  " J X B = Es DOT [J X B]: covariant component of J X B force",/,
     8  " J * B = J DOT B * SQRT(g')",/,
     9  " p' = dp/dV: pressure gradient based on volume derivative",//)

      lz = nzeta*ntheta2
      allocate (bdotj(ns,nznt), bsubuv(ns,nznt),
     1          bsubvu(ns,nznt), jperpu(nznt), jperpv(nznt),
     2          sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq),jp2(nznt),
     3          bsubua(nznt,0:1), bsubva(nznt,0:1), jxb(nznt), 
     4          jxb2(nznt), bsupu1(nznt), bsupv1(nznt), bsubu1(nznt), 
     5          bsubv1(nznt), bsubsmn(ns,0:mnyq,-nnyq:nnyq),
     6          bsubs_s(lz,0:1), bsubs_a(lz,0:1),
     7          bsubu_s(lz,0:1), bsubu_a(lz,0:1),
     8          bsubv_s(lz,0:1), bsubv_a(lz,0:1), stat = j)
      if (j .ne. 0) stop 'allocation error in jxbforce'


!
!     NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
!
      bsubsu = 0; bsubsv = 0; bsubuv = 0; bsubvu = 0; bdotj  = 0

!
!     IT IS ASSUMED THAT MSCALE, NSCALE = [1/sqrt(2), 1, ....] HERE
!     IF NOT, MUST PUT THESE M,N DEPENDENT NORMALIZATION FACTORS INTO DNORM1
!
      radial: do js = 1, ns
         if (js.gt.1 .and. js.lt.ns) then     !!Put on full mesh
            bsubs(js,:) = p5*(bsubs(js,:) + bsubs(js+1,:))
         end if
         bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
         bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
         bsubua = 0;   bsubva = 0

         if (lasym)  then
            call fsym_fft (bsubs(js,:), bsubu(js,:,:), bsubv(js,:,:),
     1             bsubs_s, bsubu_s, bsubv_s, bsubs_a, bsubu_a, bsubv_a)
         else
            bsubs_s(:,0) = bsubs(js,:); bsubu_s = bsubu(js,:,:)
            bsubv_s = bsubv(js,:,:)
         end if

         do m = 0, mnyq
            mparity = mod(m, 2)
            do n = 0, nnyq
               dnorm1 = 4*mscale(m)*nscale(n)        !! One factor mscale*nscale cancels in trig...
               if (m.eq.mnyq) dnorm1 = p5*dnorm1
               if (n.eq.nnyq .and. n.ne.0) dnorm1 = p5*dnorm1
               bsubsmn1 = 0;  bsubsmn2 = 0
               if (lasym) bsubsmn3 = 0;  bsubsmn4 = 0
               if (m.gt.mpol1 .or. n.gt.ntor) goto 222
               bsubumn1 = 0;  bsubumn2 = 0;  bsubvmn1 = 0;  bsubvmn2 = 0
               if (lasym)
     1         bsubumn3 = 0;  bsubumn4 = 0;  bsubvmn3 = 0;  bsubvmn4 = 0

               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1
                     tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1
                     bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk,0)
                     bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk,0)
                     bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity)
                     bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity)
                     bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity)
                     bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity)

                     if (lasym) then
                     bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk,0)
                     bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk,0)
                     bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity)
                     bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity)
                     bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity)
                     bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity)
                     end if

                     lk = lk + nzeta
                  end do
               end do

               dnorm1 = one/(mscale(m)*nscale(n))

!
!              Compute on half u grid (must add symmetric, antisymmetric parts for lasym)
! 
               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tcos1 = cosmu(j,m)*cosnv(k,n)*dnorm1
                     tcos2 = sinmu(j,m)*sinnv(k,n)*dnorm1
!
!                    MUST DEALIAS BSUBU,V IF LBSUBS = TRUE
!
                     bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 +
     1                  tcos2*bsubumn2
                     bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 +
     1                  tcos2*bsubvmn2

                     tcosm1 = cosmum(j,m)*cosnv(k,n)*dnorm1
                     tcosm2 = sinmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)*dnorm1
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     bsubvu(js,lk) = bsubvu(js,lk) + dnorm1*(
     1                               sinmum(j,m)*cosnv(k,n)*bsubvmn1 +
     2                               cosmum(j,m)*sinnv(k,n)*bsubvmn2)
                     bsubuv(js,lk) = bsubuv(js,lk) + dnorm1*(
     1                               cosmu(j,m)*sinnvn(k,n)*bsubumn1 +
     2                               sinmu(j,m)*cosnvn(k,n)*bsubumn2)

                     if (lasym) then
                     tsin1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsin2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 +
     1                  tsin2*bsubumn4
                     bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 +
     1                  tsin2*bsubvmn4

                     tsinm1 = sinmum(j,m)*cosnv(k,n)*dnorm1
                     tsinm2 = cosmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                   tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)*dnorm1
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                   tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                     end if
                     lk = lk + nzeta
                  end do
               end do

 222           continue

!
!        FIX MULTIPLIER FOR M != 0, N != 0 (p5)
!
               if (m .eq. 0) then
                  bsubsmn(js,m,n) =-bsubsmn2
                  if (n.ne.0) bsubsmn(js,m,-n)= 0
               else if (n .eq. 0) then
                  bsubsmn(js,m,n) = bsubsmn1
               else
                  bsubsmn(js,m,n)  = p5*(bsubsmn1 - bsubsmn2)
                  bsubsmn(js,m,-n) = p5*(bsubsmn1 + bsubsmn2)
               end if

            end do
         end do

         if (lasym) call fsym_invfft (bsubua, bsubva, 1)
         bsubu(js,:,0) = bsubua(:,0)
         bsubv(js,:,0) = bsubva(:,0)

      end do radial

      if (lasym) call fsym_invfft (bsubsu, bsubsv, ns)

      if (.not.lbsubs) go to 1500          !!SKIPS Bsubs Correction - uses Bsubs from metric elements

!
!     Compute corrected B-sub-s (impacts currents)
!
      correct_bsubs: do js = 2, ns-1
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
         brho(js,:) = ohs*
     1   ( bsupu1(:)*(bsubu(js+1,:,0) - bsubu(js,:,0))
     2   + bsupv1(:)*(bsubv(js+1,:,0) - bsubv(js,:,0)))
     3   + (pres(js+1) - pres(js))*ohs*jxb(:)
         brho00(js) = DOT_G(nznt,brho(js,1),ns,wint(js),ns)
         brho(js,:) = brho(js,:) - signgs*jxb(:)*brho00(js)/
     1      (p5*(vp(js) + vp(js+1)))

         jxb(:) = brho(js,:)
         call getbrho (brhomn, jxb, bsupu1, bsupv1, mnyq, nnyq)
         if (lprint) then
            write (33, *) ' JS = ', js
            write (33, *) '  M    N       BSUBS(old)        BSUBS(new)'
            do m = 0, mnyq
               nmin = -nnyq
               if (m .eq. 0) nmin = 0
               do n = nmin, nnyq
                  write(33,1223) m, n, bsubsmn(js,m,n), brhomn(m,n)
               end do
            end do
         end if
 1223    format (i4,1x,i4,2(6x,1pe12.3))

!
!        RECOMPUTE bsubsu,v now using corrected bsubs
!
         itheta(js,:) = bsubsu(js,:,0)        !!Store old values here
         izeta (js,:) = bsubsv(js,:,0)
         bsubsu(js,:,:) = 0
         bsubsv(js,:,:) = 0

         do m = 0, mnyq
            do n = 0, nnyq
               dnorm1 = one/(mscale(m)*nscale(n))
               if (n .eq. 0) then
                  bsubsmn1 = brhomn(m,n)
                  bsubsmn2 = 0
               else
                  bsubsmn1 = brhomn(m,n) + brhomn(m,-n)
                  bsubsmn2 =-brhomn(m,n) + brhomn(m,-n)
               end if
               do k = 1, nzeta
                  lk = k
                  do j = 1, ntheta2
                     tcosm1 = cosmum(j,m)*cosnv(k,n)*dnorm1
                     tcosm2 = sinmum(j,m)*sinnv(k,n)*dnorm1
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)*dnorm1
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)*dnorm1
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     lk = lk + nzeta
                  end do
               end do
            end do
         end do

      end do correct_bsubs

      if (lasym) call fsym_invfft (bsubsu, bsubsv, ns)

!
!     CHECK FORCE BALANCE: sqrt(g)*(bsupu*bsubsu + bsupv*bsubsv) = brho
!
      check_fb: do js = 2, ns-1
         bsupu1(:) = p5*(bsupu(js,:)*gsqrt(js,:)
     1    + bsupu(js+1,:)*gsqrt(js+1,:))
         bsupv1(:) = p5*(bsupv(js,:)*gsqrt(js,:)
     1    + bsupv(js+1,:)*gsqrt(js+1,:))
!        jp2(:) = bsupu1(:)*bsubsu(js,:,0) + bsupv1(:)*bsubsv(js,:,0)
!        jxb(:) = bsupu1(:)*itheta(js,:) + bsupv1(:)*izeta(js,:)

!        print *,'JS = ',js
!        do lk = 1, nznt
!           write(*,1224) lk, brho(js,lk), jxb(lk), jp2(lk)
!        end do

!        pause

      end do check_fb

 1224    format ('lk = ',i4,' brho(rhs) = ', 1pe12.4,
     1   ' B dot grad Bs(old) = ', 1pe12.4, ' B dot grad Bs(new) = ',
     2   1pe12.4)

 1500 continue

      deallocate (bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a)
!
!     Now compute currents on the FULL radial mesh
!     Itheta = sqrt(g) * Jsupu, Izeta = sqrt(g) * Jsupv,
!     where Jsupx = J dot grad(x)
!     jxb = J X B   bdotj = sqrt(g)*J dot B
!     jperp-x = (B X gradp) dot grad(x) / |B|**2, x=(u,v)
!     Here, we compute |j-perp|**2 = (j-perp-u)**2 * guu + ...
!     This was compared to the alternative expression (agreed very well):
!     |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
!
!     Note: Multiply currents, pressure by 1/mu0 to get in mks units!
!     Note: IZETA <=> GVV OVERWRITES GVV!!
!
      do js = 2, ns1
         ovp = two/(vp(js+1) + vp(js))
         tjnorm = ovp*signgs
         pprime = ovp*ohs*(pres(js+1)-pres(js))/dmu0
         sqgb2(:nznt) = (gsqrt(js+1,:nznt)*(bsq(js+1,:nznt)- pres(js+1))
     1                +  gsqrt(js,:nznt)  *(bsq(js,:nznt) - pres(js)))
         jperpu(:nznt) = p5*(bsubv(js+1,:nznt,0) +
     1                       bsubv(js,:nznt,0))*pprime/sqgb2
         jperpv(:nznt) =-p5*(bsubu(js+1,:nznt,0) +
     1                       bsubu(js,:nznt,0))*pprime/sqgb2
         jp2(:nznt)=p5*(jperpu**2*(guu(js+1:nrzt:ns) + guu(js:nrzt:ns))
     1        + two*jperpu*jperpv*(guv(js+1:nrzt:ns) + guv(js:nrzt:ns))
     2        +         jperpv**2*(gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))
         itheta(js,:nznt) = bsubsv(js,:nznt,0)
     1                   - ohs*(bsubv(js+1,:nznt,0) - bsubv(js,:nznt,0))
         izeta(js,:nznt) = -bsubsu(js,:nznt,0)
     1                   + ohs*(bsubu(js+1,:nznt,0) - bsubu(js,:nznt,0))
         itheta(js,:nznt) = itheta(js,:nznt)/dmu0
         izeta(js,:nznt)  = izeta(js,:nznt)/dmu0
         jxb(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:nznt) = p5*(bsupu(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupu(js,:nznt)  *gsqrt(js,:))  / jxb(:)
         bsupv1(:nznt) = p5*(bsupv(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupv(js,:nznt)  *gsqrt(js,:))  / jxb(:)
         bsubu1(:nznt) = p5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))
         bsubv1(:nznt) = p5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))
         jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt)
     1              -      izeta (js,:nznt) * bsupu1(:nznt))
         bdotj(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) +
     1                     izeta (js,:nznt) * bsubv1(:nznt)
         pnorm = one/(abs(pprime) + epsilon(pprime))
         amaxfor = max(maxval(jxb(:nznt)-pprime)*pnorm, zero)
         aminfor = min(minval(jxb(:nznt)-pprime)*pnorm, zero)
         avforce = sum(wint(2:nrzt:ns)*(jxb(:nznt) - pprime))


!        Compute <j dot B>, <B sup v> = signgs*phip
!        jpar2 = <j||**2>, jperp2 = <j-perp**2>

         jdotb(js) = tjnorm*sum(bdotj(js,:nznt)*wint(2:nrzt:ns))
         bdotb(js) = tjnorm*sum(sqgb2(:nznt)*wint(2:nrzt:ns))
         bdotgradv(js) = p5*(phip(js) + phip(js+1))*tjnorm
         jpar2(js) = tjnorm*sum(bdotj(js,:nznt)**2*wint(2:nrzt:ns)
     1                         /sqgb2(:nznt))
         jperp2(js) = tjnorm*sum(jp2(:nznt)*wint(2:nrzt:ns)*
     1                        p5*(gsqrt(js+1,:nznt) + gsqrt(js,:nznt)))

         if (mod(js,ns_skip) .eq. 0) then
            amaxfor = min(amaxfor,9.999_dp)
            aminfor = max(aminfor,-9.999_dp)
            write (njxbout, 200) phi(js)/phi(ns), avforce, jdotb(js),
     1         bdotgradv(js), 100.0_dp*amaxfor, 100.0_dp*aminfor
            write (njxbout, 90)
            do lz = 1, nzeta, nv_skip
               write (njxbout, 100) 360._dp*(lz-1)/nzeta, lz
               do lt = 1, ntheta2, nu_skip
                  lk = lz + nzeta*(lt - 1)
                  write (njxbout, 110) lt, ovp*itheta(js,lk),
     1              ovp*izeta(js,lk), ovp*(bsubuv(js,lk) -
     2              bsubvu(js,lk))/dmu0, bsupu1(lk), bsupv1(lk),
     3              jxb(lk), pprime,
     4        (jxb(lk) - pprime), ovp*bdotj(js,lk), bsubu(js,lk,0),
     5              bsubv(js,lk,0), bsubs(js,lk)
               end do
            end do
         endif
      end do

      close (njxbout)

      izeta(1,:nznt) = two*izeta(2,:nznt) - izeta(3,:nznt)           !!For print out in wrout
      izeta(ns,:nznt)= two*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)     !!For print out in wrout
      jdotb(1) = two*jdotb(2) - jdotb(3)
      jdotb(ns) = two*jdotb(ns-1) - jdotb(ns-2)
      bdotb(1) = two*bdotb(3) - bdotb(2)
      bdotb(ns) = two*bdotb(ns-1) - bdotb(ns-2)
      bdotgradv(1) = two*bdotgradv(2) - bdotgradv(3)
      bdotgradv(ns) = two*bdotgradv(ns-1) - bdotgradv(ns-2)
      jpar2(1)   = 0; jpar2(ns)  = 0; jperp2(1)  = 0; jperp2(ns) = 0

      deallocate (jperpu, jperpv, sqgb2, jp2, brhomn, bsubsmn, bsubua,
     1    bsubva, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, stat = j)
!
!     COMPUTE MERCIER CRITERION
!
      bdotj = dmu0*bdotj
      call Mercier(gsqrt,bsq,bdotj,iotas,wint,r1,ru,rv,zu,zv,bsubu,
     1             vp,phips,pres,ns,nznt)


   90 format(/"   LU      JSUPU      JSUPV      JSUPS      BSUPU",
     1   "      BSUPV      J X B       p'    J X B - p'     J * B",
     2   "      BSUBU      BSUBV      BSUBS   "/)
  100 format( " TOROIDAL ANGLE (PER PERIOD) = ", f8.3," DEGREES",
     1        " (PLANE #", i3,")")
  110 format(i5,1p12e11.3)
  200 format(/" TOROIDAL FLUX = ",1pe12.3,3x,"<J X B - p'> = ",
     1   1pe12.3,3x,"<J DOT B> = ",1pe12.3,3x,
     2   "<B DOT GRAD(V)> = ",1pe12.3,/,
     3   " MAXIMUM FORCE DEVIATIONS (RELATIVE TO p'): ",sp,0pf7.2,"%",
     4     3x,f7.2,"%")

      deallocate (bdotj, bsubuv, bsubvu, stat = j)

      end subroutine jxbforce


      subroutine fsym_fft (bs, bu, bv, bs_s, bu_s, bv_s, 
     1    bs_a, bu_a, bv_a)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nzeta,ntheta3), intent(in) :: bs
      real(rprec), dimension(nzeta,ntheta3,0:1), intent(in) :: bu, bv
      real(rprec), dimension(nzeta,ntheta2,0:1), intent(out) :: 
     1   bs_s, bu_s, bv_s, bs_a, bu_a, bv_a
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ir, i, jk, jka, mpar, kz, kr
C-----------------------------------------------
 
!
!       SYMMETRIZE CONTRAVARIANT COMPONENTS OF B
!       SO COS,SIN INTEGRALS CAN BE PERFORMED ON HALF-THETA INTERVAL
!
!       bs_s(v,u) = .5*( bs_s(v,u) - bs_s(-v,-u) )     ! * sin(mu - nv)
!       bs_a(v,u) = .5*( bs_s(v,u) + bs_s(-v,-u) )     ! * cos(mu - nv)
!
!     
      do mpar = 0, 1
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            if (i == 1) ir = 1
            do kz = 1, nzeta
               kr = ireflect(ns*kz)/ns           !-zeta
               bs_a(kz,i,mpar) = cp5*(bs(kz,i)+bs(kr,ir))
               bs_s(kz,i,mpar) = cp5*(bs(kz,i)-bs(kr,ir))
               bu_a(kz,i,mpar) = cp5*(bu(kz,i,mpar)-bu(kr,ir,mpar))
               bu_s(kz,i,mpar) = cp5*(bu(kz,i,mpar)+bu(kr,ir,mpar))
               bv_a(kz,i,mpar) = cp5*(bv(kz,i,mpar)-bv(kr,ir,mpar))
               bv_s(kz,i,mpar) = cp5*(bv(kz,i,mpar)+bv(kr,ir,mpar))
            end do
         end do
      end do
     
      end subroutine fsym_fft


      subroutine fsym_invfft (bsubsu, bsubsv, ns)
      use vmec_main, only: rprec, nzeta, ntheta1, ntheta2, 
     1    ntheta3, ireflect
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns
      real(rprec), dimension(ns*nzeta,ntheta3,0:1),
     1   intent(inout) :: bsubsu, bsubsv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ir, i, jkz, jkr
C-----------------------------------------------

      do i = 1 + ntheta2, ntheta1
         ir = ntheta1 + 2 - i                 !-theta
         do jkz= 1, ns*nzeta
            jkr = ireflect(jkz)               !-zeta
            bsubsu(jkz,i,0) = bsubsu(jkr,ir,0) - bsubsu(jkr,ir,1)
            bsubsv(jkz,i,0) = bsubsv(jkr,ir,0) - bsubsv(jkr,ir,1)
         end do
      end do

      bsubsu(:,:ntheta2,0)=bsubsu(:,:ntheta2,0) + bsubsu(:,:ntheta2,1)
      bsubsv(:,:ntheta2,0)=bsubsv(:,:ntheta2,0) + bsubsv(:,:ntheta2,1)

      end subroutine fsym_invfft


      subroutine getbrho (brhomn, frho, bsupu, bsupv, mmax, nmax)
      use kind_spec
      use vmec_input, only: nfp, nzeta
      use vmec_dim, only: ntheta1, ntheta2, ntheta3
      use vmec_persistent, only: cosmu, sinmu, cosnv, sinnv
      use vmec_params, only: mscale, nscale
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: mmax, nmax
      real(rprec), intent(out) :: brhomn(0:mmax, -nmax:nmax)
      real(rprec), dimension(nzeta, ntheta3), intent(in) ::
     1    bsupu, bsupv, frho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, m, n, nmax1, itotal, ijtot, mntot
      real(rprec) :: amn, dnorm, ccmn, ssmn, dm, dn, termsc, termcs
      real(rprec), allocatable :: amatrix(:,:), save_matrix(:,:),
     1   brhs(:)
      logical :: lpior0
C-----------------------------------------------
      external solver
C-----------------------------------------------
!
!     Solves the radial force balance B dot Brho = Fs for Brho in Fourier space
!     Here, Fs = frho(mn) is the Fourier transform sqrt(g)*F (part of radial force
!     balance sans the B dot Brho term)
!

      nmax1 = max(0,nmax-1)
      itotal = ntheta2*nzeta - 2*nmax1
      allocate (amatrix(itotal, itotal),
     1      brhs(itotal), save_matrix(itotal, itotal), stat=m)
      if (m .ne. 0) stop 'Allocation error in getbrho'
      

      amatrix = 0

!
!     BRHO = BSC(M,N)*SIN(MU)COS(NV) + BCS(M,N)*COS(MU)SIN(NV) 
!          + BCC(M,N)*COS(MU)COS(NV) + BSS(M,N)*SIN(MU)SIN(NV)   (ONLY IF ASYMMETRIC MODE ALLOWED)
!

      ijtot = 0
      brhs = 0
      do i = 1, ntheta2
         do j = 1, nzeta
!           IGNORE u=0,pi POINTS FOR v > pi: REFLECTIONAL SYMMETRY         
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            if (lpior0) cycle
            ijtot = ijtot + 1
            brhs(ijtot) = frho(j,i)
            mntot = 0
            do m = 0, mmax
               do n = 0, nmax
                  dnorm = 1._dp/(mscale(m)*nscale(n))
                  ccmn = cosmu(i,m)*cosnv(j,n)*dnorm
                  ssmn = sinmu(i,m)*sinnv(j,n)*dnorm
                  dm = m * bsupu(j,i)
                  dn = n * bsupv(j,i) * nfp
                  termsc = dm*ccmn - dn*ssmn
                  termcs =-dm*ssmn + dn*ccmn
                  if (n.eq.0 .or. n.eq.nmax) then      
                     mntot = mntot + 1
                     if (m .gt. 0) then
                        amatrix(ijtot,mntot) = termsc     !!only bsc != 0 for n=0, nmax1
                     else if (n .eq. 0) then
                        amatrix(ijtot,mntot) = bsupv(j,i) !!pedestal for m=0,n=0 mode, which should = 0
                     else 
                        amatrix(ijtot,mntot) = termcs     !!bcs(m=0,n=nmax)
                     end if
                  else if (m.eq.0 .or. m.eq.mmax) then
                     mntot = mntot + 1
                     amatrix(ijtot,mntot) = termcs        !!only bcs != 0 for m=0,mmax
                  else 
                     amatrix(ijtot,mntot+1) = termsc
                     amatrix(ijtot,mntot+2) = termcs
                     mntot = mntot + 2
                  end if   
               end do
            end do   
         end do
      end do   

      save_matrix = amatrix
      
      if (ijtot .ne. itotal .or. mntot .ne. itotal) then
         print *,' itotal = ', itotal,' ijtot = ', ijtot,
     1   ' mntot = ', mntot
         stop
      end if
      if (mmax+1 .ne. ntheta2) stop 'Error 1 in getbrho'
      if (nmax   .ne. nzeta/2) stop 'Error 2 in getbrho'

      call solver (amatrix, brhs, itotal)
      
!
!     CHECK SOLUTION FROM SOLVER
!
      
      ijtot = 0
      do i = 1, ntheta2
         do j = 1, nzeta
            lpior0 = ((i.eq.1 .or. i.eq.ntheta2) .and. (j.gt.nzeta/2+1))
            if (lpior0) cycle
            ijtot = ijtot + 1
            amn = sum(save_matrix(ijtot,:)*brhs(:))
            if (abs(frho(j,i) - amn) .gt. 1.e-8_dp*abs(amn))
     1      print *,' i = ',i,' j = ',j,' Original force = ',
     2      frho(j,i),' Final force = ', amn
         end do
      end do   
      
!
!     CONVERT TO BS*SIN(MU - NV) REPRESENTATION
!
      mntot = 0
      brhomn = 0
      do m = 0, mmax
         do n = 0, nmax
            if (n.eq.0 .or. n.eq.nmax) then     
               mntot = mntot + 1
               if (m .gt. 0) then
                  if (n .eq. 0) then
                     brhomn(m,n) = brhs(mntot)         !!bcs(m,0), m > 0
                  else   
                     brhomn(m,n) = p5*brhs(mntot)      !!bsc(m,nmax), m > 0
                     brhomn(m,-n)= p5*brhs(mntot)
                  end if   
               else if (n .ne. 0) then
                  brhomn(m,n) = -brhs(mntot)            !!bcs(0,nmax): needed for derivative
               end if   
            else if (m.eq.0 .or. m.eq.mmax) then
               mntot = mntot + 1
               if (m .eq. 0) then
                  brhomn(m,n)  = -brhs(mntot)          !!bcs(0,n) for n !=0, nmax
               else 
                  brhomn(m,n)  = -p5*brhs(mntot)       !!bcs(mmax,n)
                  brhomn(m,-n) =  p5*brhs(mntot)
               end if
            else
               brhomn(m,n)  = p5*(brhs(mntot+1) - brhs(mntot+2))
               brhomn(m,-n) = p5*(brhs(mntot+1) + brhs(mntot+2))
               mntot = mntot + 2
            end if
         end do
      end do      

      if (mntot .ne. ijtot) stop 'mntot != ijtot at end of getbrho'


      deallocate (amatrix, save_matrix, brhs)

      end subroutine getbrho
      

      subroutine mercier(gsqrt, bsq, bdotj, iotas, wint,
     1           r1, rt, rz, zt, zz, bsubu, vp, phips, pres, ns, nznt)
      use safe_open_mod
      use vmercier
      use vmec_input, ONLY: input_extension
      use vparams, only: one, zero, twopi, nmercier0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns, nznt
      real(rprec), dimension(ns,nznt), intent(in) ::
     1     gsqrt, bsq
      real(rprec), dimension(ns,nznt), intent(inout) :: bdotj
      real(rprec), dimension(ns*nznt), intent(in) :: wint, bsubu
      real(rprec), dimension(ns,nznt,0:1), intent(in) ::
     1     r1, rt, rz, zt, zz
      real(rprec), dimension(ns), intent(in) ::
     1     iotas, vp, phips, pres
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ns1, i, imercier0, nmerc = nmercier0
      real(rprec) :: sign_jac, hs, sqs, denom
      real(rprec), dimension(:,:), allocatable :: 
     1     gpp, gsqrt_full, b2
      real(rprec), dimension(nznt) :: gtt, rtf, ztf,
     1   rzf, zzf, r1f, jdotb, ob2, b2i
      real(rprec), dimension(ns) :: vp_real, phip_real, 
     1   shear, vpp, presp, torcur, ip, sj, tpp, tjb, tbb, tjj
      character*120 mercier_file
      real(rprec), external :: dot_g
C-----------------------------------------------

      mercier_file = 'mercier.'//trim(input_extension)
      call safe_open (nmerc, imercier0, mercier_file, 'replace',
     1   'formatted')
      if (imercier0 .ne. 0) return


      allocate (gpp(ns,nznt), gsqrt_full(ns,nznt), b2(ns,nznt), stat=i)
      if (i .ne. 0) stop 'allocation error in Mercier'

!
!     SCALE VP, PHIPS TO REAL UNITS (VOLUME, TOROIDAL FLUX DERIVATIVES)
!     AND PUT GSQRT IN ABS UNITS (SIGNGS MAY BE NEGATIVE)      
!     NOTE: VP has (coming into this routine) the sign of the jacobian multiplied out
!          i.e., vp = signgs*<gsqrt>
!     THE SHEAR TERM MUST BE MULTIPLIED BY THE SIGN OF THE JACOBIAN
!     (OR A BETTER SOLUTION IS TO RETAIN THE JACOBIAN SIGN IN ALL TERMS, INCLUDING
!      VP, THAT DEPEND EXPLICITLY ON THE JACOBIAN. WE CHOOSE THIS LATTER METHOD...)
!
!     COMING INTO THIS ROUTINE, THE JACOBIAN(gsqrt) = 1./(grad-s . grad-theta X grad-zeta)
!     WE CONVERT THIS FROM grad-s to grad-phi DEPENDENCE BY DIVIDING gsqrt by PHIP_REAL
!
!     NOTE: WE ARE USING 0 < s < 1 AS THE FLUX VARIABLE, BEING CAREFUL
!     TO KEEP d(phi)/ds == PHIP_REAL FACTORS WHERE REQUIRED
!     THE V'' TERM IS d2V/d(PHI)**2, PHI IS REAL TOROIDAL FLUX
!
!     SHEAR = d(iota)/d(phi)   :  FULL MESH
!     VPP   = d(vp)/d(phi)     :  FULL MESH
!     PRESP = d(pres)/d(phi)   :  FULL MESH  (PRES IS REAL PRES*mu0)
!     IP    = d(Itor)/d(phi)   :  FULL MESH
!
!     ON ENTRY, BDOTJ = Jacobian * J*B  ON THE FULL RADIAL GRID
!               BSQ = 0.5*|B**2| + p IS ON THE HALF RADIAL GRID
!

      ns1 = ns - 1 
      if (ns1 .le. 0) return
      hs = one/ns1
      sign_jac = zero
      if (gsqrt(ns,1) .ne. zero)
     1    sign_jac = abs(gsqrt(ns,1))/gsqrt(ns,1)

      if (sign_jac .eq. zero) return
      phip_real = twopi * phips * sign_jac
!
!     NOTE: phip_real should be > 0 to get the correct physical sign of real-space gradients
!     For example, grad-p, grad-Ip, etc. However, with phip_real defined this way,
!     Mercier will be correct
!
      vp_real(2:ns) = sign_jac*(twopi*twopi)*vp(2:ns)/phip_real(2:ns)  !!dV/d(PHI) on half mesh
      
!
!     COMPUTE INTEGRATED TOROIDAL CURRENT
!
      do i = 2,ns
         torcur(i)=sign_jac*twopi*dot_g(nznt,bsubu(i),ns,wint(i),ns)
      end do

!
!     COMPUTE SURFACE AVERAGE VARIABLES ON FULL RADIAL MESH
!
      do i = 2,ns1
        phip_real(i) = p5*(phip_real(i+1) + phip_real(i))
        denom     = one/(hs*phip_real(i))
        shear(i)  = (iotas(i+1) - iotas(i))*denom       !!d(iota)/d(PHI)
        vpp(i)    = (vp_real(i+1) - vp_real(i))*denom   !!d(VP)/d(PHI)
        presp(i)  = (pres(i+1) - pres(i))*denom         !!d(p)/d(PHI)
        ip(i)     = (torcur(i+1) - torcur(i))*denom     !!d(Itor)/d(PHI)
      end do  
      
!
!     COMPUTE GPP == |grad-phi|**2 = PHIP**2*|grad-s|**2           (on full mesh)
!             GSQRT_FULL = JACOBIAN/PHIP == jacobian based on flux (on full mesh)
!

      do i = 2, ns1
        gsqrt_full(i,:) = p5*(gsqrt(i,:) + gsqrt(i+1,:))
        bdotj(i,:) = bdotj(i,:)/gsqrt_full(i,:)
        gsqrt_full(i,:) = gsqrt_full(i,:)/phip_real(i)
        sj(i) = hs*(i-1)
        sqs = sqrt(sj(i))
        rtf(:) = rt(i,:,0) + sqs*rt(i,:,1)
        ztf(:) = zt(i,:,0) + sqs*zt(i,:,1)
        gtt(:) = rtf(:)*rtf(:) + ztf(:)*ztf(:)
        rzf(:) = rz(i,:,0) + sqs*rz(i,:,1)
        zzf(:) = zz(i,:,0) + sqs*zz(i,:,1)
        r1f(:) = r1(i,:,0) + sqs*r1(i,:,1)
        gpp(i,:) = gsqrt_full(i,:)**2/(gtt(:)*r1f(:)**2 + 
     1             (rtf(:)*zzf(:) - rzf(:)*ztf(:))**2)     !!1/gpp
      end do
      
!
!     COMPUTE SURFACE AVERAGES OVER dS/|grad-PHI|**3 => |Jac| du dv / |grad-PHI|**2      
!     WHERE Jac = gsqrt/phip_real
!
      do i = 2,ns
        b2(i,:) = two*(bsq(i,:) - pres(i))
      end do
        
      do i = 2,ns1
        b2i(:) = p5*(b2(i+1,:) + b2(i,:))
        ob2(:) = gsqrt_full(i,:)/b2i(:)
        tpp(i) = DOT_G(nznt,ob2,1,wint(i),ns)                !<1/B**2>
        ob2(:) = b2i(:) * gsqrt_full(i,:) * gpp(i,:)
        tbb(i) = DOT_G(nznt,ob2,1,wint(i),ns)                !<b*b/|grad-phi|**3>
        jdotb(:) = bdotj(i,:) * gpp(i,:) * gsqrt_full(i,:)
        tjb(i) = DOT_G(nznt,jdotb,1,wint(i),ns)              !<j*b/|grad-phi|**3>
        jdotb(:) = jdotb(:) * bdotj(i,:) / b2i(:)
        tjj(i) = DOT_G(nznt,jdotb,1,wint(i),ns)              !<(j*b)2/b**2*|grad-phi|**3>
      end do
            
      deallocate (gpp, gsqrt_full, b2, stat=i)

!
!     REFERENCE: BAUER, BETANCOURT, GARABEDIAN, MHD Equilibrium and Stability of Stellarators
!     We break up the Omega-subs into a positive shear term (Dshear) and a net current term, Dcurr
!     Omega_subw == Dwell and Omega-subd == Dgeod (geodesic curvature, Pfirsch-Schluter term)
!
!     Include (eventually) Suydam for reference (cylindrical limit)
!
       
      write(nmerc,90)
 90   format(6x,'S',10x,'PHI',9x,'IOTA',8x,'SHEAR',7x,' VP ',8x,'WELL',
     1       8x,'ITOR',7x,'ITOR''',7x,'PRES',7x,'PRES''',/,120('-'))
      
      do i = 2,ns1
         sqs = p5*(vp_real(i) + vp_real(i+1))*sign_jac
         if (sqs .eq. zero) cycle
         write(nmerc,100) sj(i), hs*sum(phip_real(2:i)),
     1   p5*(iotas(i+1)+iotas(i)), shear(i)/sqs,
     2   sqs, -vpp(i)*sign_jac, 
     3   p5*(torcur(i) + torcur(i+1)), ip(i)/sqs,
     4   p5*(pres(i) + pres(i+1)), presp(i)/sqs
      end do

 100  format(1p10e12.4)      
      
      write(nmerc,190)
 190  format(/,6x,'S',8x,'DMerc',8x,'DShear',7x,'DCurr',7x,'DWell',
     1     7x,'Dgeod',/,100('-'))      
     
      do i = 2,ns1
         tpp(i) = (twopi*twopi)*tpp(i)
         tjb(i) = (twopi*twopi)*tjb(i)
         tbb(i) = (twopi*twopi)*tbb(i)
         tjj(i) = (twopi*twopi)*tjj(i)
         Dshear(i) = shear(i) * shear(i)/4
         Dcurr(i)  =-shear(i) * (tjb(i) - ip(i) *tbb(i))
         Dwell(i)  = presp(i) * (vpp(i) - presp(i) *tpp(i))*tbb(i) 
         Dgeod(i)  = tjb(i) *tjb(i)  - tbb(i) *tjj(i) 
         DMerc(i)  = Dshear(i) + Dcurr(i) + Dwell(i) + Dgeod(i)
         write(nmerc,100) sj(i), Dmerc(i), Dshear(i),
     1         Dcurr(i), Dwell(i), Dgeod(i)
      end do

      close (nmerc)
     
      end subroutine mercier


      subroutine residue(gcr, gcz, gcl, facmul)
      use vmec_main
      use vmec_params, only: meven, modd, ntmax, signgs
      use vsvd
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax) :: 
     1  gcr, gcz, gcl, facmul
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: n0 = 0, m0 = 0, m1 = 1
      real(rprec), parameter :: two = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, nsfix
      real(rprec) :: r1, r2
      real(rprec) :: tnorm, t1, t2, fac, bz1
C-----------------------------------------------
 
!
!       IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!       INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!       ( ZCS = RSS, ZCC = RSC ARE THE CORRECT POLAR RELATIONS; HERE WE USE
!         THE SIMPLER RELATION ZCS = 0 INSTEAD)
!
 
      if (lthreed .and. fsqz.lt.c1pm8 .and. iter1.ne.1) then
         gcz(2:ns,1:ntor,m1,1) = zero
      endif

      if (lasym) then
      gcr(2:,:,m1,4) = gcr(2:,:,m1,4) + gcz(2:,:,m1,3)
      gcz(2:,:,m1,3) = 0
      do n = 0, ntor
         t1 = azd(ns,2) + bzd(ns,2) + (n*nfp)**2*crd(ns)
         t2 = ard(ns,2) + brd(ns,2) + (n*nfp)**2*crd(ns)
         t1 = t2/(t1 + t2)
         gcr(2:,n,m1,4) = t1*gcr(2:,n,m1,4) 
      end do

      n = max (3, mpol1)
      gcr(2:, :, n:,4) = 0              !!NEED THIS TO IMPROVE CONVERGENCE....
      gcz(2:, :, 4:,3) = 0

      end if
!
!       PUT FORCES INTO PHIFSAVE UNITS THAT PRECONDITIONERS, FNORM ARE IN
!
      if (phifac .eq. zero) then
         stop 'phifac = 0 in residue'
      else
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      end if   
 
 
      if (lrecon) then
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS
 
        r1 = sum(gcr(:ns,n0,m0,1))
        if (r0scale .eq. zero) stop 'r0scale = 0'
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        if (iopt_raxis.gt.0 .and. iresidue.eq.2 
     1     .and. fsq.lt.fopt_axis) iresidue = 3
        if (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      else
!
!     ADJUST PHIEDGE
!
         if (imovephi .gt. 0) call movephi1 (gphifac)
      endif
      gc(neqs1) = gphifac
!
!       CONSTRUCT INVARIANT RESIDUALS
!
 
      call getfsq (gcr, gcz, fsqr, fsqz, fnorm, meven)
 
      r2 = sum(gcr(ns,:ntor,:mpol1,:)*gcr(ns,:ntor,:mpol1,:)
     1   +     gcz(ns,:ntor,:mpol1,:)*gcz(ns,:ntor,:mpol1,:))
      fedge = fnorm*r2
!
!       PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      call scalfor (gcr, arm, brm, ard, brd, crd)
      call scalfor (gcz, azm, bzm, azd, bzd, crd)
!
!     APPLY ZCC = RSC CONSTRAINT FOR M=1 MODES
!
      if (lasym) gcz(2:,:,m1,3) = gcr(2:,:,m1,4)

      call getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, modd)
      fac = one/(one + (fsqr + fsqz))
      gcr = fac*gcr
      gcz = fac*gcz

!
!     CONSTRUCT INVARIANT AND PRECONDITIONED LAMBDA FORCES
!     NOTE: IF PHIP**2 USED IN BCOVAR FOR BU, BV COMPONENTS (RATHER THAN PHIP)
!     MUST MULTIPLY RBTOR IN BZ1 BY ADDITIONAL HS*SUM(ABS(PHIPS(2:NS))) FACTOR
!
      bz1 = two*hs/(rbtor*tnorm)**2
      fsql = bz1*sum(gcl*gcl)
      gcl = facmul*gcl
      fsql1 = hs*sum(gcl*gcl)
 
      if (fsqr.gt.1.e+4_dp .or. fsqz.gt.1.e+4_dp) then
         if (iter2.eq.iter1) gcl = 2.e-1_dp*gcl
         if (iter2.eq.1) irst = 4
      end if   
          
      end subroutine residue
      

      subroutine scalfor(gcx, axm, bxm, axd, bxd, cx)
      use vmec_main
      use vmec_params, only: jmin2, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax), 
     1  intent(inout) :: gcx
      real(rprec), dimension(ns + 1,2), intent(in) :: 
     1  axm, bxm, axd, bxd
      real(rprec), dimension(ns), intent(in) :: cx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: ftol_edge = 1.e-9_dp
      integer :: m , mp, n, js, mmin, nmin, jmax, jmin4(0:mnsize)
      real(rprec), dimension(:,:,:), allocatable :: ax, bx, dx
      logical, save :: ledge
C-----------------------------------------------
      allocate (ax(ns,0:ntor,0:mpol1), bx(ns,0:ntor,0:mpol1),
     1    dx(ns,0:ntor,0:mpol1))

      jmax = ns
      if (ivac .lt. 1) jmax = ns1
!
!     ACCELERATE (IMPROVE) CONVERGENCE OF FREE BOUNDARY. THIS WAS ADDED
!     TO DEAL WITH CASES WHICH MIGHT OTHERWISE DIVERGE. BY DECREASING THE
!     FSQ TOLERANCE LEVEL WHERE THIS KICKS IN (FTOL_EDGE), THE USER CAN
!     TURN-OFF THIS FEATURE
!
      if (fsqr + fsqz .lt. ftol_edge) ledge = .true.
      if (iter2.lt.400 .or. ivac.lt.1) ledge = .false.

      do m = 0, mpol1
         mp = mod(m,2) + 1
         do n = 0, ntor
            do js = jmin2(m), jmax
               ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
               bx(js,n,m) = axm(js,mp) + bxm(js,mp)*m**2
               dx(js,n,m) = axd(js,mp) + bxd(js,mp)*m**2
     1                    + cx(js)*(n*nfp)**2
            end do

            if (m .eq. 1) dx(2,n,m) = dx(2,n,m) + bx(2,n,m)

            if (jmax .lt. ns) then
               dx(ns,n,m) = zero
            else if (m .le. 1) then
               dx(ns,n,m) = 1.05_dp*dx(ns,n,m)   !May need to raise 1.05 for 3D
            else
               dx(ns,n,m) = 1.10_dp*dx(ns,n,m)
            endif
         end do
      end do
      
!     FOR DATA MATCHING MODE (0 <= IRESIDUE < 3),
!     MAGNETIC AXIS IS FIXED SO JMIN3(0) => 2 FOR M=0,N=0
 
      jmin4 = jmin3       
      if (iresidue.ge.0 .and. iresidue.lt.3) jmin4(0) = 2
      
!     DIAGONALIZE (DX DOMINANT) AND REDUCE FORCE (DX ENHANCED) AT EDGE TO IMPROVE CONVERGENCE
 
      if (ledge) then
         bx(ns,:,:) = 0*bx(ns,:,:)
         dx(ns,:,:) = 3*dx(ns,:,:)
      end if
      
      call tridslv (ax, dx, bx, gcx, jmin4, 
     1      jmax, mnsize - 1, ns, ntmax)

      deallocate (ax, bx, dx)

      end subroutine scalfor

      subroutine symforce(ars, brs, crs, azs, bzs, czs, bls, cls, rcs, 
     1   zcs, ara, bra, cra, aza, bza, cza, bla, cla, rca, zca)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), 
     1   intent(inout) :: ars, brs, crs, azs, bzs, czs, 
     2   bls, cls, rcs, zcs
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(out) :: 
     1   ara, bra, cra, aza, bza, cza, bla, cla, rca, zca
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mpar, ir, i, jk, jka
C-----------------------------------------------
 
!
!       SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
!       SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
!
!       ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * cos(mu - nv)
!       ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * sin(mu - nv)
!
!
      do mpar = 0, 1
         do i = 1, ntheta2
            ir = ntheta1 + 2 - i                 !-theta
            if (i == 1) ir = 1
            do jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               ara(jk,i,mpar) = cp5*(ars(jk,i,mpar)-ars(jka,ir,mpar))
               ars(jk,i,mpar) = cp5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
               bra(jk,i,mpar) = cp5*(brs(jk,i,mpar)+brs(jka,ir,mpar))
               brs(jk,i,mpar) = cp5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
               aza(jk,i,mpar) = cp5*(azs(jk,i,mpar)+azs(jka,ir,mpar))
               azs(jk,i,mpar) = cp5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
               bza(jk,i,mpar) = cp5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar))
               bzs(jk,i,mpar) = cp5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
               bla(jk,i,mpar) = cp5*(bls(jk,i,mpar)-bls(jka,ir,mpar))
               bls(jk,i,mpar) = cp5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
               rca(jk,i,mpar) = cp5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar))
               rcs(jk,i,mpar) = cp5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
               zca(jk,i,mpar) = cp5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar))
               zcs(jk,i,mpar) = cp5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  cra(jk,i,mpar) = cp5*(crs(jk,i,mpar)+crs(jka,ir,mpar))
                  crs(jk,i,mpar) = cp5*(crs(jk,i,mpar)-crs(jka,ir,mpar))
                  cza(jk,i,mpar) = cp5*(czs(jk,i,mpar)-czs(jka,ir,mpar))
                  czs(jk,i,mpar) = cp5*(czs(jk,i,mpar)+czs(jka,ir,mpar))
                  cla(jk,i,mpar) = cp5*(cls(jk,i,mpar)-cls(jka,ir,mpar))
                  cls(jk,i,mpar) = cp5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
               end do
            endif
         end do
      end do
 
      end subroutine symforce

        
      subroutine symrzl(r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, 
     1   zcons, r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(inout) :: 
     1   r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, zcons
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(in) :: 
     1   r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mpar, ir, i, jk, jka, n2
C-----------------------------------------------
 
!
!     FIRST SUM SYMMETRIC, ANTISYMMETRIC PIECES ON EXTENDED INTERVAL, THETA = [PI,2*PI]
!
      do mpar = 0, 1
         do i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i                 !-theta
            do jk = 1, ns*nzeta
               jka = ireflect(jk)                !-zeta
               r1s(jk,i,mpar) = r1s(jka,ir,mpar) - r1a(jka,ir,mpar)
               rus(jk,i,mpar) =-rus(jka,ir,mpar) + rua(jka,ir,mpar)
               z1s(jk,i,mpar) =-z1s(jka,ir,mpar) + z1a(jka,ir,mpar)
               zus(jk,i,mpar) = zus(jka,ir,mpar) - zua(jka,ir,mpar)
               lus(jk,i,mpar) = lus(jka,ir,mpar) - lua(jka,ir,mpar)
               rcons(jk,i,mpar)= rcons(jka,ir,mpar)-rcona(jka,ir,mpar)
               zcons(jk,i,mpar)=-zcons(jka,ir,mpar)+zcona(jka,ir,mpar)
            end do
            if (lthreed) then
               do jk = 1, ns*nzeta
                  jka = ireflect(jk)
                  rvs(jk,i,mpar)=(-rvs(jka,ir,mpar))+rva(jka,ir,mpar)
                  zvs(jk,i,mpar) = zvs(jka,ir,mpar) - zva(jka,ir,mpar)
                  lvs(jk,i,mpar) = lvs(jka,ir,mpar) - lva(jka,ir,mpar)
               end do
            endif
         end do
!
!        NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES FOR THETA = [0,PI]
!
         n2 = ntheta2
         r1s(:,:n2,mpar) = r1s(:,:n2,mpar) + r1a(:,:n2,mpar)
         rus(:,:n2,mpar) = rus(:,:n2,mpar) + rua(:,:n2,mpar)
         z1s(:,:n2,mpar) = z1s(:,:n2,mpar) + z1a(:,:n2,mpar)
         zus(:,:n2,mpar) = zus(:,:n2,mpar) + zua(:,:n2,mpar)
         lus(:,:n2,mpar) = lus(:,:n2,mpar) + lua(:,:n2,mpar)
         rcons(:,:n2,mpar) = rcons(:,:n2,mpar) + rcona(:,:n2,mpar)
         zcons(:,:n2,mpar) = zcons(:,:n2,mpar) + zcona(:,:n2,mpar)
         if (lthreed) then
            rvs(:,:n2,mpar) = rvs(:,:n2,mpar) + rva(:,:n2,mpar)
            zvs(:,:n2,mpar) = zvs(:,:n2,mpar) + zva(:,:n2,mpar)
            lvs(:,:n2,mpar) = lvs(:,:n2,mpar) + lva(:,:n2,mpar)
         endif
      end do
 
      end subroutine symrzl


      subroutine totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1, 
     1   lv1, rcn1, zcn1)
      use vmec_main
      use vmec_params, only: jmin1, jlam, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax), 
     1   intent(inout) :: rzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1), intent(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncs = 3, rmnsc = 4
      integer, parameter :: m0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncc, zmnss, lmncc, lmnss
      integer :: m, n, mparity, k, js, i, jk, j1, l
      real(rprec), dimension(:,:), allocatable :: work1
      real(rprec) :: cosmux, sinmux
C-----------------------------------------------

      allocate (work1(ns*nzeta,12))

      zmncc = rmncs + ntmax
      zmnss = rmnsc + ntmax
      lmncc = rmncs + 2*ntmax
      lmnss = rmnsc + 2*ntmax

!
!                INITIALIZATION BLOCK
!
      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0
  
      if (jlam(m0) .gt. 1) then
      rzl_array(1,1:ntor,m0,lmncc) = 2*rzl_array(2,1:ntor,m0,lmncc)
     1                             -   rzl_array(3,1:ntor,m0,lmncc)
      end if

      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = 0
         j1 = jmin1(m)
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = j1, ns
                  work1(js+l,1) = work1(js+l,1) + 
     1               rzl_array(js,n,m,rmnsc)*cosnv(k,n)
                  work1(js+l,6) = work1(js+l,6) + 
     1               rzl_array(js,n,m,zmncc)*cosnv(k,n)
                  work1(js+l,10) = work1(js+l,10) + 
     1               rzl_array(js,n,m,lmncc)*cosnv(k,n)
               end do
               if (lthreed) then
                  do js = j1, ns
                     work1(js+l,2) = work1(js+l,2) + 
     1                  rzl_array(js,n,m,rmncs)*sinnv(k,n)
                     work1(js+l,3) = work1(js+l,3) + 
     1                  rzl_array(js,n,m,rmnsc)*sinnvn(k,n)
                     work1(js+l,4) = work1(js+l,4) + 
     1                  rzl_array(js,n,m,rmncs)*cosnvn(k,n)
                     work1(js+l,5) = work1(js+l,5) + 
     1                  rzl_array(js,n,m,zmnss)*sinnv(k,n)
                     work1(js+l,7) = work1(js+l,7) + 
     1                  rzl_array(js,n,m,zmnss)*cosnvn(k,n)
                     work1(js+l,8) = work1(js+l,8) + 
     1                  rzl_array(js,n,m,zmncc)*sinnvn(k,n)
                     work1(js+l,9) = work1(js+l,9) + 
     1                  rzl_array(js,n,m,lmnss)*sinnv(k,n)
                     work1(js+l,11) = work1(js+l,11) + 
     1                  rzl_array(js,n,m,lmnss)*cosnvn(k,n)
                     work1(js+l,12) = work1(js+l,12) + 
     1                  rzl_array(js,n,m,lmncc)*sinnvn(k,n)
                  end do
               endif
            end do
         end do
 
!
!        INVERSE TRANSFORM IN M-THETA
!
         do i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux
 
            if (lthreed) then
               r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               cosmu(i,m)
               ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               sinmum(i,m)
               z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               sinmu(i,m)
               zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               cosmum(i,m)
               lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               cosmum(i,m)
               rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               cosmux
               zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               sinmux
               rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               sinmu(i,m) + work1(:,4)*cosmu(i,m)
               zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               sinmu(i,m) + work1(:,8)*cosmu(i,m)
               lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               sinmu(i,m)+work1(:,12)*cosmu(i,m))
            endif
 
         end do
      end do
 
      deallocate (work1)

      end subroutine totzspa


      subroutine totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1, 
     1   lv1, rcn1, zcn1)
      use vmec_main
      use vmec_params, only: jmin1, jlam, ntmax
      use vmec_persistent
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,3*ntmax),
     1   intent(inout) :: rzl_array
      real(rprec), dimension(ns*nzeta,ntheta3,0:1),
     1   intent(out) :: r11, ru1, 
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: rmncc = 1, rmnss = 2
      integer, parameter :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: zmncs, zmnsc, lmncs, lmnsc
      integer :: n, m, mparity, k, js, i, j1, l
      real(rprec), dimension(:,:), allocatable :: work1
      real(rprec) :: cosmux, sinmux
C-----------------------------------------------

      zmncs = rmncc + ntmax
      zmnsc = rmnss + ntmax
      lmncs = rmncc + 2*ntmax
      lmnsc = rmnss + 2*ntmax

      allocate (work1(ns*nzeta,12))

      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0; zu1 = 0
      zv1 = 0;  rcn1 = 0; zcn1 = 0
      lu1(:,:,0) = 1;  lu1(:,:,1) = 0;  lv1 = 0

!
!     EXTRAPOLATION AT JS=1 FOR M=ODD MODES AND M=0 MODES OF LAMBDA
!
      rzl_array(1,:,m1,:) = 2*rzl_array(2,:,m1,:)
     1                    -   rzl_array(3,:,m1,:)

      if (jlam(m0) .gt. 1) then
      rzl_array(1,1:ntor,m0,lmncs) = 2*rzl_array(2,1:ntor,m0,lmncs)
     1                               - rzl_array(3,1:ntor,m0,lmncs)
      end if

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     BEGIN INVERSE TRANSFORM IN N-ZETA
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv    
!
      do m = 0, mpol1
         mparity = mod(m,2)
         work1 = 0
         j1 = jmin1(m)
         do n = 0, ntor
            do k = 1, nzeta
               l = ns*(k-1)
               do js = j1,ns
               work1(js+l,1) = work1(js+l,1) + 
     1            rzl_array(js,n,m,rmncc)*cosnv(k,n)
               work1(js+l,6) = work1(js+l,6) + 
     1            rzl_array(js,n,m,zmnsc)*cosnv(k,n)
               work1(js+l,10) = work1(js+l,10) +
     1            rzl_array(js,n,m,lmnsc)*cosnv(k,n)
               end do
               if (lthreed) then
                  do js = j1,ns
                     work1(js+l,2) = work1(js+l,2) + 
     1                  rzl_array(js,n,m,rmnss)*sinnv(k,n)
                     work1(js+l,4) = work1(js+l,4) + 
     1                  rzl_array(js,n,m,rmnss)*cosnvn(k,n)
                     work1(js+l,9) = work1(js+l,9) + 
     1                  rzl_array(js,n,m,lmncs)*sinnv(k,n)
                     work1(js+l,11) = work1(js+l,11) + 
     1                  rzl_array(js,n,m,lmncs)*cosnvn(k,n)
                  end do
                  do js = j1,ns
                     work1(js+l,5) = work1(js+l,5) + 
     1                  rzl_array(js,n,m,zmncs)*sinnv(k,n)
                     work1(js+l,7) = work1(js+l,7) + 
     1                  rzl_array(js,n,m,zmncs)*cosnvn(k,n)
                     work1(js+l,3) = work1(js+l,3) + 
     1                  rzl_array(js,n,m,rmncc)*sinnvn(k,n)
                     work1(js+l,8) = work1(js+l,8) + 
     1                  rzl_array(js,n,m,zmnsc)*sinnvn(k,n)
                     work1(js+l,12) = work1(js+l,12) + 
     1                  rzl_array(js,n,m,lmnsc)*sinnvn(k,n)
                  end do
               endif
            end do
         end do
!
!        INVERSE TRANSFORM IN M-THETA
!
         do i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            cosmux
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            sinmux
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            cosmum(i,m)
 
            if (lthreed) then
               r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               sinmu(i,m)
               ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               cosmum(i,m)
               rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               sinmux
               z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               cosmu(i,m)
               zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               sinmum(i,m)
               zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               cosmux
               lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               sinmum(i,m)
               rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               cosmu(i,m) + work1(:,4)*sinmu(i,m)
               zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               cosmu(i,m) + work1(:,8)*sinmu(i,m)
               lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               cosmu(i,m) + work1(:,12)*sinmu(i,m))
            endif
 
         end do
      end do
 
      deallocate (work1)

      if (rzl_array(1,n0,m1,rmncc) .eq. zero)
     1   stop 'r01(0) = 0 in totzsps'      
!     dkappa = rzl_array(1,n0,m1,zmnsc)/rzl_array(1,n0,m1,rmncc)
!     r01 = rzl_array(ns,n0,m1,rmncc)

      end subroutine totzsps


      subroutine wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu,
     1   izeta, ier_flag)
      use vmec_main
      use vmec_params, only: mscale, nscale, signgs, version_
      use vmercier
      use vmec_persistent
      use vsvd
      use vspline
      use xstuff
      use vmec_io
      use realspace
      use safe_open_mod
      use vacmod, ONLY: potvac,xmpot,xnpot,mnpd   !added for diagno, J.Geiger
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ier_flag
      real(rprec), dimension(ns,nznt), intent(in) ::
     1   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu, izeta
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: two = 2, p5 = 0.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iasym, ireconstruct, i, j, js, lj, mn, lk, isgn,
     1   m, nmin0, n, k, iwout0, l, ntotal, n1, nwout, nfort
      real(rprec), dimension(nznt) :: bmod
      real(rprec), dimension(mnmax) :: rmnc, zmns, lmns,
     1   rmns, zmnc, lmnc, bmodmn, bmodmn1, lmns_half, lmnc_half
      real(rprec) :: dmult, gmn, bmn, currvmn, bsubumn, bsubvmn,
     1   bsubsmn, bsupumn, bsupvmn, tcosi, tsini, presfactor, zeta,
     2   lmnsh, lmnch
!     4  ,potsin, potcos                            !added for diagno, J.Geiger
      character*120 :: wout_file, fort_file
C-----------------------------------------------

      rmns = 0
      zmnc = 0
      lmnc = 0
!
!     THIS SUBROUTINE CREATES THE TEXT FILE WOUT WHICH
!     CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (FULL MESH), LMN (HALF MESH - CONVERTED FROM
!              INTERNAL FULL MESH REPRESENTATION)
!
!     IZETA (FULL), BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF)
!

      wout_file = 'wout.'//input_extension
      nwout = nwout0
      call safe_open(nwout, iwout0, wout_file, 'replace', 'formatted')
      if (iwout0 .ne. 0) stop 'Error writing WOUT file in VMEC WROUT'

      if (loldout) then
         nfort = nfort8
         fort_file = 'fort.8.'//input_extension
         call safe_open(nfort,iwout0,fort_file,'replace','formatted')
         if (iwout0 .ne. 0) then
             print *,'Problem opening fort.8 file'
             print *,'iwout0=',iwout0, 'for ',input_extension
             stop 'Error writing fort.8. file in VMEC WROUT '
         endif

         nfort = nfort18
         fort_file = 'fort.18.'//input_extension
         call safe_open(nfort,iwout0,fort_file,'replace','formatted')
         if (iwout0 .ne. 0) 
     1       stop 'Error writing fort.18. file in VMEC WROUT'
      endif

      if (lasym) then
         iasym = 1                                  ! asymmetric mode
      else
         iasym = 0
      end if
      if (lrecon) then
         ireconstruct = 1
      else
         itse = 0
         imse2 = 0
         ireconstruct = 0
      end if

!
!     Insert version information into wout file. This will be parsed in
!     read_wout_file to return the real value version_ to check the version number.
!
      write (nwout, '(a15,a)') 'VMEC VERSION = ', version_
      write (nwout, *) wb, wp, gamma, pfac,
     1  rmax_surf, rmin_surf, zmax_surf

      write (nwout, *) nfp, ns, mpol, ntor, mnmax, itfsq, niter,
     1  iasym, ireconstruct, ier_flag

      write (nwout, *) imse2 - 1, itse, nbsets, nobd, nextcur,
     1   nstore_seq
      if (nbsets .gt. 0) write (nwout, *) (nbfld(i),i=1,nbsets)
      write (nwout, '(a)') mgrid_file

!-----------------------------------------------
!     Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) then
        write (nfort8, 710) gamma, nfp, ns, mpol, ntor, mnmax, itfsq,
     1         niter/100+1
!       the variables helsym (4th) and kpres (last) were set to 0.
        write (nfort18, 40)voli,gamma,1.0/nfp,0.,mnmax,ns,mpol,ntor,
     1              ntor+1,1,itfsq,niter/100+1,0
      end if
 40   format(1x,1p,e22.12,1x,3e12.5,8i4,i1)
 710  format(f10.3,7i6)

      if (ier_flag.ne.0 .and. ier_flag.ne.4) goto 1000

      pres(1) = pres(2)
      ntotal = 0

      do js = 1, ns
         lj = (js - 1)*mnmax

         call convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, js)

         mn = 0
         bmod = sqrt(two*abs(bsq(js,:nznt)-pres(js)))
         do m = 0, mpol1
            nmin0 = -ntor
            if (m .eq. 0) nmin0 = 0
            do n = nmin0, ntor
               mn = mn + 1
               ntotal = max(mn,ntotal)
               dmult = two/(mscale(m)*nscale(abs(n)))
               if (m .eq. 0 .and. n .eq. 0) dmult = p5*dmult
               n1 = abs(n)
               isgn = sign(1, n)
               gmn = 0
               bmn = 0
               currvmn = 0
               bsubumn = 0
               bsubvmn = 0
               bsubsmn = 0
               bsupumn = 0
               bsupvmn = 0
               do j = 1, ntheta2
                  do k = 1, nzeta
                     lk = k + nzeta*(j - 1)
                     tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     1                         isgn*sinmui(j,m)*sinnv(k,n1))
                     tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     1                         isgn*cosmui(j,m)*sinnv(k,n1))
                     bmn = bmn + tcosi*bmod(lk)
                     gmn = gmn + tcosi*gsqrt(js,lk)
                     currvmn = currvmn + tcosi*izeta(js,lk)
                     bsubumn = bsubumn + tcosi*bsubu(js,lk)
                     bsubvmn = bsubvmn + tcosi*bsubv(js,lk)
                     bsubsmn = bsubsmn + tsini*bsubs(js,lk)
                     bsupumn = bsupumn + tcosi*bsupu(js,lk)
                     bsupvmn = bsupvmn + tcosi*bsupv(js,lk)
                  end do
               end do
               if (js .eq. ns/2) bmodmn(mn) = bmn
               if (js .eq. ns) bmodmn1(mn) = bmn
               if(js .eq. 1) then
                  if (m .ne. nint(xm(mn))) stop 'm != nint(xm)'
                  if (n*nfp .ne. nint(xn(mn))) stop ' n = nint(xn)'
                  write (nwout, *) m, n*nfp
                  raxis(0:ntor,1) = rmnc(1:ntor+1)
                  zaxis(0:ntor,1) = zmns(1:ntor+1)
!                 ZERO HALF-MESH QUANTITIES
                  bsubumn = 0
                  bsubvmn = 0
                  bsubsmn = 0
                  bsupumn = 0
                  bsupvmn = 0
                  bmn = 0
                  gmn = 0
               end if
!
!       PUT LMNS, LMNC ONTO HALF-MESH FOR BACKWARDS CONSISTENCY
!
               if (js .eq. 1) then
                  lmnsh = 0
                  lmnch = 0
               else if (mod(m,2) .eq. 0) then
                  lmnsh = cp5*(lmns(mn) + lmns_half(mn))
                  lmnch = cp5*(lmnc(mn) + lmnc_half(mn))
               else
                  if (js.eq.2 .and. m.eq.1) then
                     lmns_half(mn) = lmns(mn)
                     lmnc_half(mn) = lmnc(mn)
                  end if
                  lmnsh = cp5*(sm(js)*lmns(mn) + sp(js-1)*lmns_half(mn))
                  lmnch = cp5*(sm(js)*lmnc(mn) + sp(js-1)*lmnc_half(mn))
               end if

               write (nwout, *) rmnc(mn), zmns(mn), lmnsh,
     1            bmn, gmn, bsubumn,
     2            bsubvmn, bsubsmn, bsupumn, bsupvmn, currvmn
               if (lasym) then
                  write (nwout, *) rmns(mn), zmnc(mn), lmnch
               endif

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
               if (loldout) then
                 if (js .eq. 1)
     1              write (nfort8, 721) nint(xm(mn)),nint(xn(mn))
                 write (nfort18,50)  xm(mn),xn(mn),rmnc(mn),zmns(mn),gmn
                 write (nfort8, 731) rmnc(mn), zmns(mn), lmns(mn), bmn,
     1            gmn, bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
               end if
            end do
         end do
         lmns_half = lmns        !!Store previous full point values for averaging
         lmnc_half = lmnc
      end do
 50   format(1x,1p,5e14.6)
 721  format(2i10)
 731  format(5e20.13)

!
!     HALF-MESH QUANTITIES (except phi, jcuru, jcurv which are FULL MESH - computed in eqfor)
!     NOTE: jcuru are local current densities, NOT integrated in s
!     NOTE: In version <= 6.00, mass, press are written out in INTERNAL units
!     and should be multiplied by 1/mu0 to transform to pascals. In version > 6.00,
!     the pressure, mass are in correct (physical) units
!
      write (nwout, *) (iotas(js), mass(js)/dmu0, pres(js)/dmu0,
     1   beta_vol(js), phip(js), buco(js), bvco(js), phi(js), vp(js),
     2   overr(js), jcuru(js)/dmu0, jcurv(js)/dmu0, specw(js),js=2,ns)
!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) then
         write (nfort8, 731) (iotas(js),mass(js),pres(js),phips(js),
     1          buco(js),bvco(js),phi(js),vp(js),jcuru(js)/dmu0,
     2          jcurv(js)/dmu0,specw(js),js=2,ns)
        write (nfort18,50)(iotas(js),mass(js),pres(js),-phips(js),
     1          vp(js),js=1,ns)
      end if
!-----------------------------------------------

      write (nwout, *) aspect, betatot, betapol, betator, betaxis, b0

!-----------------------------------------------
!     New output added to version 6.20
!-----------------------------------------------
      write (nwout, *) nint(signgs)
      write (nwout, '(a)') input_extension
      write (nwout, *) IonLarmor, VolAvgB, rbtor0, rbtor, ctor/dmu0,
     1  Aminor_p, Rmajor_p, volume_plasma
!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      write (nwout, *) (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     1       Dgeod(js), equif(js), js=2,ns-1)

      if (nextcur .gt. 0) then
         write (nwout, *) (extcur(i),i=1,nextcur)
         write (nwout, *) (curlabel(i),i=1,nextcur)
      endif

!-----------------------------------------------
!     NOTE: jdotb is in units of A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
!-----------------------------------------------
      write (nwout, *) (fsqt(i),wdot(i),i=1,nstore_seq)
      write (nwout, *) (jdotb(js),bdotgradv(js),js=1,ns)

!-----------------------------------------------
!   Modification to obtain old fort.8 file
!-----------------------------------------------
      if (loldout) write (nfort8, 731) (fsqt(i),wdot(i),i=1,100)

!-----------------------------------------------
!   Write diagno file (J.Geiger)
!-----------------------------------------------
      IF(ldiagno)THEN         !added for diagno, J.Geiger start
         IF(lfreeb .and. .not.lasym)THEN
            nfort = nfort21
            fort_file = 'diagno_in.'//input_extension
            call safe_open(nfort,iwout0,fort_file,'replace',
     1                     'formatted')
            if (iwout0 .ne. 0) 
     1          stop 'Error writing diagno_in. file in VMEC WROUT'

            write(nfort21,'(a)') "vmec2000"
            write(nfort21,*) "nfp  mpol  ntor"
            write(nfort21,*) nfp, mpol, ntor

            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc, ns)
            write(nfort21,*) "rmnc"
            write(nfort21,'(1p,e24.16)') (rmnc(mn),mn=1,mnmax)
            write(nfort21,*) "zmns"
            write(nfort21,'(1p,e24.16)') (zmns(mn),mn=1,mnmax)

            write(nfort21,*) "potsin"
            DO i = 1, mnpd
               write(nfort21,'(1p,e24.16)') potvac(i)
            END DO
            write(nfort21,*) "phiedge"
            write(nfort21,*) phiedge
            write(nfort21,*) "nextcur"
            write(nfort21,*) nextcur
            write(nfort21,*) "external currents"
            write(nfort21,*) extcur(1:nextcur)

            write(nfort21,*) "plasma current"
            write(nfort21,*) ctor
            CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, xc,1)
            write(nfort21,*) "plasma current filament fc R"
            write(nfort21,*) rmnc(1:ntor+1)
            write(nfort21,*) "plasma current filament fc z"
            write(nfort21,*) zmns(1:ntor+1)

            close(unit=nfort21)
         ELSE
            write(6,*)"Diagno-file request not completed!"
            write(6,*)"VMEC2000 not running in free-boundary mode!"
            write(6,*)"-or- LASYM = .true. !"
            write(6,*)"Check mgrid-file and namelist!"
         ENDIF
      ENDIF                   !added for diagno, J.Geiger end

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      if (.not.lrecon) goto 900

      if (imse2 - 1.gt.2 .or. itse.gt.0) then
         write (nwout, *) tswgt, msewgt
         call smoothdata(nwout)

!       These knot values are on sqrt(s) grid
         presfactor = dmu0*pthommax             !!*pfac moved to getthom
         write (nwout, *) isnodes, (sknots(i),ystark(i),y2stark(i),i=
     1      1,isnodes)
         write (nwout, *) ipnodes, (pknots(i),presfactor*ythom(i),
     1      presfactor*y2thom(i),i=1,ipnodes)
         write (nwout, *)(datamse(i),rmid(i),qmid(i),shear(i),
     1      presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
         write (nwout, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
         write (nwout, *)(rthom(i),datathom(i),i=1,itse)
      endif
      if (nobd .gt. 0) then
         write (nwout, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
         write (nwout, *) flmwgt
      endif
      if (nbfldn .gt. 0) then
         do n = 1, nbsets
            write (nwout, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1         i=1,nbfld(n))
         end do
         write (nwout, *) bcwgt
      endif

      write (nwout, *) phidiam, delphid
!
!     Write Limiter & Prout plotting specs
!
      write (nwout, *) nsets, nparts, nlim
      write (nwout, *) (nsetsn(i),i=1,nsets)
      write (nwout, *) (((pfcspec(i,j,k),i=1,nparts),j=1,nsetsn(k)),
     1   k=1,nsets)
      write (nwout, *) (limitr(i), i=1,nlim)
      write (nwout, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      write (nwout, *) nrgrid, nzgrid
      write (nwout, *) tokid
      write (nwout, *) rx1, rx2, zy1, zy2, condif
      write (nwout, *) imatch_phiedge
      if (imatch_phiedge .eq. 2) call wroutlim(nwout)

 900  continue
!-----------------------------------------------
!     FREE BOUNDARY DATA
!-----------------------------------------------
      if (ier_flag .eq. 0) call freeb_data(rmnc, zmns, bmodmn, bmodmn1)


 1000 continue

      close (nwout)
      if( loldout) then
         close(nfort8)
         close(nfort18)
      endif

      end subroutine wrout
            
      
      subroutine freeb_data (rmnc, zmns, bmodmn, bmodmn1)
      use vmec_main
      use vacmod
      use realspace, only: r1, z1
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(mnmax) :: rmnc, zmns, bmodmn, bmodmn1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iprint, nzskip, l, k, lk, mnlim, mnres, mn,
     1           mn0, n, nedge, nedge0 = 99, iu, iv, nl, lkr
      real(rprec), allocatable, dimension(:) :: rb, phib, zb
C-----------------------------------------------
!
!     WRITE OUT EDGE VALUES OF FIELDS TO FORT.NEDGE0 (INCLUDE REFLECTED POINT)
!
!     NOTE: BR, BPHI, BZ WERE COMPUTED IN BSS, CALLED FROM EQFOR
!
      if (.not.lfreeb .and. .not.ledge_dump) return

      allocate (rb(2*nznt), phib(2*nznt), zb(2*nznt), stat=l)
      if (l .ne. 0) stop 'allocation error in freeb_data'

      nedge = 0
      lkr = nznt
      do iv = 1,nzeta
         zeta = (twopi*(iv-1))/(nzeta*nfp)
         do iu = 1,ntheta3
            lk = iv + nzeta*(iu-1)
            nl = ns*lk
            nedge = nedge+1
            rb(lk)   = r1(nl,0) + r1(nl,1)
            phib(lk) = zeta
            zb(lk)   = z1(nl,0) + z1(nl,1)
!
!      INCLUDE -u,-v POINTS HERE BY STELLARATOR SYMMETRY
!
            if (.not.lasym .and. (iu.ne.1 .and. iu.ne.ntheta2)) then
              lkr = lkr + 1
              nedge = nedge+1
              rb(lkr)   = rb(lk)
              phib(lkr) =-phib(lk)
              zb(lkr)   =-zb(lk)  
              bredge(lkr) = -bredge(lk)
              bpedge(lkr) =  bpedge(lk)         
              bzedge(lkr) =  bzedge(lk)
            endif  
         end do
      end do  

      if (ledge_dump) then
        write(NEDGE0,*) 'INPUT FILE = ',arg1
        write(NEDGE0,*) 'NEDGE = ',nedge
        write(NEDGE0,*) 'RB = ',  (rb(i), i=1,nedge)
        write(NEDGE0,*) 'PHIB = ',(phib(i), i=1,nedge)
        write(NEDGE0,*) 'ZB = ',  (zb(i), i=1,nedge)
        write(NEDGE0,*) 'BREDGE = ', (bredge(i), i=1,nedge)
        write(NEDGE0,*) 'BPEDGE = ', (bpedge(i), i=1,nedge)
        write(NEDGE0,*) 'BZEDGE = ', (bzedge(i), i=1,nedge)
      end if

!
!     WRITE OUT (TO THREED1 FILE) VACUUM INFORMATION
!

      if (.not.lfreeb) then
         deallocate (rb, phib, zb, stat=l)
         return
      end if
      
      do iprint = 1, 2
         if (iprint .eq. 1) write (nthreed, 750)
         if (iprint .eq. 2) write (nthreed, 760)
         nzskip = 1 + nzeta/6
         do l = 1, nzeta, nzskip
            zeta = (360.0_dp*(l - 1))/nzeta
            if (iprint .eq. 1) then
               do k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  write (nthreed, 770) zeta, rb(lk),
     1            zb(lk), (bsqsav(lk,n),n=1,3), bsqvac(lk)
               end do
            else
               do k = 1, ntheta2
                  lk = l + nzeta*(k - 1)
                  write (nthreed, 780) zeta, rb(lk), zb(lk), 
     1            bredge(lk), bpedge(lk), bzedge(lk),
     2            brv(lk), bphiv(lk), bzv(lk)
               end do
            endif
         end do
      end do

      deallocate (rb, phib, zb, bredge, bpedge, bzedge, stat=l)

      write (nthreed, 800)
      mnlim = (mnmax/2)*2
      do mn = 1, mnlim, 2
         write (nthreed, 810) nint(xn(mn)/nfp), nint(xm(mn)), rmnc(mn),
     1      zmns(mn), bmodmn(mn), bmodmn1(mn), nint(xn(mn+1)/nfp),
     2      nint(xm(mn+1)), rmnc(mn+1), zmns(mn+1), bmodmn(mn+1),
     3      bmodmn1(mn+1)
      end do

      write (nthreed, 820)

      mnlim = mnpd/3
      do mn = 1, mnlim
         mn0 = 3*mn - 2
         write (nthreed, 830) nint(xnpot(mn0)/nfp), nint(xmpot(mn0)),
     1      potvac(mn0), nint(xnpot(mn0+1)/nfp), nint(xmpot(mn0+1)),
     2      potvac(mn0+1), nint(xnpot(mn0+2)/nfp), 
     3      nint(xmpot(mn0+2)), potvac(mn0+2)
      end do
      mnres = mnpd - mnlim*3
      mn = mnlim*3 + 1
      if (mnres .eq. 2) then
         write (nthreed, 830) nint(xnpot(mn)/nfp), nint(xmpot(mn)),
     1      potvac(mn), nint(xnpot(mn+1)/nfp), nint(xmpot(mn+1)),
     2      potvac(mn+1)
      else if (mnres .eq. 1) then
         write (nthreed, 830) nint(xnpot(mn)/nfp), nint(xmpot(mn)),
     1      potvac(mn)
      endif

  750 format(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BSQMHDI',5x,'BSQVACI'
     1   ,5x,'BSQMHDF',5x,'BSQVACF',/)
  760 format(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',6x,'BR',8x,'BPHI',6x,'BZ'
     1   ,8x,'BRv',7x,'BPHIv',5x,'BZv',/)
  770 format(1pe10.2,1p6e12.4)
  780 format(1pe10.2,1p2e12.4,1p6e10.2)
  790 format(i5,/,(1p3e12.4))
  800 format(//,3x,'nb',2x,'mb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',3x,
     1   '|B|(s=1.)',6x,'nb',2x,'mb',9x,'rbc',9x,'zbs',3x,'|B|(s=.5)',3x
     2   ,'|B|(s=1.)'/)
  810 format(i5,i4,1p4e12.4,3x,i5,i4,1p4e12.4)
  820 format(/,3x,'nf',2x,'mf',5x,'potvacs',6x,'nf',2x,'mf',5x,'potvacs'
     1   ,6x,'nf',2x,'mf',5x,'potvacs'/)
  830 format(i5,i4,1pe12.4,3x,i5,i4,1pe12.4,3x,i5,i4,1pe12.4)

      end subroutine freeb_data
      

      subroutine lamfull(luef, luof, lue, luo)
      use vmec_main
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,*), intent(out) :: luef, luof
      real(rprec), dimension(ns,nzeta,*), intent(in) :: lue, luo
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lt, lt1
C-----------------------------------------------

C
C     COMPUTES LAMBDA ON FULL RADIAL MESH AT THETA=0, PI
C 
      do lt = 1, 2
         lt1 = 1
         if (lt .eq. 2) lt1 = ntheta2
         luof(1,lt1) = p5*(c1p5*luo(2,1,1)-p5*luo(3,1,1)+c1p5*luo(2,1,
     1      ntheta2)-p5*luo(3,1,ntheta2))
         luef(1,lt1) = p5*(c1p5*lue(2,1,1)-p5*lue(3,1,1)+c1p5*lue(2,1,
     1      ntheta2)-p5*lue(3,1,ntheta2))
         luof(ns,lt1) = c1p5*luo(ns,1,lt1) - p5*luo(ns1,1,lt1)
         luef(ns,lt1) = c1p5*lue(ns,1,lt1) - p5*lue(ns1,1,lt1)
         luof(2:ns1,lt1) = p5*(luo(3:ns1+1,1,lt1)+luo(2:ns1,1,lt1))
         luef(2:ns1,lt1) = p5*(lue(3:ns1+1,1,lt1)+lue(2:ns1,1,lt1))
      end do
 
      end subroutine lamfull

      
      subroutine magnetics_data
      use vmec_main
      use vacmod
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iobset, j, i, n, ibfld
      real(rprec) :: denwgt, dsitotal,     
     1   chierr, wght0, abs_flux, chiwgt, raddeg, denbwgt, btotal, 
     2   relerr_b, wghtb, chisqx, chibwgt
C-----------------------------------------------
!
!       ALGEBRAIC SUM OF FLUXES FOR GROUPS OF OBSERVER POSITIONS
!
      if (nobd .ne. 0) then
         flmwgt = 0
         denwgt = 0
         write (nthreed, 400)
  400    format(/' External Poloidal Flux Loop Measurements (Wb)',/,
     1      '   I     Connection          Coils       Plasma    ',
     2      '    Total     Measured   Chi-Sq Err        Sigma',
     3      '   Designation',/,1x,3('-'),1x,16('-'),7(1x,12('-')),'-')
 
         do iobset = 1, nobd
            dsitotal = dsiext(iobset) + plflux(iobset)
            chierr = dsiobt(iobset) - dsitotal
            wght0 = cbig + one
            if (any(iobset.eq.indxflx(:nflxs)))
     1          wght0 = sigma_flux(iobset)
            flmwgt = flmwgt + (chierr/wght0)**2
            denwgt = denwgt + (dsitotal/wght0)**2
            chierr = (chierr/wght0)**2
            if (wght0 .lt. cbig) then
               write (nthreed, 410) iobset, (iconnect(j,iobset),j=1,4), 
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, wght0, dsilabel(iobset)
            else
               write (nthreed, 420) iobset, (iconnect(j,iobset),j=1,4), 
     1            dsiext(iobset), plflux(iobset), dsitotal, dsiobt(
     2            iobset), chierr, dsilabel(iobset)
            endif
         end do
 
         write (nthreed, 600)
         do i = 1, nobser
            abs_flux = plflux(nobd+i)
            write (nthreed, 610) i, xobser(i), zobser(i), psiext(i), 
     1         abs_flux, psiext(i) + abs_flux
         end do
  600    format(//,' Poloidal Flux at Individual Coil Locations',/,
     1      '   I         R [M]        Z [M]        Coils       Plasma',
     2      '        Total',/,2x,3('-'),5(1x,12('-')))
  610    format(i4,1x,5f13.4)
      endif
 
      nchisaddle = nflxs
      total_saddle_chi = 0.
      if (nflxs .GT. 0) then
         chiwgt = flmwgt/(nflxs)
         flmwgt = sqrt(flmwgt/denwgt)
         total_saddle_chi = (nflxs)*chiwgt
         write (nthreed, 430) chiwgt, flmwgt
      endif
  410 format(i4,1x,4i4,6f13.4,7x,a8)
  420 format(i4,1x,4i4,5f13.4,5x,'No Match',7x,a8)
  430 format(/' MEAN CHI-SQ ERROR IN SADDLE LOOP MATCH: ',1pe10.3,/,
     1   ' RMS ERROR IN SADDLE LOOP MATCH: ',1pe10.3)
 
!
!       BR, BZ FIELD COMPONENTS FOR MATCHING AT OBSERVER POSITIONS
!
      total_b_chi = 0.
      if (nbfldn .ne. 0) then
  500    format(//' BFIELD GROUP: ',a,' // ',
     1      'External B Loop Magnetic Field Measurements (T)',/,
     2      '  I  Rcoil[M]  Zcoil[M]    Degree    B[coil]',
     3      '  B[plasma]   B[Total]  ','  B[Meas]    Sigma  Chi-Sq Err',
     4      /,1x,2('-'),3(1x,9('-')),4(1x,10('-')),1(1x,8('-')),1(1x,11(
     5      '-')))
 
         raddeg = 360.0_dp/twopi
         do n = 1, nbsets                        ! over all allowed sets
            write (nthreed, 500) bloopnames(n)
            bcwgt = 0.
            denbwgt = 0.
            b_chi(n) = 0.
            ibfld = 0
            do j = 1, nbcoils(n)
               ibfld = ibfld + 1
               btotal = bcoil(j,n) + plbfld(j,n)
               relerr_b = btotal - bbc(j,n)
               wghtb = cbig + 1.0_dp
c                                                ! over namelist sets
               if (any(j.eq.indxbfld(:nbfld(n),n))) wghtb = sigma_b(j,n)
               bcwgt = bcwgt + (relerr_b/wghtb)**2
               denbwgt = denbwgt + (btotal/wghtb)**2
               chisqx = (relerr_b/wghtb)**2
               if (wghtb .LT. cbig) then
                  write (nthreed, 520) ibfld, rbcoil(j,n), zbcoil(j,n), 
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n), wghtb, chisqx
               else
                  write (nthreed, 530) ibfld, rbcoil(j,n), zbcoil(j,n), 
     1               abcoil(j,n)*raddeg, bcoil(j,n), plbfld(j,n), btotal
     2               , bbc(j,n)
               endif
            end do
            if (nbfld(n) .gt. 0) then
               chibwgt = bcwgt/(nbfld(n))
               bcwgt = sqrt(bcwgt/denbwgt)
               b_chi(n) = (nbfld(n))*chibwgt
               write (nthreed, 540) bloopnames(n), chibwgt, 
     1            bloopnames(n), bcwgt
               total_b_chi = total_b_chi + b_chi(n)
            endif
         end do
 
      endif
  520 format(i3,3f10.3,4f11.3,f9.4,f12.4)
  530 format(i3,3f10.3,4f11.3,3x,'No Match')
  540 format(/' MEAN CHI-SQ ERROR IN B-COIL GROUP ',a,': ',1pe10.3,/,
     1   ' RMS ERROR IN B-COIL GROUP ',a,': ',1pe10.3)
 
      end subroutine magnetics_data


      subroutine movephi1 (gphifacx)   !wiecode
      use vmec_main
      use vsvd
      implicit none 
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: gphifacx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: phifac0
C-----------------------------------------------
!
!     UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!     FDOT/F = -2*(1 - Ipexp/Ipcode), WHERE F = PHISCALE
!
      if (ctor .eq. zero) stop 'ctor = 0 in movephi1'
      phifac0 = phifac*(currv/ctor)
      gphifacx = rsfac*c1pm2*(phifac0 - phifac)
 
      end subroutine movephi1


      subroutine newphi(phipog)
      use kind_spec
      use vmec_main
      use realspace
      use vsvd
      use xstuff
c-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nrzt) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: phifold, tnorm
C-----------------------------------------------

!
!     UPDATE PHI SCALE FACTOR (PHIFAC)
!
      phifold = phifac
      phifac = xc(neqs1)

      if (phifold .eq. zero) stop 'phifac = 0 in newphi'
      tnorm = phifac/phifold
      phip = tnorm * phip
      phipog = tnorm * phipog
 
      end subroutine newphi


      subroutine printout(i0, delt0, w0, lscreen)
      use vmec_main
      use realspace
      use vsvd
!     use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer i0
      real(rprec) :: delt0, w0
      logical lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: betav, w, avm, den
C-----------------------------------------------
      betav = wp/wb
      w = w0*twopi*twopi
      den = zero
      specw(1) = one
      call spectrum (xc(:irzloff), xc(1+irzloff:2*irzloff))
      den = sum(phip(2:ns))
      avm = dot_product(phip(2:ns),specw(2:ns)+specw(:ns-1))
      avm = 0.5_dp*avm/den
      if (ivac .ge. 1 .and. iter2.gt.1) delbsq =
     1     sum(dbsq(:nznt)*wint(2:nrzt:ns))/
     1     sum(bsqsav(:nznt,3)*wint(2:nrzt:ns))
      if (i0.eq.1 .and. lfreeb) then
         if (lscreen) print 20
         if (imatch_phiedge .eq. 1) then
            write (nthreed, 15)
         else
            write (nthreed, 16)
         endif
      else if (i0.eq.1 .and. .not.lfreeb) then
         if (lscreen) print 30
         write (nthreed, 25)
      endif
   15 format(/,' ITER    FSQR      FSQZ       FSQL      fsqr      ',
     1   'fsqz      DELT    RAX(v=0)       WMHD      BETA',
     2   '      <M>   DEL-BSQ   FEDGE',/)
   16 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz      DELT      RAX(v=0)     WMHD       BETA',
     2   '     PHIEDGE  DEL-BSQ  FEDGE',/)
   20 format(/,' ITER    FSQR      FSQZ      FSQL    ',
     1   'RAX(v=0)      WMHD      DEL-BSQ',/)
   25 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz      DELT     RAX(v=0)      WMHD',
     2   '       BETA     <M>        ',/)
   30 format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     1   'fsqz     RAX(v=0)     WMHD',/)
      if (.not.lfreeb) then
         if (lscreen) print 45, i0, fsqr, fsqz, fsql, fsqr1, fsqz1,r00,w
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, avm
         return
      endif
      if (lscreen) print 50, i0, fsqr, fsqz, fsql, r00, w, delbsq
      if (imatch_phiedge .eq. 1) then
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, avm, delbsq, fedge
      else
         write (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, delt0, 
     1      r00, w, betav, abs(phiedge*phifac), delbsq, fedge
      endif
   40 format(i5,1p6e10.2,1pe11.3,1pe12.4,1pe11.3,1pe9.2,0pf7.3,1p2e9.2)
   45 format(i5,1p5e10.2,1pe11.3,1pe12.4)
   50 format(i5,1p3e10.2,1pe10.3,1pe11.3,1pe12.4,1pe11.3)
!     if (lrecon) then
!        if (ystark(1)*ystark(isnodes).lt.zero .or. 
!    1       ystark(2)*ystark(isnodes-1).lt.zero) then
!           if (lscreen) print 200
!           write (nthreed, 200)
!        endif
! 200    format(' CHECK THE SIGN OF MSE DATA: INCONSISTENT WITH ',
!    1      'CURRENT, TOROIDAL FIELD!')
!       if (pressum0 .ne. zero .and. apres .ne. zero .and. lscreen)
!    1     print 60, errsvd, pfac, phifac, 
!    2               abs(fsqsum0/pressum0), aminor/apres
!  60    format(2x,'chisq-err   = ',f8.3,3x,'pres-scale = ',f8.3,3x,
!    1      'phifac    = ',f8.3,/,2x,'<F00> force = ',f8.3,3x,
!    2      '<a>/apres  = ',f8.3)
!        if (imatch_phiedge.eq.2 .and. ivac>1) then
!           if (lscreen) print 65, dlim_min, rlim_min, zlim_min
!  65       format(2x,'dlim (min)  = ',f8.3,3x,'at rlim    = ',f8.3,
!    1         '   zlim = ',f8.3/)
!        else
!           if (lscreen) print *, ' '
!        endif
!     endif
      
      end subroutine printout


      subroutine spectrum(rmn, zmn)
      use vmec_main
      use vmec_params, only: mscale, nscale, ntmax
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(ns,0:ntor,0:mpol1,ntmax) :: 
     1   rmn, zmn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, ntype, n, m
      real(rprec), dimension(ns) :: t1, dnumer, denom
      real(rprec) :: scale
C-----------------------------------------------
 
      dnumer(2:ns) = zero
      denom(2:ns) = zero
      do ntype = 1,ntmax
        do n = 0,ntor
          do m = 1,mpol1
             scale = (mscale(m)*nscale(n))**2
             do js = 2,ns
                t1(js) =(rmn(js,n,m,ntype)**2 + zmn(js,n,m,ntype)**2)
     2               *scale
             end do
             dnumer(2:ns) = dnumer(2:ns) + t1(2:ns)*xmpq(m,3)
             denom (2:ns) = denom (2:ns) + t1(2:ns)*xmpq(m,2)
          end do
        end do
      enddo

      specw(2:ns) = dnumer(2:ns)/denom(2:ns)

      end subroutine spectrum


      subroutine readin(input_file, iseq_count, lfirst, ier_flag,
     1      lscreen)
      use vmec_main
      use vmec_params, only: ntmax, signgs
      use vacmod
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq_count, ier_flag
      logical lfirst, lscreen
      character*(*) :: input_file
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ns_default = 31
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: iexit, ipoint, n, iunit, ier_flag_init,
     1   i, n1, m, nsmin, igrid, m1, isgn
      real(rprec), dimension(:,:), pointer ::
     1  rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      real(rprec) rtest, ztest, treadon, treadoff, tzc, trc, delta
      character*(80) :: line
C-----------------------------------------------
!
!       LOCAL VARIABLES
!
!       rbcc,rbss,rbcs,rbsc
!                boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
!       zbcc,zbss,zbcs,zbsc
!                boundary Fourier coefficient arrays for Z

!
!       STACKING ORDER (DEPENDS ON LASYM):
!       XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!               LASYM = T                             LASYM = F
!          1      ,  mns => rmncc                1      ,  mns => rmncc
!          1+  mns,2*mns => rmnss                1+  mns,2*mns => rmnss
!          1+2*mns,3*mns => rmncs
!          1+3*mns,4*mns => rmnsc
!          1+4*mns,5*mns => zmncs                1+2*mns,3*mns => zmncs
!          1+5*mns,6*mns => zmnsc                1+3*mns,4*mns => zmnsc
!          1+6*mns,7*mns => zmncc
!          1+7*mns,8*mns => zmnss
!          1+8*mns,9*mns => lmncs                1+4*mns,5*mns => lmncs
!          1+9*mns,10*mns=> lmnsc                1+5*mns,neqs  => lmnsc
!          1+10*mns,11*mns=> lmncc
!          1+11*mns,neqs  => lmnss
!          neqs1          => phifac
!          neqs2          => MSE axis correction
!


!
!                STANDARD INPUT DATA AND RECOMMENDED VALUES
!
!   Plasma parameters (MKS units)
!          ai:   iota (ncurr=0) or toroidal current density (=ac, ncurr=1)
!                expansion coefficients (series in s)
!          am:   mass or pressure (gamma=0) expansion coefficients (series in s)
!                in MKS units (NWT/M**2)
!      curtor:   value of toroidal current (A). Used if ncurr = 1 to specify
!                current profile, or if in data reconstruction mode.
!      extcur:   array of currents in each external current group. Used to
!                multiply Green''s function for fields and loops read in from
!                MGRID file. Should use real current units (A).
!       gamma:   value of compressibility index (gamma=0 => pressure prescribed)
!         nfp:   number of toroidal field periods ( =1 for Tokamak)
!         rbc:   boundary coefficients of cos(m*theta-n*zeta) for R
!         zbs:   boundary coefficients of sin(m*theta-n*zeta) for Z
!         rbs:   boundary coefficients of sin(m*theta-n*zeta) for R
!         zbc:   boundary coefficients of cos(m*theta-n*zeta) for Z
!
!
!   Numerical parameters
!  mgrid_file:   full path for vacuum green''s function data
!       ncurr:   flux conserving (=0) or prescribed toroidal current (=1)
!    ns_array:   array of radial mesh sizes to be used in multigrid sequence
!
!   Convergence control parameters
!  ftol_array:   array of value of residual(s) at which each multigrid
!                iteration ends
!       niter:   number of iterations (used to terminate run)
!       nstep:   number of timesteps between printouts on screen
!    nvacskip:   iterations skipped between full update of vacuum solution
!       tcon0:   weight factor for constraint force (=1 by default)
!
!   Equilibrium reconstruction parameters
!      phifac:   factor scaling toroidal flux to match apres or limiter
!   datastark:   pitch angle data from stark measurement
!    datathom:   pressure data from Thompson, CHEERS (Pa)
!     imatch_         = 1 (default),match value of PHIEDGE in input file
!     phiedge:   = 0, use pressure profile width to determine PHIEDGE
!                = 2, use LIMPOS data (in mgrid file) to find PHIEDGE
!                = 3, use Ip to find PHIEDGE (fixed-boundary only)
!        imse:   number of Motional Stark effect data points
!                >0, use mse data to find iota; <=0, fixed iota profile ai
!        itse:   number of pressure profile data points
!                = 0, no thompson scattering data to read
!     isnodes:   number of iota spline points (computed internally unless specified explicitly)
!     ipnodes:   number of pressure spline points (computed internally unless specified explicitly)
!       lpofr:   logical variable. =.true. if pressure data are
!                prescribed in real space. =.false. if data in flux space.
!      pknots:   array of pressure knot values in sqrt(s) space
!      sknots:   array of iota knot values in sqrt(s) space
!       tensp:   spline tension for pressure profile
!
!       tensi:   spline tension for iota
!      tensi2:   vbl spline tension for iota
!      fpolyi:   vbl spline tension form factor (note: if tensi!=tensi2
!               then tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
!               - - - - - - - - - - - - - - - - - -
!    mseangle_   uniform experimental offset of MSE data
!     offset:    (calibration offset) ... PLUS ...
!    mseangle_   multiplier on mseprof offset array
!     offsetM:   (calibration offset)
!     mseprof:   offset array from namelist MSEPROFIL
!                so that the total offset on the i-th MSE data point is
!                taken to be
!                = mseangle_offset+mseangle_offsetM*mseprof(i)
!               - - - - - - - - - - - - - - - - - -
! pres_offset:   uniform arbitrary  radial offset of pressure data
!     presfac:   number by which Thomson scattering data is scaled
!                to get actual pressure
!     phidiam:   diamagnetic toroidal flux (Wb)
!      dsiobt:   measured flux loop signals corresponding to the
!                combination of signals in iconnect array
!     indxflx:   array giving index of flux measurement in iconnect array
!    indxbfld:   array giving index of bfield measurement used in matching
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!    bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
!                the orientation br*cos(abcoil) + bz*sin(abcoil)
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!                of the m-th coil in the n-th set from mgrid file.
!       nflxs:   number of flux loop measurements used in matching
!    nbfld(n):   number of selected external bfield measurements in set n from nml file
!      nbfldn:   total number of external bfield measurements used in matching
!               - - - - - - - - - - - - - - - - - -
!             NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
!             AS PERCENT OF RESPECTIVE MEASUREMENT
!  sigma_thom:   standard deviation (Pa) for pressure profile data
! sigma_stark:   standard deviation (degrees) in MSE data
!  sigma_flux:   standard deviaton (Wb) for external poloidal flux data
!     sigma_b:   standard deviation (T) for external magnetic field data
!sigma_current:  standard deviation (A) in toroidal current
!sigma_delphid:  standard deviation (Wb) for diamagnetic match
!
!                - - - - - - - - - - - - - - - - - -
!               INPUTS IN MAKEGRID CODE WHICH ARE EXPORTED
!               IN MGRID FILE
!                - - - - - - - - - - - - - - - - - -
!               FOR SPECIFYING BFIELD GREEN''S FUNCTIONS
!        rmin:   minimum R position for green''s function data
!        rmax:   maximum R position for green''s function data
!        zmin:   minimum Z position for green''s function data
!        zmax:   maximum Z position for green''s function data
!          ir:   parameter specifying number of points in R mesh
!                for green''s function box
!          jz:   parameter specifying number of points in Z mesh
!                for green''s function box
!          kp:   parameter specifying number of toroidal planes
!                within one field period used for green''s
!                function box. MUST agree with NV parameter in
!                VMEC for now
!     nextcur:         no. of external current groups (eg., TF, PF, helical)
!    curlabel:   array of labels describing each current group
!                     included in green''s function BFIELD response
!
!               - - - - - - - - - - - - - - - - - -
!               FOR DIAGNOSTICS AND DATA ANALYSIS
!               (HERE,COILS ARE FOR MEASURING FIELDS, FLUXES)
!    iconnect:   two-dimensional array describing electrical
!                      connection of up to four flux loops. Specifies
!                     the sign and flux loop number of (up to) four
!                     connected individual loops (indexing based on
!                     xobser,zobser arrays).
!      nobser:   no. of distinct flux loops
!        nobd:   no. of distinct connections of flux loops
!      nbsets:   no. of b-field coil sets
!     nbcoils:   array specifying no. of b-field coils in each set
!  bloopnames:   array of labels describing b-field sets
!    dsilabel:   array of labels describing connected flux loops
!      xobser:   array of flux loop R-positions
!      zobser:   array of flux loop Z-positions
!      rbcoil:   two-dimensional array of b-field coil R-positions
!                (no. of coil groups, no. coils in a group)
!      zbcoil:   two-dimensional array of b-field coil Z-positions
!      abcoil:   two-dimensional array of angles (in degrees) of
!                       each b-field coil
!
!
!
!       THE (ABSOLUTE) CHI-SQ ERROR IS DEFINED AS FOLLOWS:
!
!          2
!       CHI      =     SUM [ EQ(K,IOTA,PRESSURE)  -  DATA(K) ] ** 2
!                     (K) -----------------------------------
!                                   SIGMA(K)**2
!
!       HERE, SIGMA IS THE STANDARD DEVIATION OF THE MEASURED DATA, AND
!       EQ(IOTA,PRESSURE) IS THE EQUILIBRIUM EXPRESSION FOR THE DATA TO BE
!       MATCHED:
!
!       EQ(I)   =    SUM [ W(I,J)*X(J) ]
!                   (J)
!
!       WHERE W(I,J) ARE THE (LINEAR) MATRIX ELEMENTS AND X(J) REPRESENT
!       THE KNOT VALUES OF IOTA (AND/OR PRESSURE). THE RESULTING LEAST-SQUARES
!       MATRIX ELEMENTS AND DATA ARRAY CAN BE EXPRESSED AS FOLLOWS:
!
!       ALSQ(I,J) = SUM [ W(K,I) * W(K,J) / SIGMA(K) ** 2]
!                   (K)
!
!       BLSQ(I)   = SUM [ W(K,I) * DATA(K)/ SIGMA(K) ** 2]
!                   (K)
!
!       THEREFORE, INTERNALLY IT IS CONVENIENT TO WORK WITH THE 'SCALED'
!       W'(K,I) = W(K,I)/SIGMA(K) AND DATA'(K) = DATA(K)/SIGMA(K)
!
!       ****!   I - M - P - O - R - T - A - N - T     N - O - T - E   *****
!
!       THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
!       SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
!       TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
!       AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
!       THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
!       SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.
!
      ier_flag_init = ier_flag
      ier_flag = 0
      call second0(treadon)
      if (ier_flag_init .eq. 4) goto 1000

!
!     READ IN DATA FROM INDATA FILE
!
      call read_indata(input_file, iunit, ier_flag)
      if (ier_flag .ne. 0) return


      if (tensi2 .eq. zero ) tensi2 = tensi

!
!     open output files here, print out heading to threed1 file
!
      call heading(input_extension, time_slice,
     1      iseq_count, lmac, lscreen, lfirst)

!
!     READ IN COMMENTS DEMARKED BY "!"
!
      rewind (iunit)
      iexit = 0
      do while( iexit.eq.0 )
         read (iunit, '(a)') line
         iexit = index(line,'INDATA') + index(line,'indata')
         ipoint = index(line,'!')
         if (ipoint .eq. 1) write (nthreed,*) line
      enddo
      close (iunit)

!
!     READ IN AND STORE (FOR SEQUENTIAL RUNNING) MAGNETIC FIELD DATA
!     FROM MGRID_FILE FIRST TIME (lfirst = T) ONLY
!     SET LOGICAL FLAGS FOR ALL SUBSEQUENT RUNS
!
      if (lfirst .and. lfreeb) call read_mgrid (lscreen, ier_flag)
      if (ier_flag .ne. 0) return

!
!     PARSE NS_ARRAY
!
      nsin = max (3, nsin)
      multi_ns_grid = 1
      if (ns_array(1) .eq. 0) then                   !!Old input style
          ns_array(1) = min(nsin,nsd)
          multi_ns_grid = 2
          ns_array(multi_ns_grid) = ns_default        !!Run on 31-point mesh
      else
          nsmin = 1
          do while (ns_array(multi_ns_grid) .gt. nsmin .and.
     1             multi_ns_grid .lt. 100)
             nsmin = max(nsmin, ns_array(multi_ns_grid))
            if (nsmin.le.nsd) then
               multi_ns_grid = multi_ns_grid + 1
            else                                      !!Optimizer, Boozer code overflows otherwise
               ns_array(multi_ns_grid) = nsd
               nsmin = nsd
               print *,' NS_ARRAY ELEMENTS CANNOT EXCEED ',nsd
               print *,' CHANGING NS_ARRAY(',multi_ns_grid,') to ', nsd
            end if
          end do
          multi_ns_grid = multi_ns_grid - 1
      endif
      if (ftol_array(1) .eq. zero) then
         ftol_array(1) = 1.e-8_dp
         if (multi_ns_grid .eq. 1) ftol_array(1) = ftol
         do igrid = 2, multi_ns_grid
            ftol_array(igrid) = 1.e-8_dp * (1.e8_dp * ftol)**
     1        ( real(igrid-1,rprec)/(multi_ns_grid-1) )
         end do
      endif

!
!     WRITE OUT DATA TO THREED1 FILE
!
      write (nthreed,100)
     1  ns_array(multi_ns_grid),ntheta1,nzeta,mpol,ntor,nfp,
     2  gamma,spres_ped,phiedge,curtor
 100  format(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'/,
     1  1x,45('-'),/,
     2  '     ns     nu     nv     mu     mv',/,
     3  5i7,//,' CONFIGURATION PARAMETERS:',/,1x,25('-'),/,
     4  '    nfp      gamma      spres_ped    phiedge(wb)     curtor(A)'
     5  ,/,i7,1pe11.3,1p2e15.3,1pe14.3,/)

      if (nvacskip.le.0) nvacskip = nfp
      write (nthreed,110) ncurr,niter,ns_array(1),nstep,nvacskip,
     1  ftol_array(multi_ns_grid),tcon0,lasym
 110    format(' RUN CONTROL PARAMETERS:',/,1x,23('-'),/,
     1  '  ncurr  niter   nsin  nstep  nvacskip      ftol     tcon0',
     2  '   lasym',/, 4i7,i10,1p2e10.2,L8/)
      if (nextcur.ne.0) then
        write(nthreed,120)
        do i = 1,nextcur,5
          write (nthreed, 125) trim(curlabel(i)), trim(curlabel(i+1)),
     1    trim(curlabel(i+2)), trim(curlabel(i+3)), trim(curlabel(i+4)),
     2    extcur(i), extcur(i+1), extcur(i+2), extcur(i+3), extcur(i+4)
        enddo
      endif
 120  format(' EXTERNAL CURRENTS',/,1x,17('-'))
 125  format(5a12,/,1p5e12.4,/)

      if (bloat .ne. one) then
          write (nthreed,'(" Profile Bloat Factor: ",f7.5)') bloat
          phiedge = phiedge*bloat
      endif

      write(nthreed,130)
 130  format(' MASS PROFILE COEFFICIENTS am - newton/m**2',
     1  ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
      write(nthreed,135)(am(i-1),i=1, size(am))
      if (ncurr.eq.0) then
          write(nthreed,140)
          write(nthreed,135)(ai(i-1),i=1, size(ai))
      else
          select case(ac_form)
          case (1)
             write(nthreed,1145)
          case default
             write(nthreed,145)
          end select
             
          write(nthreed,135)(ac(i-1),i=1, size(ac))
      endif
!     write(nthreed,150)
!     write(nthreed,135)(aphi(i-1),i=1, size(aphi))    !!Temporarily disabled

 135  format(1p6e12.3)
 140  format(/' IOTA PROFILE COEFFICIENTS ai',
     1   ' (EXPANSION IN TOROIDAL FLUX):',/,1x,35('-'))
 145  format(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1        ' ac (EXPANSION IN TOROIDAL FLUX):',/,1x,38('-'))
 1145 format(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS',
     1        ' ac (EXPANSION IN TOROIDAL FLUX w/ SQRT Term):',
     2          /,1x,38('-'))
 150  format(/' NORMALIZED TOROIDAL FLUX COEFFICIENTS aphi',
     1   ' (EXPANSION IN S):',/,1x,35('-'))
      write(nthreed,180)
 180  format(/,' R-Z FOURIER BOUNDARY COEFFICIENTS',/,
     1  ' R = RBC*cos(m*u - n*v) + RBS*sin(m*u - n*v),',
     2  ' Z = ZBC*cos(m*u - n*v) + ZBS*sin(m*u-n*v)'/1x,86('-'),
     3  /,'   nb  mb     rbc         rbs         zbc         zbs   ',
     4   '     raxis(cc)   raxis(cs)   zaxis(cc)   zaxis(cs)')

 1000  continue

      if (.not.lasym) then
!
!       CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
!

      delta = atan( (rbs(0,1) - zbc(0,1))/(rbc(0,1) + zbs(0,1)) )
      if (delta .ne. zero) then
        do m = 0,mpol1
          do n = -ntor,ntor
            trc = rbc(n,m)*cos(m*delta) + rbs(n,m)*sin(m*delta)
            rbs(n,m) = rbs(n,m)*cos(m*delta) - rbc(n,m)*sin(m*delta)
            rbc(n,m) = trc
            tzc = zbc(n,m)*cos(m*delta) + zbs(n,m)*sin(m*delta)
            zbs(n,m) = zbs(n,m)*cos(m*delta) - zbc(n,m)*sin(m*delta)
            zbc(n,m) = tzc
          enddo
        enddo
      endif

      endif

!
!     ALLOCATE MEMORY FOR NU, NV, MPOL, NTOR SIZED ARRAYS
!
      call allocate_nunv

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!         + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!         + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!

!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)

      rbcc => rmn_bdy(:,:,1)
      rbss => rmn_bdy(:,:,2)
      zbcs => zmn_bdy(:,:,1)
      zbsc => zmn_bdy(:,:,2)

      if (lasym) then
      rbcs => rmn_bdy(:,:,ntmax-1)
      rbsc => rmn_bdy(:,:,ntmax)
      zbcc => zmn_bdy(:,:,ntmax-1)
      zbss => zmn_bdy(:,:,ntmax)
      endif

      rmn_bdy(0:ntor,0:mpol1,:ntmax) = zero
      zmn_bdy(0:ntor,0:mpol1,:ntmax) = zero

      do 190 m=0,mpol1
         m1 = m+1
         do 190 n=-ntor,ntor
            n1 = abs(n) + 1
            if (n .eq. 0) then
               isgn = 0
            else if (n .gt. 0) then
               isgn = 1
            else
               isgn = -1
            end if
            rbcc(n1,m1) = rbcc(n1,m1) + rbc(n,m)
            rbss(n1,m1) = rbss(n1,m1) + isgn*rbc(n,m)
            zbcs(n1,m1) = zbcs(n1,m1) - isgn*zbs(n,m)
            zbsc(n1,m1) = zbsc(n1,m1) + zbs(n,m)

            if (lasym) then
            rbcs(n1,m1) = rbcs(n1,m1) - isgn*rbs(n,m)
            rbsc(n1,m1) = rbsc(n1,m1) + rbs(n,m)
            zbcc(n1,m1) = zbcc(n1,m1) + zbc(n,m)
            zbss(n1,m1) = zbss(n1,m1) + isgn*zbc(n,m)
            if (m .eq. 0) zbss(n1,m1) = zero
            if (m .eq. 0) rbsc(n1,m1) = zero
            end if

            if (ier_flag_init .ne. 0) cycle
            trc = abs(rbc(n,m)) + abs(rbs(n,m))
     1          + abs(zbc(n,m)) + abs(zbs(n,m))
            if (m .eq. 0) then
               rbss(n1,m1) = zero
               zbsc(n1,m1) = zero
               if (n .lt. 0) cycle
               if (trc.eq.zero .and. abs(raxis(n,1)).eq.zero .and.
     1             abs(zaxis(n,1)).eq.zero) cycle
               write (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m), raxis(n,1), raxis(n,2),
     2                   zaxis(n,2), zaxis(n,1)
            else
               if (trc .eq. zero) cycle
               write (nthreed,195) n, m, rbc(n,m), rbs(n,m),
     1                   zbc(n,m), zbs(n,m)
            end if
 190  continue
 195  format(i5,i4,1p8e12.4)


!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!
      m = 1
      m1 = m+1
      rtest = sum(rbcc(1:ntor1,m1))
      ztest = sum(zbsc(1:ntor1,m1))
      signgs = one
      if (rtest*ztest .gt. zero) signgs = -one


      iresidue = -1
      if (lrecon) then
!
!       DETERMINE CURRENT-FLUX CONSISTENCY CHECK
!
        signiota = one
        if (signgs*curtor*phiedge .lt. zero)signiota = -one
        if (sigma_current .eq. zero) then
          write (*,*) 'Sigma_current cannot be zero!'
          ier_flag = 1
          return
        end if

!
!       SET UP RECONSTRUCTION FIXED PROFILES
!
        dcon = atan(one)/45
        call readrecon                   !Setup for reconstruction mode
        call fixrecon(ier_flag)          !Fixed arrays for reconstruction
        if (ier_flag .ne. 0) return
      end if

      currv = dmu0*curtor              !Convert to Internal units

      call second0(treadoff)
      timer(2) = timer(2) + (treadoff-treadon)

      end subroutine readin


      subroutine read_indata(in_file, iunit, ier_flag)
      use vmec_main
      use vmec_input, only: bloat, ncurr
      use vmec_params, only: ntmax
      use vacmod
      use read_namelist_mod
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ier_flag, iunit
      character*(*) :: in_file
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ireadseq, lt, iosnml
      real(rprec) :: denom
C-----------------------------------------------
!     INDATA0 MUST BE CLOSED EXTERNAL TO THIS RETURN

      iunit = indata0
      call safe_open (iunit, ireadseq, in_file, 'old', 'formatted')
      if (ireadseq .ne. 0) then
         print *, ' In VMEC, error opening input file: ',
     1   trim(in_file), '. Iostat = ', ireadseq
         ier_flag = 5
         return
      endif


      call read_namelist (iunit, iosnml, 'indata')
      if (iosnml .ne. 0) then
         print *, ' In VMEC, indata namelist error: iostat = ', iosnml
         ier_flag = 5
         return
      endif

!
!     This is temporarily disabled in torflux, torflux_deriv routines...
!
      aphi(0) = zero
      denom = sum (aphi)
      if (denom .ne. zero) then
         aphi = aphi/denom
      else
         aphi = ( / zero, one, (zero, lt=2,10) / )
         denom = sum (aphi)
         aphi = aphi/denom
      end if   
 
      call read_namelist (iunit, iosnml, 'mseprofile')

      if (lrecon .and. itse.le.0 .and. imse.le.0) lrecon = .false.
      if (lfreeb .and. mgrid_file.eq.'NONE') lfreeb = .false.

      if (bloat .eq. zero) bloat = one
      if ((bloat.ne.one) .and. (ncurr.ne.1)) then
         ier_flag = 3
         return 
      endif
!
!     COMPUTE NTHETA, NZETA VALUES
!
      mpol = abs(mpol)
      ntor = abs(ntor)
      if (mpol .gt. mpold) stop 'mpol>mpold: lower mpol'
      if (ntor .gt. ntord) stop 'ntor>ntord: lower ntor'
      mpol1 = mpol - 1
      ntor1 = ntor + 1
      if (ntheta .le. 0) ntheta = 2*mpol + 6    !number of theta grid points (>=2*mpol+6)
      ntheta1 = 2*(ntheta/2)
      ntheta2 = 1 + ntheta1/2                   !u = pi
      if (ntor .eq. 0) lthreed = .false.
      if (ntor .gt. 0) lthreed = .true.

      if (nzeta .le. 0) nzeta = 2*ntor + 4      !number of zeta grid points (=1 if ntor=0)
      if (ntor .eq. 0) nzeta = 1
      mnmax = ntor1 + mpol1*(1 + 2*ntor)        !size of rmnc,  rmns,  ...
      mnsize = mpol*ntor1                       !size of rmncc, rmnss, ...

      mf = mpol + 1
      nf = ntor
      nu = ntheta1
      nv = nzeta
      mf1 = 1 + mf
      nf1 = 2*nf + 1
      mnpd = mf1*nf1
      nfper = nfp
!     NEED FOR INTEGRATION IN BELICU FOR AXISYMMETRIC PLASMA
      if (nf.eq.0 .and. nv.eq.1) nfper = 64

      if (lasym) then
         ntmax = 4
         ntheta3 = ntheta1
         mnpd2 = 2*mnpd
      else
         ntmax = 2
         ntheta3 = ntheta2
         mnpd2 = mnpd
      end if
      
      nvp = nv*nfper
      nuv = nu*nv
      nu2 = nu/2 + 1
      nu3 = ntheta3
      nznt = nzeta*ntheta3
      nuv2 = nznt
!     if (nuv2 < mnpd) then
!        print *, ' nuv2 < mnpd: not enough integration points'
!        stop 11
!     endif

      if (ncurr.eq.1 .and. all(ac.eq.cbig)) ac = ai            !!Old format: may not be reading in ac
      where (ac .eq. cbig) ac = zero
      
      end subroutine read_indata

 
      subroutine read_mgrid(lscreen, ier_flag)
      use vmec_main
      use vacmod
      use vsvd
      use vspline
      use mgrid_mod
      use safe_open_mod
      use system_mod
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      integer :: ier_flag
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
#ifdef VMS
      character*(*), parameter :: mgrid_defarea = 'vmec$:[makegrid]'
#else
      character*(*), parameter :: mgrid_defarea = '$HOME/vmec/MAKEGRID'
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ipath, ierr, ierr1, ierr2, istat0, istat, ii, 
     1   i, j, n, n1, m, nsets_max, k, iunit = nbfld0
      character*200 :: home_dir
      logical :: lgrid_exist
C-----------------------------------------------
      mgrid_path = trim(mgrid_file)
      inquire (file=mgrid_path,exist=lgrid_exist)
      if (.not.lgrid_exist) then
          if (lscreen) print *,' MGRID FILE NOT FOUND IN SPECIFIED ',
     1       'PATH: TRYING TO FIND MGRID-FILE IN DEFAULT AREA'
          ii = index(mgrid_file,'/',back=.true.)
          istat = index(mgrid_defarea, '$HOME')
          if (istat .ne. 0) then
             call getenv('HOME', home_dir)
             if (istat .gt. 1) then
                home_dir = mgrid_defarea(1:istat-1) // trim(home_dir)
     1                 // trim(mgrid_defarea(istat+5:))
             else
                home_dir = trim(home_dir) // 
     1                     trim(mgrid_defarea(istat+5:))
             end if
          else
             home_dir = mgrid_defarea
          end if      
          mgrid_path = trim(home_dir) // mgrid_file(ii+1:)
          inquire (file=mgrid_path,exist=lgrid_exist)
      end if

      if (lgrid_exist) then
          call safe_open(iunit, ipath, mgrid_path, 'old', 'unformatted')
          if (ipath .ne. 0)then
             if (lscreen) then
                print *,' ',trim(mgrid_file(ii+1:)),' NOT FOUND IN ', 
     1                      trim(home_dir) 
                print *, ' MUST SUPPLY VACUUM BFIELD ON GRID TO ',
     1                   'RUN FREE-BOUNDARY'
                print *,' WILL PROCEED TO RUN VMEC IN',
     1                  ' FIXED BOUNDARY MODE'
             end if
             lfreeb = .false.
          endif
      else
         lfreeb = .false.
      end if

      if (.not.lfreeb) then
        lrecon = .false.
        return
      end if

! 
!     PARTITION B-LOOPS INTO SETS
!
      if (lscreen)
     1    print '(2x,2a)','Opening vacuum field file: ',mgrid_file

      read(iunit,iostat=ierr) nr0b, nz0b, np0b, nfper0, nextcur
      if (ierr .ne. 0) ier_flag = 9
      read(iunit,iostat=ierr) rminb, zminb, rmaxb, zmaxb
      if (ierr .ne. 0) ier_flag = 9

      if (np0b .ne. nv) then
        print *,' NV=',nv,' NOT EQUAL TO NP0B=',np0b,' IN MGRID FILE'
        ier_flag = 9
      else if ((nf.ne.0) .and. (nfper0.ne.nfp)) then
        print *,' NFP(read in) = ',nfp,' DOES NOT AGREE WITH ',
     1    'NFPER (in vacuum field file) = ',nfper0
        ier_flag = 9
      else if (nextcur .le. 0) then
        print *,' NEXTCUR = ',nextcur,' IN READING MGRID FILE'
        ier_flag = 9
      end if

      if (ier_flag .eq. 10) return
      
      allocate (curlabel(5*(nextcur/5+1)), stat=istat0)    !min of 5 for printing
      curlabel = " "
      read(iunit,iostat=istat) (curlabel(n),n=1,nextcur)
      if (istat.ne.0 .or. istat0.ne.0) then
         print *,' reading mgrid-files (curlabel) failed'
         ier_flag = 9
         return
      end if

!
!       NOTE: IF THE BTEMP ARRAY GETS TOO LARGE, USER CAN
!           ADD UP BVAC DIRECTLY FOR EACH EXTERNAL CURRENT GROUP
!           IN LOOP 50 WITHOUT STORING BTEMP(...,II)

      nbvac = nr0b*nz0b*nv
      allocate (btemp(nbvac,3,nextcur), stat=istat)
      if (istat .ne. 0) then
        print *,' allocation for b-vector storage failed'
        ier_flag = 9
        return
      end if

      do ii = 1,nextcur
         read(iunit) (btemp(i,1,ii),btemp(i,2,ii),btemp(i,3,ii),
     1              i=1,nbvac)
      enddo
      
!
!       READ IN EXTERNAL POLOIDAL FLUX, FIELD MEASURMENT
!       LOOP COORDINATES AND LABELS
!
      read(iunit,iostat=ierr1) nobser, nobd, nbsets
      if (ierr1.ne.0) then
         nobser = 0
         nobd   = 0
         nbsets = 0
         if (lscreen) print *,' No observation data in mgrid data'
         go to 900
      end if

      nbfldn = sum(nbfld(:nbsets))
      allocate (nbcoils(nbsets), stat=istat0)
      read(iunit) (nbcoils(n),n=1,nbsets)

      nbcoil_max = maxval(nbcoils(:nbsets))

      allocate (xobser(nobser), zobser(nobser), dsilabel(nobd),
     1       iconnect(4,nobser+nobd), unpsiext(nobser,nextcur),
     2       xobsqr(nobser), needflx(nobser), plflux(nobser+nobd),
     3       dsiext(nobd), psiext(nobser), bloopnames(nbsets),
     4       needbfld(nbcoil_max,nbsets), plbfld(nbcoil_max,nbsets),
     5       rbcoil(nbcoil_max,nbsets), zbcoil(nbcoil_max,nbsets),
     6       abcoil(nbcoil_max,nbsets), bcoil(nbcoil_max,nbsets),
     7       rbcoilsqr(nbcoil_max,nbsets), b_chi(nbsets),
     8       dbcoil(nbcoil_max,nbsets,nextcur), stat = istat)
      if (istat .ne. 0) then
          if (lscreen) 
     1       print *,' allocation error for xobser: istat = ',istat
          ier_flag = 9
          return
      end if
     
      if (nobser .gt. nfloops) then
         print *, 'NOBSER>NFLOOPS'
         ier_flag = 9
      end if   
      if (nobd .gt. nfloops) then
         print *, 'NOBD>NFLOOPS'
         ier_flag = 9
      end if
      if (nflxs .gt. nfloops) then
         print *, 'NFLXS>NFLOOPS'
         ier_flag = 9
      end if   
      if (nbfldn .gt. nbctotp) then
         print *, 'NBFLDN>NBCTOTP'
         ier_flag = 9
      end if   
      if (nbcoil_max .gt. nbcoilsp) then
         print *, 'NBCOIL_MAX>NBCOILSP'
         ier_flag = 9
      end if   

      if (ier_flag .eq. 10) return

      if (nobser+nobd .gt. 0) iconnect(:4,:nobser+nobd) = 0

      read(iunit) (xobser(n), zobser(n),n=1,nobser)
      read(iunit) (dsilabel(n),n=1,nobd)
      read(iunit) ((iconnect(j,n),j=1,4),n=1,nobd)

      if (nbcoil_max.gt.0 .and. nbsets.gt.0) then
         rbcoil(:nbcoil_max,:nbsets) = 0
         zbcoil(:nbcoil_max,:nbsets) = 0
         abcoil(:nbcoil_max,:nbsets) = 0

         do n=1,nbsets
           if (nbcoils(n).gt.0) then
           read(iunit) n1,bloopnames(n1)
           read(iunit)(rbcoil(m,n),zbcoil(m,n),abcoil(m,n),
     1             m=1,nbcoils(n))
           endif
         enddo

         dbcoil(:nbcoil_max,:nbsets,:nextcur) = 0
      end if
      do ii = 1,nextcur
        !un-connected coil fluxes
         read(iunit) (unpsiext(n,ii),n=1,nobser)
         do n = 1,nbsets
            read(iunit) (dbcoil(m,n,ii),m=1,nbcoils(n))
         enddo
      enddo

!
!     READ LIMITER & PROUT PLOTTING SPECS
!
      read(iunit,iostat=ierr2) nlim,(limitr(i),i=1,nlim)
      if (ierr2 .ne. 0)then
        nlim = 0
        if (lscreen) print *,' No limiter data in mgrid file'
        go to 900
      end if

      nlim_max = maxval(limitr)

      if (nlim .gt. nlimset) then
         print *, 'nlim>nlimset'
         ier_flag = 9
         return
      end if   
 
      allocate( rlim(nlim_max,nlim),   zlim(nlim_max,nlim),
     1          reslim(nlim_max,nlim) ,seplim(nlim_max,nlim),
     2          stat=istat)
      if (istat .ne. 0) then
         print *, 'rlim istat!=0'
         ier_flag = 9
         return
      end if   

      read(iunit, iostat=ierr2)
     1   ((rlim(i,j),zlim(i,j),i=1,limitr(j)),j=1,nlim)
      read(iunit, iostat=ierr2) nsets,(nsetsn(i), i=1,nsets)
      
      if (nsets .gt. nigroup) then
         print *, 'nsets>nigroup'
         ier_flag = 9
         return
      else if (ierr2 .ne. 0) then
         ier_flag = 9
         return
      end if

      nsets_max = maxval(nsetsn)
 
      if (nsets_max .gt. npfcoil) then
         print *, 'nsetsn>npfcoil'
         ier_flag = 9
         return
      end if   

      allocate (pfcspec(nparts,nsets_max,nsets), stat=istat)

!     NOTE TO RMW: SHOULD READ IN NPARTS HERE (PUT INTO MGRID FILE)

      read(iunit, iostat=ierr2) (((pfcspec(i,j,k),i=1,nparts),
     1        j=1,nsetsn(k)), k=1,nsets)
      read(iunit, iostat=ierr2) rx1,rx2,zy1,zy2,condif,
     1  nrgrid,nzgrid,tokid

      if (ierr2.ne.0 .or. istat.ne.0) then
         ier_flag = 9
         return
      end if
      
      if (nobser .gt. 0) xobsqr(:nobser) = sqrt(xobser(:nobser))
! 
!       PARTITION MGRID B-LOOPS INTO SETS
!
      nbcoilsn = sum(nbcoils(:nbsets))

      do n = 1,nbsets
        rbcoilsqr(:nbcoils(n),n) = sqrt(rbcoil(:nbcoils(n),n))
      enddo

 900  continue

      close (iunit)

      delrb = (rmaxb-rminb)/(nr0b-1)
      delzb = (zmaxb-zminb)/(nz0b-1)

!
!     SUM UP CONTRIBUTIONS FROM INDIVIDUAL COIL GROUPS
!
      if (lfreeb) then
         write (nthreed,20) nr0b, nz0b, np0b, rminb, rmaxb,
     1        zminb, zmaxb, trim(mgrid_file)
 20      format(//,' VACUUM FIELD PARAMETERS:',/,1x,24('-'),/,
     1  '  nr-grid  nz-grid  np-grid      rmin      rmax      zmin',
     2  '      zmax     input-file',/,3i9,4f10.3,5x,a)

         istat = 0
         if (.not.allocated(bvac)) allocate (bvac(nbvac,3), stat=istat)
         if (istat.ne.0) then
           print *,' bvac allocation failed'
           ier_flag = 9
           return
         end if

         bvac(:nbvac,:3) = 0
         if (nobser .gt. 0) psiext(:nobser) = 0
         if (nbcoil_max.gt.0 .and. nbsets.gt.0)
     1       bcoil(:nbcoil_max, :nbsets) = 0
 
         do ii = 1,nextcur
            do i = 1,3
               bvac(:nbvac, i) = bvac(:nbvac, i) + 
     1                           extcur(ii)*btemp(:nbvac,i,ii)
            enddo
            if (nobser .gt. 0)
     1      psiext(:nobser) = psiext(:nobser) +
     2                        extcur(ii)*unpsiext(:nobser,ii)
            do n=1,nbsets
               n1 = nbcoils(n)
               bcoil(:n1,n) = bcoil(:n1,n) + 
     1                        extcur(ii)*dbcoil(:n1,n,ii)
            enddo
         enddo
      endif                   !!IF LFREEB
      
      end subroutine read_mgrid


      subroutine restart(time_step)
      use vmec_main
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: time_step
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p03 = 1.03_dp, cp90 = 0.90_dp
c-----------------------------------------------
 
      select case (irst) 
      case default
         xstore(:neqs2) = xc(:neqs2)
         return 
      case (2:3) 
         xcdot(:neqs2) = zero
         xc(:neqs2) = xstore(:neqs2)
         time_step = time_step*((irst-2)/c1p03 + cp90*(3-irst))
         if (irst .eq. 2) ijacob = ijacob + 1
         irst = 1
         return 
      end select

      end subroutine restart
      

      subroutine runvmec(input_file, iseq_count, lreseta, ier_flag, 
     1   lfirst_call, lscreen, reset_file_name)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer iseq_count, ier_flag
      logical lreseta, lfirst_call, lscreen
      character*(*) :: input_file, reset_file_name
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: bad_jacobian_flag = 6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: igrid0, igrid, nsval, index_end, index_dat, 
     1   ns_min, ier_init, jacob_off
      real(rprec) :: timeon, timeoff
      logical :: interp_flag, lreset
C-----------------------------------------------
 
!
!
!     INDEX OF LOCAL VARIABLES
!
!     ier_flag   specifies error condition if nonzero
!     interp_flag 
!                = FALSE,compute xc array from scaled boundary data
!                = TRUE, compute xc array by interpolation
!     lfirst_call
!                = T, initial call to runvmec
!                = F, runvmec has been previously called
!
      call second0 (timeon)

!
!     PARSE input_file into path/input.ext
!
      index_dat = index(input_file,'input.')
      index_end = len_trim(input_file)
      if (index_dat .gt. 0) then
         input_extension  = input_file(index_dat+6:index_end)
      else
         input_extension = input_file(1:index_end)
         input_file = 'input.'//trim(input_extension)
      end if   

!
!     INITIALIZE PARAMETERS
!
      lreset = .false.
      ier_init = ier_flag
      if (lfirst_call .or. lreseta) lreset = .true.
      if (ier_init.ne.4 .and. ier_init.ne.bad_jacobian_flag) then
         call vsetup (lreset, iseq_count)
      else
         iequi = 0
         if (lfreeb) ivac = 1    !!Must restart vacuum calculations if free boundary
      end if   

!
!     READ INPUT FILES INDATA, MGRID_FILE
!     USE ISEQ_COUNT TO AVOID REREADING MGRID FILE
!
      call readin 
     1    (input_file, iseq_count, lfirst_call, ier_flag, lscreen)
      if (ier_flag .ne. 0) goto 1000      
 
!
!     COMPUTE INVARIANT ARRAYS 
!
      call fixaray
    
      if (ier_init .ne. 4) then
         write (nthreed, 230)
  230 format(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)
      end if

!
!     COMPUTE INITIAL SOLUTION ON COARSE GRID
!     IF PREVIOUS SEQUENCE DID NOT CONVERGE, DO COARSE RESTART
!

      if (lreseta) then
        igrid0 = 1
      else
        igrid0 = multi_ns_grid
      endif


      imovephi = 0
      ns_min = 0
      jacob_off = 0
      ier_flag = ier_init
      if (ier_flag .eq. bad_jacobian_flag) jacob_off = 1
      if (all(ns_array .eq. 0)) then
         ier_flag = 8
         goto 1000
      end if
         
      do igrid = igrid0, multi_ns_grid + jacob_off
         if (jacob_off.eq.1 .and. igrid.eq.igrid0) then
!           TRY TO GET NON-ZERO JACOBIAN ON A 3 PT RADIAL MESH         
            nsval = 3
            ftolv = 1.e-7_dp
         else
            nsval = ns_array(igrid-jacob_off)
            if (nsval .le. ns_min) cycle
            ns_min = nsval
            ftolv = ftol_array(igrid-jacob_off)
         end if
         if (igrid .eq. igrid0) then
            interp_flag = .false.
         else
            interp_flag = .true.
         endif
         call eqsolve (nsval, interp_flag, ier_flag, lreseta, 
     1                 lfirst_call, lscreen, reset_file_name)
         if (imatch_phiedge .eq. 3) imovephi = 1
         if (ier_flag .ne. 0 .and. ier_flag .ne. 4) exit
      end do
 
  100 continue
 
      call second0 (timeoff)
      timer(0) = timer(0) + timeoff - timeon
 
!
!     WRITE OUTPUT TO THREED1, WOUT FILES; FREE MEMORY ALLOCATED GLOBALLY
!
 1000 call fileout (iseq_count, ier_flag, lscreen) 
      
      end subroutine runvmec
      
      subroutine vsetup (lreset, iseq_count)
      use vmec_main
      use vacmod
      use realspace
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iseq_count
      logical :: lreset
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: mpol_default = 6
      integer, parameter :: ntor_default = 0
C----------------------------------------------- 
!
!     Default values
!
      lasym = .false.
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0
      nzeta  = 0
      lfreeb = .true.              !reset later after input file read
      lrecon = .true.
      loldout = .false.
      ldiagno = .false.
      lmac = .false.
      ledge_dump = .false.

      z00 = zero
      mgrid_file = 'NONE'

!
!     ZERO ARRAYS WHICH MAY BE READ IN
!
      rbc = zero
      rbs = zero
      zbc = zero
      zbs = zero
 
      am = zero
      ai = zero
      ac = cbig
      ac_form = 0
      ac_form = 0
      aphi = zero
      extcur = zero
      bloat = one

      ns_array = 0
      nsin = 0
      ftol_array = zero
      raxis = zero
      zaxis = zero
      gamma = zero
      spres_ped = one
       
      iequi = 0
      ivac  = -1
      delbsq = one
      delt = 1.1
      tcon0 = one
      curtor = 1.e30_dp
      time_slice = zero
      if (lreset) then
         pfac   = one
         phifac = one
         timer = zero
      endif
      fsqr = one
      fsqz = one
      ftolv = fsqr
!
!     Reconstruction stuff
!
      lpofr = .true.
      lpprof = .true.
      icurrout = 0
      total_chi_square_n = one
      total_chisq_n0 = total_chi_square_n
      nchistp = 0
      imse = -1
      itse = 0
      isnodes = 0
      ipnodes = 0
      iopt_raxis = 1
      imatch_phiedge = 1
      nflxs = 0
      nbfld = 0
      mseangle_offset = zero
      mseangle_offsetm = zero
      pres_offset = zero
      sigma_current = 1.e30_dp
      sigma_delphid = 1.e30_dp
      tensi = one
      tensp = one
      tensi2 = zero
      fpolyi = one
      presfac = one
      phidiam = 1.e30_dp
 
      mseprof = one
      indxflx = 0
      indxbfld = 0
      sigma_stark = 1.1*cbig
      sigma_thom = 1.1*cbig
      sigma_flux = 1.1*cbig
      sigma_b = 1.1*cbig
!
!     FREE-BOUNDARY STUFF, READ IN FIRST TIME ONLY
!
      if (iseq_count .eq. 0) then
        nbcoil_max = 0
        nlim_max = 0
      end if

      end subroutine vsetup

      
      subroutine wroutlim(nwout)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nwout
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer l1, l2, ltouch, limpts, n, i
      integer, dimension(nlim_max*nlim) :: l1min,l2min
      real(rprec) :: sepmin
      character, dimension(nlim_max,nlim) :: closest
C-----------------------------------------------
 
 
      if (nlim .le. 0) return 
 
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            closest(l2,l1) = ' '
            seplim(l2,l1) = sqrt(abs(seplim(l2,l1)))
         end do
      end do

      sepmin = seplim(1,1)
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            if (seplim(l2,l1) .lt. sepmin)
     1      sepmin = seplim(l2,l1)
         end do
      end do

!     Do it this way to catch multiple mins
      ltouch = 0
      do l1 = 1, nlim
         do l2 = 1, limitr(l1)
            if (abs(seplim(l2,l1)-sepmin) .lt. 1.e-5_dp) then
               closest(l2,l1) = '*'
               ltouch = ltouch + 1
               l2min(ltouch) = l2
               l1min(ltouch) = l1
            endif
         end do
      end do
      write (nwout, 702) ltouch
      write (nwout, 704) (l1min(i),l2min(i),rlim(l2min(i),l1min(i)),
     1   zlim(l2min(i),l1min(i)),seplim(l2min(i),l1min(i)),i=1,ltouch)
  702 format(8i10)
  704 format(2i5,3e20.13)
      if (lrecon) then
      write (nthreed, 705)
  705 format(/,' PLASMA BOUNDARY-LIMITER PROXIMITY'/
     1   ' POINTS OUTSIDE PLASMA HAVE RESIDUE < 0.5')
      write (nthreed, 707)
  707 format(/,13x,' nlim    n       R         Z        Residue',
     1   '          min |d|        nearest'/13x,' ----',4x,'-',7x,'-',9x
     2   ,'-',8x,'-------',10x,'-------',8x,'-------')
 
      do n = 1, nlim
         limpts = limitr(n)
         do i = 1, limpts
            write (nthreed, 708) n, i, rlim(i,n), zlim(i,n), reslim(i,n)
     1         , seplim(i,n), closest(i,n)
         end do
      end do
  708 format(13x,i5,i5,2f10.3,1pe15.4,5x,1pe12.4,8x,a)
      end if

      end subroutine wroutlim

