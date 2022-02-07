#!/bin/sh
# LAB eliminated splines for dapha/dlambda calculation, changed to numerical
# derivatives, and used "packed"  (94+40) scheme for evaluation w of lambda
# integral.  Meshes for trapped fraction and W integrals are now seperate.
# 41 mesh points was sufficient for ftrapped convergenec to better than 10-3.
# This coding is "hardwired" and must be changed explicitly if so desired. 
# The integrals subroutine was subsumed into woflam.
# Output files were reduced to two: answers.ext with details, and
# jBsb.ext for values of JdotB
# LAB made analytic combination of two pairs of terms in W
# LAB fixed error with extra <B**2>/Bmax**2 in "other" terms.  Increased
# their value by almost a factor of two.
# SPH finished most of allocation modifications
# LAB  changed calculation of H2 to reflect proper flux surface average
# LAB  shifted aipsi and gpsi (BOOZER I and G) over one mesh point--
# original copy was 2:nsb to 2:nsb rather than 1:irup
# LAB  fixed trapped particle fraction calculation--1/B**2 jacobian was not 
# in flux surface average of v||/b.  This impacts other parts of the calculations 
# that us this quantity.
# LAB fixed output mixup for bsout, cleaned up formating, added al31s to write
# SPH started array allocation rewrite
# LAB added JB to calculation, removed extra output files, added jBbs.ext output. 
# LAB fixed formatting problem for last part of "answers" output files
# LAB/SPH added sign_jacobian to bs current
# LAB improved tokamak trapped fraction by using aspect ratio from fields
# instead of scaling outermost surface
# LAB  6/25/00 version 7 has allocated angular increments--based on 2*mbuse + 1 and
# 2*nbuse + 1  (tested on tj-2).  In addition, a fine grained search for
# bmax was implemented--this eliminated jumps in current at radii were there
#  were relatively large changes in the location of max due to flat fields.
# damp was changed to damp**2 to better relflect "collision frequency"
# damp variable name was changed to damp_bs to reflect change
# demo, read me files added for users
# end version 7
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
      module vmec0 
      use kind_spec
      character*5, parameter :: version_ = '7.00' 
      end module vmec0
      
      module trig
      use kind_spec
      real(rprec), allocatable, dimension(:) :: trigsu
      real(rprec), allocatable, dimension(:) :: trigsv
      end module trig 

      module parambs 
      use vmec0
      use bootsj_input
      real(rprec), parameter :: pi = 3.141592653589793238462643_dp
      real(rprec), parameter :: dmu0 = 4.0e-7_dp*pi
      integer, parameter :: nlambda = 101       !this now only governs the trapped calcualtion
      integer :: nthetah, nzetah                !poloidal, toridal grid parameters
      integer, parameter :: nthetahm = 32       !poloidal, toroidal  grid refinement 
      integer, parameter :: nzetahm = 16        !for bmax   
      integer, parameter :: iotasign = 1 
      real(rprec), parameter :: zetasign = -1
      integer :: irdim, irup  
      integer, dimension(:), allocatable :: jlist, jlist_idx
      real(rprec)  delta_rho_1
      real(rprec), dimension(:), allocatable :: flux,
     1    aiogar, aipsi, gpsi, pres1, 
     2    betar, dense, densi, tempe1, tempi1 
      integer, dimension(:), allocatable :: idx
      real(rprec), dimension(:,:), allocatable ::
     1  bfield, gsqrt_b,  b2obm, omb32, bfieldm,    !  arrays to store quantities on the
     2  sinmi, sinnj, cosmi, cosnj,                 !  theta phi grid
     2  sinmim, sinnjm, cosmim, cosnjm              !  quantities that are on the theata

      real(rprec), dimension(:), allocatable ::
     1    dibs, aibs, dibst, aibst, phip,
     2    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti, 
     3    bstempte, bstempti, qsafety, capr, caps, h2, 
     4    ftrapped, fpassing, epsttok, fttok, b2avg, 
     5    gbsnorm, aiterm1, other1, ajBbs,
     6    rhoar, bsnorm, fptok, amain, d_rho,  
     7    bmax1, thetamax, zetahmax  
      real(rprec)   alphae, alphai , psimax
      real(rprec)   temperho1, tempirho1, densrho1
      real(rprec), dimension(:,:,:), allocatable :: amnfit 
      real(rprec)   periods  
      real(rprec), dimension(:), allocatable :: theta, zetah
      real(rprec), dimension(nthetahm) :: thetam 
      real(rprec), dimension(nzetahm) :: zetahm 
      complex(rprec), dimension(:,:), allocatable ::
     1    dmn, fmn, alpha1mn
      real(rprec), dimension(:,:), allocatable :: rfmn
      real(rprec)   avgbobm2, sum_gsqrt_b
      real(rprec), dimension(136) :: w1, alamb_h
CCCCfix needed--put nlambda variables into the trapped fraction calculation
      real(rprec)   dlambda, drho
      real(rprec)  sign_jacobian
      logical, dimension(:), allocatable :: lsurf
      logical l_boot_all, lscreen
      end module parambs 
EOF
cat > bootsj.f << "EOF"
      program driver
      use kind_spec
      use parambs, only: lscreen
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, i, numargs, iunit_in = 10
      real(rprec) :: t1, t2
      character*120 :: extension
      character*120 :: arg1, arg2
      real(rprec) :: curtor
C-----------------------------------------------
!
!     driver: reads from command line file the wout file extension and surfaces
!     (half-radial) at which the bootstrap current is  required
!     writes the bootstrap current, JdotB to a file, jBbs.extension
!           
!     call this as follows:   xbootjs input.boots (T or F)
!
!     where input.boz contains the wout file extension and the jrad values (as a
!     blank-delimited list, not necessarily on a single line):
!
!     FILE_EXTENSION
!     2  3   5   10  12
!
!     The surface numbers are relative to half-mesh vmec quantities.
!     thus, 2 is the first half mesh point.
!
!     The optional (T) or (F) argument allows (default) or suppresses output to screen.
!
      lscreen = .true.
      
      call getcarg(1, arg1, numargs)
      if (numargs .gt. 1) call getcarg(2, arg2, numargs)

      if (numargs .lt. 1 .or.
     1   (arg1 .eq. '-h' .or. arg1 .eq. '/h')) then
         print *,
     1   ' ENTER INPUT FILE NAME ON COMMAND LINE'
         print *,' For example: xbootsj in_bootsj.tftr'
         print *
         print *,
     1   ' Optional command line argument to suppress screen output:'
         print *,' xbootsj input_file <(T or F)>'
         print *
         print *,' where F will suppress screen output'
         stop
      else if (numargs .gt. 1) then
         if (arg2(1:1).eq.'f' .or. arg2(1:1).eq.'F') lscreen = .false.   
      endif


!      INPUT FILE
!      1st line:   extension WOUT file
!      2nd line:   surfaces to be computed 


      call safe_open(iunit_in, istat, trim(arg1), 'old', 'formatted')
      if (istat .ne. 0) stop 'error opening input file in bootsj'
      
      read (iunit_in, *,iostat=istat) extension      
      
      call bootsj(curtor, trim(extension), iunit_in)
       
      end program driver


      subroutine bootsj(aibstot, extension, iunit_in)
c
!  The code BOOTSJ calculates the bootstrap current for 3D configurations.
!  The main part of the code, calculation of the geometrical factor GBSNORM,
!  was written by Johnny S. Tolliver of ORNL on the basis of
!  Ker-Chung Shaing's formulas. Other parts of the code, as well as modifications
!  to the Tolliver's part, were written by Paul Moroz of UW-Madison.
c
!  References:
c
!  1. K.C. Shaing, B.A. Carreras, N. Dominguez, V.E. Lynch, J.S. Tolliver
!    "Bootstrap current control in stellarators", Phys. Fluids B1, 1663 (1989).
!  2. K.C. Shaing, E.C. Crume, Jr., J.S. Tolliver, S.P. Hirshman, W.I. van Rij
!     "Bootstrap current and parallel viscosity in the low collisionality
!     regime in toroidal plasmas", Phys. Fluids B1, 148 (1989).
!  3. K.C. Shaing, S.P. Hirshman, J.S. Tolliver "Parallel viscosity-driven
!     neoclassical fluxes in the banana regime in nonsymmetri! toroidal
!     plasmas", Phys. Fluids 29, 2548 (1986).
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use safe_open_mod
      use trig
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer     :: iunit_in
      real(rprec) :: aibstot
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: nfax = 13
      integer, parameter :: indata0 = 7 
      integer, parameter :: jbs_file=59, ans_file=18, ans_dat_file=19
      real(rprec) :: one = 1, p5 = .5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(nfax) :: ifaxu, ifaxv
      integer :: ntrigu, ntrigv
      integer :: irho, irho1, ierr, iunit, ijbs, ians, ians_plot
      real(rprec), dimension(:), allocatable :: cputimes
      real(rprec) :: time1, timecpu, unit, file, status, err,
     1   time2, r, x, al31t, gradbs1, gradbs2, 
     2   gradbs3, gradbs4,  al31s      
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: al31
C-----------------------------------------------
 
c
c  define the starting CPU time
c
      call second0 (time1)
      timecpu = time1
      if (lscreen) write (*, 4) version_ 
    4 format(/,' Start BOOTSJ: Version ', a)

c  open files


      iunit = indata0
      call safe_open(iunit, ierr, 'input.' // trim(extension), 'old',
     1     'formatted')
      if (ierr .ne. 0) then
         print *,' Error opening input file: input.', trim(extension)
         return
      end if   

      ijbs = jbs_file
      call safe_open(ijbs, ierr, 'jBbs.'//trim(extension), 'replace',
     1     'formatted')

      ians = ans_file
      call safe_open(ians, ierr, 'answers.'//trim(extension), 'replace',
     1     'formatted')

      ians_plot = ans_dat_file
      call safe_open(ians_plot, ierr, 'answers_plot.'//trim(extension), 
     1     'replace', 'formatted')

c  read and initialize data

      call datain(trim(extension), iunit_in, iunit, ijbs, ians)
      close (iunit)

      ntrigu = 3*nthetah/2 + 1
      ntrigv = 2*nzetah
      allocate (cputimes(irup))
      allocate (
     1  dmn(-mbuse:mbuse,0:nbuse), fmn(-mbuse:mbuse,0:nbuse),
     2  rfmn(-mbuse:mbuse,0:nbuse),alpha1mn(-mbuse:mbuse,0:nbuse),
     3  trigsv(ntrigv),trigsu(ntrigu), stat=irho)

      if (irho .ne. 0) stop 'allocation error in bootsj main'
c     convert jlist to idx form, then "and" the two lists
c     SPH: Corrected error here, writing off end of jlist error when jlist <= 1
c

      jlist_idx = 0
      do irho = 1, irup
        if (jlist(irho) .gt. 1) jlist_idx(jlist(irho)-1) = 1
        idx(irho) = idx(irho)*jlist_idx(irho)
      enddo 
 
 
      l_boot_all = .true.
            
!  if any of the boozer evaluation surfaces are missing, or not requested
!  then set l_boot_all=false, total current can not be calculated

      do irho=1, irup
         if(idx(irho) .eq. 0) l_boot_all = .false.
      enddo
      if(.not.l_boot_all) then
         if (lscreen) write (*,*) 'partial surface evaluation'
         write (ians,*) 'partial surface evaluation'
      endif
      

      call fftfax (nthetah, ifaxu, trigsu)        
      call cftfax (nzetah, ifaxv, trigsv)     
 

c start main radial loop

      do irho = 1, irup

c  initialize timing for radius step
        
         call second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu

c  if there is no boozer information available, skip radial point

         if(idx(irho) .eq. 0)  cycle
         irho1 = irho - 1
         r = sqrt(rhoar(irho) + 1.E-36_dp)

c  initialize  angle grids grid for first radial evaluation.  For this 
c  and all subsequent radial evaluate B and related and quantites as well
c  as plasma derivatives and the tokamak trapped fraction.

         call bongrid(irho, ians)

c  calculate bootstrap current for equivalent tokamak

         x = fttok(irho)/(fptok(irho)+1.E-36_dp)
         
         al31t = al31(x,zeff1,alphae,alphai)

c  calculate gradient factors overall normalization inlcuding q and 
c  boozer g.

         call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)

         bsdenste(irho) = gradbs1*al31t          !due to dens gradient
         bsdensti(irho) = gradbs2*al31t          !due to dens gradient
         bstempte(irho) = gradbs3*al31t          !due to temp gradient
         bstempti(irho) = gradbs4*al31t          !due to temp gradient
         
         dibst(irho) = bsdenste(irho) + bsdensti(irho) + bstempte(irho)
     1       + bstempti(irho)                    !total Jbst

         if (l_boot_all) then
            if (irho .eq. 1) then
               aibst(1) = dibst(1)*d_rho(1)
            else
               aibst(irho) = aibst(irho1)+dibst(irho)*d_rho(irho)
            endif
         end if

c  Now start the general evaluation.
c  Find coefficients d(n,m) and evaluate the fraction trapped.

         call denmf (trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate R, S, and H2.

         call caprsh2(irho)

c  Evaluate the integral term.

         call woflam (trigsu, trigsv, ifaxu, ifaxv, irho)
 
c  Evaluate the summation term in W(lambda) that does not depend on lambda.
c  note that while the paper shows an itegral, it is the fpassing integral
c  that cancels the fpassing in the over all multiplier
 
         call othersums(irho)

c  Calculate the final answer
 
         amain(irho) = p5*(one - aiogar(irho)/qsafety(irho)) + p5*(one + 
     1      aiogar(irho)/qsafety(irho))*h2(irho)
     
         gbsnorm(irho) = amain(irho) + other1(irho) +  aiterm1(irho)
 
c- derivative of the enclosed bootstrap current over the normalized
c  toroidal flux, dIbs/ds, and the enclosed Ibs (in MA)
 
         x = ftrapped(irho)/(fpassing(irho)+1.E-36_dp)
         al31s = al31(x,zeff1,alphae,alphai)
         call grad (gradbs1, gradbs2, gradbs3, gradbs4, irho)
c
c                                                !due to dens gradient
         bsdense(irho) = gbsnorm(irho)*gradbs1*al31s
c                                                !due to dens gradient
         bsdensi(irho) = gbsnorm(irho)*gradbs2*al31s
c                                                !due to temp gradient
         bstempe(irho) = gbsnorm(irho)*gradbs3*al31s
c                                                !due to temp gradient
         bstempi(irho) = gbsnorm(irho)*gradbs4*al31s
         
         dibs(irho) = bsdense(irho) + bsdensi(irho) + bstempe(irho) + 
     1      bstempi(irho)                        !total Jbst (dI/ds)
         
c   convert to j dot B.  2*dmu0 from beta, 10**6 from MA, psimax is 1/dpsi/ds
c   and sign_jacobian takes out the sign previously used for dpsi to dA
CCCCC  flux was changed to real flux in boot vmec so that the psimax now 
CCCCC  needs an additional sign_jacobian so that the two will cancel 
             ajBbs(irho) = (2.0e6_dp)*dmu0*dibs(irho)*
     1      (pres1(irho)/betar(irho))/psimax

         if (l_boot_all) then
            if (irho .eq. 1) then
               aibs(1) = dibs(1)*d_rho(1)
            else
               aibs(irho) = aibs(irho1) + dibs(irho)*d_rho(irho)
            endif
         end if

c- the ratio of bootstrap current to that in ESD (equivalent symmetric device):

         bsnorm(irho) = dibs(irho)/(dibst(irho)+1.E-36_dp)

c  get time at end of loop if completed

         call second0 (time2)
         timecpu = time2 - time1
         cputimes(irho) = timecpu
c 

      end do
c
c- Output answers for BOOTSJ
c
      call output (cputimes, aibstot, ijbs, ians, ians_plot)
      close(ians)
      close(ians_plot)
      close(ijbs)

      call deallocate_all
      if (lscreen) write (*,400) (cputimes(irup)-cputimes(1))
  400 format(1x,'Finished BOOTSJ, time =  ', f8.3, '  sec')   

      deallocate (cputimes, trigsu, trigsv)
      deallocate (dmn, fmn, rfmn, alpha1mn)

      end subroutine bootsj
      
      
      subroutine allocate_angles
      use parambs
      use vmec0
      implicit none
      integer :: istat
      allocate(bfield(nthetah,nzetah), gsqrt_b(nthetah,nzetah),
     1  sinmi(-mbuse:mbuse,nthetah), sinnj(0:nbuse,nzetah),
     2  cosmi(-mbuse:mbuse,nthetah), cosnj(0:nbuse,nzetah),
     3  theta(nthetah), zetah(nzetah),
     4  sinmim(-mbuse:mbuse,nthetahm), sinnjm(0:nbuse,nzetahm),
     5  cosmim(-mbuse:mbuse,nthetahm), cosnjm(0:nbuse,nzetahm),
     6  b2obm(nthetah,nzetah), bfieldm(nthetahm,nzetahm), stat=istat)

      if (istat .ne. 0) stop 'allocation error in allocate_angles'
      
      end subroutine allocate_angles
      

      subroutine allocate_radial
      use parambs
      use vmec0
      use read_boozer_mod
      
      implicit none
      integer :: istat
      
      allocate (flux(irdim), qsafety(irdim), aiogar(irdim), 
     1    idx(irup), aipsi(irdim), gpsi(irdim), pres1(irdim), 
     2    betar(irdim), dense(irdim), densi(irdim), 
     3    tempe1(irdim), tempi1(irdim), lsurf(irdim), 
     4    jlist(irdim), jlist_idx(irdim), stat=istat)
      if (istat .ne. 0) return

      allocate (amnfit(irup,-mboz_b:mboz_b,0:nboz_b), stat=istat)
      if (istat .ne. 0) return
      
      allocate (dibs(irup), aibs(irup), dibst(irup), aibst(irup), 
     1    bsdense(irup), bsdensi(irup), bstempe(irup), bstempi(irup), 
     2    bsdenste(irup), bsdensti(irup), bstempte(irup), 
     3    bstempti(irup), capr(irup), 
     4    caps(irup), h2(irup), ftrapped(irup), fpassing(irup), 
     5    epsttok(irup), fttok(irup), gbsnorm(irup), aiterm1(irup), 
     6    other1(irup),  
     7    rhoar(irup), bsnorm(irup), fptok(irup), amain(irup),
     8    bmax1(irup), thetamax(irup), zetahmax(irup), 
     9    ajBbs(irup), phip(irup), d_rho(irdim),  
     A    b2avg(irup), stat = istat) 

      if (istat .ne. 0) stop 'allocation error in allocate_radial'
            
      end subroutine allocate_radial
      

      subroutine deallocate_all
      use parambs
      use vmec0
      implicit none
      
      if (allocated(amnfit)) deallocate (amnfit)
      if (allocated(dibs)) deallocate (dibs, aibs, dibst, aibst, 
     1    bsdense, bsdensi, bstempe, bstempi, bsdenste, bsdensti, 
     2    bstempte, bstempti, qsafety, capr, caps, h2, 
     3    ftrapped, fpassing, epsttok, fttok, gbsnorm, aiterm1,
     4    other1, rhoar, bsnorm, fptok, amain,  
     5    bmax1, thetamax, zetahmax, d_rho, b2avg,
     6    ajBbs, phip ,theta, zetah) 
      if (allocated(flux)) deallocate (flux, aiogar, aipsi, 
     1    gpsi, pres1, betar, dense, densi, tempe1, tempi1, lsurf,
     2    jlist, jlist_idx)

      deallocate(idx)
      deallocate(bfield, b2obm, gsqrt_b, sinmi, sinnj, cosmi, cosnj,
     1  bfieldm, sinmim, sinnjm, cosmim, cosnjm)
      
      end subroutine deallocate_all


      subroutine grad(gradbs1, gradbs2, gradbs3, gradbs4, irho)
c
c     calculate gradient factors, gradbs1 and gradbs2
c     gradbs1 - due to electron density gradient
c     gradbs2 - due to ion density gradient
c     gradbs3 - due to electron temperature gradient
c     gradbs4 - due to ion temperature gradient
c
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      real(rprec) :: gradbs1, gradbs2, gradbs3, gradbs4
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: gradbs, de, de1, di, di1, p
C-----------------------------------------------

      gradbs =-2.5_dp*gpsi(irho)*qsafety(irho)*betar(irho)*sign_jacobian
      de = dense(irho)                               !electron density
      de1 = densrho1
      di = dense(irho)/zeff1                         !ion density
      di1 = densrho1/zeff1
      p = de*tempe1(irho) + di*tempi1(irho) + 1.E-36_dp     !plasma pressure
      gradbs1 = gradbs*de1*tempe1(irho)/p  !due to electron density gradient
      gradbs2 = gradbs*di1*tempi1(irho)/p  !due to ion density gradient
c                                     !due to electron temperature gradi
      gradbs3 = alphae*gradbs*temperho1*de/p
c                                       !due to ion temperature gradient
      gradbs4 = alphai*gradbs*tempirho1*di/p
      
      end subroutine grad


      subroutine datain(extension, iunit_in, iunit, ijbs, ians)
c--
c  read and initialize data for BOOTSJ
c--
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      use read_boozer_mod
      use read_namelist_mod
      use bootsj_input
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit, ijbs, ians, iunit_in
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, n, m, j, i1, ntheta_min, nzeta_min,
     1   mn, ir, idum
      real(rprec) :: temp
      real(rprec), dimension(:), allocatable :: work
      real(rprec) :: unit, file, status, err,
     1    tempe0, tempi0, pres10, pres0, a, a1, dum
C-----------------------------------------------
!
!     Read in boozer |B|, iota, phip, pres AND allocate needed arrays
!
      call read_boozer(extension)

      jlist = 0
      read (iunit_in, *,iostat=idum) jlist
      if (idum .gt. 0) stop ' Error reading bootsj input file in DATAIN'
      
      close (iunit_in)


!
!
!- Give default values of control data otherwise read from the file 'bootin'
!
!  nrho is an obsolete variable that is no longer used.  It is set to ns-1
!  where ns is the number of whole mesh points used in vemc.
!
!  nrho is an obsolete variable that is no longer used.  It is set to ns-1
!  where ns is the number of whole mesh points used in vemc.
!
      nrho = 30                !number of rho values to use
      mbuse = 6                !number of m (poloidal) terms in B field.
      nbuse = 4                !number of nzetah (toroidal) terms in B field.
      zeff1 = 1.0_dp           !effective ion charge
      dens0 = 0.3_dp           !central electron density in 10**20 m-3
      teti = 2.0_dp            !ratio of Te/Ti for electron to ion 
                               !temperature profiles
      tempres = -one           !tempe1(i)=pres(ir)**tempres
                               !if(tempres.eq.-1.0_dp) then 
                               !tempe1(i)=sqrt(pres(i))
      damp = -0.01_dp          !superceded by damp_bs
      damp_bs = -0.01_dp       !damping factor for resonances
      isymm0 = 0               !if ne.0 then force a symmetric-device calculation,
      tempe1 = one             !1 keV default temperature
      tempi1 = zero
      ate    = zero
      ati    = zero
      ate(0) = -one
      ati(0) = -one
!
!- Read control data
!
      call read_namelist (iunit, i, 'bootin')
      if (i .ne. 0 ) stop 'Error reading bootin namelist'

!  If the new damping variable is not read in, set to to the value this gives
!  an unchanged damping factor, namely damp_bs**2 = damp.  If, in addition,
!  damp is not read in, then set damp_bs = 0.001
 
      if(damp_bs .lt. zero) then !in this case no damp_bs was read in
        if(damp .gt. zero) then
           damp_bs = sqrt(damp)
        else
           damp_bs = 0.001_dp      !in this case no damp was read in
        endif 
      endif
      teti = abs(teti)

!
!     CHECK DIMENSION SIZE CONSISTENCY
!            
      if(nboz_b .lt. nbuse) then
         if (lscreen) write(*,*) 'nbuse > nbos_b, nbuse = nboz_b'
         nbuse = nboz_b
      endif
      if(mboz_b .lt. mbuse) then
         if (lscreen) write(*,*) 'mbuse > mbos_b, mbuse = mboz_b'
         mbuse = mboz_b
      endif

      nzeta_min = 2*nbuse + 1
      ntheta_min = 2*mbuse + 1
      
      do i = 0, 6
         nzetah = 4*2**i
         if(nzetah .gt. nzeta_min) exit
         nzetah = 2*2**i * 3
         if(nzetah .gt. nzeta_min) exit
      enddo
      
      do i = 0, 6
         nthetah = 4*2**i
         if(nthetah .gt. ntheta_min) exit
         nthetah = 2*2**i * 3
         if(nthetah .gt. ntheta_min) exit
      enddo

      if(lscreen) print *, 'mbuse = ',mbuse,'nbuse = ',
     1    nbuse,'nthetah = ',nthetah, ' nzetah = ',nzetah

   90 format(5e16.8)

!     convert bmn's to amnfit's (positive m's only to positive n's only)

      amnfit = zero

      lsurf = .false.
      status = tiny(a1)
      do ir = 1, irup
         do mn = 1,mnboz_b
            m = ixm_b(mn)
            n = ixn_b(mn)/nfp_b
            if (m .gt. mboz_b) stop 'boozmn indexing conflict, m'
            if (abs(n) .gt. nboz_b) stop 'boozmn indexing conflict, n'
            if (n.lt.0) then
               m = -m
               n = -n
            end if   
            if (m.eq.0 .and. n.eq.0 .and. bmn_b(mn,ir).gt.status)
     1         lsurf(ir) = .true.
            amnfit(ir,m,n) = bmn_b(mn,ir+1)     !!2nd half grid == 1st grid pt. here 
         end do  
      end do
      
      call read_boozer_deallocate   

      zeff1 = max(one,zeff1)
c  setup s grid.  A nonuniform mesh is allowed.

      psimax = maxval(abs(flux))
      
c  we need to keep the sign of psimax    
      if(flux(irup) .lt. zero) psimax = -psimax
      
c  first normalize the mesh to 1
c  and then shift from the full to half mesh with the first point on 1, not 2
c  also need to calulate the deltas for differencing on the vmec grid 
c  special values at the ends will be calculated in bongrid

      do ir = 1, irup
        rhoar(ir) = 0.5_dp*(flux(ir) + flux(ir+1))/psimax
        d_rho(ir) = (flux(ir+1) - flux(ir))/psimax
      end do
 
c      write (ijbs, 130) irup, psimax
c  130 format(' Last flux surface used is number ',i2,' with PSIMAX =',1p
c     1   e11.4,' WB')
c
c- Switch sign of q, if desired.
c
      if (iotasign .lt. 0) then
         qsafety(:irup) = iotasign*qsafety(:irup)
      endif

c
c                                                !Boozer I/g values
      aiogar(:irup) = aipsi(:irup)/(gpsi(:irup)+1.0e-36_dp)
c
      call positiv (pres1, irup, 2) !to be sure that arrays are positive
      call positiv (betar, irup, 2)
c  evaluate electron and ion temperature profiles on vmec mesh
      if (any(ate .ne. zero)) then
         do ir = 1, irup
            tempe1(ir) = temp(rhoar(ir), ate)
         end do   
      end if   
      if (any(ati .ne. zero)) then
         do ir = 1, irup
            tempi1(ir) = temp(rhoar(ir), ati)
         end do   
      end if   
      
      tempe0 = tempe1(1)            !central electron temperature in keV
      tempi0 = tempi1(1)                 !central ion temperature in keV
      pres10 = pres1(1)                    !central total pressure in Pa
      pres0 = 1.6022E4_DP
c  if the leading coefficient on either te or ti is negative, the intent is
c  to assume a maximum density, and then calcuate the profiles using a
c  power law given by abs(tempres).  this coefficient must originally be 
c  be negative and ca not be greater than one.  That is to say, the profiles
c  are determined by ne(0) and tempres.  The ratio of Te to Ti is given by
c  teti.    
      if (tempe0.le.zero .or. tempi0.le.zero) tempres = abs(tempres)
      tempres = min(one,tempres)
c
      if (tempres .ge. zero) then   !in that case, calculate the temperature
         teti = teti + 1.E-36_dp
         a = one + one/(zeff1*teti)
c                                   !central electron temperature in keV
         tempe0 = pres10/(a*pres0*dens0)
         tempi0 = tempe0/teti            !central ion temperature in keV
         tempe1(:irup) = pres1(:irup)**tempres
c                             !suggested Te and Ti profiles are the same
         tempi1(:irup) = tempe1(:irup)
         a = tempe0/tempe1(1)
         a1 = tempi0/tempi1(1)
         tempe1(:irup) = tempe1(:irup)*a
         tempi1(:irup) = tempi1(:irup)*a1
      endif
c
      call positiv (tempe1, irup, 2)
      call positiv (tempi1, irup, 2)

      dense(:irup) = pres1(:irup)/(pres0*(tempe1(:irup)+tempi1(:irup)/
     1   zeff1)+1.E-36_dp)
      allocate(work(irdim))
      call smooth1 (dense, 1, irup, work, zero)
      call positiv (dense, irup, 2)
      i1 = irup - 1
      a = tempe1(irup) + tempi1(irup)/zeff1
      a1 = tempe1(i1) + tempi1(i1)/zeff1
      dense(irup) = dense(i1)*a1*betar(irup)/(a*betar(i1)+1.E-36_dp)
c  the above game is apparently to control the spline derivatives at the
c  edge.  
      densi(:irup) = dense(:irup)/zeff1
      dens0 = dense(1)                         !central electron density
 

c
c- Echo input data to output file.
c
      write (ians, bootin)
      
      deallocate (work)
 
      end subroutine datain


      subroutine read_boozer(extension)
      use read_boozer_mod
      use parambs
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
C-----------------------------------------------
!
!     NOTE: read_boozer_deallocate is called later, after bmn_b are assigned
! 
      call read_boozer_file(extension, istat)
      if (istat .ne. 0) then
         print *,'Error reading boozer file in BOOTSJ, istat = ', istat 
         stop
      end if

      irdim   = ns_b
      irup    = ns_b - 1                !No. points in radial (half-mesh) profiles
      periods = nfp_b                   !No. field periods
      
c  now that we know the number of radial points,  radial quantites
c  can be allocated

      call allocate_radial

c  LAB--change the indexing on aipsi and gpsi to relect half mesh status
      aipsi(1:irup)    = buco_b(2:ns_b)                 !Boozer I 
      gpsi (1:irup)    = bvco_b(2:ns_b)                 !Boozer g
      qsafety(1:irup)  = one/(iota_b(2:ns_b) +
     1                        sign(1.0e-14_dp,iota_b(2:ns_b)))
      pres1(1:irup)    = pres_b(2:ns_b)
      betar(1:irup)    = beta_b(2:ns_b)
      idx(1:irup)      = idx_b(2:ns_b)
      flux(2:ns_b)     = phi_b(2:ns_b)
      phip(1:irup)     = phip_b(2:ns_b)
      
      sign_jacobian = one   !version 6.1 phi_b has the sign of the physical flux
c     phip_b retains the internal vmec convention.
      if( gpsi(irup)*phip_b(ns_b) <= zero) sign_jacobian = -one 
      flux(1) = zero
  
      end subroutine read_boozer


      subroutine positiv(ya, n, ivar)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ivar
      real(rprec), dimension(*) :: ya
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       real(rprec), parameter :: zero = 0, D36 = 1.E-36_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
C-----------------------------------------------
c--
c     make corrections to ensure a positive array
c--
c
      if (ivar .eq. 1) then
         ya(:n) = abs(ya(:n)) + D36
      else
         where (ya(:n) .le. zero) ya(:n) = D36
      endif
 
      end subroutine positiv


      subroutine smooth1(ya, n1, n2, wk, frac)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n1, n2
      real(rprec) frac
      real(rprec), dimension(*) :: ya, wk
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
       real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n11, n21, n, i
      real(rprec) :: as, a, a1, a2, a3
C-----------------------------------------------
c--
c       smoothing real array ya(*) for i: n1.le.i.le.n2
c     frac - defines how strong is smoothing;
c     if frac=0 then only smoothes the peaks
c--
      n11 = n1 + 1
      n21 = n2 - 1
      n = n21 - n11
      if (n .le. 0) return                    !too little number of points
c
c- check the edges
c
      as = sum(abs(ya(n11:n21-1)-ya(n11+1:n21)))
      as = as/n
c
      if (as .le. zero) then
         ya(n1) = ya(n11)
         ya(n2) = ya(n21)
         return 
      else
         a = 3*(as + abs(ya(n11)-ya(n11+1)))
         a1 = abs(ya(n1)-ya(n11))
         if (a1 > a) ya(n1) = 2*ya(n11) - ya(n11+1)
         a = 3*(as + abs(ya(n21)-ya(n21-1)))
         a1 = abs(ya(n2)-ya(n21))
         if (a1 > a) ya(n2) = 2*ya(n21) - ya(n21-1)
      endif
c
c- work array
c
      wk(n1:n2) = ya(n1:n2)
c
c- check for strong peaks and remove
c
      do i = n11, n21
         a1 = .5_dp*(ya(i+1)+ya(i-1))
         a2 = 3*abs(ya(i+1)-ya(i-1))
         a3 = abs(ya(i)-a1)
         if (a3 > a2) ya(i) = a1
      end do
c
c- smoothing with a factor frac
c
      if (frac .le. zero) return 
      ya(n11:n21) =(ya(n11:n21)+.5_dp*(wk(n11+1:n21+1)+wk(n11-1:n21-1))*
     1   frac)/(one + frac)

      end subroutine smooth1


      subroutine bongrid(irho, ians)

c  First time through, form the THETA-ZETAH grid 
c  Every time through, evaluate B and QSAFETY on the grid for this RHO.
c  Here, m is the poloidal mode number and n (and nh, for n-hat) is the
c  toroidal mode number.

C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho, ians
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, nh, m, mbuse1, imax, jmax
      integer, save :: ihere = 0
      integer :: ij_max(2)
      real(rprec), save :: twopi, dth, dzetah, dthm, dzetahm
      real(rprec) :: d, dbmintok, sbmaxtok,
     1 b, del
C-----------------------------------------------


      if (ihere .eq. 0) then
      
c  allocate all of the variables needed on the theta-zeta grid

         call allocate_angles
     
c-----------------------------------------------------------------------     
c  Form the THETA-ZETAH grid.

         twopi = 8*atan(one)
         dth = twopi/nthetah
         dzetah = twopi/nzetah
  
         do i = 1, nthetah
            theta(i) = (i - 1)*dth
         end do
         do j = 1, nzetah
            zetah(j) = (j - 1)*dzetah
         end do
         
!   Form a finer mesh to evaluate fine grained bmax.

         dthm = dth/(nthetahm-1)
         dzetahm = dzetah/(nzetahm-1)

!  load sin and cosine arrays         
         
         do j = 1, nzetah
            do nh = 0, nbuse
               sinnj(nh,j) = sin(zetasign*nh*zetah(j))
               cosnj(nh,j) = cos(zetasign*nh*zetah(j))
            enddo
         enddo
         do i = 1, nthetah         
            do m = -mbuse, mbuse
               sinmi(m,i) = sin(m*theta(i))
               cosmi(m,i) = cos(m*theta(i))
            end do
         end do
         ihere = 1
      endif          !end of ihere initial evaluation
         
c  Now evaluate B on the theta-zetah grid for this RHO using the epsilons just
c  found.  Loop over theta and phihat, summing up all the (m,nh) terms in B.

      do j = 1, nzetah
         do i = 1, nthetah
            b = zero
            do m = -mbuse, mbuse
               do nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1               (cosmi(m,i)*cosnj(nh,j)-sinmi(m,i)*sinnj(nh,j))
               end do
            end do
            bfield(i,j) = abs(b)
         end do
      end do
      
c   find max of b on global mesh      

      ij_max = maxloc(bfield)
      imax = ij_max(1)
      jmax = ij_max(2)

c  use the theta and zeta from this search as the center of a finer search

c  first form the grid

      thetam(1) = theta(imax) - dth/2
      zetahm(1) = zetah(jmax) - dzetah/2
      
      do i = 2, nthetahm
         thetam(i) = thetam(i-1) + dthm
      enddo
      do j = 2, nzetahm
         zetahm(j) = zetahm(j-1) + dzetahm
      enddo

c  load the sines and cosines on the finer mesh
         
      do j = 1, nzetahm
         do nh = 0, nbuse
            sinnjm(nh,j) = sin(zetasign*nh*zetahm(j))
            cosnjm(nh,j) = cos(zetasign*nh*zetahm(j))
         enddo
      enddo
      do i = 1, nthetahm         
         do m = -mbuse, mbuse
            sinmim(m,i) = sin(m*thetam(i))
            cosmim(m,i) = cos(m*thetam(i))
         end do
      end do

c  evaluate b on the finer mesh
      
      do j = 1, nzetahm
         do i = 1, nthetahm
            b = zero
            do m = -mbuse, mbuse
               do nh = 0, nbuse
                  b = b + amnfit(irho,m,nh)*
     1              (cosmim(m,i)*cosnjm(nh,j)-sinmim(m,i)*sinnjm(nh,j))
               end do
            end do
            bfieldm(i,j) = abs(b)
         end do
      end do
      

c- evaluate bmax1(irho), thetamax(irho), b2avg(irho), and zetahmax(irho)
c  based on finer mesh evaluation

      ij_max = maxloc(bfieldm)
      imax = ij_max(1)
      jmax = ij_max(2)
      thetamax(irho) = thetam(imax)
      zetahmax(irho) = zetahm(jmax)
      bmax1(irho) = bfieldm(imax,jmax)

c  evaluate jacobian.  Leave off flux surface quantites.  Evaluate
c  the sum for later use.  Note that the value is not scaled to
c  Bmax.
      
      gsqrt_b(:nthetah,:nzetah) = one/bfield(:nthetah,:nzetah)**2
      sum_gsqrt_b = sum(gsqrt_b)
     
c  find b2avg LAB--boozer paper uses both bzero2 (1/leading coefficient of 1/b**2 expansion
c  and <b**2>.  They are the same (in boozer coordinates).  Both result from 
c  flux surface averages. I will use <b2> everywhere
 
      b2avg(irho) = sum(bfield**2 * gsqrt_b)/sum_gsqrt_b
      
c  Scale the array BFIELD so that it contains B/Bmax instead of B.

      bfield(:nthetah,:nzetah) = bfield(:nthetah,:nzetah)/bmax1(irho)
      where(bfield .gt. one) bfield = one
      b2obm = bfield**2
               
c   pressure related derivatives are needed
c   first calculate difference denominators for pressure related derivatives
c   difference array has differences on the full mesh

c   use a parabolic fit near rho = 0 to get derivatives at 1st half mesh

      if(irho .eq. 1) then
         drho = (rhoar(2)**2 - rhoar(1)**2)/(2*rhoar(1))

c   use slope of last two points at outer rho point

      elseif(irho .eq. irup) then
         drho = 0.5_dp*(d_rho(irho)+d_rho(irho-1))
            
c  all other points

      else
         drho = d_rho(irho) + 0.5_dp*(d_rho(irho+1)+d_rho(irho-1))
      endif
             
c  evaluate Electron temperature gradients in Kev

      if (irho .ne. 1 .and. irho .ne. irup) temperho1 = 
     1   (tempe1(irho+1)-tempe1(irho-1))/drho
      if (irho .eq. 1) temperho1 = (tempe1(irho+1)-tempe1(irho))/drho
      if (irho .eq. irup) temperho1=(tempe1(irho)-tempe1(irho-1))/drho
         
c evaluate Ion temperature gradients in Kev

      if (irho .ne. 1 .and. irho .ne. irup) tempirho1 =
     1  (tempi1(irho+1)-tempi1(irho-1))/drho
      if (irho .eq. 1)tempirho1 = (tempi1(irho+1)-tempi1(irho))/drho
      if (irho .eq. irup)temperho1=(tempi1(irho)-tempi1(irho-1))/drho
         
c  evaluate electron density gradients in 10**20 m-3

      if (irho .ne. 1 .and. irho .ne. irup) densrho1 = 
     1  (dense(irho+1)-dense(irho-1))/drho
      if (irho .eq. 1) densrho1 = (dense(irho+1)-dense(irho))/drho
      if (irho .eq. irup) densrho1 = (dense(irho)-dense(irho-1))/drho
                   
c  write out the lower order Boozer coefficients

      mbuse1 = mbuse
      if (mbuse > 5) mbuse1 = 5            !mbuse1 > 5 will not fit on page

      write (ians, 400)
  400 format(/' nh ',$)
      do m = -mbuse1, mbuse1
         write (ians, 402) m
      end do
  402 format('   m=',i2,'   ',$)
      write (ians, '(a1)') ' '
c
      do nh = 0, nbuse
         write (ians, 406) nh, (amnfit(irho,m,nh),m=(-mbuse1),mbuse1)
      end do
  406 format(1x,i2,1p13e10.3)

c  Calculate the fraction trapped and fraction passing for the "equivalent"

      call tok_fraction(fttok(irho),irho)  

      fptok(irho) = one - fttok(irho)

      end subroutine bongrid
      
      
      subroutine tok_fraction(ft_tok,irho)
c--
c  This calculates the fraction passing and the fraction trapped
c  for a tokamak using the Lin-Liu Miller approximation with Houlberg
c  weighting factors.  (Phys. Plas. 2, 1966 (1995).
c--
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho
      real(rprec) ft_tok
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer i
      real(rprec) fub, flb, bavg, b2av, sum_gsqrt_tok
      real(rprec) bmax
C-----------------------------------------------

      real(rprec), dimension(:), allocatable ::
     1    bfield_tok,  gsqrt_tok, b2obm_tok, one_b2
      allocate(bfield_tok(nthetah),  gsqrt_tok(nthetah),
     1   b2obm_tok(nthetah), one_b2(nthetah), stat = i)
      if(i .ne. 0) stop 'allocation error of tokamak fields'
     

c  make symetric copies of field quantites--average over zeta to extract
c  symmetric component.

      do i =1,nthetah
         bfield_tok(i) = sum(bfield(i,:nzetah))/nzetah
      enddo
      
c  renormalize to max tok field equivalent

      bmax = maxval(bfield_tok)
      bfield_tok = bfield_tok/bmax
      where(bfield_tok .gt. one) bfield_tok = one

c  calculate 1/b**2, b**2

      one_b2 = one/bfield_tok**2
      b2obm_tok = bfield_tok**2
      
c  jacobian only includes 1/b**2 component and is normalized to bmax tok
c  integrate over jacobian to obtain normalization for averages      
      sum_gsqrt_tok= sum(one_b2)
      
      
c  find average b, b**2
     
      bavg = sum(bfield_tok*one_b2)/sum_gsqrt_tok
      b2av = sum(b2obm_tok*one_b2)/sum_gsqrt_tok
      
c  find upper bound

      fub = one-(one-sqrt(one-bavg)*(one+0.5_dp*bavg))*b2av/bavg**2
      
c  find lower bound

      
c  minus <1/b**2>

      flb = - sum(one/b2obm_tok**2)/sum_gsqrt_tok
      
c  plus <sqrt(1-b)/b**2>

      flb = flb + sum(sqrt(one-bfield_tok)
     1    /b2obm_tok**2)/sum_gsqrt_tok
     
c  plus <sqrt(1-b)/b)>/2

      flb = flb + 0.5_dp*sum(sqrt(one-bfield_tok)
     1    /bfield_tok/b2obm_tok)/sum_gsqrt_tok
     
c  1+(previous stuff)*<b**2>

            
      flb = one + flb*b2av

c finally combine upper and lower limits with Houlberg factors      
      ft_tok = 0.25_dp*flb + 0.75_dp*fub
      
      deallocate(bfield_tok, gsqrt_tok, b2obm_tok, one_b2)
      
      return
      end subroutine tok_fraction


      function al31 (x, zeff, alphae, alphai)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, zeff, alphae, alphai
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: z, z2, d, a, al31
C-----------------------------------------------
c-
c calculates L31 transport coefficient, according to
c S. Hirshman, Phys. Fluids 31, 3150 (1988)
c x - ratio of trapped to circulated particles
c-
      z = zeff                              !effective ion charge number
      z2 = z**2
      d=1.414_dp*z+z2+x*(0.754_dp+2.657_dp*z+2.0_dp*z2)+
     1   x*x*(0.348_dp+1.243_dp*z+z2)
      a = 0.754_dp + 2.21_dp*z + z2 + x*(0.348_dp + 1.243_dp*z + z2)
      al31 = x*a/d
      alphae = 1 - (0.884_dp + 2.074_dp*z)/a
      alphai = 1 - 1.172_dp/(1 + 0.462_dp*x)

      end function al31


      subroutine denmf(trigsu, trigsv, ifaxu, ifaxv, irho)
c  Evaluate the coefficients d(m,n) using CRAY fft991, cfft99
c  vectorized 1D FFT routines.  Also evaluate fraction trapped and fraction
c  passing.  
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      integer, dimension(13) :: ifaxu, ifaxv
      real(rprec), dimension(3*nthetah/2 + 1) :: trigsu
      real(rprec), dimension(2*nzetah) :: trigsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, D18 = 1.0e-18_dp,
     1 one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, i
      real(rprec), dimension(nthetah + 2,nzetah) :: a11
C-----------------------------------------------
c
c- First evaluate the complex coefficients d(n,m).
c
c- Load the the array a11 with (Bm/B)**2 * (1-B/Bm)**1.5.
c  Remember that the arrays BFIELD now contains B/Bm.
c  this is a flux surface average

      if (isymm0 .eq. 0) then
         a11(:nthetah,:nzetah) = 
     1      (abs(one - bfield(:nthetah,:nzetah)) + D18)**1.5_dp
     2      *b2avg(irho)/(bmax1(irho)*bfield(:nthetah,:nzetah))**2
         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero
 
         call do_fft (a11, dmn, trigsu, trigsv, ifaxu, ifaxv, nthetah, 
     1      nzetah, mbuse, nbuse)
    
      endif
   
      avgbobm2 = b2avg(irho)/bmax1(irho)**2

c  Now calculate the fraction passing and the fraction trapped.

      call fraction(irho)
 
      end subroutine denmf


      subroutine reorganz(coefs, mbuse, nbuse, factor, a1, 
     1   ntheta, nzeta)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mbuse, nbuse, ntheta, nzeta
      real(rprec) factor
      real(rprec), dimension(ntheta + 2,nzeta) :: a1
      complex(rprec), dimension(-mbuse:mbuse,0:nbuse) :: coefs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, m, i, n, ia
C-----------------------------------------------
c--
c  Version reorganz creates a COMPLEX coefs array.
c  Subroutine to reorganize and scale the output from fft991 and cft99 into complex
c  coefficients with index (m,n).  The coefficients are scaled by the
c  factor FACTOR.  This factor should be 1/(number of points in the
c  forward complex transform direction).  The original coefficients appear
c  in array a1.  The scaled coefficients will be
c  copied to array COEFS in the calling list.  Here, m is the poloidal
c  mode number and n is the toroidal mode number/periods. Because we use two
c  FORWARD transforms, we must flip the sign of m to get the desired nu
c  argument (u=theta, v=zeta)
c--
c
c- Because of (anti)symmetry, only nonnegative values of n are needed.  We also
c  only fill to m = mbuse and n = nbuse since this is all that will be
c  used in the sums over m and n.
c  Therefore, only i = 1 to mbuse+1 (for m = 0 to mbuse) and i = NTH+1-mbuse
c  to nth (for m = -mbuse to -1) and only j = 1 to nbuse+1 (for n = 0 to nbuse)
c  are needed.
c
      do j = 1, nbuse + 1
         n = j - 1
         do i = 1, mbuse + 1
            m = i - 1
            ia = 2*i - 1
            coefs(-m,n) = factor*cmplx(a1(ia,j),a1(ia+1,j))
         end do
      end do
 
      do j = nzeta + 1 - nbuse, nzeta
         n = -(nzeta + 1) + j
         do i = 1, mbuse + 1
            m = i - 1
            ia = 2*i - 1
            coefs(m,-n) = factor*cmplx(a1(ia,j),(-a1(ia+1,j)))
         end do
      end do
 
      coefs(1:mbuse,0) = coefs(-1:-mbuse:-1,0)
 
      end subroutine reorganz


      subroutine do_fft(a11, answer_mn, trigsu, trigsv, ifaxu, ifaxv, 
     1   ntheta, nzeta, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ntheta, nzeta, mbuse, nbuse
      integer, dimension(*) :: ifaxu, ifaxv
      real(rprec), dimension(ntheta + 2,nzeta) :: a11
      real(rprec), dimension(3*ntheta/2 + 1) :: trigsu
      real(rprec), dimension(2*nzeta) :: trigsv
      complex(rprec), dimension(-mbuse:mbuse,0:nbuse) :: answer_mn
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: inc, jump, isign, jumpv, incv
      real(rprec), dimension(:), allocatable :: work11
      real(rprec) :: factor
C-----------------------------------------------
 
      allocate (work11(nzeta*(ntheta+2)))
c
c- Forward real to complex transform in the theta direction with index n.
c  i.e., integral of exp(-i*n*theta) * f(theta,zetah).
c
      inc = 1
      jump = ntheta + 2
      isign = -1
 
      call fft991(a11,work11,trigsu,ifaxu,inc,jump,ntheta,nzeta,isign)
 
c
c- now forward transform in the zetah direction with index m.
c  i.e., integral of exp(-i*m*zetah) * [theta transform of f(theta,zetah)]
c
      jumpv = 1
      incv = jump/2
      call cfft99(a11,work11,trigsv,ifaxv,incv,jumpv,nzeta,incv,isign)
      
      deallocate (work11)

c
c- Now reorganize and scale these to get the complex d(n,m)
c  FACTOR = 1 / (number of points used in the forward transform direction).
c  Because of (anti)symmetry, only nonnegative m values are needed.  We also
c  only fill to n = mbuse and m = nbuse since this is all that will be
c  used in the sums over n and m.
c
      factor = one/nzeta
c  store a11 in answer_mn array
      call reorganz (answer_mn, mbuse, nbuse, factor, a11, ntheta, 
     1   nzeta)
 
      end subroutine do_fft


      subroutine fraction(irho)
c--
c  This calculates the fraction passing and the fraction trapped.
c--
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: one = 1, zero = 0
      integer, parameter :: n_lambda = 41
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l
      real(rprec), dimension(n_lambda) ::
     1   alambda, avgvpov, antgrand, answer      
      
C-----------------------------------------------

      dlambda = one/(n_lambda - 1)

c  Fill lambda mesh, evaluate V||/V

      do l = 1, n_lambda
         alambda(l) = (l - 1)*dlambda
         avgvpov(l) = sum(sqrt(abs(one - alambda(l)*bfield))*gsqrt_b)
      end do


c  Form integrand

      where (avgvpov .gt. zero) antgrand = alambda/avgvpov

c  Do integral for the passing fraction

      call simpun (alambda, antgrand, n_lambda, answer)

      fpassing(irho) = .75_dp*avgbobm2*answer(n_lambda)*sum_gsqrt_b
      ftrapped(irho) = one - fpassing(irho)

      end subroutine fraction


      subroutine caprsh2(irho)
c--
c  Evaluate CAPR, CAPS, and H2.  Here, m is the poloidal mode number
c  and n is the toroidal mode number/periods.
c--LAB--changed calculation of H2 to use proper 1/b**2 Jacobian for
c  flux surface average
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: nh, m, i, j
      real(rprec) :: h2top_th, h2top_phi, den, qn, e2
      real(rprec) ::  b2, h2top
      real(rprec) :: h2top_tmp, den_tmp, sin_eps
C-----------------------------------------------------------------------

c  the following code was added to evaluate H2 LAB
      h2top = zero
      den   = zero
      do i=1, nthetah
         do j=1, nzetah
            h2top_th = zero
            h2top_phi = zero
            den_tmp = zero
            do nh = 0, nbuse            !  evaluate i,j terms in num and denom of H2
               qn = qsafety(irho)*periods*nh*zetasign
               do m = -mbuse, mbuse
                  sin_eps = amnfit(irho,m,nh)*
     1            (sinmi(m,i)*cosnj(nh,j)+cosmi(m,i)*sinnj(nh,j))                  
                  h2top_th = h2top_th - m*sin_eps
                  h2top_phi = h2top_phi - qn*sin_eps
                  den_tmp = den_tmp - (m + qn)*sin_eps
               enddo
            enddo
            b2 = bfield(i,j)**2
            h2top = h2top + (h2top_th**2 - h2top_phi**2)/b2
            den  = den + den_tmp**2/b2
         enddo
      enddo               
          
      if (den .eq. zero) stop 'den = 0 in caprsh2'
      h2(irho) = h2top/den
      capr(irho) = (one - h2(irho))/(2*qsafety(irho))
      caps(irho) = (one + h2(irho))/2

      end subroutine caprsh2


      subroutine woflam(trigsu, trigsv, ifaxu, ifaxv, irho)

c  Evaluate the lambda-dependent part of the functions W1(lambda) 
c  for all LAMBDA.
c  Calculates the sum over n and m of (mR+nS)/(m-nq) * dalphamn/dlambda
c      * <V||/V (mn)>/<V||/V>*lambda*(-2).
c  Here, m is the poloidal mode number and n is the toroidal mode
c  number/periods.
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
      integer, dimension(13) :: ifaxu, ifaxv
      real(rprec), dimension(3*nthetah/2 + 1) :: trigsu
      real(rprec), dimension(2*nzetah) :: trigsv
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1,
     1   d18 = 1.0e-18_dp, xlam = 0.96_dp
      integer, parameter :: n_lam_coarse = 97, n_lam = 137

c  as discussed below, the mesh is split with 96 + 40 intervals

C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, j, i, n, m
      real(rprec), dimension(:,:), allocatable :: a11
      complex(rprec), dimension(:,:), allocatable ::
     1   alphamn, vmn
      real(rprec) :: qn, numer, avg_vpov, denom
      real(rprec), dimension(n_lam) :: xlam_f
      real(rprec), dimension(n_lam-1) :: xlam_h
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: sumit
C-----------------------------------------------
c
c
      allocate (a11(nthetah+2,nzetah), 
     5   alphamn(-mbuse:mbuse,0:nbuse),
     8   vmn(-mbuse:mbuse,0:nbuse), stat=l)

      if (l .ne. 0) stop 'allocation error in woflam'
      
c-----------------------------------------------------------------------
c///////////////////////////////////////////////////////////////////////
c-----------------------------------------------------------------------

c  Meshing for lamda integral.  The lambda integral has a near singular
c  value at lambda = 1.  To evaluate, construct dalpha/dlamda on full 
c  mesh.  Evaluate the rest on the full mesh.  Becasue of the (almost)
c  singularity, split the mesh into coarse (0 to 0.96) and fine (.96 to
c  1.0) components.  This was tested on both qas and qos configuratons
c  and gives results that are a percent or better on the integral and
c  ten to twenty times better for the total current.  To save storage,
c  there is only one loop on lambda.  All grids are local.  

c-----------------------------------------------------------------------
c///////////////////////////////////////////////////////////////////////
c-----------------------------------------------------------------------

 
      if (isymm0 .ne. 0) then
         aiterm1(irho) = zero
         return 
      endif

c  First contruct the meshes.  This scheme is hardwired.


      xlam_f(1) = zero
      xlam_h(1) = 0.005_dp

      do l = 2, n_lam_coarse-1
         xlam_f(l) = xlam_f(l-1) + 0.01_dp
         xlam_h(l) = xlam_h(l-1) + 0.01_dp
      enddo

      xlam_f(n_lam_coarse) = 0.96_dp               !!End coarse mesh/Begin fine mesh
      xlam_h(n_lam_coarse) = 0.9605_dp             !!Begin of fine half-mesh

      do l = n_lam_coarse+1, n_lam-1
         xlam_f(l) = xlam_f(l-1) + 0.001_dp
         xlam_h(l) = xlam_h(l-1) + 0.001_dp
      enddo
         
      xlam_f(n_lam) = one
         
c  Loop over LAMBDA = 0 to 1, calculating alpha on the full mesh, and
c  the rest of the integral on the half mesh.  The jacobian for alpha
c  puts in a b2avg but one of the b values in the denominator cancels the
c  b/bmax in alpha.  This form of the expression paralles that in the
c  hamada paper but with flux surface averages for the alpha.  This results
c  in a form that does not have the beta found in the boozer paper
c  see refereces at beginning.

      do l = 1, n_lam
         a11(:nthetah,:nzetah) =
     2      sqrt(abs(one - xlam_f(l)*bfield(:nthetah,:nzetah))+D18)
     3      *b2avg(irho)/bmax1(irho)**2/bfield(:nthetah,:nzetah)
         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero
         call do_fft (a11, fmn, trigsu, trigsv, ifaxu, ifaxv, nthetah, 
     1      nzetah, mbuse, nbuse)
         if(l .eq. 1) then
            alphamn = fmn
            cycle                       !need to calculate to alphas and difference
         end if                         !before a term for the integral can be evaluated.
         
     
c  save alpha at lambda = 1

         if(l .eq. n_lam) alpha1mn = fmn 
           
c  form d_alphalmn  at this point we have the value of alpha for l in fmn
c  and the previous value of alpha (l-1) in alphamn.  Store the difference in
c  fmn using vmn as a temp variable to hold alpha, then update alpha for
c  the next cycle.

         vmn = fmn
         fmn = fmn - alphamn
         alphamn = vmn
            
c- Now do the fft to get <exp(i(m*theta-n*zeta)V||/V> on the half mesh

         a11(:nthetah,:nzetah) = 
     1      sqrt(abs(one - xlam_h(l-1)*bfield(:nthetah,:nzetah)) + D18)
     2      *(b2avg(irho)/bfield(:nthetah,:nzetah)*bmax1(irho))**2

         a11(nthetah+1,:nzetah) = zero
         a11(nthetah+2,:nzetah) = zero
 
         call do_fft (a11, vmn, trigsu, trigsv, ifaxu, ifaxv, nthetah, 
     1      nzetah, mbuse, nbuse)
     
c  This needs to be divided by <V||/V>.  But this is exactly the real part of the
c  0,0 term of this transform.
         avg_vpov = real(vmn(0,0))
 
 
c- for only those harmonics that are going to be used in the sum.

         qn = periods*qsafety(irho)*zetasign
         rfmn(0,0) = zero
         do m = -mbuse, mbuse
            do n = 0, nbuse

               denom = m + n*qn
               if (n.ne.0 .or. m.ne.0) then

                  numer = denom/(denom**2 + (damp_bs*m)**2)
 
                  rfmn(m,n) = (m*capr(irho) + n*periods*caps(irho))*
     1             real(fmn(m,n)*vmn(m,n))*(-2*numer)
               endif
            end do
         end do
c
c- The sum in W is the SUM over all n and m (excluding (0,0)) of rfmn.
c  But only the nonegative m have been calculated and stored.  Therefore
c  the sum becomes
c    2 * SUM over m = -mbuse to mbuse and n = 1 to nbuse of rfmn
c     plus SUM over n = -mbuse to mbuse of fn0
c     minus f(0,0)
 
         w1(l-1) = xlam_h(l-1)*sumit(rfmn,mbuse,nbuse)/avg_vpov
         
c- End loop over lambda.

      end do
c  evaluate integral

      aiterm1(irho) = sum(w1)
      aiterm1(irho) = -0.75_dp*aiterm1(irho)*qsafety(irho)
     1    /ftrapped(irho) * (one + aiogar(irho)/qsafety(irho))
     
      deallocate (a11, alphamn, vmn, stat=l)
      
      end subroutine woflam


      function sumit (f, mbuse, nbuse)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mbuse, nbuse
      real(rprec), dimension(-mbuse:mbuse,0:nbuse) :: f
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, m
      real(rprec) :: sumit
C-----------------------------------------------
c-
c specific sum of terms in f(i,j) table
c-
c
      sumit = 0
      do m = -mbuse, mbuse
         sumit = sumit + sum(f(m,1:nbuse))
      end do
      sumit = 2*sumit
      sumit = sumit + sum(f(-mbuse:mbuse,0))
      sumit = sumit - f(0,0)

      end function sumit


      subroutine othersums(irho)
c--
c     calculates other1 
c     of integral of W(lambda) which is not dependant on lambda
c--
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: irho
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, m
      real(rprec) :: denom, qn
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rprec) , EXTERNAL :: sumit
C-----------------------------------------------
c
c-  m - poloidal, n - toroidal
c
      if (isymm0 .ne. 0) then
         other1(irho) = zero
         return 
      endif
c
c- Form the "other" sums.  The first uses alpha1(m,n), calculated in WOFLAM,
c  and d(m,n) from DENM.
c
c- Load  rfmn with
c
c  (m*R+n*periods*S)/(m-n*periods*q) *
c         exp(m*thetamax-n*zetahmax) * (1.5*alpha1(m,n)+d(m,n))
c
c- for only those harmonics that are going to be used in the sum.
c
      qn = periods*qsafety(irho)*zetasign
      do m = -mbuse, mbuse
         do n = 0, nbuse
            denom = m + n*qn
            if (n.ne.0 .or. m.ne.0) then
               denom = denom/(denom**2 + (damp_bs*m)**2)
               rfmn(m,n) = (m*capr(irho)+n*periods*caps(irho))*denom*(
     1            cos(m*thetamax(irho)-n*zetahmax(irho))*real(1.5_dp*
     2            alpha1mn(m,n)+dmn(m,n))-sin(m*thetamax(irho)-n*
     3            zetahmax(irho))*aimag(1.5_dp*alpha1mn(m,n)+dmn(m,n)))
            else
               rfmn(m,n) = zero
            endif
         end do
      end do
      
c  First sum,

      other1(irho) = sumit(rfmn,mbuse,nbuse)
c
c- Then multiply by all the stuff out front.
c
      other1(irho) = -other1(irho)*qsafety(irho)/ftrapped(irho)*(one + 
     1   aiogar(irho)/qsafety(irho))


      end subroutine othersums
      

      subroutine output(cputimes, aibstot, ijbs, ians, ians_plot)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use parambs
      use vmec0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ijbs, ians, ians_plot
      real(rprec) :: aibstot
      real(rprec), dimension(*) :: cputimes
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1,
     1  D18 = 1.E-18_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, ir, j, m
      real(rprec) :: aibsttot, st, x, d,
     1   err, a, b, hs
C-----------------------------------------------
      write (ijbs, 87)
      do i = 1, irup
         if(idx(i) .eq. 1) then
            write (ijbs, *) (i+1), ajBbs(i), bsnorm(i)
         endif 
      enddo
   87 format(1x,'surface #  <Jbs_dot_B>    <Jbs>/<Jbs-tok>')

      write (ians, 90)
      do i=1,irup
         if(idx(i) .eq. 1) then 
            write (ians, 100) rhoar(i),capr(i),caps(i),ftrapped(i),
     1         h2(i), amain(i),aiterm1(i),
     2         other1(i), gbsnorm(i)
         endif
      enddo
   90 format(/,1x,
     1   '  s        R           S        ftrapped       H2     ',
     2   '    amain      lam int     other      gnorm')
  100 format(1x,0pf5.3,1p8e12.4)
  
        
      write (ians, 101)
      do i=1,irup
         if(idx(i) .eq. 1) then
            write (ians, 102) rhoar(i),qsafety(i),thetamax(i),
     1          zetahmax(i),bmax1(i),ftrapped(i)/fttok(i),
     2          fptok(i)/fpassing(i),cputimes(i)
         endif
      enddo      
  101 format(/,1x,'  s        Q        thetamax     zetamax   ',
     1   '    Bmax    ft/fttok    fptok/fp    cpu secs')
  102 format(1x,0pf5.3,1p8e12.4)

      write(ians,*) '   s    gnorm       jbsnorm  ',
     1   '   dI/ds      I(s)(hm)  '
     2    ,'  j_grad_ne   j_grad_ni   j_grad_Te   j_grad_Ti '
      do i=1,irup
         if(idx(i) .eq. 1) then
         write(ians,32)rhoar(i), gbsnorm(i),bsnorm(i), dibs(i),
     1   aibs(i),  bsdense(i), bsdensi(i), bstempe(i), 
     1   bstempi(i)
         endif
      enddo
   32 format(1x,0pf5.3,1p9e12.4) 
      
        
      write (ians, 104)      
      do i=1,irup
         if(idx(i) .eq. 1) then
            write (ians, 105)i,tempe1(i),tempi1(i),dense(i),
     1      dense(i)/zeff1,   betar(i),ajBbs(i)
         end if
      enddo
  104 format(/,1x,'     Te          Ti          Ne         Ni
     1Beta          jB')     
  105 format(1x,i3,1p6e12.4)


      if(l_boot_all .and. lscreen) then
         aibstot = aibs(irup)
         write (*, '(a,f12.7,a,i3,a,i3)') ' Total bootstrap current =',
     1      aibstot , ' MA'
         aibsttot = aibst(irup)
         aibstot = aibstot*1.e6_dp                     !in Amperes
      endif

  107 format('Ibs = ',f7.4,' MA')
  175 format(2x,1p5e14.6)
  
      do i = 1, irup
         if(idx(i) .eq. 1) then
           write(ians_plot, *) rhoar(i),  gbsnorm(i), amain(i), 
     1      aiterm1(i), other1(i),
     2      dibs(i), bsdense(i), bsdensi(i), bstempe(i), bstempi(i),
     3      qsafety(i), ftrapped(i), bsnorm(i),
     4      tempe1(i),tempi1(i),dense(i), dense(i)/zeff1, betar(i),
     5      ajBbs(i)
         endif
      enddo
 
      end subroutine output


      function temp (x, T)
      use kind_spec
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, temp, T(0:10)
C-----------------------------------------------
      temp = (T(0) + x*(T(1) + x*(T(2) + x*(T(3) + x*(T(4) + 
     1     x*(T(5) + x*(T(6) + x*(T(7) + x*(T(8) + x*(T(9) + x*T(10)
     2   ))))))))))
      end function temp
      

!========================================================================
!  The following are supplemental subroutines from various libraries
!========================================================================


      SUBROUTINE SIMPUN(XX, FX, NX, AX)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      use kind_spec
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NX
      REAL(rprec), DIMENSION(NX) :: XX, FX, AX
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, two=2, six=6
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: II, IX, IC
      REAL(rprec) :: D1, D2, D3, A2, A3
C-----------------------------------------------
c---
c  II>0 then ax(i) = Integral(f(x)dx) from x=xx(1) to x=xx(i). Here f(x)=fx(i)
c  II<0 then ax(i) = Integral(f(x)dx) from x=xx(nx) to x=xx(i)
c  PROGRAM AUTHOR      J. BARISH,
c  COMPUTING TECHNOLOGY CENTER, UNION CARBIDE CORP., NUCLEAR DIV.,
c  OAK RIDGE, TENN.
c---
c
      II = 1
      IF (II < 0) GO TO 30
      AX(1) = ZERO
      DO IX = 2, NX, 2
         D1 = XX(IX) - XX(IX-1)
         AX(IX) = AX(IX-1) + D1/TWO*(FX(IX)+FX(IX-1))
         IF (NX .eq. IX) EXIT 
         D2 = XX(IX+1) - XX(IX-1)
         D3 = D2/D1
         A2 = D3/SIX*D2**2/(XX(IX+1)-XX(IX))
         A3 = D2/TWO - A2/D3
         AX(IX+1)=AX(IX-1)+(D2-A2-A3)*FX(IX-1)+A2*FX(IX)+A3*FX(IX+1)
      END DO
   20 CONTINUE
      RETURN 
   30 CONTINUE
      AX(NX) = ZERO
      DO IX = 2, NX, 2
         IC = NX + 1 - IX
         D1 = XX(IC+1) - XX(IC)
         AX(IC) = AX(IC+1) + D1/TWO*(FX(IC+1)+FX(IC))
         IF (NX .eq. IX) GO TO 20
         D2 = XX(IC+1) - XX(IC-1)
         D3 = D2/(XX(IC)-XX(IC-1))
         A2 = D3/SIX*D2**2/D1
         A3 = D2/TWO - A2/D3
         AX(IC-1)=AX(IC+1)+(D2-A2-A3)*FX(IC-1)+A2*FX(IC)+A3*FX(IC+1)
      END DO

      END SUBROUTINE SIMPUN

EOF
EOC
