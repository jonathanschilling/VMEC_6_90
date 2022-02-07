      subroutine dkes_input_prepare (arg, dkes_input_file)
c
c   This code prepares an input file for DKES based on data
c   in the boozmn.* file which is periodically written out
c   by the VMEC optimizer.  After compiling (see instructions
c   given below) and running a file called input.dkes is written.
c   DKES can then be run by typing "xdkes input.dkes".
c   This input file is based on a fixed set of parameters
c   specified below that is intended to provide
c   a rapidly evaluated (high collisionality) optimization target.
c   These can be modified to consider other parameter ranges and
c   a more complete mode spectrum as the need arises.  Parameters
c   which may be of interest to vary are:
c
c   The variables in the call to dkes_input_prepare are:
c
c     booz_file_name (=arg(1))
c             = name of boozmn file containing |B| data at various surfaces
c
c     nsurf (=arg(2))
c             = flux surface index (in boozer file) where
c                   DKES is to evaluate transport coefficients
c
c     cmul (arg(3))
c             = DKES collisionality parameter: nu/v  (in meter**-1)
c
c     efield (arg(4))
c             = DKES electric field parameter: E_s/v
c
c     lscreen (arg(5))
c             = logical variable which controls screen output
c            (= .true. screen output on, = .false. screen output off)
c
c
c
c     max_bmns = Number of Bmns which are retained for the
c                  DKES spectrum (this also influences the
c                  Fourier spectrum used for the distribution
c                  function in DKES)
c
c     legendre_modes = number of Legendre polynomials used in
c                  DKES to represent the pitch angle variation
c                  of the distribution function
c
c     coupling_order = Parameter for mode coupling order.
c                  This is the number of iterations used to
c                  generate the optimal m,n Fourier spectrum for
c                  the DKES distribution function.  At high
c                  collisionalities (cmul > 0.1), 2 iterations seems
c                  to be adequate (here adequate means that the upper
c                  and lower bounds coming out of DKES are close to
c                  to each other for the transport coefficient of
c                  interest).  As one goes to lower collisionalities,
c                  progressively more iterations must be used to get
c                  convergence of the upper/lower bounds.  This will
c                  imply increasingly longer run times for DKES.
c
c
c  Authors:  W.I. van Rij  --->  wrote original version of SPECFIL code
c            H. Maassberg  --->  Adapted SPECFIL to work within an interactive
c                                code to prepare input for DKES
c            D. Spong      --->  converted to f90, made to read new
c                                    form of Bmn file, different spectrum
c                                    generation rules tried, new output files,
c                                    interactive input turned off.  Adapted
c                                    for running with VMEC optimizer. Added
c                                    comments to the spectrum optimization.
c
c
      use kind_spec
      use date_and_computer
      use read_boozer_mod
      use read_wout_mod, rmnc_w=>rmnc, zmns_w=>zmns, lmns_w=>lmns,
     &    xm_w=>xm, xn_w=>xn, phip_w=>phip, mpol_w=>mpol,
     &    ntor_w=>ntor, nfp_w=>nfp, ns_w=>ns, mnmax_w=>mnmax
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: arg(5), dkes_input_file
c-----------------------------------------------
c   L o c a l   P a r a m e t e r s
c-----------------------------------------------
c     real(rprec), parameter :: radial_position = 0.5_dp
c     real(rprec), parameter :: cmul = 0.05_dp
c     real(rprec), parameter :: efield = 0._dp
      integer, parameter :: max_bmns = 20
      integer, parameter :: legendre_modes = 100
      integer, parameter :: coupling_order = 3
c
      integer, parameter :: mhi = 50
      integer, parameter :: ibbi = 1
      integer, parameter :: nhi = 50
      integer, parameter :: nfcd = 100
      integer, parameter :: mtd = 15
      integer, parameter :: nzd = 15
      integer, parameter :: mtd1 = mtd + 1
      integer, parameter :: nzd1 = 2*nzd + 1
      real(rprec), parameter :: aval1 = 1.e-4_dp
      real(rprec), parameter :: aval2 = 1.e-4_dp
      real(rprec), parameter :: one = 1, two = 2,
     1  half = 0.5_dp, zero = 0
      character*(50), parameter ::
     1   banner = 'THIS IS THE DKES input prep CODE Version 1.0'
c
c    aval1 = lower limit for the |B| Fourier coefficients
c    aval2 = lower limit for dominant zeta dependent spectrum
c         (note: either one or the other of these is active, but not
c          both.- DAS)
c    The parameters nhi and mhi are described in comments to follow.
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      real(rprec), dimension(:), allocatable :: bmn_local,
     >   abs_blocal, blocal_ordered
      integer, dimension(:), allocatable :: m_ordered, n_ordered
      integer, dimension(:), allocatable :: nsurf
      integer :: iunit, js, iloc, n_norm, i, j, k, ierr, istat
      real(rprec) btheta, bzeta, check, iota_local, psip,
     >  chip, TWOPI
      real(rprec) phi_edge, cmul, efield, bmn_test
      integer, dimension(mhi,nhi,3) :: nplus, nmins
      integer, dimension(2*nhi,mhi) :: mns_out = 0
      integer, dimension(2*nhi) :: n_out = 0
      integer, dimension(nfcd) :: mtheta, ntheta, mzeta, nzeta, ivec
      integer, dimension(nhi) :: nvalsb
      integer, dimension(nfcd) :: intor, impol
      integer m, n, mp, nt, mp0, nt0, mpu, ntu, ntorb, mpolb,
     >  nszt, nsth, kk, jj, ii, ntor, mpol, ij, ik, im, in, iout,
     >  ier, nord, nmin, nmax, mmin, mmax, numargs
      real(rprec), dimension(nfcd) :: afc
      real(rprec), dimension(mtd1,nzd1) :: alf
      real(rprec) a, alfm, alfu, eps, au, aa, alfz
      real(rprec) hs, vnorm, vol_inner, vol_outer,
     >  rminor_i, rminor_o
      integer index_dat, index_end
      real(rprec) time_begin, time_end
      character*(120) :: booz_input_file, wout_input_file
      character*(10) :: date0, time0, zone0
      character*(40) :: dateloc
!     character*120 dkes_input_file
      integer :: imon, index_surf
      logical :: lscreen = .true.
c-----------------------------------------------
      TWOPI = 8*atan(one)

      call second0(time_begin)
c
      read(arg(2),'(i20)') js
      read(arg(3),'(f10.5)') cmul
      read(arg(4),'(f10.5)') efield
      if (arg(5)(1:1).eq.'f' .or. arg(5)(1:1).eq.'F') lscreen = .false.
c
      if (lscreen) then
         print 48
         call date_and_time(date0,time0,zone0)
         read (date0(5:6),'(i2)') imon
         write (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1         time0(1:2),time0(3:4),time0(5:6)
         write (*,'(1x,a,/,1x,a)') banner, dateloc
         print *
         write(*,'(" js = ",i4,", ",i4," Bmns")') js,max_bmns
         write(*,'(" cmul = ",f10.6,", efield = ",f10.6)')cmul,efield
         write(*,'(" coupling order = ",i2,
     1    ", ",i4," Legendre modes")') coupling_order,legendre_modes
         print *
      end if
 100  format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
 120  format(/,' TIME IN DKES input prep CODE:',1pe12.2,' SEC')
  48  format('====================================================')
c
      index_dat = index(arg(1),'.')
      index_end = len_trim(arg(1))
      booz_input_file  = arg(1)(index_dat+1:index_end)
      call read_boozer_file (booz_input_file, k)
      if (k .ne. 0) stop 'Error reading boozmn file in DKES_INPUT_PREP'
      allocate(bmn_local(mnboz_b), abs_blocal(mnboz_b),
     1    blocal_ordered(mnboz_b), stat=ierr)
      allocate(m_ordered(mnboz_b), n_ordered(mnboz_b), stat=ierr)
      allocate(nsurf(ns_b), stat=ierr)
c
c     Read data from the wout file and allocate storage:
c
      wout_input_file = 'wout.' // booz_input_file
      call readw_and_open(wout_input_file,k)
c
c     Find out which surfaces are not included in the boozmn file
c         (i.e., by checking if Bmn's are zero)
c
      index_surf = 0
      do j = 1, ns_b
         bmn_test = zero
         do i = 1,mnboz_b
            abs_blocal(i) = abs(bmn_b(i,j))
         end do
         bmn_test = sum(abs_blocal)
         if(bmn_test .gt. 1.e-8_dp) then
            index_surf = index_surf + 1
            nsurf(index_surf) = j
         endif
      end do


      do i = 1,mnboz_b
         bmn_local(i) = bmn_b(i,js)
         abs_blocal(i) = abs(bmn_local(i))
      end do
      bmn_test = sum(abs_blocal)
      if(bmn_test .lt. 1.e-8_dp .and. lscreen) then
         write(*,'(" This surface is not in the boozmn file")')
         write(*,'("  - need to pick a different surface")')
         write(*,'(" The following surfaces are available:")')
         write(*,'(6(3x,i2))') (nsurf(j), j=1,index_surf-1)
c
         call second0(time_end)
         if (lscreen) then
            write (*,120) time_end - time_begin
            write (*,48)
         end if
c
         stop 23
      endif
c
      iota_local = iota_b(js)
c      btheta = buco_b(js)
c      bzeta = bvco_b(js)
      btheta = half*(buco(js) + buco(js-1))
      bzeta = half*(bvco(js) + bvco(js-1))
      phi_edge = abs(phi_b(ns_b))
c
      hs = one/real((ns_w - 1), rprec)
      vnorm = TWOPI*TWOPI*hs
      vol_inner = sum(vp(2:js-1))
      vol_outer = vol_inner + vp(js)
      vol_inner = vol_inner*vnorm
      vol_outer = vol_outer*vnorm
      rminor_i = sqrt(two*vol_inner/(TWOPI*TWOPI*Rmajor))
      rminor_o = sqrt(two*vol_outer/(TWOPI*TWOPI*Rmajor))
      psip = -phip_w(js)*hs/(rminor_o - rminor_i)
      chip = psip*iota_local
c       write(*,'(" ns_w = ",i3)') ns_w

c
c     Sort the Bmn's (along with associated m's and n's in order
c     of increasing abs(Bmn):
c
      do i=1,mnboz_b
         check = maxval(abs_blocal)
c    Find location of this value of check in the abs_blocal array
         do j=1,mnboz_b
            if (abs(check - abs_blocal(j)) .lt. 1.e-6_dp*abs_blocal(j))
     1      iloc = j
         end do
         blocal_ordered(i) = bmn_local(iloc)
         m_ordered(i) = ixm_b(iloc)
         n_ordered(i) = ixn_b(iloc)
         abs_blocal(iloc) = zero
      end do
!
!     Find min/max m and n out of the first max_bmns of Bmn
!     modes. Compute ntorb and mpolb for DKES.
!
      nmin = 100
      mmin = 100
      nmax = -100
      mmax = -100
      do i=1,max_bmns
         n_norm = n_ordered(i)/nfp_b
         if(m_ordered(i) .lt. mmin) mmin = m_ordered(i)
         if(m_ordered(i) .gt. mmax) mmax = m_ordered(i)
         if(n_norm .lt. nmin) nmin = n_norm
         if(n_norm .gt. nmax) nmax = n_norm
      end do
      ntorb = nmax - nmin + 1
      mpolb = mmax - mmin + 1
!
!     Write out the initial part of the DKES input file.
!
      dkes_input_file = 'input_dkes.' // booz_input_file
      iunit = 15
      call safe_open(iunit, istat, dkes_input_file, 'replace', 
     1    'formatted')
      write (iunit,'(1x,"&datain")')
      write (iunit,'(1x,"nzperiod= ",i2,",")') nfp_b
      write (iunit,'(1x,"lalpha= ",i3,", nrun = 1,")') legendre_modes
      write (iunit,'(1x,"cmul = ",e12.4,",")') cmul
      write (iunit,'(1x,"efield = ",f7.4,",")') efield
      write (iunit,'(1x,"mpolb = ",i2,",",2x,"ntorb = ",
     1     i2,",",2x,"ibbi = 1,")') mpolb, ntorb
      write (iunit,'(1x,"chip = ",f7.4,",","  psip = ",f7.4,",")') 
     1     chip, psip
      write (iunit,'(1x,"btheta = ",f7.4,","," bzeta = ",f7.4,",")')
     1     btheta, bzeta
c
c     Write out only the top max_bmns largest Bmn's:
c
      do i=1,max_bmns
         n_norm = n_ordered(i)/nfp_b
         if(n_norm .ge. 0 .and. m_ordered(i) .le. 9) then
            write (iunit,'(" borbi(",i1,",",i1,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if (n_norm .lt. 0 .and. m_ordered(i) .le. 9) then
            write (iunit,'(" borbi(",i2,",",i1,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if(n_norm .ge. 0 .and. m_ordered(i) .gt. 9) then
            write (iunit,'(" borbi(",i1,",",i2,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         else if(n_norm .lt. 0 .and. m_ordered(i) .gt. 9) then
            write (iunit,'(" borbi(",i2,",",i2,")= ",1x,e11.5,",")')
     1      n_norm,m_ordered(i),blocal_ordered(i)
         end if
      end do
c
c     Translate m, n, and Bmn arrays into arrays used in the
c       GIDKES2 code:
c
      nord = coupling_order
      eps = max(zero,min(aval1,0.01_dp))
      alfm = zero
      do i=1,max_bmns
         mp = m_ordered(i)
         nt = n_ordered(i)/nfp_b
         a = blocal_ordered(i)
         impol(i) = mp
         intor(i) = nt
         afc(i) = a
         if (mp<=mtd .and. iabs(nt)<=nzd) alf(mp+1,nzd+1+nt) = a
         if (mp.eq.0 .and. nt.eq.0) then
            mp0 = mp
            nt0 = nt
         else if (abs(a) > alfm) then
            alfm = abs(a)
            alfu = a
            mpu = mp
            ntu = nt
         endif
      end do
      alfz = aval2
c     alfz = min(half*alfm,max(0.2_dp*alfm,10*eps))
c
c     The following loop assigns modes to S(zeta) (non-zero n modes through
c     the mzeta and nzeta arrays) and S(theta) (non-zero m modes through
c     the mtheta and ntheta arrays).  There is some flexibility as how this
c     is done through the alfz parameter defined above.  If alfz is
c     relatively large (i.e., if above alfz is used) not too many modes
c     make it into the mzeta and nzeta arrays, more go into the mtheta,
c     ntheta array.  If alfz is small (comment out the above alfz line)
c     then more modes go into the mzeta and nzeta arrays and not so many
c     into the mtheta and ntheta arrays.
c
c
c-----------------------------------------------------------------------
c
c     Description of the SPECFIL code:      (W.I. van Rij)
c
c     Purpose:
c     Estimate the Fourier spectrum for the DKES2 distribution function
c
c     Let S denote the fourier spectrum of 1/B [except for (0,0) term]
c     Let S(zeta) denote the dominant zeta-dependent spectrum of 1/B
c     Define S(theta) = S - S(zeta)
c
c     Let nszt (<=nhi) be the number of Fourier modes in S(zeta)
c     Let nsth (<=mhi) be the number of Fourier modes in S(theta)
c              (nhi, mhi are defined in parameter statement)
c
c     For i-th mode of S(zeta) load m value in mzeta(i) and n value in
c     nzeta(i) [i=1,2,....,nszt]
c     For i-th mode of S(theta) load m value in mtheta(i) and n value in
c     ntheta(i) [i=1,2,....,nsth]
c
c     n values are expressed in units of nzperiod
c
c     n values must satisfy   -nhi <= n <= nhi-1
c     m values must satisfy      0 <= m <= mhi-1
c     These bounds may be violated depending on the value of the
c     coupling order NORD (error message is printed on terminal).
c     In this case, you are prompted for a reduced value of NORD.
c
c-----------------------------------------------------------------------
c     subroutines required:
c     FILIT, SORTI
c-----------------------------------------------------------------------
c
      nszt = 0
      nsth = 0
      ntorb = 0
      do m = 0, mtd
         au = zero
         l12: do n = -nzd, nzd
            a = alf(m+1,nzd+1+n)
            aa = abs(a)
            if (aa > eps) then
               if (n.ne.0 .and. aa>alfz) then
                  nszt = nszt + 1
                  mzeta(nszt) = m
                  nzeta(nszt) = n
               else if (m.ne.0 .or. n.ne.0) then
                  nsth = nsth + 1
                  mtheta(nsth) = m
                  ntheta(nsth) = n
               endif
c
c       Keep track of the number of modes (mpolb, ntorb); also set
c       up the array of n's (nvalsb) as used in the old dkes input
c       file format (DAS).
c
               mpolb = m + 1
               if (ntorb .eq. 0) then
                  ntorb = 1
                  nvalsb(1) = n
               else
                  do i = 1, ntorb
                     if (nvalsb(i) .eq. n) cycle  l12
                  end do
                  ntorb = ntorb + 1
                  nvalsb(ntorb) = n
               endif
            endif
         end do l12
      end do
c
      call sorti (nvalsb, ntorb)
c
c     estimate the Fourier spectrum for the distribution function
c     -----   original SPECFIL version   -----
c
c     The nplus and nmins arrays are used to indicate where m,n
c     modes are used for the distribution function. nplus holds the
c     n .ge. 0 modes while nmins holds the n .lt. 0 modes. These
c     arrays contain the current 3 interates of S through the third
c     index.
c
      nord = max(2,nord)

      nplus(:mhi,:nhi,:) = 0
      nmins(:mhi,:nhi,:) = 0

      nplus(1,1,1) = 1
c
c     Generate S(0) = S(theta) + S(zeta)
c
      do i = 1, nsth
         m = mtheta(i)
         n = ntheta(i)
         call filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      end do
      if (ier .ne.0) stop

      do i = 1, nszt
         m = mzeta(i)
         n = nzeta(i)
         call filit (m, n, 1, nplus, nmins, mhi, nhi, ier)
      end do
      if (ier .ne.0) stop

      nplus(:mhi,:nhi,2) = nplus(:mhi,:nhi,1)
      nmins(:mhi,:nhi,2) = nmins(:mhi,:nhi,1)
c
c     Generate S(1) = S(0)*S(zeta)
c

      do i = 1, nszt
         do in = 1, nhi
            do im = 1, mhi
               if (nplus(im,in,1) .eq. 1) then
                  m = mzeta(i) + im - 1
                  n = nzeta(i) + in - 1
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
                  m = mzeta(i) - im + 1
                  n = nzeta(i) - in + 1
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
               endif
               if (nmins(im,in,1) .eq. 1) then
                  m = mzeta(i) + im - 1
                  n = nzeta(i) - in
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
                  m = mzeta(i) - im + 1
                  n = nzeta(i) + in
                  call filit (m, n, 2, nplus, nmins, mhi, nhi, ier)
                  if (ier .ne.0) stop
               endif
            end do
         end do
      end do

      kk = 2
      if (nord .ne.2) then
c
c
c     Generate S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2):
c
c     Note: For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-2) recurrence
c      relation (this orders eps_tor << eps_hel) set ii = 0 and jj = 1 initially.
c      For the S(n) = S(zeta)*S(n-1) + S(theta)*S(n-1) recurrence
c      relation set ii = 1 and jj = 0 initially.
c
         ii = 1
         jj = 0
         do k = 3, nord
            ii = ii + 1
            if (ii .eq. 4) ii = 1
            jj = jj + 1
            if (jj .eq. 4) jj = 1
            kk = kk + 1
            if (kk .eq. 4) kk = 1
            nplus(:mhi,:nhi,kk) = nplus(:mhi,:nhi,jj)
            nmins(:mhi,:nhi,kk) = nmins(:mhi,:nhi,jj)
c
c        Generate the S(zeta)S(n-1) convolution:
c
            do i = 1, nszt
               do in = 1, nhi
                  do im = 1, mhi
                     if (nplus(im,in,jj) .eq. 1) then
                        m = mzeta(i) + im - 1
                        n = nzeta(i) + in - 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) - in + 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                     if (nmins(im,in,jj) .eq. 1) then
                        m = mzeta(i) + im - 1
                        n = nzeta(i) - in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mzeta(i) - im + 1
                        n = nzeta(i) + in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                  end do
               end do
            end do
c
c        Generate the S(theta)S(n-2) convolution:
c

            do i = 1, nsth
               do in = 1, nhi
                  do im = 1, mhi
                     if (nplus(im,in,ii) .eq. 1) then
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in - 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in + 1
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                     if (nmins(im,in,ii) .eq. 1) then
                        m = mtheta(i) + im - 1
                        n = ntheta(i) + in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                        m = mtheta(i) - im + 1
                        n = ntheta(i) - in
                        call filit(m,n,kk,nplus,nmins,mhi,nhi,ier)
                        if (ier .ne.0) go to 72
                     endif
                  end do
               end do
            end do
         end do
      endif

      ntor = 0
      mpol = 0
      ii = 0
      do n = 1, nhi
         jj = 0
         do m = 1, mhi
            if (nplus(m,n,kk) .eq. 1) then
               jj = jj + 1
               ivec(jj) = m - 1
            endif
         end do
         if (jj .ne.0) then
            ii = ii + jj
            ntor = ntor + 1
            mpol = max(mpol,jj)
c         Capture the n's and m's into arrays n_out and mns_out
c           for later printout:
            n_out(n+nhi) = n-1
            do iout=1,jj
              mns_out(n+nhi,iout) = ivec(iout)
            end do
c
            ij = min(13,jj)
            if (ntor <= 9 .and. ntor >= -9) then
               write (iunit, 101) ntor, n - 1, jj, (ivec(i),i=1,ij)
            else
               write (iunit, 201) ntor, n - 1, jj, (ivec(i),i=1,ij)
            endif
            if (jj > 13) then
               ik = -4
   61          continue
               ik = ik + 18
               ij = min(jj,ik + 17)
               write (iunit, 202) (ivec(i),i=ik,ij)
               if (ij < jj) go to 61
            endif
c            endif
         endif
      end do

      do n = 1, nhi
         jj = 0
         do m = 1, mhi
            if (nmins(m,n,kk) .eq. 1) then
               jj = jj + 1
               ivec(jj) = m - 1
            endif
         end do
         if (jj .ne.0) then
            ii = ii + jj
            ntor = ntor + 1
            mpol = max(mpol,jj)
c         Capture the n's and m's into arrays n_out and mns_out
c           for later printout:
            n_out(nhi-n+1) = -n
            do iout=1,jj
              mns_out(nhi-n+1,iout) = ivec(iout)
            end do
c
            ij = min(13,jj)
            if (ntor <= 9 .and. ntor >= -9) then
               write (iunit, 101) ntor, (-n), jj, (ivec(i),i=1,ij)
            else
               write (iunit, 201) ntor, (-n), jj, (ivec(i),i=1,ij)
            endif
            if (jj > 13) then
               ik = -4
   64          continue
               ik = ik + 18
               ij = min(jj,ik + 17)
               write (iunit, 202) (ivec(i),i=ik,ij)
               if (ij < jj) go to 64
            endif
         endif
      end do
c
      write (iunit, 203) mpol, ntor
  101 format(' mmnn(1,',i1,')=',15(i3,','))
  201 format(' mmnn(1,',i2,')=',15(i3,','))
  202 format(1x,18(i3,','))
  203 format(' mpol=',i2,', ntor=',i3,' /')

      go to 73
  72  write(*,'("error in call to filit")')
  73  continue
      close(unit= iunit)
c
       call second0(time_end)
       if (lscreen) then
          write (*,120) time_end - time_begin
          write (*,48)
       end if
c

      deallocate(bmn_local, abs_blocal, blocal_ordered,
     >   stat=istat)
      deallocate(m_ordered, n_ordered, stat=istat)
      deallocate(nsurf, stat=istat)

      call read_boozer_deallocate
      call read_wout_deallocate

      end subroutine dkes_input_prepare


      subroutine filit(m, n, k, nplus, nmins, mhi, nhi, ier)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, k, mhi, nhi, ier
      integer, dimension(mhi,nhi,3) :: nplus, nmins
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mp1, np1
C-----------------------------------------------

c-----------------------------------------------------------------------
c     auxilary routine for program GIDKES2
c     Author:     W.I. van Rij
c     modified:   HNM   STOP replaced by error parameter IER and RETURN
c-----------------------------------------------------------------------


      ier = 0
      if (m < 0) then
         m = -m
         n = -n
      endif

      mp1 = m + 1
      if (mp1 > mhi) then
         write (6, '(a,i3)') ' error in FILIT:   m+1 exceeds mhi = ',
     1      mhi
         ier = 1
         return
      endif

      if (n<0 .and. m.eq.0) n = -n

      if (n >= 0) then
         np1 = n + 1
         if (np1 > nhi) then
            write (6, '(a,i3)') ' error in FILIT:   n+1 exceeds nhi = '
     1         , nhi
            ier = 2
            return
         endif
         nplus(mp1,np1,k) = 1
         return

      else

         if ((-n) > nhi) then
            write (6, '(a,i3)') ' error in FILIT:   -n  exceeds nhi = '
     1         , nhi
            ier = 2
            return
         endif
         nmins(mp1,(-n),k) = 1
         return

      endif

      end subroutine filit


      SUBROUTINE SORTI(IX, NM)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NM
      INTEGER, DIMENSION(*) :: IX
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ML, MG, IIL, IIG, ML1, MG1, NL, NG, N, II, IX0
C-----------------------------------------------

Ccom*    Sorting of an Integer Array (size)
C----------------------------------------------------------------------
C     Sorting of an integer array. On exit, the elements of IX increase
C     with index i.
C
C     IX      integer array
C     NM      no. of elements
C----------------------------------------------------------------------
Cend
C
      IF (NM <= 1) RETURN
      ML = 1
      MG = NM
      IF (NM .ne.2) THEN
C
    1    CONTINUE
         IIL = IX(ML)
         IIG = IIL
         ML1 = ML + 1
         MG1 = MG - 1
         NL = ML
         NG = ML
C
         DO N = ML1, MG
            II = IX(N)
            IF (II < IIL) THEN
               NL = N
               IIL = II
            ENDIF
            IF (II > IIG) THEN
               NG = N
               IIG = II
            ENDIF
         END DO
C
         IF (NL .ne.ML) THEN
            IX0 = IX(NL)
            IX(NL) = IX(ML)
            IX(ML) = IX0
         ENDIF
         IF (NG .ne.MG) THEN
            IF (NG .eq. ML) NG = NL
            IX0 = IX(NG)
            IX(NG) = IX(MG)
            IX(MG) = IX0
         ENDIF
         ML = ML1
         MG = MG1
         IF (MG - ML >= 2) GO TO 1
C
         IF (MG <= ML) RETURN
      ENDIF
      IF (IX(MG) >= IX(ML)) RETURN
      IX0 = IX(MG)
      IX(MG) = IX(ML)
      IX(ML) = IX0
      RETURN
      END SUBROUTINE SORTI
