#!/bin/sh
cat > temp.c << "EOC"
cat > vmodules.f << "EOF"
       module BJdata
       use kind_spec
       implicit none
       integer, parameter :: max_bmns = 120
       real(rprec), dimension(max_bmns) :: xm_rdcd, xn_rdcd
       integer, dimension(max_bmns) :: m_rdcd, n_rdcd
       real(rprec), dimension(max_bmns) :: blocal_rdcd
       real(rprec) theta0, phi0, iota_local, ep_mu
       real(rprec) TWOPI,PI
       real(rprec) length_factor
       integer istop, J_star_opt
       character*3 device
       end module BJdata
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
       module B_and_J_Library
       use BJdata
       
       contains
       
       subroutine b_eval(theta,zeta,bfield)
       implicit none
       real(rprec), parameter :: zero = 0._dp
       real(rprec) :: theta, zeta, bfield
       real(rprec), dimension(max_bmns) :: temp
       temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
       bfield = sum(temp)

       end subroutine b_eval
       
       function b_along_fld_line(arg)  result(bb)
       implicit none
       real(rprec) :: zeta, theta, bb, arg
       real(rprec), dimension(max_bmns) :: temp
       select case (device)
        case ("qos")
         zeta = arg
         theta = theta0 + iota_local*zeta
        case ("qas")
         theta = arg
         zeta = theta/iota_local + phi0
        case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
        end select
        temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
        bb = sum(temp)
c        call b_eval(theta,zeta,bb)
       end function b_along_fld_line
c
       function invb_along_fld_line(arg)  result(bi)
       implicit none
       real(rprec) :: zeta, theta, bb, bi, arg
       real(rprec), dimension(max_bmns) :: temp
       select case (device)
        case ("qos")
         zeta = arg
         theta = theta0 + iota_local*zeta
        case ("qas")
         theta = arg
         zeta = theta/iota_local + phi0
        case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
        end select
        temp(:) = blocal_rdcd(:)*cos(theta*xm_rdcd(:)
     >   - zeta*xn_rdcd(:))
        bb = sum(temp)
c        call b_eval(theta,zeta,bb)
        bi = 1./bb
       end function invb_along_fld_line
c

       function b_along_const_theta(zeta)  result(bb)
       implicit none
       real(rprec) :: zeta, theta, bb
         theta = theta0
        call b_eval(theta,zeta,bb)
       end function b_along_const_theta
c

       function invb_along_const_theta(zeta)  result(bi)
       implicit none
       real(rprec) :: zeta, theta, bb, bi
        theta = theta0
        call b_eval(theta,zeta,bb)
        bi = 1./bb
       end function invb_along_const_theta
c       
       end module B_and_J_Library
EOF

cat > jinvar.f << "EOF"
      program J_invariant
c ******************************************************************************
c  J_INVARIANT
c
c  originally written by D. Spong, ORNL
c  modified by M. Zarnstorff, PPPL, Aug. 2002 to allow explicit
c       specification of the range of pitch values, for use by chisq_jconf
c
c  modified by M. Zarnstorff, PPPL, Oct. 2002 to directly implement the JCONF
c       calculation in a single call (as opposed to a succession of calls,
c       as done previously)
c
c
c  J_Invariant can now be invoked in three different ways, depending on the 
c  number of command line arguments
c
c  1) the standard case, as in original code, calculates on a single surface 
c     (index js), automatically choose the ep/mu values
c  xj_invariant  ext js nep_mu numJstar lscreen
c
c  2) similar to (1), but specify fixed range of ep/mu values
c  xj_invariant  ext js nep_mu numJstar lscreen epmu_min epmu_max
c
c  3) calculate fractions of J-inv range for each pitch on surface js that 
c     intersects the surface ks
c  xj_invariant  ext js nep_mu numJstar lscreen ks
c
c
c ******************************************************************************
      use read_boozer_mod
      use B_and_J_Library
      use date_and_computer, only: months
      use safe_open_mod
      implicit none
      external rhs
      external jac
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: j_invar = 108
      character*(50), parameter ::
     1   banner = ' THIS IS THE J-Invariant CODE Version 1.0.2z'
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      real(rprec), parameter :: pm4 = 1.e-4_dp, zero = 0._dp,
     >   one = 1._dp      
      integer js, i, j, k, ierr, istat, nfp, iloc, numargs, iunit,
     >        iunit2, is, ks
      real(rprec) :: btheta, bzeta, check, psip, chip, an
      real(rprec) ::phi_edge, phi_center, phi_lower, phi_upper
      integer :: neqn, itol, istate, itask, jt, istep, iopt
      integer :: liw, lrw, n_upper, jrad, istep_max
      real(rprec) :: thet, phi, bf, bf1, time_begin, 
     >  time_end,time_all
      real(rprec), dimension(2) :: y, f
      real(rprec) :: delta_phi, rtol, atol, phi_in, phi_out, 
     >  maxbmin, minbmax, minbmin, maxbmax, theta_min
      real(rprec) :: J_avg, J_min, J_max, width_epmu
      real(rprec) :: min_epmu, max_epmu, avg
      integer index_dat, index_end, num_ep_mu, NumJstar, i_ep_mu
      character*120 arg1, booz_input_file, output_file, output2_file
      character*120 arg2, arg3, arg4, arg5, arg6, arg7
      character*(10) :: date0, time0, zone0
      character*(40) :: dateloc
      real(rprec), dimension(:), allocatable :: epmu, 
     1        jinv_min, jinv_max, jsrc_min, jsrc_max
      real(rprec), dimension(:,:), allocatable :: J_inv_all
      integer :: imon
      logical :: lscreen, fix_pitch, ljconf

      TWOPI = 8*atan(1._dp)
      ep_mu = 1._dp
      device = "qas"
      J_star_opt = 0
      if (device .eq. "qas") J_star_opt = 0
      PI = TWOPI/2
      lscreen = .true.
      fix_pitch = .false.
      ljconf = .false.
c
c     Read data from the boozmn file and allocate storage:
c
      call second0(time_begin)
      call getcarg(1, arg1, numargs)
      if( numargs .eq. 7) then
         fix_pitch = .true.
      else if( numargs .eq. 6) then
         ljconf = .true.
      else if( numargs .ne. 5) then
       write(*,'("Error: 5 - 7 command line arguments are required")')
       stop 20
      endif

      call getcarg(2, arg2, numargs)
      call getcarg(3, arg3, numargs)
      call getcarg(4, arg4, numargs)
      call getcarg(5, arg5, numargs)
c      
      read(arg2,'(i20)') js
      read(arg3,'(i20)') num_ep_mu
      read(arg4,'(i20)') NumJstar
      if (arg5(1:1).eq.'f' .or. arg5(1:1).eq.'F') lscreen = .false.

      if( ljconf) then
         call getcarg(6, arg6, numargs)
         read(arg6,*) ks

      else if( fix_pitch) then
         call getcarg(6, arg6, numargs)
         call getcarg(7, arg7, numargs)

         read(arg6,*) min_epmu
         read(arg7,*) max_epmu
      endif

      allocate (epmu(num_ep_mu), jinv_min(num_ep_mu), 
     1   jinv_max(num_ep_mu),
     2   jsrc_min(num_ep_mu), jsrc_max(num_ep_mu) )
      allocate (J_inv_all(NumJstar, num_ep_mu))

      if (lscreen) then
          write(*,48)
          call date_and_time(date0,time0,zone0)
          read (date0(5:6),'(i2)') imon
          write (dateloc,100) months(imon),date0(7:8),date0(1:4),
     1      time0(1:2),time0(3:4),time0(5:6)
          write (*,'(1x,a,/,2x,a)') banner, dateloc 
          write(*,*) ' '
          write(*,110)
       end if
 100   format('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)
 110   format('  js   ep/mu    J_min    J_avg    J_max')
 120   format(/,'TIME IN J-Invariant CODE:',1pe12.2,' SEC')
  48   format('====================================================')
c      
c      
      index_dat = index(arg1,'.')
      index_end = len_trim(arg1)
      booz_input_file  = arg1(index_dat+1:index_end)
      output_file = "j_invar_out."//booz_input_file
      output2_file = "j_invar_sum."//booz_input_file

      call read_boozer_file (booz_input_file, k)
      if (k .ne. 0) stop 'Error reading boozmn file in J_INVARIANT'

      iunit = j_invar
      iunit2 = iunit + 1

      call safe_open(iunit2, istat, trim(output2_file), 'unknown', 
     1    'formatted')

      if( .not. ljconf) then

        call safe_open(iunit, istat, trim(output_file), 'unknown', 
     1    'formatted')

        call j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1     max_epmu, J_inv_all, epmu)

        do i=1, num_ep_mu

          J_min = minval(J_inv_all(:,i))
          J_max = maxval(J_inv_all(:,i))
          J_avg = sum(J_inv_all(:,i))/real(NumJstar)

          if(lscreen) then
            write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          endif

          write(iunit,*) (J_inv_all(j,i), j=1,NumJstar)
          write(iunit2,*) ep_mu, J_min, J_avg, J_max
        enddo

        close(unit=iunit)

      else  ! ljconf
!     Get the source surface
        fix_pitch = .false.   ! make sure...
        call j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1     max_epmu, J_inv_all, epmu)

        do i=1, num_ep_mu

          J_min = minval(J_inv_all(:,i))
          J_max = maxval(J_inv_all(:,i))
          J_avg = sum(J_inv_all(:,i))/real(NumJstar)

          jinv_min(i) = j_min
          jinv_max(i) = j_max

          if(lscreen) then
            write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          endif
        enddo

        fix_pitch = .true.
        min_epmu = epmu(1)
        max_epmu = epmu(num_ep_mu)
        jsrc_min = jinv_min
        jsrc_max = jinv_max

        if(lscreen) then
          write(*,*)
          write(*,110)
        endif

!      Now work out to the outer surface
c        do is = js+1, ks
         is = ks
          call j_inv_calc(is, num_ep_mu, NumJstar, fix_pitch, min_epmu,
     1      max_epmu, J_inv_all, epmu)

          do i=1, num_ep_mu

            J_min = minval(J_inv_all(:,i))
            J_max = maxval(J_inv_all(:,i))
            J_avg = sum(J_inv_all(:,i))/real(NumJstar)

            jinv_min(i) = max(jinv_min(i), j_min)
            jinv_max(i) = min(jinv_max(i), j_max)

            if(lscreen .and. is == ks)
     1       write(*,'(1x,i3,4(2x,f7.4))') js,epmu(i),J_min,J_avg,J_max
          enddo
c        enddo

        jsrc_min = max(jinv_max-jinv_min, 0._dp) / (jsrc_max-jsrc_min)
        avg = sum(jsrc_min)/num_ep_mu

        if( lscreen) then
          print *,' '

          print *,' J-contour confinement, nu=',NumJstar
          print *,'   E/mu      fraction lost'

          do i=1, num_ep_mu
            print '(2x,f6.3,4x,f6.3)', epmu(i), jsrc_min(i)
          enddo

          print '(1x,a7,4x,f6.3)', 'Average', avg
        endif

        do i=1, num_ep_mu
          write(iunit2,*) jsrc_min(i)
        enddo

        write(iunit2,*) ' '
        write(iunit2,*) avg
      endif   ! ljconf

      close(unit=iunit2)


      call second0(time_end)
      if (lscreen) then
          write (*,120) time_end - time_begin
          write (*,48)
      end if   

      deallocate( epmu, jinv_min, jinv_max, jsrc_min, jsrc_max, 
     1            J_inv_all)
      call read_boozer_deallocate

      end program J_invariant



      subroutine j_inv_calc( js, num_ep_mu, NumJstar, fix_pitch, 
     1     min_epmu, max_epmu, J_inv_all, epmu)
c ******************************************************************************
      use read_boozer_mod
      use B_and_J_Library
      implicit none
      external rhs
      external jac
      real(rprec), external :: sfmin
!-----------------------------------------------
!   Arguments
!-----------------------------------------------
      integer :: js, num_ep_mu, NumJstar
      logical :: fix_pitch
      real(rprec) :: min_epmu, max_epmu
      real(rprec), dimension(num_ep_mu) :: epmu
      real(rprec), dimension(NumJstar, num_ep_mu) :: J_inv_all

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), parameter :: pm4 = 1.e-4_dp, zero = 0._dp,
     >   one = 1._dp      
      real(rprec), dimension(mnboz_b) :: bmn_local, 
     >  abs_blocal, blocal_ordered
      integer, dimension(mnboz_b) :: m_ordered, n_ordered
      integer i, j, k, ierr, istat, nfp, iloc, numargs, iunit,
     >        iunit2, ks
      real(rprec) :: btheta, bzeta, psip, chip, an
      real(rprec) :: phi_edge, phi_center, phi_lower, phi_upper
      real(rprec), dimension(1) :: check
      real(rprec), dimension(NumJstar) :: phimin, phimax, 
     >  bmin, bmax
      real(rprec), dimension(NumJstar) :: theta_J, thmn, thmx
      real(rprec), dimension(NumJstar) :: J_plus, J_minus, J_inv
      real(rprec), dimension(NumJstar) :: L_plus, L_minus, L
      integer :: neqn, itol, istate, itask, jt, istep, iopt
      integer :: liw, lrw, n_upper, jrad, istep_max
      real(rprec) :: thet, phi, bf, bf1, time_begin, 
     >  time_end,time_all
      real(rprec), dimension(2) :: y, f
      real(rprec), dimension(55) :: rwork
      integer, dimension(22) :: iwork
      integer, dimension(ns_b,NumJstar) :: int_steps
      integer, dimension(NumJstar) :: int_steps_plus, int_steps_minus
      real(rprec) :: delta_phi, rtol, atol, phi_in, phi_out, 
     >  maxbmin, minbmax, minbmin, maxbmax, theta_min
      real(rprec) :: theta_center, theta_lower, theta_upper, 
     >  bmn_test
      real(rprec) :: J_avg, J_min, J_max, width_epmu
      integer index_dat, index_end, i_ep_mu
      integer :: imon
      logical :: lscreen, ljconf
 
      ep_mu = 1._dp
 
c
c   Get local (at flux surface js) value of quantities
c     needed for DKES run
c
      bmn_test = zero
      do i = 1,mnboz_b
         bmn_local(i) = bmn_b(i,js)
         abs_blocal(i) = abs(bmn_local(i))
         bmn_test = bmn_test + abs_blocal(i)
      end do
      if (bmn_test .eq. zero) then
        write(*,'("   Requested surface has not been mapped")')
        write(*,'("(i.e., sum of |Bmns| at this surface is zero)")')
        stop 30
      endif
      iota_local = iota_b(js)
      length_factor = sqrt(one + iota_local*iota_local)
      btheta = buco_b(js)
      bzeta = bvco_b(js)
      phi_edge = abs(phi_b(ns_b))
      psip = phi_edge/TWOPI
      chip = psip*iota_local
      nfp = nfp_b
      an = nfp

c
c     Sort the Bmn's (along with associated m's and n's in order
c     of increasing abs(Bmn):
c
      do i=1,mnboz_b

c      check = maxval(abs_blocal)
c    Find location of this value of check in the abs_blocal array
c        do j=1,mnboz_b
c         if(abs(check - abs_blocal(j)) .lt. 1.d-6*abs_blocal(j))
c    >      iloc = j
c        end do

       check = maxloc(abs_blocal)
       iloc = check(1)
       blocal_ordered(i) = bmn_local(iloc)
       m_ordered(i) = ixm_b(iloc)
       n_ordered(i) = ixn_b(iloc)
       abs_blocal(iloc) = zero
      end do
c
c    Keep max_bmns components for use in module B_and_J_Library
c
      do i=1, max_bmns
        m_rdcd(i) = m_ordered(i)
        xm_rdcd(i) = m_ordered(i)
        n_rdcd(i) = n_ordered(i)
        xn_rdcd(i) = n_ordered(i)
        blocal_rdcd(i) = blocal_ordered(i)
c        write(*,'(i5,2x,i5,2x,e15.7)') m_rdcd(i), n_rdcd(i),
c     >     blocal_rdcd(i)
      end do
c
c    Find bmin and bmax's along adjacent field lines:
c
      do i=1,NumJstar
c
      select case (device)
c
      case("qos")
       if(J_star_opt .eq. 0) then
         theta0 = PI*(one - (iota_local/an))*
     >      (i-1)/real(NumJstar-1,rprec)
         phi_center = theta0/(an - iota_local)
c
       else if(J_star_opt .ne. 0) then
         theta0 = PI*(i-1)/real(NumJstar-1,rprec)
         phi_center = theta0/an
       endif
       if(i .eq. 1) then
        phi_lower = phi_center - (PI/(5*an))
        phi_upper = phi_center + (PI/(5*an))
       else if(i .gt. 1) then
        phi_lower = phi_center - (PI/(2*an))
        phi_upper = phi_center + (PI/(2*an))
       endif
       if(J_star_opt .eq. 0) then
        phimin(i) = sfmin(phi_lower,phi_upper,
     >      b_along_fld_line,pm4)
        bmin(i) = b_along_fld_line(phimin(i))
c
       else if(J_star_opt .ne. 0) then
        phimin(i) = sfmin(phi_lower,phi_upper,
     >      b_along_const_theta,pm4)
        bmin(i) = b_along_const_theta(phimin(i))
       endif
c
       phi_center = phi_center + (PI/an)
       phi_lower = phi_lower + (PI/an)
       phi_upper = phi_upper + (PI/an)
       if(J_star_opt .eq. 0) then
        phimax(i) = sfmin(phi_lower,phi_upper,
     >       invb_along_fld_line,pm4)
        bmax(i) = b_along_fld_line(phimax(i))
        thmn(i) = theta0 + iota_local*phimin(i)
        thmx(i) = theta0 + iota_local*phimax(i)
c
       else if(J_star_opt .ne. 0) then
        phimax(i) = sfmin(phi_lower,phi_upper,
     >      invb_along_const_theta,pm4)
        bmax(i) = b_along_const_theta(phimax(i))
        thmn(i) = theta0
        thmx(i) = theta0        
       endif
c
       case("qas")
c
        phi0 = TWOPI*(i-1)/(an*(NumJstar-1))
        theta_center = zero
        theta_lower = theta_center - (PI/4)
        theta_upper = theta_center + (PI/4)

        thmn(i) = sfmin(theta_lower,theta_upper,
     >      b_along_fld_line,pm4)
        bmin(i) = b_along_fld_line(thmn(i))
        phimin(i) = thmn(i)/iota_local + phi0
c
        theta_center = PI
        theta_lower = theta_center - (PI/4)
        theta_upper = theta_center + (PI/4)
        thmx(i) = sfmin(theta_lower,theta_upper,
     >      invb_along_fld_line,pm4)
        bmax(i) = b_along_fld_line(thmx(i))
c
       case default
         write(*,'("Need to select either qo or qa device")')
         stop 23
       end select
c
      end do         ! do i=1,NumJstar
c
      maxbmin = maxval(bmin)
      minbmax = minval(bmax)
      minbmin = minval(bmin)
      maxbmax = maxval(bmax)
      width_epmu = minbmax - maxbmin

      if( .not. fix_pitch) then
         if( num_ep_mu <= 10) then
            min_epmu = maxbmin + 0.05_dp * width_epmu
            width_epmu = 0.9_dp * width_epmu
         else
            min_epmu = maxbmin + width_epmu/num_ep_mu/2
            width_epmu = width_epmu*(num_ep_mu-1)/num_ep_mu
         endif
      else
         width_epmu = max_epmu - min_epmu
      endif
c
      do i_ep_mu = 1,num_ep_mu
       if(num_ep_mu .ne. 1) then
c        ep_mu = 1.05_dp*maxbmin + ((0.95_dp*minbmax - 1.05_dp*maxbmin)
c     >    *(i_ep_mu - 1)/real(num_ep_mu - 1,rprec))
c        ep_mu = maxbmin + 0.05_dp*width_epmu + (0.9_dp*width_epmu
c    >    *(i_ep_mu - 1)/real(num_ep_mu - 1,rprec))
         ep_mu = min_epmu + width_epmu
     >    * (i_ep_mu - 1)/real(num_ep_mu - 1,rprec)

       else if(num_ep_mu .eq. 1) then
         ep_mu = min_epmu + 0.5_dp * width_epmu
       endif

c      write(*,'("ep/mu = ",f8.4)') ep_mu
c
      do i=1,NumJstar
       J_plus(i)=zero; L_plus(i)=zero; J_inv(i)=zero
       J_minus(i)=zero; L_minus(i)=zero; L(i) = zero
      end do
c       
      do i=1,NumJstar
c
      if (ep_mu .ge. maxbmax) cycle     !Passing orbit
      if (ep_mu .le. minbmin) cycle     !Forbidden (i.e., v||**2 < 0) orbit
c
      select case (device)
c
      case ("qos")
       if(J_star_opt .eq. 0) then
        theta0 = PI*(one - (iota_local/an))*(i-1)
     >               /real(NumJstar-1,rprec)
        theta_min = theta0 + iota_local*phimin(i)
       else if(J_star_opt .ne. 0) then
        theta0 = PI*(i-1)/real(NumJstar-1,rprec)
       endif
        delta_phi = TWOPI/(10*an)
        istep_max = 10
      case ("qas")
       phi0 = TWOPI*(i-1)/(an*(NumJstar-1))
       delta_phi = TWOPI/(10*an)
       istep_max = int(20/iota_local)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
       y(1)=zero; y(2)=zero; istop=0; neqn=2; itol=1; iopt=0
       istate=1; itask=1; jt=10; rtol=1.d-6; atol=1.d-7
       lrw = 55; liw = 22
       rwork(5) = delta_phi
       do istep = 1,istep_max
c
      select case (device)
c
      case ("qos")
        phi_in = phimin(i) + delta_phi*(istep - 1)
        phi_out = phimin(i) + delta_phi*istep
      case ("qas")
        phi_in = phi0 + delta_phi*(istep - 1)
        phi_out = phi0 + delta_phi*istep
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
        call lsode(rhs,neqn,y,phi_in,phi_out,itol,rtol,
     >  atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
        if(istop .eq. 1) exit
       end do
       int_steps_plus(i) = istep
       if(istep .lt. istep_max) J_plus(i) = abs(y(1))
       if(istep .ge. istep_max) J_plus(i) = 0
       L_plus(i) = abs(y(2))
c
      select case (device)
c
      case ("qos")
       delta_phi = TWOPI/(10*an)
       istep_max = 10
      case ("qas")
       delta_phi = TWOPI/(10*an)
       istep_max = int(20/iota_local)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
       y(1)=zero; y(2)=zero; istop=0; neqn=2; itol=1;iopt=0
       istate=1; itask=1; jt=10; rtol=1.d-6; atol=1.d-7
       lrw = 55; liw = 22
       rwork(5) = delta_phi
       do istep = 1,istep_max
c
      select case (device)
      case ("qos")
        phi_in = phimin(i) - delta_phi*(istep - 1)
        phi_out = phimin(i) - delta_phi*(istep)
      case ("qas")
        phi_in = phi0 - delta_phi*(istep - 1)
        phi_out = phi0 - delta_phi*(istep)
      case default
       write(*,'("Need to select either qo or qa device")')
       stop 23
      end select
c
        call lsode(rhs,neqn,y,phi_in,phi_out,itol,rtol,
     >  atol,itask,istate,iopt,rwork,lrw,iwork,liw,jac,jt)
        if(istop .eq. 1) exit
       end do
       int_steps_minus(i) = istep
       if(istep .lt. istep_max) J_minus(i) = abs(y(1))
       if(istep .ge. istep_max) J_minus(i) = zero
       L_minus(i) = abs(y(2))
       J_inv(i) = J_plus(i) + J_minus(i)
       L(i) = L_plus(i) + L_minus(i)
       int_steps(js,i) = int_steps_plus(i) + int_steps_minus(i)
c
      end do                                         !do i=1,NumJstar
 
      j_inv_all(:,i_ep_mu) = j_inv(:)
      epmu(i_ep_mu) = ep_mu
      end do                                !do i_ep_mu = 1,num_ep_mu

      end subroutine j_inv_calc

 
      function sfmin(ax, bx, f, tol)
c ******************************************************************************
      use kind_spec
      implicit none
      real(rprec), parameter :: zero = 0._dp, one = 1._dp
      real(rprec) :: ax, bx, f, tol, sfmin
c
c      an approximation  x  to the point where  f  attains a minimum  on
c  the interval  (ax,bx)  is determined.
c
c  input..
c
c  ax    left endpoint of initial interval
c  bx    right endpoint of initial interval
c  f     function subprogram which evaluates  f(x)  for any  x
c        in the interval  (ax,bx)
c  tol   desired length of the interval of uncertainty of the final
c        result (.ge.0.)
c
c  output..
c
c  fmin  abcissa approximating the point where  f  attains a
c        minimum
c
c      the method used is a combination of  golden  section  search  and
c  successive parabolic interpolation.  convergence is never much slower
c  than  that  for  a  fibonacci search.  if  f  has a continuous second
c  derivative which is positive at the minimum (which is not  at  ax  or
c  bx),  then  convergence  is  superlinear, and usually of the order of
c  about  1.324....
c      the function  f  is never evaluated at two points closer together
c  than  eps*abs(sfmin)+(tol/3), where eps is  approximately  the  square
c  root  of  the  relative  machine  precision.   if   f   is a unimodal
c  function and the computed values of   f   are  always  unimodal  when
c  separated  by  at least  eps*abs(x)+(tol/3), then  sfmin  approximates
c  the abcissa of the global minimum of  f  on the interval  ax,bx  with
c  an error less than  3*eps*abs(sfmin)+tol.  if   f   is  not  unimodal,
c  then sfmin may approximate a local, but perhaps non-global, minimum to
c  the same accuracy.
c      this function subprogram is a slightly modified  version  of  the
c  algol  60 procedure  localmin  given in richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
c
      real(rprec) :: a,b,c,d,e,eps,xm,p,q,r,tol1,t2,u,v,w,fu,fv,fw,
     2    fx,x,tol3
c
c  c is the squared inverse of the golden ratio
      c=0.5_dp*(3.0_dp- sqrt(5.0_dp))
c
c  eps is approximately the square root of the relative machine
c  precision.
c
c   10 eps=r1mach(4)
   10 eps=epsilon(eps)    !Use f90 routine to get machine eps rather than r1mach
      tol1=eps+one
      eps= sqrt(eps)
c
      a=ax
      b=bx
      v=a+c*(b-a)
      w=v
      x=v
      e=zero
      fx=f(x)
      fv=fx
      fw=fx
      tol3=tol/3
c
c  main loop starts here
c
   20 xm=0.5_dp*(a+b)
      tol1=eps* abs(x)+tol3
      t2=2*tol1
c
c  check stopping criterion
c
      if ( abs(x-xm).le.(t2-0.5_dp*(b-a))) go to 190
      p=zero
      q=zero
      r=zero
      if ( abs(e).le.tol1) go to 50
c
c  fit parabola
c
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2*(q-r)
      if (q .le. zero) go to 30
      p=-p
      go to 40
   30 q=-q
   40 r=e
      e=d
   50 if (( abs(p).ge. abs(0.5_dp*q*r)).or.(p.le.q*(a-x))
     2          .or.(p.ge.q*(b-x))) go to 60
c
c  a parabolic-interpolation step
c
      d=p/q
      u=x+d
c
c  f must not be evaluated too close to ax or bx
c
      if (((u-a).ge.t2).and.((b-u).ge.t2)) go to 90
      d=tol1
      if (x.ge.xm) d=-d
      go to 90
c
c  a golden-section step
c
   60 if (x.ge.xm) go to 70
      e=b-x
      go to 80
   70 e=a-x
   80 d=c*e
c
c  f must not be evaluated too close to x
c
   90 if ( abs(d).lt.tol1) go to 100
      u=x+d
      go to 120
  100 if (d .le. zero) go to 110
      u=x+tol1
      go to 120
  110 u=x-tol1
  120 fu=f(u)
c
c  update  a, b, v, w, and x
c
      if (fx.gt.fu) go to 140
      if (u.ge.x) go to 130
      a=u
      go to 140
  130 b=u
  140 if (fu.gt.fx) go to 170
      if (u.ge.x) go to 150
      b=x
      go to 160
  150 a=x
  160 v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
      go to 20
  170 if ((fu.gt.fw).and.(w.ne.x)) go to 180
      v=w
      fv=fw
      w=u
      fw=fu
      go to 20
  180 if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) go to 20
      v=u
      fv=fu
      go to 20
c
c  end of main loop
c
  190 sfmin=x

      end function sfmin


      subroutine rhs(neq, phi_loc, y, f)
      use B_and_J_Library
      integer neq
      real(rprec), parameter :: zero = 0._dp, one = 1._dp
      real(rprec), dimension(neq) :: y, f
      real(rprec) phi_loc, theta_loc, bf
      
      select case (device)
       case ("qos")
        if(J_star_opt .ne. 0) call b_eval(theta0,phi_loc,bf)
        if(J_star_opt .eq. 0) bf = b_along_fld_line(phi_loc)
       case ("qas")
        theta_loc = (phi_loc - phi0)*iota_local
        bf = b_along_fld_line(theta_loc)
      end select
c
      if (ep_mu .le. bf .or. istop .eq. 1) then
        f(1) = zero; f(2) = zero; istop = 1
      else if (ep_mu .gt. bf) then
        f(1) = length_factor*sqrt(one - (bf/ep_mu))
        f(2) = length_factor
        istop = 0
      end if

!     end subroutine rhs
      end
      

      subroutine jac(neq,t,y,ml,mu,pd,nrowd)
      use kind_spec
c
c     dummy jacobian subroutine (not used) for LSODE
c
      real(rprec), dimension(neq) :: y(neq)
      real(rprec), dimension(nrowd,neq) :: pd
      real(rprec) t
      integer neq, ml, mu, nrowd

      end subroutine jac
EOF
EOC
