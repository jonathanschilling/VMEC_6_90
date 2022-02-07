      program diagno
!***********************************************************************
!***********************************************************************
!! The program calculates coil responses using results from the vmec code
!! @author H.Gardener, Joachim Geiger
!! @version 1.0
!***********************************************************************
!***********************************************************************
!
!    /* when called by vmec this calculates the response of
!   the flux loops and poloidal beta coils on wviias */
!
      use constants
      use diagno_input
      use coilsystem
      use surface_values

      implicit none
      integer numargs
      character(150) :: arg1


!     character(10) :: machine_string = 'W7-AS'
!     character(10) :: machine_string = 'W7-X'


      id_string = ''
      CALL GETCARG(1, arg1, numargs)                        ! MCZ
      if (numargs .gt. 0) then
         id_string = trim(arg1)                             ! MCZ

      else
         write(6,*)'Usage: xdiagno data_extension'
         write(6,*)'Data set identification string :'
         read(5,*) id_string
      endif

! read control information

      call read_control_input

! read vmec/nemec output

      call read_input

!           /* read in external coil coordinates */
      call readcoil(machine_string)

      if( ltrace_progress) write(6,*)'Calling initialize_surface_values!'

      call initialize_surface_values

      if( ltrace_progress) write(6,*)'Before diagnostics!'

      if (idiag(4)==1 .or. idiag(5)==1 .or. idiag(7)==1 .or.      &
     &    idiag(9)==1) call bthdiag

      if (idiag(6)==1) call fluxdiag

!     if (idiag(3)==1) call barrow

      end program diagno


      subroutine becoil (x ,y ,z ,bx,by,bz)
!***********************************************************************
!! The subroutine returns the cartesian components of the magnetic field
!! at a point given its cartesian coordinates. The calcuation is done
!! using Biot-Savart or other supplied methods.
!***********************************************************************

      use constants
      use coilsystem

      implicit none

      real(rprec),intent(in)  :: x ,y ,z   !cartesian coordinates
      real(rprec),intent(out) :: bx,by,bz  !cartesian components of B

      real(rprec) :: vertbx,vertby,vertbz,brhos,bzets,brho,bzet,rho
      integer :: i

      bx=0;by=0;bz=0
      call biot_savart(x,y,z,bx,by,bz)
      if(ivertb.eq.0)then
        vertbx = 0
        vertby = 0
        vertbz = deltbz
      elseif (ivertb.eq.1)then
        brhos = 0
        bzets = 0
        do i=1,nc_vert
          rho = sqrt(x**2+y**2)
          call bvert(brho,bzet,rho,z,r_v(i),z_v(i))
          brhos = brhos + brho*f_v(i)
          bzets = bzets + bzet*f_v(i)
        enddo
        vertbx = current(4)*brhos*x/sqrt(x**2+y**2)
        vertby = current(4)*brhos*y/sqrt(x**2+y**2)
        vertbz = current(4)*bzets
      else
        write(6,*)'ivertb = ',ivertb,'! value undefined!'
        stop
      endif

      bx = bx+vertbx
      by = by+vertby
      bz = bz+vertbz

      contains

      subroutine biot_savart(x,y,z,bx,by,bz)
!! Evaluates Biot-Savart.

      use constants
      implicit none

      real(rprec),intent(in)  :: x,y,z     !cartesian coordinates
      real(rprec),intent(out) :: bx,by,bz  !cartesian components of B
      real(rprec), dimension(ntopo) :: x1,y1,z1,rw,fa
      real(rprec) :: ax,ay,az
      integer :: i_alloc
      real(rprec), dimension(:), allocatable, save :: vx,vy,vz,dx,dy,dz
      real(rprec), save :: fac
      integer, save :: i_call = 0

      if(i_call == 0) then
         fac   =  1.e-7_rprec        !if currents in A then B is in T.
         allocate(vx(nall),vy(nall),vz(nall),dx(nall),dy(nall),dz(nall),&
     &           stat = i_alloc)
         if(i_alloc /= 0) stop 'Allocation of arrays biot_savart failed!'

         dx=0;dy=0;dz=0;vx=0;vy=0;vz=0

         dx(1:nall-1) = (xw(2:nall)-xw(1:nall-1))*curre(1:nall-1)
         dy(1:nall-1) = (yw(2:nall)-yw(1:nall-1))*curre(1:nall-1)
         dz(1:nall-1) = (zw(2:nall)-zw(1:nall-1))*curre(1:nall-1)
         vx(1:nall-1) = &
     &         yw(1:nall-1)*dz(1:nall-1)-zw(1:nall-1)*dy(1:nall-1)
         vy(1:nall-1) = &
     &         zw(1:nall-1)*dx(1:nall-1)-xw(1:nall-1)*dz(1:nall-1)
         vz(1:nall-1) = &
     &         xw(1:nall-1)*dy(1:nall-1)-yw(1:nall-1)*dx(1:nall-1)
      endif

      x1(1:nall)  = x - xw(1:nall)
      y1(1:nall)  = y - yw(1:nall)
      z1(1:nall)  = z - zw(1:nall)
      rw(1:nall)  = &
     &       sqrt(x1(1:nall)**2+y1(1:nall)**2+z1(1:nall)**2)

      fa(1:nall-1)  = (rw(2:nall)+rw(1:nall-1))/ &
     &           (rw(2:nall)*rw(1:nall-1)*( &
     &             rw(2:nall)*rw(1:nall-1) +x1(2:nall)*x1(1:nall-1) &
     &            +y1(2:nall)*y1(1:nall-1) +z1(2:nall)*z1(1:nall-1)))

      ax = sum(fa(1:nall-1)*dx(1:nall-1))
      ay = sum(fa(1:nall-1)*dy(1:nall-1))
      az = sum(fa(1:nall-1)*dz(1:nall-1))

      bx = fac*(sum(fa(1:nall-1)*vx(1:nall-1)) -y*az+z*ay)
      by = fac*(sum(fa(1:nall-1)*vy(1:nall-1)) -z*ax+x*az)
      bz = fac*(sum(fa(1:nall-1)*vz(1:nall-1)) -x*ay+y*ax)

!     if(i_call == 0) write(6,*) &
!    &  '      x           y           z     ' &
!    & ,'     bx          by          bz     '
      if(i_call < 10) then
!       write(6,'(1p,6e12.4)') x,y,z,bx ,by ,bz
        i_call=i_call+1
      endif

      end subroutine biot_savart

      subroutine bvert(brho,bzet,rho,zp,ra,z0)
!***********************************************************************
! calculation of vertical field by elliptical integrals
! using the polynomial approximations given in abramowitz pp.591 ,
! which has an error of less than 2*e-08 .
! the returned values are those of a field in tesla generated
! by a current of magnitude muo/4*pi in a.
!***********************************************************************

      use constants
      implicit none
      real(rprec), intent(in) :: rho,zp,ra,z0
      real(rprec), intent(out) :: brho,bzet
      real(rprec) :: arz,amrz,xm,are,ark
      integer ifail

      ifail = 0
      arz = (rho+ra)**2 +(zp-z0)**2
      amrz = (rho-ra)**2 +(zp-z0)**2
      xm = 4*ra*rho/arz

      are = eie(xm)
      ark = eik(xm)

      brho = (zp-z0)*2.e-7_rprec/rho/sqrt(arz)  !gives T for a current of 1 A
      brho = brho*(-ark+ &
     &       (rho**2+ra**2+(zp-z0)**2)*are/amrz)

      bzet =  2.e-7_rprec/sqrt(arz)             !gives T for a current of 1 A
      bzet = bzet*(ark+ &
     &       (ra**2-rho**2-(zp-z0)**2)*are/amrz)

      return
      end subroutine bvert

      real(rprec) function eik(x)
!***********************************************************************
!! Approxiamtion of complete Elliptical Integral K
!! Reference: Abramowitz
!***********************************************************************

      use constants
      implicit none
      real(rprec), intent(in) :: x
      real(rprec) :: omx

      real(rprec), dimension(5) :: &
     &     a=(/1.38629436112_rprec, .09666344259_rprec, .03590092383_rprec, &
     &         .03742563713_rprec, .01451196212_rprec/), &
     &     b=(/ .50000000000_rprec, .12498593597_rprec, .06880248576_rprec, &
     &         .03328355346_rprec , .00441787012_rprec/)

      if(x.ge.1) then
        print * ,'eik : x.ge.1 :',x
        stop
      endif
      if(x.lt.0) then
        print * ,'eik : x.lt.1 :',x
        stop
      endif

      omx = 1 - x

      eik = b(1)+omx*(b(2)+omx*(b(3)+omx*(b(4)+omx*b(5))))
      eik = eik*log(1/omx)
      eik = eik +a(1)+omx*(a(2)+omx*(a(3)+omx*(a(4)+omx*a(5))))

      return
      end function eik

      real(rprec) function eie(x)
!***********************************************************************
!! Approxiamtion of complete Elliptical Integral E
!! Reference: Abramowitz
!***********************************************************************

      use constants
      implicit none
      real(rprec), intent(in) :: x
      real(rprec) :: omx

      real(rprec), dimension(5) :: &
     &     a=(/1.00000000000_rprec, .44325141463_rprec, .06260601220_rprec, &
     &         .04757383546_rprec, .01736506451_rprec/), &
     &     b=(/ .0000000000_rprec, .24998368310_rprec, .09200180037_rprec, &
     &         .04069697526_rprec , .00526449639_rprec/)

      if(x.ge.1) then
        print * ,'eie : x.ge.1 :',x
        stop
      endif
      if(x.lt.0) then
        print * ,'eie : x.lt.1 :',x
        stop
      endif

      omx = 1 - x

      eie = b(1)+omx*(b(2)+omx*(b(3)+omx*(b(4)+omx*b(5))))
      eie = eie*log(1/omx)
      eie = eie +a(1)+omx*(a(2)+omx*(a(3)+omx*(a(4)+omx*a(5))))

      return
      end function eie

      end subroutine becoil

      subroutine belicu (xp,yp,zp,bx,by,bz)
!! Returns the cartesian components of the magnetic field at a
!! given point specified by its cartesian coordinates. The field
!! originates from a line current specified by subroutin tolicu.

      use constants
      use diagno_input
      implicit none

      real(rprec), intent(in)  :: xp,yp,zp
      real(rprec), intent(out) :: bx,by,bz
      real(rprec), dimension(:), allocatable, save :: &
     &     xcw,ycw,zcw,dx,dy,dz,vx,vy,vz
      integer :: i_alloc
      integer, save :: ibelicu = 0 
      real(rprec), dimension(nvp+1) :: x1,y1,z1,rw,fa,r12
      real(rprec) :: ax,ay,az
      real(rprec), save :: fac

      if(ibelicu == 0) then
        fac = 1.e-7_rprec        !if current is in A, B will be in T
        call tolicu
        allocate(dx(nvp+1),dy(nvp+1),dz(nvp+1) &
     &          ,vx(nvp+1),vy(nvp+1),vz(nvp+1) &
     &          ,stat = i_alloc)
        if(i_alloc /= 0) stop 'Allocation in tolicu failed!'

        dx=0;dy=0;dz=0;vx=0;vy=0;vz=0

        dx(1:nvp) = (xcw(2:nvp+1)-xcw(1:nvp))
        dy(1:nvp) = (ycw(2:nvp+1)-ycw(1:nvp))
        dz(1:nvp) = (zcw(2:nvp+1)-zcw(1:nvp))
        vx(1:nvp) =  ycw(1:nvp)*dz(1:nvp)-zcw(1:nvp)*dy(1:nvp)
        vy(1:nvp) =  zcw(1:nvp)*dx(1:nvp)-xcw(1:nvp)*dz(1:nvp)
        vz(1:nvp) =  xcw(1:nvp)*dy(1:nvp)-ycw(1:nvp)*dx(1:nvp)
        ibelicu = 1
      endif
      x1  = xp - xcw
      y1  = yp - ycw
      z1  = zp - zcw
      rw  = sqrt(x1**2+y1**2+z1**2)

      r12(1:nvp) = rw(2:nvp+1)*rw(1:nvp)
      fa(1:nvp) = (rw(2:nvp+1)+rw(1:nvp))/ &
     &   (r12(1:nvp)*(r12(1:nvp) + &
     &           x1(2:nvp+1)*x1(1:nvp)+y1(2:nvp+1)*y1(1:nvp) &
     &          +z1(2:nvp+1)*z1(1:nvp)))

      ax = sum(fa(1:nvp)*dx(1:nvp))
      ay = sum(fa(1:nvp)*dy(1:nvp))
      az = sum(fa(1:nvp)*dz(1:nvp))

      bx = fac*(sum(fa(1:nvp)*vx(1:nvp)) -yp*az+zp*ay)
      by = fac*(sum(fa(1:nvp)*vy(1:nvp)) -zp*ax+xp*az)
      bz = fac*(sum(fa(1:nvp)*vz(1:nvp)) -xp*ay+yp*ax)

!----------------------------------------------------------
      contains

      subroutine tolicu
!! specifies the geometry of the line current carrying the net plasma
!! current.

      use surface_values

      implicit none

      integer :: nuh,i,l,i_alloc

      allocate(xcw(nvp+1),ycw(nvp+1),zcw(nvp+1) &
     &      ,stat = i_alloc)
      if(i_alloc /= 0) stop 'Allocation in tolicu failed!'

      select case (input_form)
      case ("vmec2000")
        xcw = xax
        ycw = yax
        zcw = zax
      case default
        nuh    = nu/2
        i      = 0
        do l=1,nuvp,nu
          i      = i + 1
          xcw(i)  = (xg(l)+xg(l+nuh))/2
          ycw(i)  = (yg(l)+yg(l+nuh))/2
          zcw(i)  = (zg(l)+zg(l+nuh))/2
!         write(19,'(1p,9e14.5)') xcw(i),ycw(i),zcw(i) &
!    &      ,xg(l),xg(l+nuh),yg(l),yg(l+nuh),zg(l),zg(l+nuh)
        enddo
        xcw(nvp+1) = xcw(1)
        ycw(nvp+1) = ycw(1)
        zcw(nvp+1) = zcw(1)
      end select

      end subroutine tolicu

      end subroutine belicu

      subroutine bpoly(xp,yp,zp,bx,by,bz)
!***********************************************************************
!! returns the cartesian components of the magnetic field at a
!! point specified by its cartesian coordinates by performing
!! the surface integrals and adding the field due to a net current
!! (calculated by subroutine belicu).
!***********************************************************************


      use constants
      use surface_values
      use diagno_input

      implicit none

      real(rprec), intent(in)  :: xp,yp,zp
      real(rprec), intent(out) :: bx,by,bz
      real(rprec) :: fac,bxt,byt,bzt
      real(rprec), allocatable, dimension(:) :: an1,an2,an3,  &
     &             temp,temp1
!      real(rprec), dimension(nuvp) :: an1,an2,an3,temp,temp1

      allocate (an1(nuvp), an2(nuvp), an3(nuvp), temp(nuvp), temp1(nuvp))

      fac    = 1/(2*twopi*nuv)
      bx     = 0
      by     = 0
      bz     = 0
!
!     !* integrate over all periods of the toroidal surface */
!
      an1 = 1/sqrt((xp-xg)**2+(yp-yg)**2+(zp-zg)**2)
      an2 = an1**2
      an3 = an1**3
      temp = (xp-xg)*snx+(yp-yg)*sny+(zp-zg)*snz
      temp1= bexn+3*an2*temp

      bx = sum(((xp-xg)*temp1-snx)*an3)
      by = sum(((yp-yg)*temp1-sny)*an3)
      bz = sum(((zp-zg)*temp1-snz)*an3)

      deallocate (an1, an2, an3, temp, temp1)

      bx = -bx*fac
      by = -by*fac
      bz = -bz*fac

      call belicu(xp,yp,zp,bxt,byt,bzt)

      bxt = tor_cur*bxt
      byt = tor_cur*byt
      bzt = tor_cur*bzt

      bx = bx + bxt
      by = by + byt
      bz = bz + bzt

      end subroutine bpoly

      subroutine bthdiag
!***********************************************************************
!! Calculates the magnetic field response due to the plasma based on the
!! precalculated surface functions. Performs different tasks at the
!! moment and is currently a spaghetti-like code.
!***********************************************************************

      use constants
      use diagno_input
      use surface_values

      implicit none

      integer :: i,j,k,ncoils,i_alloc,nseg,npts, istat
      real(rprec), dimension(:), allocatable :: x_c,y_c,z_c,th_inc &
     &       ,phi_inc,eff_area,b_c,b_n,flux
      real(rprec), dimension(:), allocatable :: r_c2th,z_c2th,phi_c2th &
     &      ,x_c2th,y_c2th,wa_c2th,flux_c2th,bint
      integer :: nseg_c2th,npts_c2th
      real(rprec) :: xp,yp,zp,bx,by,bz,br,bphi,dx,dy,dz,wa
      character(80) :: seg_coil_name
      character(5) :: action
!----- temporary -------------------------
      real(rprec) :: bl,dl
!-----------------------------------------

      if( ltrace_progress) write(6,*) ' Start of bthdiag!'
!
!    /* read in br, br and bphi field coils */
!
      if( idiag(4)==1 ) then
        if( ltrace_progress) write(6,*)'Read in dataset of field coils'
        open(unit=11,file=trim(bprobes_file),status='old' &
     &      ,action='read',iostat=istat)
        if (istat .eq. 0) then
          read(11,*)ncoils
          allocate(x_c(ncoils),y_c(ncoils),z_c(ncoils) &
     &       ,th_inc(ncoils),phi_inc(ncoils),eff_area(ncoils) &
     &       ,b_c(ncoils), b_n(ncoils), flux(ncoils) &
     &       , stat = i_alloc)
          if(i_alloc /= 0) stop 'Allocation 1 in bthdiag failed!'
          do j=1,ncoils
            read(11,*)x_c(j),y_c(j),z_c(j) &
     &            ,th_inc(j),phi_inc(j),eff_area(j)
          enddo
          th_inc = th_inc*onerad
          phi_inc = phi_inc*onerad
          close(unit=11)
          if( ltrace_progress) write(6,*)&
     &      '#nc       x/m           y/m           z/m    ' &
     &     ,'       B/T          flux/Wb '
          do j=1,ncoils
            xp = x_c(j)
            yp = y_c(j)
            zp = z_c(j)
            call bpoly(xp,yp,zp,bx,by,bz)
            br =(xp*bx+yp*by)/sqrt(xp*xp+yp*yp)
            bphi =(-yp*bx+xp*by)/sqrt(xp*xp+yp*yp)
            b_c(j) = br*cos(th_inc(j)) + bz*sin(th_inc(j))
            b_n(j) = br*sin(th_inc(j)) - bz*cos(th_inc(j))
            flux(j) = eff_area(j)*b_c(j)
            if( ltrace_progress) write(6,'(i3,1p,10e14.4)') &
     &          j,x_c(j),y_c(j),z_c(j) &
     &         ,b_c(j),flux(j)
          enddo
          open(unit=11,file='diagno_bth.'//id_string,status='unknown' &
     &        ,action='write')
          write(11,'(i6,1p,6e14.6)')(i,x_c(i),y_c(i),z_c(i) &
     &          ,b_c(i),b_n(i),flux(i),i=1,ncoils)
          close(unit=11)
          deallocate(x_c,y_c,z_c,th_inc,phi_inc,eff_area,b_c,flux &
     &             , stat = i_alloc)
          if(i_alloc /= 0) stop 'Deallocation 1 in bthdiag failed!'

        else
          stop 'Failed to open B-probes definition file'
        endif
      endif

      if( idiag(9) == 1) then
        if( ltrace_progress) write(6,*)'Read in dataset of field coils in MIR-1'
        open(unit=11,file=trim(mir1_file),status='old' &
     &      ,action='read',iostat=istat)
        if (istat .eq. 0) then
          read(11,*)ncoils
          allocate(x_c(ncoils),y_c(ncoils),z_c(ncoils) &
     &         ,th_inc(ncoils),phi_inc(ncoils),eff_area(ncoils) &
     &         ,b_c(ncoils), flux(ncoils) &
     &         , stat = i_alloc)
          if(i_alloc /= 0) stop 'Allocation 1 in bthdiag failed!'
          do j=1,ncoils
            read(11,*)x_c(j),y_c(j),z_c(j) &
     &            ,th_inc(j),phi_inc(j),eff_area(j)
          enddo
          th_inc = th_inc*onerad
          phi_inc = phi_inc*onerad
          close(unit=11)
          if( ltrace_progress) write(6,*)&
     &      '#nc       x/m           y/m           z/m    ' &
     &     ,'       B/T          flux/Wb '
          do j=1,ncoils
            xp = x_c(j)
            yp = y_c(j)
            zp = z_c(j)
            call bpoly(xp,yp,zp,bx,by,bz)
!       call becoil(xp,yp,zp,bx,by,bz)
!       call belicu(xp,yp,zp,bx,by,bz)       !is called in bpoly
            br =(xp*bx+yp*by)/sqrt(xp*xp+yp*yp)
            bphi =(-yp*bx+xp*by)/sqrt(xp*xp+yp*yp)
          b_c(j) = br*cos(th_inc(j)) + bz*sin(th_inc(j))
!       b_c(j) = br*cos(th_inc(j))*cos(phi_inc(j)) &
!    &         + bphi*cos(th_inc(j))*sin(phi_inc(j)) &
!    &         + bz*sin(th_inc(j))
            flux(j) = eff_area(j)*b_c(j)
            if( ltrace_progress) write(6,'(i3,1p,10e14.4)') &
     &          j,x_c(j),y_c(j),z_c(j) &
     &         ,b_c(j),flux(j)
          enddo
          open(unit=11,file='diagno_mir1.'//id_string,status='unknown' &
     &      ,action='write')
          write(11,'(1p,30e14.5)')(flux(i),i=1,ncoils)
          close(unit=11)
          deallocate(x_c,y_c,z_c,th_inc,phi_inc,eff_area,b_c,flux &
     &             , stat = i_alloc)
          if(i_alloc /= 0) stop 'Deallocation 1 in bthdiag failed!'
        else
          stop 'Failed to open Mir1 definition file'
        endif
      endif

      if( idiag(5)==1 ) then
        seg_coil_name = "cos2theta"
        action = "read"
        call seg_coil(seg_coil_name, action, istat)

        if( istat == 0 ) then
          nseg=nseg_c2th
          npts=npts_c2th
          if( ltrace_progress) write(6,*)nseg,' Segments for segmented coil!'
          if( ltrace_progress) write(6,*)npts,' Points for each coil segment!'
          allocate(bint(nseg), stat = i_alloc)
          if(i_alloc /= 0) stop 'Allocation of bint in bthdiag failed!'
          bint = 0
          do i=1,nseg
            k=(i-1)*npts+1
            do j=1,npts-1
              xp=(x_c2th(k+1)+x_c2th(k))/2
              yp=(y_c2th(k+1)+y_c2th(k))/2
              zp=(z_c2th(k+1)+z_c2th(k))/2
              wa=(wa_c2th(k+1)+wa_c2th(k))/2
              dx=    (x_c2th(k+1)-x_c2th(k))*wa
              dy=    (y_c2th(k+1)-y_c2th(k))*wa
              dz=    (z_c2th(k+1)-z_c2th(k))*wa
              call bpoly(xp,yp,zp,bx,by,bz)
!         call becoil(xp,yp,zp,bx,by,bz)
!         call belicu(xp,yp,zp,bx,by,bz)       !is called in bpoly
              dl = sqrt(dx**2+dy**2+dz**2)
              bl = (bx*dx + by*dy + bz*dz)
              bint(i) = bint(i)+bx*dx + by*dy + bz*dz
!!            write(6,'(a45,i5,1p,4e14.5)') &
!!     &           'k , vecb*vecdl/dl , |vecdl| , wa , vecb*vecdl' &
!!     &           ,k,bl/dl,dl/wa,wa,bl
              k=k+1
            enddo
          enddo

          if( ltrace_progress) then
             write(6,*)'Cos2theta Coil:'
             write(6,*)'#Flux in Wb of Coil  1 to ',nseg,' :'
             write(6,'(1p,30e14.5)')(bint(i),i=1,nseg)
             write(6,*)'Sum of fluxes : ',sum(bint),' Wb'
          endif

          if(idiag(51)==1) call write_Cos2th_database_w7as(nseg,bint)

          bint(:nseg) = bint(:nseg) * seg_rog_turns(:nseg)     ! scale appropriately

          open(unit=11,file='diagno_seg.'//id_string,status='unknown' &
     &        ,action='write')
          write(11,'(1p,30e14.5)')(bint(i),i=1,nseg)
          do i=1, nseg
             if( i < 10) then
                write(11,'(a,i1)') ' seg',i
             else if (i<100) then
                write(11,'(a,i2)') ' seg',i
             else 
                write(11,'(a,i3)') ' seg',i
             endif          
          enddo

          close(unit=11)
          call seg_coil(seg_coil_name,"close", istat)

          deallocate(bint, stat = i_alloc)
          if(i_alloc /= 0) write(6,*) 'Deallocation of bint in bthdiag failed!'

        else
          stop  'Failure reading segmented-rogowski definition file'
        endif
      endif

      if(idiag(7)==1)then
        if( ltrace_progress) write(6,*)'Read in dataset for testpoints'
        open(unit=11,file=trim(bfield_points_file),status='old' &
     &    ,action='read', iostat=istat)
        if( istat /= 0) stop 'Failure reading b-field point definnition file'

        read(11,*)ncoils
        if( ltrace_progress) then
           write(6,*)ncoils,' test points from file bfield_points.data'
           write(6,'(a3,6a14)') ' No','      x       ','      y       ' &
     &      ,'      z       ','      Br      ' &
     &      ,'      Bphi    ','      Bz      '
        endif
        do j=1,ncoils
          read(11,*)xp,yp,zp
          call bpoly(xp,yp,zp,bx,by,bz)
!         call becoil(xp,yp,zp,bx,by,bz)
!         call belicu(xp,yp,zp,bx,by,bz)       !is called in bpoly
          br =(xp*bx+yp*by)/sqrt(xp*xp+yp*yp)
          bphi =(-yp*bx+xp*by)/sqrt(xp*xp+yp*yp)
          if( ltrace_progress) write(6,'(i3,1p,6e14.5)') &
     &        j,xp,yp,zp,br,bphi,bz
        enddo
        close(unit=11)
      endif

      contains

      subroutine seg_coil(coil_name,action, istat)
!! Handle different segmented Rogowski coils.

      implicit none
      character(80), intent(in) :: coil_name
      character(5), intent(in) :: action
      integer istat

      select case(coil_name)
        case ("cos2theta")
          if( ltrace_progress) write(6,*)'Read data of Coil:',coil_name
          call cos2theta_coord(action, istat)

        case default
          write(6,*)'Coil name ',coil_name,' unknown!'
          stop 'Check coil_name in bthdiag!'
      end select

      end subroutine seg_coil

      subroutine cos2theta_coord(action, istat)
!! Read the coordinates of the Cosine-2-Theta Coil

      implicit none
      character(5), intent(in) :: action
      integer :: j,n, istat

      select case(action)
        case("read")
          open(unit=10,file=trim(seg_rog_file)         &
     &    ,status='unknown',action='read',iostat=istat)
          if( istat /= 0) return

          read(10,*)nseg_c2th,npts_c2th
          n = nseg_c2th*npts_c2th
          allocate(r_c2th(n),z_c2th(n) &
     &      ,phi_c2th(n),wa_c2th(n) &
     &      ,x_c2th(n),y_c2th(n) &
     &      ,flux_c2th(n) &
     &    ,stat=i_alloc)
          if(i_alloc /= 0) &
     &        stop 'Allocation failed in read_cos2theta_coord!'
          read(10,*)(r_c2th(j),z_c2th(j),phi_c2th(j) &
     &              ,wa_c2th(j),j=1,n)
          close(unit=10)
          x_c2th=r_c2th*cos(onerad*phi_c2th)
          y_c2th=r_c2th*sin(onerad*phi_c2th)
        case("close")
          deallocate(r_c2th,z_c2th,phi_c2th,wa_c2th &
     &              ,x_c2th,y_c2th,flux_c2th &
     &              ,stat=i_alloc)
          if(i_alloc /= 0) &
     &        write(6,*) 'Deallocation failed in read_cos2theta_coord!'
        case default
          write(6,*)'Action ',action,' unknown!'
          stop 'Check requested action in bthdiag!'
      end select

      end subroutine cos2theta_coord
!
! Internal subroutine
!
      subroutine write_Cos2th_database_w7as(nseg,bint)
!! Writes results of the Cosine-2-Theta Calculations to a database
!! file : Cos2th.database_new
      use constants
      implicit none
      integer,intent(in) :: nseg
      real(rprec), dimension(nseg) :: bint
      integer :: i
      real(rprec) :: cratio,c0t,c1t,c2t,s1t

      open(unit=11,file=trim(fdb_SegRogCoils),status='old' &
     &    ,action='write',position='append')
      
      write(11,'(60a)',advance='no') id_string(1:len_trim(id_string))
      write(11,'(1pe14.5)',advance='no') current(1)
      cratio = current(2)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(3)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(4)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = tor_cur/current(1)    !plasma current/I_mod
      write(11,'(1pe14.5)',advance='no') cratio
      write(11,'(1pe14.5)',advance='no') rav
      write(11,'(1p,6e14.5)',advance='no') &
     &               r00min,r00max,z00max,rmidmin,rmidmax,zmidmax
      write(11,'(1p,8e14.5)',advance='no') (bint(i),i=1,nseg)
      c0t=bint(1)+bint(2)+bint(3)+bint(4)
      c1t=bint(1)-bint(2)-bint(3)+bint(4)
      c2t=bint(1)-bint(2)+bint(3)-bint(4)
      s1t=bint(1)+bint(2)-bint(3)-bint(4)
      write(11,'(4(1pe14.5))',advance='yes') c0t,c1t,c2t,s1t
      close(unit=11)

      end subroutine write_Cos2th_database_w7as

      end subroutine bthdiag

      subroutine fluxdiag
!***********************************************************************
!* rewritten by J.Geiger, based on fluxdiag written by H.Gardener 1988**
!
!    /* calculates flux due to plasma through external diagnostic
!       coils. coil coordinates are read in from unit 10.
!          nb: loops may be closed or "open" (i.e. toroidal but
!              input over one period only). only the first end-point
!              is to be input so that the no. of points = no. of
!              segments = nseg.        h.j.g.      8.4.88   */
!
!***********************************************************************
      use constants
      use diagno_input
      use surface_values

      implicit none

      real(rprec) :: xp,yp,zp,cop,sip,sfac
!      real(rprec), dimension(nuvp) :: disinv
      real(rprec), dimension(:), allocatable :: disinv
      character(48), dimension(:), allocatable :: title
      integer, dimension(:), allocatable :: nseg_ar
      real(rprec), dimension(:), allocatable :: xfl,yfl,zfl,dx,dy,dz &
     &                                  ,xmid,ymid,zmid,flux 
      integer :: i,n,nfl,nseg,iflflg,idia,i_alloc,nsegmx

!     write(3,*)'phiedge = ',phiedge,'[Wb]'
      sfac = 1/(2*twopi*nuv)      !factor for surface integral

!
! Read number of flux loops and number of segments.
!
      open(unit=10,file=trim(flux_diag_file),status='old',action='read')

      read(10,*) nfl
      if( ltrace_progress) write(6,*) nfl,' flux loops!'

      allocate(nseg_ar(nfl), title(nfl+2), flux(nfl+2), stat=i_alloc)
      if(i_alloc /= 0) stop &
     &    'Allocation of int.array nseg_ar fluxdiag failed!'

      nseg_ar = 0
      do n = 1,nfl
!       read(10,*) &
        read(10,'(3i5,a48)') nseg_ar(n),iflflg,idia,title(n)
        if( ltrace_progress) write(6,*)nseg_ar(n),' Segments for flux loop ',n

        do i = 1,nseg_ar(n)
          read(10,*) xp,yp,zp
        enddo
      enddo

      nsegmx=maxval(nseg_ar)+1
      if( ltrace_progress) write(6,*)'Maximum number of segments+1 = ',nsegmx

      allocate(xfl(0:nsegmx),yfl(0:nsegmx),zfl(0:nsegmx), &
     &         dx(0:nsegmx),dy(0:nsegmx),dz(0:nsegmx), &
     &         xmid(nsegmx),ymid(nsegmx),zmid(nsegmx)  &
     &        ,stat=i_alloc)
      if(i_alloc /= 0) stop &
     &    'Allocation of pos.arrays xlf,yfl,zfl in fluxdiag failed!'
!
!
!
! Reread input unit to get flux loops and calculate fluxes
!
      rewind(unit=10)
      read(10,*) nfl

      sip = sin(twopi/float(np))
      cop = cos(twopi/float(np))
      idia = 0 ; iflflg = 0
      flux = 0

      do n = 1,nfl
        read(10,'(3i5,a48)') nseg,iflflg,idia,title(n)
        do i=1,nseg
          read(10,*) xfl(i),yfl(i),zfl(i)
        enddo
        xfl=xfl/100 ; yfl=yfl/100 ; zfl=zfl/100  ![cm] -> [m]
!
!        /* calculate zeroth and second end-point */
!
        if (iflflg.eq.0) then
          xfl(0) = xfl(nseg)
          yfl(0) = yfl(nseg)
          zfl(0) = zfl(nseg)
          xfl(nseg+1) = xfl(1)
          yfl(nseg+1) = yfl(1)
          zfl(nseg+1) = zfl(1)
        else
!
!           /* toroidal loop input over one period */
!
          xfl(0) = xfl(nseg)*cop + yfl(nseg)*sip
          yfl(0) = yfl(nseg)*cop - xfl(nseg)*sip
          zfl(0) = zfl(nseg)
          xfl(nseg+1) = xfl(1)*cop - yfl(1)*sip
          yfl(nseg+1) = yfl(1)*cop + xfl(1)*sip
          zfl(nseg+1) = zfl(1)
        endif
!
!        /* calculate distance of points */
!
        dx(1:nseg) = xfl(2:nseg+1) - xfl(1:nseg)
        dy(1:nseg) = yfl(2:nseg+1) - yfl(1:nseg)
        dz(1:nseg) = zfl(2:nseg+1) - zfl(1:nseg)
        if (iflflg.eq.0) then   !for poloidally closed loops
          dx(0) = dx(nseg)
          dy(0) = dy(nseg)
          dz(0) = dz(nseg)
          dx(nseg+1) = dx(1)
          dy(nseg+1) = dy(1)
          dz(nseg+1) = dz(1)
        else   !for toroidally closed loops
!
!           /* toroidal loop input over one period */
!
          dx(0) = dx(nseg)*cop + dy(nseg)*sip
          dy(0) = dy(nseg)*cop - dx(nseg)*sip
          dz(0) = dz(nseg)
          dx(nseg+1) = dx(1)*cop - dy(1)*sip
          dy(nseg+1) = dy(1)*cop + dx(1)*sip
          dz(nseg+1) = dz(1)
        endif
        xmid(1:nseg) = (xfl(1:nseg) + xfl(2:nseg+1))/2
        ymid(1:nseg) = (yfl(1:nseg) + yfl(2:nseg+1))/2
        zmid(1:nseg) = (zfl(1:nseg) + zfl(2:nseg+1))/2
!
!        /* perform line integrals */
!
        allocate (disinv(nuvp), stat = i_alloc)
        if( ltrace_progress) write(6,*)'Performing line integral using closed simpson formula ...'
        do i = 1,nseg
          xp = xfl(i)
          yp = yfl(i)
          zp = zfl(i)
          disinv(1:nuvp) = 1/sqrt((xp - xg(1:nuvp))**2 &
     &                            +(yp - yg(1:nuvp))**2 &
     &                            +(zp - zg(1:nuvp))**2)
          flux(n) = flux(n) &
     &       + (sum(atopx(1:nuvp)*disinv(1:nuvp))*(dx(i-1)+dx(i)) &
     &        + sum(atopy(1:nuvp)*disinv(1:nuvp))*(dy(i-1)+dy(i)) &
     &        + sum(atopz(1:nuvp)*disinv(1:nuvp))*(dz(i-1)+dz(i)) )/2
          xp = xmid(i)
          yp = ymid(i)
          zp = zmid(i)
          disinv(1:nuvp) = 1/sqrt((xp - xg(1:nuvp))**2 &
     &                            +(yp - yg(1:nuvp))**2 &
     &                            +(zp - zg(1:nuvp))**2)
          flux(n) = flux(n) &
     &          + 2*(sum(atopx(1:nuvp)*disinv(1:nuvp))*dx(i) &
     &             + sum(atopy(1:nuvp)*disinv(1:nuvp))*dy(i) &
     &             + sum(atopz(1:nuvp)*disinv(1:nuvp))*dz(i) )
        enddo

        deallocate (disinv)

        flux(n)=flux(n)/3
        flux(n) = sfac*flux(n)         !finish surface integral
        if (iflflg.eq.1) flux(n)=flux(n)*np

!       if (idia.eq.1) flux(n) = -flux(n) + phiedge !up to 21.2.1999
        if (idia.eq.1)then
          if( ltrace_progress) write(6,*)'flux(n),before = ',flux(n),'  phiedge = ',phiedge
          flux(n) = flux(n) - phiedge !changed sign:22.2.1999
        endif
        if( ltrace_progress) then
           write(6,'(" coil number ",i5,10x,a)') n, trim(title(n))
           write(6,'(a7,i4,8(a7,i2,a3,1p,e18.10,a3))') &
     &          ' nseg= ',nseg,'  flux(',n,')= ',flux(n),' Wb'
           write(6,'(a33,1p,e14.5)')'  Flux normalized to phiedge/2 : ' &
     &                             ,2*flux(n)/phiedge
           write(6,*) ' '
        endif
        write(3,'(" coil number ",i5,10x,6a8)') n, trim(title(n))
        write(3,'(a7,i4,8(a7,i2,a3,1p,e18.10,a3))') &
     &    ' nseg= ',nseg,'  flux(',n,')= ',flux(n),' Wb'
        write(3,'(a33,1p,e14.5)')'  Flux normalized to phiedge/2 : ' &
     &            ,2*flux(n)/phiedge
        write(3,*) ' '
      enddo
      close(unit=10)
!
! write database file, if requested
!
      if(idiag(61)==1) call write_flux_database_w7as

!
! scale fluxes by number of turns in coils
!

      flux(:nfl) = flux(:nfl) * flux_turns(:nfl)

!
! compensate diamagnetic, if requested
!
      if( idiag(60)==1 ) then
         flux(nfl+1) = sum( flux(:nfl) * dia_comp_coefs(:nfl) )
         title(nfl+1) = 'comp. diamag. loop'
         flux_comb_coefs(nfl+1) = 0
         nfl = nfl + 1
      endif
!
! make linear combination of flux loops, if requested
!
      if( idiag(62)==1 ) then
         flux(nfl+1) = sum( flux(:nfl) * flux_comb_coefs(:nfl) )
         title(nfl+1) = 'combined flux loops'
         nfl = nfl + 1
      endif

!
! write fluxdiag.result
!
      open(unit=11,file='diagno_flux.'//id_string,status='unknown' &
     &    ,action='write')
      write(11,'(1p,60e18.10)')(flux(n), n=1,nfl)

      do n=1, nfl
         write(11, '(a)') title(n)
      enddo

      close(unit=11)
!
! do deallocations
!
      deallocate(xfl,yfl,zfl,dx,dy,dz,xmid,ymid,zmid,flux, &
     &    nseg_ar, title, stat=i_alloc)
      if(i_alloc /= 0) stop &
     &    'Deallocation of arrays xlf etc. in fluxdiag failed!'


      contains
!
! Internal subroutine
!
      subroutine write_flux_database_w7as
!! Writes results of flux calculations into a database file,
!! currently: flux.database_highbeta2001
!! for W7-AS

      use constants
      implicit none
      real(rprec) :: cratio, olddia1, newdia1, sc_equi

      open(unit=11,file=trim(fdb_fluxloops),status='old' &
     &    ,action='write',position='append')
      
      write(11,'(60a)',advance='no') id_string(1:len_trim(id_string))
      write(11,'(1pe14.5)',advance='no') current(1)
      cratio = current(2)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(3)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(4)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = tor_cur/current(1)    !plasma current/I_mod
      write(11,'(1pe14.5)',advance='no') cratio
      write(11,'(1pe14.5)',advance='no') rav
      write(11,'(1p,6e14.5)',advance='no') &
     &               r00min,r00max,z00max,rmidmin,rmidmax,zmidmax
      olddia1 = 14*flux(1)+60*5*flux(2)
      newdia1 = 14*flux(4)+3.56*60*flux(5)+2.66*60*flux(6)
      write(11,'(1p,8e14.5)',advance='no') &
     &            14*flux(1),60*flux(2),60*flux(3),olddia1 &
     &           ,14*flux(4),60*flux(5),60*flux(6),newdia1
      write(11,'(1pe14.5)',advance='no') 9*flux(7)       !diamagn.loop 2
      sc_equi = flux(9)+flux(13)  ! excentric saddle loops
      write(11,'(1p,4e14.5)',advance='no') &
     &               flux(9),flux(13),sc_equi,flux(10)
      write(11,'(1p,2e14.5)',advance='yes') &
     &               flux(14),flux(15)
      close(unit=11)

      end subroutine write_flux_database_w7as

      subroutine write_flux_database_w7x
!! Writes results of flux calculations into a database file,
!! currently: flux.database_w7x
!! for W7-X

      use constants
      implicit none
      real(rprec) :: cratio

      open(unit=11,file=trim(fdb_fluxloops),status='old' &
     &    ,action='write',position='append')
      
      write(11,'(60a)',advance='no') id_string(1:len_trim(id_string))
      write(11,'(1pe14.5)',advance='no') current(1)
      cratio = current(2)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(3)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(4)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(5)/current(1)
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(6)/current(1)    !curA/I_1
      write(11,'(1pe14.5)',advance='no') cratio
      cratio = current(7)/current(1)    !curB/I_1
      write(11,'(1pe14.5)',advance='no') cratio
      write(11,'(1pe14.5)',advance='no') rav
      write(11,'(1p,6e14.5)',advance='no') &
     &               r00min,r00max,z00max,rmidmin,rmidmax,zmidmax
      write(11,'(1p,e14.5)',advance='no') flux
      close(unit=11)

      end subroutine write_flux_database_w7x

      end subroutine fluxdiag

