!----------------------------------------------------------------------
!     DRIVER ROUTINE FOR BNORM CODE
!----------------------------------------------------------------------
      program bfield_n 
      use kind_spec
      use safe_open_mod
      implicit none
      
      integer, parameter :: nu=64, nv=64, mf=10, nf=10, md=20, nd=20
      integer, parameter :: ibnorm=15
      integer :: iunit = ibnorm, m, n, istat
      real(rprec), dimension(0:mf,-nf:nf) :: bnfou
      character*(200) :: extension, separation
c----------------------------------------------------------------- 
c     bnormal reads WOUT file (from VMEC) and delivers the
c     Fourier harmonics of B dot n :  bnfou, writing these out to the bnorm file 
c----------------------------------------------------------------- 
!
!     CALL BNORM AS FOLLOWS FROM COMMAND LINE:   
!
!     xbnorm wout.extension coil_separation(optional)
!
!     For example:
!
!     xbnorm wout.extension 0.30     
!
!     corresponds to the using the wout.ipp file for vmec input, and
!     a coil separation of 0.30 (in the correct units, of course).
!
!     This codes produces two files needed by NESCOIL (and possibly other codes):
!
!     bnorm.extension           normal b-field coefficients
!     nescin.extension          input file read for NESCOIL
!
!
!     NOTE: nescin has been modified from original format: np added to front end
!      

!
!     READ IN ARGUMENT FROM COMMAND LINE
!
      call bn_read_commandline(extension, separation)

      call bnormal(nu, nv, mf, nf, md, nd, bnfou, extension)

!
!     WRITE OUT BNFOU (FOURIER COEFFICIENTS OF BNORMAL) TO BNORM.EXTENSION FILE
!
      call safe_open(iunit, istat, 'bnorm.' // extension, 'replace',
     1      'formatted')

      do m = 0, mf
         do n = -nf,nf
            write(iunit, 1100) m, n, bnfou(m,n)
         end do
      end do

 1100 format(1x,2i5,1pe24.16)
      
      close (iunit)
         

!
!     Write out input file read by NESCOIL (SPH: 9/98)
!
      call bn_write_nescoil_input(extension)

!     print *, ' Finished writing bnorm and nescin files'

      end program bfield_n


!----------------------------------------------------------------------
!     NESCOIL ROUTINES
!----------------------------------------------------------------------
      subroutine bn_read_commandline(extension, separation)
      use neswrite, only: coil_separation, dp
      implicit none
      character(len=*) :: extension, separation
      integer :: numargs, ia, istat
c ----------------------------------------------------------------------
c                                                            11.09.99
c     purpose:  read command line arguments (SPH)
c
c
c ----------------------------------------------------------------------

      call getcarg(1, extension, numargs)
      if (numargs .eq. 2) call getcarg(2, separation, numargs)
      coil_separation = 0
      if (numargs .lt.1 .or. extension .eq. '-h' .or. 
     1                       extension .eq. '/h') then
         print *,
     1   ' Syntax:  xbnorm wout.extension coil_separation(optional)'
         print *, ' '
         print *, 
     3   ' where the coil separation (in same units as R00) is optional'
         stop 'Invalid command line'
      else if (numargs .eq. 2) then
         read (separation, *, iostat=istat) coil_separation
         if (istat .ne. 0) stop 'BNORM Error reading coil_separation'
      endif
      
      ia = index(extension, 'wout.')
      if (ia .gt. 0) then
         extension = extension(ia+5:len_trim(extension))   
      end if   
      
      end subroutine bn_read_commandline
      

      subroutine bn_alloc
      use meshes
      use bnvariables
      implicit none
      integer ierr

      allocate ( bu(nuv),                  bv(nuv) 
     .          ,sinu(nuv,  0:md),conu(nuv,  0:md)
     .          ,sinv(nuv,-nd:nd),conv(nuv,-nd:nd)
     .          ,tanu(-nu:nu),      tanv(-nvp:nvp)
     .          ,indu(nuvp),            indv(nuvp)       
     .          ,    x(nuvp),   y(nuvp),   z(nuvp)
     .          ,  djx(nuvp), djy(nuvp), djz(nuvp)    
     .          ,    xu(nuv),   yu(nuv),   zu(nuv)

     .          ,    xv(nuv),   yv(nuv),   zv(nuv)         
     .          , guu(nuvh1),guv(nuvh1),gvv(nuvh1) 
     .          , dju(nuvh1),djv(nuvh1),sqf(nuvh1)          
     .          , au(nuvh1) , av(nuvh1),               stat=ierr)
     
      if (ierr .ne. 0) stop 'Allocation in BNORM subroutine bn_alloc'

      end subroutine bn_alloc


! ----------------------------------------------------------------------
      subroutine bn_bfield_parallel 
! ----------------------------------------------------------------------
c                                                            11.09.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n
      real(rprec) :: cop
c ----------------------------------------------------------------------
      curpol = alp*bsubv(0,0)
c     write(6,3000) curpol
 3000 format('   curpol= ',1pe12.4)
      do m=  0,md
      do n=-nd,nd
        bsubu(m,n) = pi2*bsubu(m,n)/curpol
        bsubv(m,n) = alp*bsubv(m,n)/curpol
      enddo
      enddo

      do  i = 1 , nuv                                               
         bu (i)   = 0._dp         
         bv (i)   = 0._dp         
      enddo

      do  m = 0,mb                                                   
       do  n = -nb,nb                                                   
        do  i  = 1,nuv                                               
         cop      = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)
         bu (i)   = bu (i) + bsubu(m,n)*cop 
         bv (i)   = bv (i) + bsubv(m,n)*cop 
        enddo
       enddo
      enddo 

      end subroutine   bn_bfield_parallel


!----------------------------------------------------------------------
      subroutine bnormal(n_u,n_v,mfou,nfou,mdim,ndim,bnfou,extension)
!----------------------------------------------------------------------
c                                                          15.08.98
c    purpose:
c

c    compuce
c                   / 
c               1   |   [ j', x-x'] . n 
c     B.n   =  ---- |  ----------------- df ,  x,x' on the boundary
c              4 pi |     |x-x'|^(3/2)
c                   /
c 
c          
c       j   = B_subu . xv - B_subv . xu       
c                             
c the result ist normalized so that the toroidal integral 
c                
c              / 
c              |   
c     I_p  =   |  B. ds  per field period   is I_p  =1.
c              |     
c              /
c ----------------------------------------------------------------------
      use meshes
      use kind_spec
      implicit none

      integer :: n_u, n_v, mfou, nfou, mdim, ndim
      real(rprec), dimension(0:mfou,-nfou:nfou) :: bnfou
      character*(*) :: extension
         nu     = n_u
         nv     = n_v
         nuv    = nu*nv
         nuvh1  = nu*nv/2+nu
         mf     = mfou
         nf     = nfou
         md     = mdim
         nd     = ndim
      call bn_read_vmecf90(extension)
      call bn_alloc
      call bn_precal
      call bn_bfield_parallel
      call bn_surface
      call bn_vecpot
      call bn_fouri(bnfou)

      end subroutine bnormal


      subroutine bn_fouri(bnfou)
c ----------------------------------------------------------------------
c     purpose:                                               06.04.00   
c                 
c  
c ---------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n, nlo
c ---------------------------------------------------------------------
      real(rprec), dimension(0:mf,-nf:nf)   :: bnfou
      real(rprec), dimension(0:md,-nd:nd)   :: aufou,avfou
      real(rprec), dimension(nuvh1)         :: ampl   
      real(rprec), dimension(nuv)           :: bn   
      real(rprec)                           :: faz
c ----------------------------------------------------------------------
      do i=1,nu
         ampl(i) = 1
         ampl(i+nuv/2) = 1
      enddo
      do i=nu+1,nuv/2
         ampl(i)   = 2
      enddo
      do  m  = 0,mf
         faz   = 2*fnuv
      if(m.eq.0) faz = fnuv 
       do  n  = -nf,nf
         aufou(m,n) = 0
         avfou(m,n) = 0
        do  i = 1,nuvh1
         aufou(m,n) = aufou(m,n) +au(i)*ampl(i)*
     .             (conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n))*faz
         avfou(m,n) = avfou(m,n) +av(i)*ampl(i)*
     .             (conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n))*faz
        enddo
       enddo
      enddo
        bn(:nuv) = 0
      do  m  = 0,mf
       do  n  = -nf,nf
        do  i = 1,nuv
        bn(i) = bn(i) +pi2*(m*avfou(m,n)-n*aufou(m,n))*
     .             (sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n))
        enddo
       enddo
      enddo
      do  m  = 0,mf
         faz   = 2*fnuv
       if(m.eq.0) faz = fnuv 
       do  n  = -nf,nf

          bnfou(m,n) = 0
        do  i = 1,nuvh1
         bnfou(m,n) = bnfou(m,n) +bn(i)/sqf(i)*ampl(i)*
     .             (sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n))*faz
        enddo
       enddo
      enddo
!     end  subroutine bn_fouri
      end

!----------------------------------------------------------------------
      subroutine bn_precal
! ----------------------------------------------------------------------
c                                                          11.09.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
      integer :: i, m, n, ku, kv
c ----------------------------------------------------------------------
         pi     = 2*asin(1._dp)
         pi2    = 2 * pi
         alp    = pi2 / (np)
         alu    = pi2 / (nu)
         alv    = pi2 / (nv)
         alvp   = pi2 / (nvp)
         fnuv   = 1._dp/(nuv)
c
      do  m=0,md
         i      = 0
       do  kv = 1,nv
        do  ku = 1,nu
         i      = i + 1
         conu(i,m) = cos(alu*m*(ku-1))
         sinu(i,m) = sin(alu*m*(ku-1))
        enddo
       enddo
      enddo
      do  n=-nd,nd 
         i      = 0
       do  kv = 1,nv
        do  ku = 1,nu
         i      = i + 1
         sinv(i,n) = sin(alv*n*(kv-1))
         conv(i,n) = cos(alv*n*(kv-1))
        enddo

       enddo
       do  i=1,nuv
         sinv(i,-n) = -sinv(i,n)
         conv(i,-n) =  conv(i,n)
       enddo
      enddo
      do kv=1,nvp
       do ku=1,nu
         i    = ku + nu*(kv-1)
         indv(i)  = kv
         indu(i)  = ku
       enddo
      enddo
      do  ku = -nu+1,nu-1
       if(iabs(ku).ne.nu/2) then
         tanu(ku) = tan(.5_dp*alu*ku)/pi
       else
         tanu(ku)=1.e+20_dp
       endif
      enddo
      do  kv = -nvp+1,nvp-1
       if(iabs(kv).ne.nvp/2) then
         tanv(kv) = tan(.5_dp*alvp*kv)/pi
        else
         tanv(kv) = 1.e+20_dp
        endif
      enddo

      end subroutine bn_precal


      subroutine bn_read_vmecf90(extension)
      use  meshes
      use  bnvariables 
      use neswrite, only: coil_separation, mnmax_in => mnmax, ixm, ixn,
     1   raxis_in => raxis, zaxis_in => zaxis, nfp_in => nfp,
     2   iota_edge, phip_edge
      use read_wout_mod                  !Error-free reading of wout file
      implicit none
c-----------------------------------------------
c   local variables
c
      integer :: ierr, iopen, m, n, mn, mpol1
      character*(*) :: extension

c-----------------------------------------------
      allocate (bsubu(0:md,-nd:nd), bsubv(0:md,-nd:nd), 
     1      cr(0:md,-nd:nd), cz(0:md,-nd:nd),cl(0:md,-nd:nd), stat=ierr)
      bsubu = 0;  bsubv = 0; cr = 0; cz = 0; cl = 0;
      if (ierr .ne. 0) stop 'Allocation error in bn_read_vmecf90'
c-----------------------------------------------
!
!     THIS MODULE SUBROUTINE LOADS UP ARRAYS, CONSTANTS READ IN FROM WOUT FILE
!
      call read_wout_file ('wout.' // extension, ierr, iopen)
      if (iopen .ne. 0) stop 'error opening wout in bn_read_vmecf90'
      if (ierr .ne. 0) stop 'error reading wout in bn_read_vmecf90'
       
      mpol1  = mpol-1

      if(mf.lt.mpol1 .or. nf.lt.ntor) then
         print *, 'increase number of poloidal and/or toroidal modes:',
     1            ' mf,nf' 
         stop
      endif
c---------------------------------------------------------------------
      allocate (ixm(mnmax), ixn(mnmax), raxis_in(0:ntor), 
     1          zaxis_in(0:ntor))

      raxis_in(0:ntor) = raxis(0:ntor,1)
      zaxis_in(0:ntor) = zaxis(0:ntor,1)
      mnmax_in = mnmax
      nfp_in = nfp

      do mn = 1, mnmax
         ixm(mn) = nint(xm(mn))
         ixn(mn) =-nint(xn(mn))/nfp   !!Flip sign: NESCOIL convention
         m = ixm(mn)
         n = ixn(mn)
         bsubu(m,n) = 1.5_dp*bsubumn(mn,ns) - 0.5_dp*bsubumn(mn,ns-1)
         bsubv(m,n) = 1.5_dp*bsubvmn(mn,ns) - 0.5_dp*bsubvmn(mn,ns-1)
         cl(m,n) = 1.5_dp*lmns(mn,ns) - 0.5_dp*lmns(mn,ns-1)
         cr(m,n) = rmnc(mn,ns)
         cz(m,n) = zmns(mn,ns)
      end do

      iota_edge = 1.5_dp*iotas(ns) - 0.5_dp*iotas(ns-1)
      phip_edge = 1.5_dp*phip (ns) - 0.5_dp*phip (ns-1)

      if (coil_separation .le. 0._dp) then
         coil_separation = abs(cr(1,0))
         print *,' A default coil-plasma separation was chosen: ',
     1   coil_separation
         print *,' You may enter this value as the 2nd arg ',
     1           'on the command line'
      end if

      np  = nfp
      nvp = nv*np
      nuvp = nu*nv*np
      mb  = mpol1
      nb  = ntor

!
!     Deallocate memory in READ_WOUT module
!
      call read_wout_deallocate

      end subroutine bn_read_vmecf90


! ----------------------------------------------------------------------
      subroutine regint(f,a,b,c,n)
! ----------------------------------------------------------------------
      use kind_spec
      implicit none
      integer :: i, n
      real(rprec), dimension(n) :: f, a, b, c
      real(rprec) :: sqp, sqm, sqa, sqc, top, tom
c ----------------------------------------------------------------------
      do  i=1,n
         sqp    = sqrt(a(i)+2._dp*b(i)+c(i))       
         sqm    = sqrt(a(i)-2._dp*b(i)+c(i))       
         sqa    = sqrt(a(i))    
         sqc    = sqrt(c(i))
         top    = log((sqc*sqp+c(i)+b(i))/(sqa*sqp-a(i)-b(i)))/sqp     
         tom    = log((sqc*sqm+c(i)-b(i))/(sqa*sqm-a(i)+b(i)))/sqm     
         f(i)   = top + tom
      enddo

      end subroutine regint
      

! ----------------------------------------------------------------------
      subroutine bn_surface
! ----------------------------------------------------------------------
c                                                            11.08.99
c     purpose:
c
c
c ----------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
c ----------------------------------------------------------------------
      integer :: i, m, n, ku, kv, k, np2
      real(rprec), dimension(:), allocatable :: r, ru, rv
      real(rprec) :: snx, sny, snz, coh, sih, co, si, cofp, sifp, 
     1    cofm, sifm, cm, cn
c ----------------------------------------------------------------------
      allocate (r(nuv), ru(nuv), rv(nuv), stat=i)
      
      do m=0,mb
        cr(m,0) = .5_dp*cr(m,0)
        cz(m,0) = .5_dp*cz(m,0)
      enddo
      do  i = 1 , nuv                                               
         r(i)    = 0._dp         
         z(i)    = 0._dp                                                    
         ru(i)   = 0._dp                                         
         rv(i)   = 0._dp                                                   
         zu(i)   = 0._dp                                                   
         zv(i)   = 0._dp                                                   
      enddo
      do  m = 0,mb                                                   
      do  n = 0,nb                                                   
         cm     = m*pi2                                                 
         cn     = n*pi2                                                 
      do  i  = 1 , nuv                                              
         cofp   = conu(i,m)*conv(i,n)-sinu(i,m)*sinv(i,n)
         cofm   = conu(i,m)*conv(i,n)+sinu(i,m)*sinv(i,n)
         sifp   = sinu(i,m)*conv(i,n)+conu(i,m)*sinv(i,n)
         sifm   = sinu(i,m)*conv(i,n)-conu(i,m)*sinv(i,n)
         r(i)   = r(i)    +       cr(m,n)*cofp + cr(m,-n)*cofm
         z(i)   = z(i)    +       cz(m,n)*sifp + cz(m,-n)*sifm
         ru(i)  = ru(i)   - cm *( cr(m,n)*sifp + cr(m,-n)*sifm )    
         rv(i)  = rv(i)   - cn *( cr(m,n)*sifp - cr(m,-n)*sifm )    
         zu(i)  = zu(i)   + cm *( cz(m,n)*cofp + cz(m,-n)*cofm )    
         zv(i)  = zv(i)   + cn *( cz(m,n)*cofp - cz(m,-n)*cofm )    
      enddo
      enddo
      enddo
      do m=0,mb
        cr(m,0) = 2._dp*cr(m,0)
        cz(m,0) = 2._dp*cz(m,0)
      enddo
c----------------------------------------------------------
      do  kv = 1, nv                                            
         coh    = cos(alvp*(kv-1))
         sih    = sin(alvp*(kv-1))
      do  ku = 1, nu                                                 
         i      = nu*(kv-1)+ku                                         
         x(i)   = coh * r(i)                                   
         y(i)   = sih * r(i)                                   
         xu(i)  = coh * ru(i)                                  
         yu(i)  = sih * ru(i)                                  
         xv(i)  = coh * rv(i) - alp*y(i)                   
         yv(i)  = sih * rv(i) + alp*x(i) 
      enddo
      enddo
      do   i = 1 , nuvh1
         snx    = yu(i)*zv(i)-zu(i)*yv(i)                          
         sny    = zu(i)*xv(i)-xu(i)*zv(i)                          
         snz    = xu(i)*yv(i)-yu(i)*xv(i)                          
         sqf(i) = sqrt(snx*snx+sny*sny+snz*snz)
         guu(i) = xu(i)*xu(i)+yu(i)*yu(i)+zu(i)*zu(i)               
         guv(i) = xu(i)*xv(i)+yu(i)*yv(i)+zu(i)*zv(i)               
         gvv(i) = xv(i)*xv(i)+yv(i)*yv(i)+zv(i)*zv(i)  
         dju(i) = bu(i)*guv(i)- bv(i)*guu(i)
         djv(i) = bu(i)*gvv(i)- bv(i)*guv(i) 
      enddo
      do i=1,nuv
         djx(i)  = bu(i)*xv(i)-bv(i)*xu(i)
         djy(i)  = bu(i)*yv(i)-bv(i)*yu(i)
         djz(i)  = bu(i)*zv(i)-bv(i)*zu(i)
      enddo

      do k=1,np-1
        co  =cos(alp*k)
        si  =sin(alp*k)
       do i=1,nuv
           x(i+k*nuv) = x(i)*co - y(i)*si
           y(i+k*nuv) = y(i)*co + x(i)*si
           z(i+k*nuv) = z(i)
         djx(i+k*nuv) = djx(i)*co-djy(i)*si
         djy(i+k*nuv) = djy(i)*co+djx(i)*si
         djz(i+k*nuv) = djz(i)
       enddo
      enddo
         np2 =  np*np
      do i=1,nuvh1
         guv(i) = np*guv(i)
         gvv(i) = np2*gvv(i)
      enddo

      deallocate (r, ru, rv, stat=i)

      end subroutine bn_surface
      

      subroutine bn_vecpot
c ----------------------------------------------------------------------
c     purpose:                                               11.08.99   

c
c
c ---------------------------------------------------------------------
      use meshes
      use bnvariables
      implicit none
c ----------------------------------------------------------------------
      integer :: i, ip
      real(rprec), dimension(:), allocatable :: ax, ay, az, analyt
      real(rprec) :: pi41, dintu, dintv      
c ----------------------------------------------------------------------
      allocate (ax(nuvh1), ay(nuvh1), az(nuvh1), analyt(nuvh1),stat=i)

      pi41  = .5_dp/pi2
      ax   = 0;      ay   = 0;      az   = 0 

      do i = 1,nuvh1
         call vecpot (ax(i), ay(i), az(i), i, 1, i-1)
         call vecpot (ax(i), ay(i), az(i), i, i+1, nuvp)
      enddo


!     FREE UP BIG ARRAYS
      deallocate (x, y, z, indu, indv, djx, djy, djz)
      
      call regint(analyt, guu, guv, gvv, nuvh1) ! compute  singular  integral                

      do  i  = 1,nuvh1
         dintu = (ax(i)*xu(i)+ay(i)*yu(i)+az(i)*zu(i))*fnuv
         dintv = (ax(i)*xv(i)+ay(i)*yv(i)+az(i)*zv(i))*fnuv
         au(i) = pi41*(dintu +dju(i)*analyt(i)*np)                      
         av(i) = pi41*(dintv +djv(i)*analyt(i)*np)                      
      enddo

      deallocate (ax, ay, az, analyt)

!     end  subroutine bn_vecpot
      end


      subroutine vecpot (ax, ay, az, i, nlo, nhi)
c ----------------------------------------------------------------------
      use kind_spec
      use bnvariables, only: indu, indv, guu, guv, gvv, djx, djy, djz,
     1    x, y, z, tanu, tanv      
      implicit none
      integer, intent(in) :: i, nlo, nhi
      integer :: ip, istat
      real(rprec), intent(inout) :: ax, ay, az
      real(rprec), target, allocatable :: dx(:), dy(:), dz(:)
      real(rprec), pointer :: sq(:), sqs(:), du(:), dv(:)
      real(rprec) :: sqsum
c ----------------------------------------------------------------------
      ip = nhi - nlo + 1
      if (ip .le. 0) return

      allocate (dx(nlo:nhi), dy(nlo:nhi), dz(nlo:nhi), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in BNORM routine vecpot'

      sq => dz
      du => dx
      dv => dy
      sqs => dy

      dx    = x(i) - x(nlo:nhi)
      dy    = y(i) - y(nlo:nhi)
      dz    = z(i) - z(nlo:nhi)
      sq    = 1._dp/sqrt(dx*dx + dy*dy + dz*dz)
      du    = tanu(indu(nlo:nhi) - indu(i))
      dv    = tanv(indv(nlo:nhi) - indv(i))
      sqs   = 1._dp/sqrt(guu(i)*du*du + 2._dp*guv(i)*du*dv 
     1                 + gvv(i)*dv*dv)
      
      sqsum = sum(sqs)
      
      ax = ax + sum(djx(nlo:nhi)*sq) - djx(i)*sqsum 
      ay = ay + sum(djy(nlo:nhi)*sq) - djy(i)*sqsum 
      az = az + sum(djz(nlo:nhi)*sq) - djz(i)*sqsum 
      
      deallocate (dx, dy, dz)
      
      end subroutine vecpot


! ----------------------------------------------------------------------
!     ROUTINES FOR WRITING NESCOIL INPUT FILE (Added 9/98: SPH)
! ----------------------------------------------------------------------
      subroutine bn_write_nescoil_input(extension)
      use meshes
      use bnvariables
      use neswrite, only: coil_separation, iota_edge, phip_edge, mnmax,
     1   ixm, ixn, nfp, ntheta, nzeta, rmnc, zmns,
     2  ixm_ws, ixn_ws, mnmax_ws, rmnc_ws, zmns_ws
      use safe_open_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      character*(*) :: extension
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ibex = 0
      integer, parameter :: nescoil0 = 15
      integer, parameter :: nmod_seg = 64, nfilaments = 10
      integer :: nu_plasma, nv_plasma
      real(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat, mn, m, n, iunit
      real(rprec) :: cut, cup
C-----------------------------------------------
      iunit = nescoil0
      call safe_open(iunit, istat, 'nescin.' // extension, 'replace',
     1      'formatted')
      if (istat .ne. 0) 
     1      stop 'Unable to open nescoil input file in bnorm'

!-----------------------------------------------
!     Write data for NESCOIL. First 3 lines (field periods) are not in
!     original NESCOIL code, but are needed here.
!-----------------------------------------------
      cut = zero         !Make these readable from input (wout file)
      cup = one          !1 for modular, 0 for saddles...
      
      write (iunit, 10) '------ Spatial dimensions ----'
      write (iunit, 10) 'nu, nv, nu1, nv1, npol, ntor'
      nu_plasma = nu; nv_plasma = nv
      write (iunit,*) nu, nv, nu_plasma, nv_plasma, nmod_seg, nfilaments
      write (iunit,*)
      
      write (iunit, 10) '------ Fourier Dimensions ----'
      write (iunit, 10) 'mf, nf, md, nd (max in surf and bnorm files)'
      write (iunit, *) mf, nf, md, nd
      write (iunit,*)

      write (iunit, 10) '------ Plasma information from VMEC ----'
      write (iunit, 10) 'np     iota_edge       phip_edge       curpol'
      write (iunit, *) nfp, iota_edge, phip_edge, curpol
      write (iunit,*)

      write (iunit, 10) '------ Current Controls ----'
      write (iunit, 10) 'cut  cup  ibex(=1,use fixed background coils)'
      write (iunit,*) cut, cup, ibex
      write (iunit,*)

      write (iunit, 10) '------ SVD controls -----'
      write (iunit, 10) 'mstrt, mstep, mkeep, mdspw, curwt, trgwt'
      write (iunit,*) 0,     0,      0,     4,    0.0_dp,   0.0_dp
      write (iunit,*)

      write (iunit, 10) '------ Output controls -----'
      write (iunit, 10) 'w_psurf w_csurf w_bnuv w_jsurf w_xerr w_svd'
      write (iunit,*) 0,       0,       0,      0,       0,      0
      write (iunit,*)


      write (iunit, 10) '------ Plasma Surface ---- '
      write (iunit, 10) 'Number of fourier modes in table'
      write (iunit,*) mnmax
      write (iunit, 10) 'Table of fourier coefficients'
      write (iunit, 10) '     m     n    cr (m,n)    cz (m,n)'
      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)        !(nfp divided out already, 
                            ! and n->-n for nescoil convention)
         write (iunit,'(x,2i6,1p2e12.4,1pe20.12)') m,
     1         n, cr(m,n), cz(m,n), cl(m,n)
      end do


      ntheta = 2*(mb+1) + 6
      ntheta = max(64, ntheta)
      nzeta  = 2*nb + 4
      nzeta  = max(32, nzeta)


      allocate (rmnc(mnmax), zmns(mnmax), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in bnorm bn_write'
      
      rmnc = zero; zmns = zero
      
      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)              !!Nescoil convention: mu + nv
         rmnc(mn) = cr(m,n)      
         zmns(mn) = cz(m,n)
      end do   
      
      call scaleup_boundary(ntheta, nzeta)

      write (iunit,*)
      write (iunit, '(a, 1pe12.4, a)')
     1  '------ Current Surface: Coil-Plasma separation = ',
     2  coil_separation,' -----'
      write (iunit, 10) 'Number of fourier modes in table'
      write (iunit,*) mnmax_ws
      write (iunit, 10) 'Table of fourier coefficients'
      write (iunit, 10) '     m     n    cr2(m,n)    cz2(m,n)'

      do mn = 1, mnmax_ws
        write (iunit,'(x,2i6,1p2e12.4)') ixm_ws(mn),
     1       ixn_ws(mn), rmnc_ws(mn), zmns_ws(mn)
      end do
      
      close (iunit)
 10   format (a)
      

      deallocate (rmnc_ws, zmns_ws, ixm_ws, ixn_ws)
      end subroutine bn_write_nescoil_input


      subroutine scaleup_boundary(nub, nvb)
      use neswrite
      use meshes, only: md, nd
      use normal_info, only: u0, v0, rb_ws, zb_ws, vb_ws
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nub, nvb
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: zero = 0, one = 1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: m, n, mn, n1, i, nuvb, iflag
      integer :: ku, kv, k, nmax, mmax
      real(rprec) :: cosnv, sinnv, alub, alvb, twopi, rtemp, ztemp
      real(rprec) :: cosmu, sinmu, theta, zeta, cosmn1, sinmn1, dnorm, 
     1   b1, c1, r1, re, ae, error, Map_v, sep_tol
      real(rprec), allocatable, dimension(:) :: rbn, zbn
C-----------------------------------------------
      external Map_v
!
!     Attempts to scale up boundary by a uniform distance coil_separation
!     On entry, rmnc_b, zmns_b contain the boundary coefficients
!     On exit,  they contain the scaled up boundary coefficients
!
      nuvb = nub * nvb
      twopi = 8*atan(one)
      alub = twopi/nub
      alvb = twopi/nvb                    !!Note: ixn has nfp factored out already...

      allocate (rmnc_ws(nuvb), zmns_ws(nuvb), ixm_ws(nuvb),
     1          ixn_ws(nuvb), stat=iflag)
      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'

      allocate (rbn(nuvb), zbn(nuvb), stat=iflag)
      if (iflag .ne. 0)stop 'Allocation error in bnorm scaleup_boundary'
 
      re = 10*epsilon(re)
      ae = re
      i = 0
      
      do ku = 1, nub
        u0 = alub*(ku - 1)
        do kv = 1, nvb
          i = i + 1
          v0 = alvb*(kv - 1)                                            !!Np*(Real toroidal angle)
          b1 = v0 - twopi/6
          c1 = v0 + twopi/6
          r1 = v0
!         Map v0 at plasma boundary to vb_ws at winding surface
          call fzero(Map_v, b1, c1, r1, re, ae, iflag)
          error = abs(Map_v (b1))
          if (iflag .gt. 2) print *,'  i = ', i,' iflag = ', iflag,
     1    ' v0 = ', v0,' vb_ws = ', vb_ws, ' v1 = ', b1, 
     2    ' Error mapping to Winding Surface in BNORM code ', error

          rbn(i) = rb_ws
          zbn(i) = zb_ws
        end do  
      end do  
       
!
!     FFT new surface (use NESCOIL convention, mu + nv)
!
      mmax = min(md, (nub-1)/2)
      nmax = min(nd, abs(nvb-1)/2)
      mnmax_ws = 0
      rmnc_ws = zero;  zmns_ws = zero
      sep_tol = 1.e-3_dp*coil_separation
      
      mloop: do m = 0, mmax
         nloop: do n = -nmax, nmax
            if (m.eq.0 .and. n.gt.0) cycle nloop
            dnorm = one/nuvb
            if (m.ne.0 .or. n.ne.0) dnorm = 2*dnorm
            i = 0
            rtemp = zero;  ztemp = zero
            do ku = 1, nub
               theta= alub*(ku-1)
               cosmu = cos(m*theta)*dnorm
               sinmu = sin(m*theta)*dnorm
               do kv = 1, nvb
                  i = i + 1
                  zeta = alvb*(kv-1)
                  cosnv = cos(n*zeta)
                  sinnv = sin(n*zeta)
                  cosmn1 = cosmu*cosnv - sinmu*sinnv          !cos(mu+nv) NESCOIL CONVENTION
                  sinmn1 = sinmu*cosnv + cosmu*sinnv          !sin(mu+nv)
                  rtemp = rtemp + rbn(i) * cosmn1
                  ztemp = ztemp + zbn(i) * sinmn1
               end do
            end do
            if (abs(rtemp).lt.sep_tol .and. abs(ztemp).lt.sep_tol)
     1      cycle nloop
            mnmax_ws = mnmax_ws+1
            rmnc_ws(mnmax_ws) = rtemp
            zmns_ws(mnmax_ws) = ztemp
            ixm_ws(mnmax_ws) = m
            ixn_ws(mnmax_ws) = n
         end do nloop
      end do mloop
       
      deallocate (rbn, zbn)


      end subroutine scaleup_boundary
            

      function Map_v(v)
      use kind_spec
      use normal_info, only: u0, v0, vb_ws
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: v, Map_v
C-----------------------------------------------

!
!     Finds phi (=vb_ws/Np) at a displaced winding surface corresponding to a
!     value of phi (=v/Np) on the plasma surface. It then computes the difference
!     compared with original phi (=v0/Np) on the plasma surface
!
      call normal_vector(u0, v)
c
      Map_v = vb_ws - v0                               !v0 is ORIGINAL Np*(Real Phi) on plasma surface

      end function Map_v


      subroutine normal_vector(u, v)
      use kind_spec
      use neswrite
      use normal_info
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: v, u
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, one = 1
      integer :: mn, m, n
      real(rprec) :: rbn, zbn, rub, rvb, zub, zvb,
     1    cosnv, sinnv, cosmn1, sinmn1, cosmu, sinmu, snx, sny
      real(rprec) :: x2, y2, v2, norm, twopi
C-----------------------------------------------
      
      twopi = 8*atan(one)

      rbn = zero;  rub = zero;   rvb = zero
      zbn = zero;  zub = zero;   zvb = zero
 
      do mn = 1, mnmax
         m = ixm(mn)
         n = ixn(mn)                          !!nfp factored out already, nescoil convention
         cosmu = cos(m*u)
         sinmu = sin(m*u)
         cosnv = cos(n*v)
         sinnv = sin(n*v)
         cosmn1 = cosmu*cosnv - sinmu*sinnv   !!cos(mu+nv), nescoil convention
         sinmn1 = sinmu*cosnv + cosmu*sinnv
         rbn = rbn + rmnc(mn) * cosmn1
         rub = rub - ixm(mn) * rmnc(mn) * sinmn1
         rvb = rvb - nfp*ixn(mn) * rmnc(mn) * sinmn1
         zbn = zbn + zmns(mn) * sinmn1
         zub = zub + ixm(mn) * zmns(mn) * cosmn1
         zvb = zvb + nfp*ixn(mn) * zmns(mn) * cosmn1
      end do
 
!
!     un-normalized r,phi,z components of outward normal vector
!
      surf_norm(1) = rbn * zub
      surf_norm(2) = rub * zvb - rvb * zub
      surf_norm(3) =-rbn * rub
 
      norm = sqrt(sum(surf_norm**2))
        
      surf_norm(:) = surf_norm(:) / norm
      
      cosnv = cos(v/nfp)
      sinnv = sin(v/nfp)
      
      snx = surf_norm(1)*cosnv - surf_norm(2)*sinnv
      sny = surf_norm(1)*sinnv + surf_norm(2)*cosnv

      x2    = rbn*cosnv + coil_separation*snx
      y2    = rbn*sinnv + coil_separation*sny
      v2    = atan2(y2,x2)                                  !REAL phi on winding surface
      if(v2 .lt. zero .and. x2 .lt. zero) v2 = v2 + twopi   !Cut at v = +- pi

      v2 = v2*nfp
      vb_ws = v2
      
      rb_ws = sqrt (x2*x2 + y2*y2)
      zb_ws = zbn + coil_separation*surf_norm(3)
      
      end subroutine normal_vector


      SUBROUTINE FZERO (F, B, C, R, RE, AE, IFLAG)
C***BEGIN PROLOGUE  FZERO
C***PURPOSE  Search for a zero of a function F(X) in a given interval
C            (B,C).  It is designed primarily for problems where F(B)
C            and F(C) have opposite signs.
C***LIBRARY   SLATEC
C***CATEGORY  F1B
C***TYPE      SINGLE PRECISION (FZERO-S, DFZERO-D)
C***KEYWORDS  BISECTION, NONLINEAR EQUATIONS, ROOTS, ZEROS
C***AUTHOR  Shampine, L. F., (SNLA)
C           Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     FZERO searches for a zero of a REAL function F(X) between the
C     given REAL values B and C until the width of the interval (B,C)
C     has collapsed to within a tolerance specified by the stopping
C     criterion,
C        ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C     The method used is an efficient combination of bisection and the
C     secant rule and is due to T. J. Dekker.
C
C     Description Of Arguments
C
C   F     :EXT   - Name of the REAL external function.  This name must
C                  be in an EXTERNAL statement in the calling program.
C                  F must be a function of one REAL argument.
C
C   B     :INOUT - One end of the REAL interval (B,C).  The value
C                  returned for B usually is the better approximation
C                  to a zero of F.
C
C   C     :INOUT - The other end of the REAL interval (B,C)
C
C   R     :IN    - A (better) REAL guess of a zero of F which could help
C                  in speeding up convergence.  If F(B) and F(R) have
C                  opposite signs, a root will be found in the interval
C                  (B,R); if not, but F(R) and F(C) have opposite signs,
C                  a root will be found in the interval (R,C);
C                  otherwise, the interval (B,C) will be searched for a
C                  possible root.  When no better guess is known, it is
C                  recommended that r be set to B or C, since if R is
C                  not interior to the interval (B,C), it will be
C                  ignored.
C
C   RE    :IN    - Relative error used for RW in the stopping criterion.
C                  If the requested RE is less than machine precision,
C                  then RW is set to approximately machine precision.
C
C   AE    :IN    - Absolute error used in the stopping criterion.  If
C                  the given interval (B,C) contains the origin, then a
C                  nonzero value should be chosen for AE.
C
C   IFLAG :OUT   - A status code.  User must check IFLAG after each
C                  call.  Control returns to the user from FZERO in all
C                  cases.
C
C                1  B is within the requested tolerance of a zero.
C                   The interval (B,C) collapsed to the requested
C                   tolerance, the function changes sign in (B,C), and
C                   F(X) decreased in magnitude as (B,C) collapsed.
C
C                2  F(B) = 0.  However, the interval (B,C) may not have
C                   collapsed to the requested tolerance.
C
C                3  B may be near a singular point of F(X).
C                   The interval (B,C) collapsed to the requested tol-
C                   erance and the function changes sign in (B,C), but
C                   F(X) increased in magnitude as (B,C) collapsed, i.e.
C                     ABS(F(B out)) .GT. MAX(ABS(F(B in)),ABS(F(C in)))
C
C                4  No change in sign of F(X) was found although the
C                   interval (B,C) collapsed to the requested tolerance.
C                   The user must examine this case and decide whether
C                   B is near a local minimum of F(X), or B is near a
C                   zero of even multiplicity, or neither of these.
C
C                5  Too many (.GT. 500) function evaluations used.
C
C***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
C                 code, Report SC-TM-70-631, Sandia Laboratories,
C                 September 1970.
C               T. J. Dekker, Finding a zero by means of successive
C                 linear interpolation, Constructive Aspects of the
C                 Fundamental Theorem of Algebra, edited by B. Dejon
C                 and P. Henrici, Wiley-Interscience, 1969.
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   700901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  FZERO
      USE KIND_SPEC
      IMPLICIT NONE
      REAL(RPREC), PARAMETER :: ZERO = 0, ONE = 1
      REAL(RPREC) :: A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R,
     +     RE,RW,T,TOL,Z,F
      INTEGER :: IC,IFLAG,KOUNT
      EXTERNAL F
C***FIRST EXECUTABLE STATEMENT  FZERO
C
C   ER is two times the computer unit roundoff value which is defined
C   here by the function EPSILON.
C
      ER = 2 * EPSILON(ER)
C
C   Initialize.
C
      Z = R
      IF (R .LE. MIN(B,C)  .OR.  R .GE. MAX(B,C)) Z = C
      RW = MAX(RE,ER)
      AW = MAX(AE,ZERO)
      IC = 0
      T = Z
      FZ = F(T)
      FC = FZ
      T = B
      FB = F(T)
      KOUNT = 2
      IF (SIGN(ONE,FZ) .EQ. SIGN(ONE,FB)) GO TO 1
      C = Z
      GO TO 2
    1 IF (Z .EQ. C) GO TO 2
      T = C
      FC = F(T)
      KOUNT = 3
      IF (SIGN(ONE,FZ) .EQ. SIGN(ONE,FC)) GO TO 2
      B = Z
      FB = FZ
    2 A = C
      FA = FC
      ACBS = ABS(B-C)
      FX = MAX(ABS(FB),ABS(FC))
C
    3 IF (ABS(FC) .GE. ABS(FB)) GO TO 4
C
C   Perform interchange.
C
      A = B
      FA = FB
      B = C
      FB = FC
      C = A
      FC = FA
C
    4 CMB = 0.5_DP*(C-B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AW
C
C   Test stopping criterion and function count.
C
      IF (ACMB .LE. TOL) GO TO 10
      IF (FB .EQ. ZERO) GO TO 11
      IF (KOUNT .GE. 500) GO TO 14
C
C   Calculate new iterate implicitly as B+P/Q, where we arrange
C   P .GE. 0.  The implicit form is used to prevent overflow.
C
      P = (B-A)*FB
      Q = FA - FB
      IF (P .GE. ZERO) GO TO 5
      P = -P
      Q = -Q
C
C   Update A and check for satisfactory reduction in the size of the
C   bracketing interval.  If not, perform bisection.
C
    5 A = B
      FA = FB
      IC = IC + 1
      IF (IC .LT. 4) GO TO 6
      IF (8*ACMB .GE. ACBS) GO TO 8
      IC = 0
      ACBS = ACMB
C
C   Test for too small a change.
C
    6 IF (P .GT. ABS(Q)*TOL) GO TO 7
C
C   Increment by TOLerance.
C
      B = B + SIGN(TOL,CMB)
      GO TO 9
C
C   Root ought to be between B and (C+B)/2.
C
    7 IF (P .GE. CMB*Q) GO TO 8
C
C   Use secant rule.
C
      B = B + P/Q
      GO TO 9
C
C   Use bisection (C+B)/2.
C
    8 B = B + CMB
C
C   Have completed computation for new iterate B.
C
    9 T = B
      FB = F(T)
      KOUNT = KOUNT + 1
C
C   Decide whether next step is interpolation or extrapolation.
C
      IF (SIGN(ONE,FB) .NE. SIGN(ONE,FC)) GO TO 3
      C = A
      FC = FA
      GO TO 3
C
C   Finished.  Process results for proper setting of IFLAG.
C
   10 IF (SIGN(ONE,FB) .EQ. SIGN(ONE,FC)) GO TO 13
      IF (ABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN
   11 IFLAG = 2
      RETURN
   12 IFLAG = 3
      RETURN
   13 IFLAG = 4
      RETURN
   14 IFLAG = 5

      END SUBROUTINE FZERO
