!     THIS ROUTINE COMPUTES .5 * B**2 ON THE VACUUM / PLASMA SURFACE
!     BASED ON THE PROGRAM BY P. MERKEL [J. Comp. Phys. 66, 83 (1986)]
!     AND MODIFIED BY W. I. VAN RIJ AND S. P. HIRSHMAN (1987)
 
!     THE USER MUST SUPPLY THE FILE << MGRID >> WHICH INCLUDES THE MAGNETIC
!     FIELD DATA TO BE READ BY THE SUBROUTINE BECOIL
!     THE "VACUUM.INC" FILE IS DEFINED IN VMEC.UNIX
!
      subroutine vacuum(rmnc, rmns, zmns, zmnc, xm, xn, rcurr, zcurr, 
     1   plascur, rbtor, wint, ns, ivac_skip, ivac, mnmax, ier_flag, 
     2   lscreen)
      use vacmod
      use vparams, only: nthreed, zero, one, dmu0
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: ns, ivac_skip, ivac, mnmax, ier_flag
      real(rprec) :: plascur, rbtor
      real(rprec), dimension(mnmax), intent(in) :: 
     1   rmnc, rmns, zmns, zmnc, xm, xn
      real(rprec), dimension(nv), intent(in) :: rcurr, zcurr
      real(rprec), dimension(*), intent(in) :: wint
      logical :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: mn, n, n1, m, i, istore_max
      real(rprec), dimension(:), pointer :: potcos, potsin
      real(rprec), allocatable :: bsubu(:), bsubv(:), potu(:), potv(:)
      real(rprec), allocatable :: amatrix(:)
      real(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv,
     1    det, bsupu, bsupv, bsubuvac, fac
C-----------------------------------------------
      if (.not.allocated(potvac)) stop 'POTVAC not allocated in VACCUM'

      allocate (amatrix(mnpd2*mnpd2), bsubu(nuv2), bsubv(nuv2),
     1    potu(nuv2), potv(nuv2), stat = i)
      if (i .ne. 0) stop 'Allocation error in vacuum'      

      potsin => potvac(1:mnpd)
      potcos => potvac(1+mnpd:)

      allocate (bexu(nuv2), bexv(nuv2), bexn(nuv2),
     1     bexni(nuv2), r1b(nuv), rub(nuv2), rvb(nuv2),
     2     z1b(nuv), zub(nuv2), zvb(nuv2), auu(nuv2), auv(nuv2), 
     3     avv(nuv2), snr(nuv2), snv(nuv2), snz(nuv2), drv(nuv2),
     4     guu_b(nuv2), guv_b(nuv2), gvv_b(nuv2), rb2(nuv),
     5     rcosuv(nuv), rsinuv(nuv), stat=i)
      if (i .ne. 0) stop 'Allocation error in vacuum'

!
!       INDEX OF LOCAL VARIABLES
!
!       rmnc,rmns,zmns,zmnc:     Surface Fourier coefficients (m,n) of R,Z
!       xm,xn:     m, n values corresponding to rc,zs array
!       rcurr,zcurr:
!                  Position of magnetic axis (toroidal current filament)
!       bsqvac:    B**2/2 at the vacuum interface
!       plascur:   net toroidal current
!       mnmax:     number of R, Z modes in Fourier series of R,Z
!       ivac_skip: regulates whether full (=0) or incremental (>0)
!                 update of matrix elements is necessary
!
!
!       compute and store mean magnetic fields (due to
!       toroidal plasma current and external tf-coils)
!       note: these are fixed for a constant current iteration
!       bfield = rbtor*grad(zeta) + plascur*grad("theta") - grad(potential)
!            where "theta" is computed using Biot-Savart law for filaments
!
!       Here, the potential term is needed to satisfy B ! dS = 0 and has the form:
!
!       potential = sum potsin*sin(mu - nv) + potcos*cos(mu - nv)
!
 
      call surface (rmnc, rmns, zmns, zmnc, xm, xn, mnmax)
      call bextern (rcurr, zcurr, plascur, rbtor, wint, ns)
!
!       Determine scalar magnetic potential
!
      istore_max = min(64,nuv2)
      call scalpot (potvac, amatrix, wint, ns, istore_max, ivac_skip)
      call solver (amatrix, potvac, mnpd2)
!
!       compute tangential covariant (sub u,v) and contravariant
!       (super u,v) magnetic field components on the plasma surface
!
      potu(:nuv2) = zero;  potv(:nuv2) = zero
 
      mn = 0
      do n = -nf, nf
         dn2 = -(n*nfper)
         n1 = abs(n)
         do m = 0, mf
            mn = mn + 1
            dm2 = (m)
            do i = 1, nuv2
               cosmn = cosu1(i,m)*cosv1(i,n1) + csign(n)*sinu1(i,m)*
     1            sinv1(i,n1)
               potu(i) = potu(i) + dm2*potsin(mn)*cosmn 
               potv(i) = potv(i) + dn2*potsin(mn)*cosmn 
            end do
            if (.not.lasym) cycle
            do i = 1, nuv2
               sinmn = sinu1(i,m)*cosv1(i,n1) - csign(n)*cosu1(i,m)*
     1            sinv1(i,n1)
               potu(i) = potu(i) - dm2*potcos(mn)*sinmn
               potv(i) = potv(i) - dn2*potcos(mn)*sinmn
            end do
         end do
      end do
      do i = 1, nuv2
         bsubu(i) = potu(i) + bexu(i)
         bsubv(i) = potv(i) + bexv(i)
         huv = 0.5_dp*guv_b(i)*(nfper)
         hvv = gvv_b(i)*(nfper*nfper)
         det = one/(guu_b(i)*hvv-huv*huv)
         bsupu = (hvv*bsubu(i)-huv*bsubv(i))*det
         bsupv = ((-huv*bsubu(i))+guu_b(i)*bsubv(i))*det
         bpolvac(i) = 0.5_dp*bsubu(i)*bsupu
         bsqvac(i) = bpolvac(i) + 0.5_dp*bsubv(i)*bsupv
         brv(i) = rub(i)*bsupu + rvb(i)*bsupv
         bphiv(i) = r1b(i)*bsupv
         bzv(i) = zub(i)*bsupu + zvb(i)*bsupv
      end do

!
!       PRINT OUT VACUUM PARAMETERS
!
      if (ivac .eq. 0) then
         ivac = ivac + 1
         if (lscreen) write (*, 200) nfper, mf, nf, nu, nv
         write (nthreed, 200) nfper, mf, nf, nu, nv
  200    format(/,2x,'In VACUUM, np =',i3,2x,'mf =',i3,2x,'nf =',i3,
     1      ' nu =',i3,2x,'nv = ',i4)
         bsubuvac = sum(bsubu(:nuv2)*wint(ns:ns*nuv2:ns))
         bsubvvac = sum(bsubv(:nuv2)*wint(ns:ns*nuv2:ns))
         fac = 1.e-6_dp/dmu0
         if (lscreen )write (*,1000) (-pi2*bsubuvac*fac),
     1       plascur*fac, bsubvvac, rbtor
         write (nthreed, 1000) (-pi2*bsubuvac*fac), plascur*fac, 
     1       bsubvvac, rbtor
 1000    format(2x,'2*pi * a * -BPOL(vac) = ',1pe10.2,
     1      ' TOROIDAL CURRENT = ',1pe10.2,/,2x,'R * BTOR(vac) = ',
     2      1pe10.2,' R-BTOR = ',1pe10.2)
         if (rbtor*bsubvvac .lt. zero) ier_flag = 7
         if (abs((plascur+pi2*bsubuvac)/rbtor) .gt. 1.e-3_dp)
     1      ier_flag = 10
      endif

      if (allocated(bexu))
     1    deallocate (bexu, bexv, bexn, bexni, r1b, rub, rvb, z1b, zub, 
     2    zvb, auu, auv, avv, snr, snv, snz, drv, guu_b, guv_b, gvv_b, 
     3    rb2, rcosuv, rsinuv, stat=i)
      if (i .ne. 0) stop 'Deallocation error in vacuum'
       
      deallocate (amatrix, bsubu, bsubv, potu, potv, stat = i)

      end subroutine vacuum

      
      subroutine analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, l, ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
      real(rprec), dimension(nuv2), intent(in) :: sl, tl
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
      real(rprec), allocatable :: sinp(:), cosp(:)
C-----------------------------------------------
      allocate (sinp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum'
      
      if (n .ne. abs(n)) stop 'error calling analysum!'
 
      sinp = (sinu1(:,m)*cosv1(:,n) - sinv1(:,n)*
     1            cosu1(:,m))*cmns(l,m,n)
      bvec(m,n,1) = bvec(m,n,1) + sum(tl*bexni*sinp)

      if (ivacskip .eq. 0) grpmn(:,m,n,1) = grpmn(:,m,n,1) + sl*sinp
      deallocate (sinp)

      if (.not.lasym) return

      allocate (sinp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum'
      cosp = (cosu1(:,m)*cosv1(:,n) + cosv1(:,n)*
     1        cosu1(:,m))*cmns(l,m,n)
      bvec(m,n,2) = bvec(m,n,2) + sum(tl*bexni*cosp)
      if (ivacskip .eq. 0) grpmn(:,m,n,2) = grpmn(:,m,n,2) + sl*cosp
      deallocate (cosp)
 
      end subroutine analysum

      
      subroutine analysum2(grpmn, bvec, slp, tlp, slm, tlm, 
     1    m, n, l, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer m, n, l, ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
      real(rprec), dimension(nuv2), intent(in) :: 
     1   slp, tlp, slm, tlm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: istat
      real(rprec), allocatable :: sinp(:), sinm(:), temp(:)
      real(rprec), allocatable :: cosp(:), cosm(:)
C-----------------------------------------------
      allocate (sinp(nuv2), sinm(nuv2), temp(nuv2), stat=istat)
      if (istat .ne. 0) stop 'Allocation error in analysum2'
      if (n .ne. abs(n)) stop 'error calling analysum2!'
 
      sinp =  sinu1(:,m)*cosv1(:,n)*cmns(l,m,n)
      temp = -sinv1(:,n)*cosu1(:,m)*cmns(l,m,n)
      sinm = sinp - temp                !sin(mu + |n|v) * cmns (l,m,|n|)
      sinp = sinp + temp                !sin(mu - |n|v) * cmns
      bvec(m,n,1)  = bvec(m,n,1)  + sum(tlp*bexni*sinp)
      bvec(m,-n,1) = bvec(m,-n,1) + sum(tlm*bexni*sinm)

      if (ivacskip .eq. 0) then 
         grpmn(:,m,n,1)  = grpmn(:,m,n,1)  + slp*sinp
         grpmn(:,m,-n,1) = grpmn(:,m,-n,1) + slm*sinm
      end if   

      deallocate (sinp, sinm, temp)

      if (.not.lasym) return
      
      allocate (cosp(nuv2), cosm(nuv2), stat=istat)  
      if (istat .ne. 0) stop 'Allocation error in analysum2'

      cosp = cosu1(:,m)*cosv1(:,n)*cmns(l,m,n)     
      temp = cosv1(:,n)*cosu1(:,m)*cmns(l,m,n)
      cosm = cosp - temp                !cos(mu + |n|v) * cmns (l,m,|n|) 
      cosp = cosp + temp                !cos(mu - |n|v) * cmns (l,m,|n|) 
      bvec(m,n,2)  = bvec(m,n,2)  + sum(tlp*bexni*cosp)
      bvec(m,-n,2) = bvec(m,-n,2) + sum(tlm*bexni*cosm)

      if (ivacskip .eq. 0) then 
         grpmn(:,m,n,2)  = grpmn(:,m,n,2)  + slp*cosp
         grpmn(:,m,-n,2) = grpmn(:,m,-n,2) + slm*cosm
      end if   

      deallocate (cosp, cosm)
 
      end subroutine analysum2


      subroutine analyt(grpmn, bvec, ivacskip)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ivacskip
      real(rprec), intent(inout) :: grpmn(nuv2,0:mf,-nf:nf,*)
      real(rprec), intent(inout) :: bvec(0:mf,-nf:nf,*)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, one=1, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, n, m
      real(rprec), dimension(:), allocatable :: 
     1   r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, tlp, tlm2,
     2    tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, slp, tlpm, slpm,
     3    delt1u, azp1u, azm1u, cma11u, sqad1u, sqad2u
      real(rprec) :: fl, sign1
C-----------------------------------------------
      allocate (r0p(nuv2), r1p(nuv2), r0m(nuv2), r1m(nuv2), 
     1          sqrtc(nuv2), sqrta(nuv2), tlp2(nuv2), tlp1(nuv2), 
     2          tlp(nuv2), tlm2(nuv2), tlm1(nuv2), tlm(nuv2), adp(nuv2),
     3          adm(nuv2), cma(nuv2), ra1p(nuv2), ra1m(nuv2), slm(nuv2), 
     4          slp(nuv2), tlpm(nuv2), slpm(nuv2), delt1u(nuv2), 
     5          azp1u(nuv2), azm1u(nuv2), cma11u(nuv2), sqad1u(nuv2), 
     6          sqad2u(nuv2), stat = l)
      if (l .ne. 0) stop 'Allocation error in analyt'

!
!     ALL EQUATIONS REFER TO THE APPENDIX OF THE PAPER BY P. MERKEL
!     IN J. COMPUT. PHYSICS 66, p83 (1986)
!
!     THE REQUIRED INTEGRALS ARE:
!
!     BVECS(m,n) = Int< sin(mu - nv) Int<BNORM`*Gan(u`,v`,u-u`,v-v`)>
!     BVECC(m,n) = Int< cos(mu - nv) Int<BNORM`*Gan(u`,v`,u-u`,v-v`)>
!
!     Where Int<...> means integration of u (theta) and v (zeta) and
!     summation over field periods. In terms of Merkels Imn integral,
!     a = guu (g theta-theta), etc., we have
!
!     BVECS(m,n) = (2*pi/nfp) * Int<BNORM *sin(mu` - nv`)[Im,-n](a,b,c)>
!     BVECC(m,n) = (2*pi/nfp) * Int<BNORM *cos(mu` - nv`)[Im,-n](a,b,c)>
!
!     Similarly, the analytic part of the matrix A(m,n;m`,n`) can be written:
!
!     A(m,n;m`,n`) = (2*pi/nfp) * Int<sin(mu` - nv`)*sin(m`u` - n`v`)
!                              [Km,-n](a`,b`,c`;A`,B`,C`)>
!
!     On exit, GRPMN(ip,m,n) = ALP * SIN(ip,m,n) * K[m,-n](ip)
!
!
!     COMPUTE ALL QUANTITIES INDEPENDENT OF THE MODE INDICES L,M,N
!
!     ADP(M): a +(-)2b + c
!     CMA:    c - a
!     DELTA:  4*(ac - b**2)
!     AZP(M): A +(-)2*B + C
!     CMA1:   C - A
!     R1P(M): Coefficient of l*Tl+(-) in eq (A17)
!     R0P(M): Coefficient of l*T(l-1)+(-) in eq (A17)
!     RA1P(M):Coefficient of Tl+(-) in eq (A17)
!
      adp  = guu_b  + guv_b  + gvv_b 
      adm  = guu_b  - guv_b  + gvv_b 
      cma  = gvv_b  - guu_b 
      sqrtc  = two*sqrt(gvv_b)
      sqrta  = two*sqrt(guu_b)
 
      if (ivacskip .eq. 0) then
         delt1u  = adp *adm  - cma *cma 
         azp1u  = auu  + auv  + avv 
         azm1u  = auu  - auv  + avv 
         cma11u  = avv  - auu 
         r1p  = (azp1u *(delt1u  - (cma *cma ))/adp 
     1        - (azm1u *adp ) + ((two*cma11u )*cma ))/delt1u 
         r1m  = (azm1u *(delt1u  - (cma *cma ))/adm 
     1        - (azp1u *adm ) + ((two*cma11u )*cma ))/delt1u 
         r0p  = ((-(azp1u *adm )*cma /adp ) - azm1u *cma 
     1        + (two*cma11u )*adm )/delt1u 
         r0m  = ((-(azm1u *adp )*cma /adm ) - azp1u *cma 
     1        + (two*cma11u )*adp )/delt1u 
         ra1p = azp1u /adp 
         ra1m = azm1u /adm 
      endif
!
!     INITIALIZE T0+ and T0-
!
!     TLP(M): TL+(-)
!     TLP(M)1:T(L-1)+(-)
!     TLP(M)2:T(L-2)+(-)
!
      sqad1u  = sqrt(adp )
      sqad2u  = sqrt(adm )
      tlp1 = zero
      tlm1 = zero
      tlp  = one/sqad1u *log((sqad1u *sqrtc  + adp 
     1     + cma )/(sqad1u *sqrta  - adp  + cma ))
      tlm  = one/sqad2u *log((sqad2u *sqrtc  + adm  
     1     + cma )/(sqad2u *sqrta  - adm  + cma ))
      tlpm = tlp  + tlm 
!
!     BEGIN L-SUM IN EQ (A14) TO COMPUTE Cmn COEFFICIENTS
!
      do l = 0, mf + nf
         fl = l
         sign1 = one - two*mod(l,2)          !!(-1)**l
!
!       COMPUTE SL+ and SL- , Eq (A17)
!       SLP(M): SL+(-)
!
         if (ivacskip .eq. 0) then
            slp = (r1p*fl + ra1p)*tlp + r0p*fl*tlp1 - (r1p + r0p)/sqrtc
     1          + sign1*(r0p - r1p)/sqrta
            slm = (r1m*fl + ra1m)*tlm + r0m*fl*tlm1 - (r1m + r0m)/sqrtc
     1          + sign1*(r0m - r1m)/sqrta
            slpm = slp + slm
         endif
!
!       BEGIN MODE NUMBER (m,n) LOOP
!
         do n = 0, nf
            do m = 0, mf
               if (cmns(l,m,n) .eq. zero) cycle 

               if (n.eq.0 .or. m.eq.0) then
!
!       1. n = 0 and  m >= 0  OR n > 0 and m = 0
!
                 call analysum (grpmn, bvec, slpm, tlpm, m, n, l, 
     1               ivacskip)
 
               else
!
!       2. n>=1  and  m>=1
!
                 call analysum2 (grpmn, bvec, slm, tlm, slp, tlp, 
     1               m, n, l, ivacskip)
               endif
            end do
         end do
!
!       UPDATE "TL's" (FOR L -> L+1) USING EQ (A15)
!
         tlp2 = tlp1
         tlm2 = tlm1
         tlp1 = tlp
         tlm1 = tlm
         tlp = ((sqrtc - (sign1*sqrta)) - (((two*fl) + one)*cma)*tlp1 - 
     1      fl*adm*tlp2)/(adp*(fl + one))
         tlm = ((sqrtc - (sign1*sqrta)) - (((two*fl) + one)*cma)*tlm1 - 
     1      fl*adp*tlm2)/(adm*(fl + one))
         tlpm = tlp + tlm
      end do
 
      deallocate (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1, 
     1          tlp, tlm2, tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm, 
     2          slp, tlpm, slpm, delt1u, azp1u, azm1u, cma11u, sqad1u,
     3          sqad2u, stat = l)

      end subroutine analyt

      
      subroutine becoil(rad, zee, br, bp, bz, brvac, bpvac, bzvac)
      use vparams, only: one, nthreed
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*), intent(in) :: rad, zee
      real(rprec), dimension(*), intent(out) :: br, bp, bz
      real(rprec), dimension(nr0b,nz0b,nv), intent(in) :: 
     1   brvac, bpvac, bzvac
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      character*(50), parameter :: warning =
     1   'Plasma Boundary exceeded Vacuum Grid Size'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, save :: icount = 0
      integer :: igrid, i, kv, ir, jz, ir1, jz1
      real(rprec) :: rad0, zee0, ri, zj, 
     1   dri, dzj, f22, f21, f12, f11
C-----------------------------------------------
!
!       THIS ROUTINE FINDS THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!       MAGNETIC FIELD AT A FIXED PHI PLANE BY 2-D INTERPOLATION
!
      igrid = 0
      icount = icount + 1

      do i = 1, nuv2
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE VACUUM GRID.  IF NOT,
!       SET THEM EQUAL TO LARGEST (OR SMALLEST) VACUUM GRID POINTS
!
         if (rad(i) .gt. rmaxb) then
            rad0 = rmaxb
            igrid = 1
         else if (rad(i) .lt. rminb) then
            igrid = 1
            rad0 = rminb
         else
            rad0 = rad(i)
         endif
         if (zee(i) .gt. zmaxb) then
            igrid = 1
            zee0 = zmaxb
         else if (zee(i) .lt. zminb) then
            igrid = 1
            zee0 = zminb
         else
            zee0 = zee(i)
         endif
!
!       DETERMINE PHI-PLANE, KV (MUST LIE IN FIRST FIELD PERIOD)
!
         kv = 1 + mod(i - 1,nv)
!
!
!       DETERMINE R, Z CORNER INDICES (IR,IZ)
!
         ir = int((rad0 - rminb)/delrb) + 1
         jz = int((zee0 - zminb)/delzb) + 1
         ir1 = min0(nr0b,ir + 1)
         jz1 = min0(nz0b,jz + 1)
!
!       COMPUTE R , Z AND DR , DZ AT MESH POINT (IR , IZ)
!
         ri = rminb + (ir - 1)*delrb
         zj = zminb + (jz - 1)*delzb
         dri = (rad0 - ri)/delrb
         dzj = (zee0 - zj)/delzb
         f22 = dri*dzj
         f21 = dri - f22
         f12 = dzj - f22
         f11 = one + f22 - (dri + dzj)
!
!       COMPUTE INTERPOLATED B FIELD
!
         br(i) = f11*brvac(ir,jz,kv) + f22*brvac(ir1,jz1,kv) + 
     1      f21*brvac(ir1,jz,kv) + f12*brvac(ir,jz1,kv)
         bz(i) = f11*bzvac(ir,jz,kv) + f22*bzvac(ir1,jz1,kv) + 
     1      f21*bzvac(ir1,jz,kv) + f12*bzvac(ir,jz1,kv)
         bp(i) = f11*bpvac(ir,jz,kv) + f22*bpvac(ir1,jz1,kv) + 
     1      f21*bpvac(ir1,jz,kv) + f12*bpvac(ir,jz1,kv)
 
      end do
 
      if (igrid.eq.1 .and. mod(icount,25).eq.0) then
         print *, warning
         write (nthreed, *) warning
      endif
 
      end subroutine becoil
      

      subroutine belicu(bx, by, bz, cos1, sin1, rp, zp, torcur)
      use vacmod
      use vacwires
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), intent(in) :: torcur
      real(rprec), dimension(*), intent(in)  :: cos1, sin1, rp, zp
      real(rprec), dimension(*), intent(out) :: bx, by, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
      integer :: i, j
      real(rprec), allocatable, dimension(:) :: x0, y0, z0, rw, fa
      real(rprec) :: sumx, sumy, sumz, ax, ay, az, curpi, 
     1   xp, yp, r12
C-----------------------------------------------
 
      curpi = torcur/pi4
      if (nvp .ne. nv*nfper) stop 'nvp != nv*nfper in belicu'
      allocate (x0(nvp+1), y0(nvp+1), z0(nvp+1), rw(nvp+1), fa(nvp+1))

      do j = 1,nuv2
        xp   = rp(j) * cos1(j)
        yp   = rp(j) * sin1(j)

        x0   = xp - xw
        y0   = yp - yw
        z0   = zp(j) - zw
        rw   = sqrt(x0*x0 + y0*y0 + z0*z0 )

        do 20 i = 1,nvp
          r12    = rw(i+1)*rw(i)
          fa(i)  = (rw(i+1)+rw(i))/
     1            (r12*(r12 +x0(i+1)*x0(i)+y0(i+1)*y0(i)+z0(i+1)*z0(i)))
 20     continue

        ax = 0;    ay = 0;      az = 0
        sumx = 0;  sumy = 0;    sumz = 0

        do 30 i = 1,nvp
          ax = ax + fa(i)*dx(i)
          ay = ay + fa(i)*dy(i)
          az = az + fa(i)*dz(i)
          sumx = sumx + fa(i)*vx(i)
          sumy = sumy + fa(i)*vy(i)
          sumz = sumz + fa(i)*vz(i)
 30     continue
        bx(j)  = curpi * (sumx - yp*az    + zp(j)*ay)
        by(j)  = curpi * (sumy - zp(j)*ax + xp*az)
        bz(j)  = curpi * (sumz - xp*ay    + yp*ax)
      enddo

      deallocate (x0, y0, z0, rw, fa)

      end subroutine belicu

      
      subroutine bextern(rcurr, zcurr, plascur, rbtor, wint, ns)
      use vacmod
      use vacwires
      use mgrid_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns
      real(rprec), intent(in) :: plascur, rbtor
      real(rprec), dimension(nv), intent(in) :: rcurr, zcurr
      real(rprec), dimension(*), intent(in) :: wint
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec), allocatable :: brad(:), bphi(:), bz(:)
C-----------------------------------------------
      if (.not.allocated(bvac)) stop 'BVAC Unallocated in bextern'
      allocate (brad(nuv2), bphi(nuv2), bz(nuv2), stat=i)
      if (i .ne. 0) stop 'allocation error in bextern'

!
!     THIS ROUTINE COMPUTES THE EXTERNAL, COIL PRODUCED B DOT DS
!     NOTE THAT BEXN = - BEX * DS IS THE EFFECTIVE SOURCE TERM
!
      call becoil(r1b,z1b,brad,bphi,bz,bvac(1,1),bvac(1,3),bvac(1,2))

!
!     COMPUTE CONTRIBUTION FROM NET TOROIDAL PLASMA CURRENT
!
      if (abs(plascur) > 1.e-8_dp*abs(rbtor)) then
         allocate (xw(nvp+1),yw(nvp+1),zw(nvp+1),vx(nvp+1),
     1     vy(nvp+1),vz(nvp+1),dx(nvp+1),dy(nvp+1),dz(nvp+1),stat=i)
         if (i .ne. 0) stop 'allocation error in bextern subroutine'
         call tolicu (rcurr, zcurr)
         call belicu (bexu, bexv, bexn, cosuv, sinuv, r1b, z1b, plascur)
         deallocate (xw, yw, zw, vx, vy, vz, dx, dy, dz)
         do i = 1, nuv2
            brad(i) = brad(i) + bexu(i)*cosuv(i) + bexv(i)*sinuv(i)
            bphi(i) = bphi(i) - bexu(i)*sinuv(i) + bexv(i)*cosuv(i)
            bz(i) = bz(i) + bexn(i)
         end do
      endif
 
      do i = 1, nuv2
        bexu(i) = rub(i)*brad(i) + zub(i)*bz(i)
        bexv(i) = rvb(i)*brad(i) + zvb(i)*bz(i) + r1b(i)*bphi(i)
        bexn(i) = -(brad(i)*snr(i)+bphi(i)*snv(i)+bz(i)*snz(i))
      end do

      deallocate (brad, bphi, bz)
      bexni(:nuv2) = wint(ns:nuv2*ns:ns)*bexn*pi2*pi2
       
      end subroutine bextern
      

      subroutine fouri(grpmn, gsource, amatrix, amatsq, bvec, wint, ns)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns
      real(rprec), dimension(nv,nu3,mnpd,*), intent(in) :: grpmn
      real(rprec), dimension(nuv), intent(in) :: gsource
      real(rprec), dimension(mnpd,mnpd,*) :: amatrix
      real(rprec), dimension(mnpd2,mnpd2), intent(out) :: amatsq
      real(rprec), dimension(0:mf,-nf:nf,*) :: bvec
      real(rprec), dimension(*) :: wint
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, i, j, n, kvi, kui, mn, m
      real(rprec), allocatable, dimension(:,:,:) :: bcos, bsin, source
      real(rprec), allocatable :: actemp(:,:,:,:), astemp(:,:,:,:)
      real(rprec) :: cosn, sinn, cosm, sinm
C-----------------------------------------------
!
!       AMATRIX(,1) = A(Sin)(Sin);  AMATRIX(,2) = A(Sin)(Cos);
!       AMATRIX(,3) = A(Cos)(Sin);  AMATRIX(,4) = A(Cos)(Cos)
!
!       SYMMETRIZE SOURCE TERMS
!
      allocate (bcos(nu2,-nf:nf,2), bsin(nu2,-nf:nf,2),
     1   actemp(nu3,-nf:nf,mnpd,2), astemp(nu3,-nf:nf,mnpd,2),
     2   source(nv,nu2,2), stat = i)
      if (i .ne. 0) stop 'allocation error in fouri'

      k = 0
      do i = 1, nu2
         do j = 1, nv
            k = k + 1
            source(j,i,1) = 0.5_dp*onp*(gsource(k)-gsource(imirr(k)))
            if (lasym)
     1      source(j,i,2) = 0.5_dp*onp*(gsource(k)+gsource(imirr(k)))
         end do
      end do
!
!       INITIALIZE SUMMED VECTORS TO ZERO
!
      bcos(:,(-nf):nf,:) = 0
      bsin(:,(-nf):nf,:) = 0
      actemp(:,(-nf):nf,:,:) = 0
      astemp(:,(-nf):nf,:,:) = 0
!
!       PERFORM KV (TOROIDAL ANGLE) TRANSFORM
!
      do n = 0, nf
         do kvi = 1, nv
            cosn = cosv(n,kvi)
            sinn = sinv(n,kvi)
            bcos(:,n,1) = bcos(:,n,1) + cosn*source(kvi,:,1)
            bsin(:,n,1) = bsin(:,n,1) + sinn*source(kvi,:,1)
            actemp(:,n,:,1) = actemp(:,n,:,1) + cosn*grpmn(kvi,:,:,1)
            astemp(:,n,:,1) = astemp(:,n,:,1) + sinn*grpmn(kvi,:,:,1)

            if (lasym) then
               bcos(:,n,2) = bcos(:,n,2) + cosn*source(kvi,:,2)
               bsin(:,n,2) = bsin(:,n,2) + sinn*source(kvi,:,2)
               actemp(:,n,:,2) = actemp(:,n,:,2) + cosn*grpmn(kvi,:,:,2)
               astemp(:,n,:,2) = astemp(:,n,:,2) + sinn*grpmn(kvi,:,:,2)
            end if
   
            if (n .ne. 0) then
               bcos(:,(-n),1) = bcos(:,n,1)
               bsin(:,(-n),1) = -bsin(:,n,1)
               actemp(:,(-n),:,1) = actemp(:,n,:,1)
               astemp(:,(-n),:,1) = -astemp(:,n,:,1)

               if (lasym) then
                  bcos(:,(-n),2) = bcos(:,n,2)
                  bsin(:,(-n),2) = -bsin(:,n,2)
               actemp(:,(-n),:,2) = actemp(:,n,:,2)
               astemp(:,(-n),:,2) = -astemp(:,n,:,2)
               end if

            endif
         end do
      end do
!
!       PERFORM KU (POLOIDAL ANGLE) TRANSFORM
!
      do m = 0, mf
         do kui = 1, nu2
            cosm = -cosui(m,kui)
            sinm = sinui(m,kui)
            bvec(m,-nf:nf,1) = bvec(m,-nf:nf,1) + 
     1       bcos(kui,-nf:nf,1)*sinm + bsin(kui,-nf:nf,1)*cosm
            if (lasym) then
            bvec(m,-nf:nf,2) = bvec(m,-nf:nf,2) - 
     1         bcos(kui,-nf:nf,2)*cosm + bsin(kui,-nf:nf,2)*sinm
            end if
         end do
!
!        NOTE: TRANSPOSE KUI, MN INDICES HERE ...
!
         do kui = 1, nu3
            cosm = -cosu(m,kui)*wint(kui*ns*nv)*pi2*pi2
            sinm = sinu(m,kui)*wint(kui*ns*nv)*pi2*pi2
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,1) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,1) + sinm*transpose(actemp(kui,-nf:nf,:,1)) + 
     2         cosm*transpose(astemp(kui,-nf:nf,:,1))
            if (lasym) then
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,2) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,2) + sinm*transpose(actemp(kui,-nf:nf,:,2)) + 
     2         cosm*transpose(astemp(kui,-nf:nf,:,2))
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,3) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,3) - cosm*transpose(actemp(kui,-nf:nf,:,1)) + 
     2         sinm*transpose(astemp(kui,-nf:nf,:,1))
            amatrix(:,m+1:nf*2*mf1+m+1:mf1,4) = amatrix(:,m+1:nf*2*mf1+m
     1         +1:mf1,4) - cosm*transpose(actemp(kui,-nf:nf,:,2)) + 
     2         sinm*transpose(astemp(kui,-nf:nf,:,2))
            end if
         end do
      end do

      deallocate (bcos, bsin, actemp, astemp, source, stat = i)

!
!       ZERO BVEC(0,n) AND AMATRIX(0,n,m`,n`) FOR n < 0
!
      bvec(0,:0,1) = 0
      if (lasym) bvec(0,:0,2) = 0
!
!     M = 0 MODES
!
      amatrix(:nf*mf1+1:mf1,:,1) = 0
!
!     ADD DIAGONAL ELEMENT TO AMATRIX
!
      do mn = 1, mnpd
         amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3
      end do

      if (lasym) then
         amatrix(:nf*mf1+1:mf1,:,2) = 0
         amatrix(:nf*mf1+1:mf1,:,3) = 0
         amatrix(:nf*mf1+1:mf1,:,4) = 0
         do mn = 1, mnpd
            amatrix(mn,mn,4) = amatrix(mn,mn,4) + pi3
         end do
      end if
 
!
!       PUT ELEMENTS INTO SQUARE MATRIX
!
      amatsq(:mnpd,:mnpd) = amatrix(:,:,1)                      !Sin-Sin
 
      if (lasym) then
         amatsq(1+mnpd:mnpd*2,:mnpd) = amatrix(:,:,2)           !Sin-Cos
         amatsq(:mnpd,1+mnpd:mnpd*2) = amatrix(:,:,3)           !Cos-Sin
         amatsq(1+mnpd:mnpd*2,1+mnpd:mnpd*2) = amatrix(:,:,4)   !Cos-Cos
      end if
 
      end subroutine fouri

      
      subroutine fourp (grpmn, grp, istore, istart, iend)
      use vacmod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: istart, iend, istore
      real(rprec), intent(in) :: grp(nuv,istore)
      real(rprec) :: grpmn(nuv2,0:mf,-nf:nf,*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, kv, ku, ip, iu, m
      real(rprec), allocatable, dimension(:,:,:,:) :: g1, g2
      real(rprec), allocatable :: kernels(:), kernelc(:), 
     1     gcos(:,:),gsin(:,:)
      real(rprec) :: cosm, sinm, cosn, sinn
C-----------------------------------------------
!
!     PERFORM KV (TOROIDAL ANGLE) TRANSFORM
!     NOTE: THE m,n INDICES HERE CORRESPOND TO THE FIRST INDEX OF AMATRIX
!     NOTE: THE .5 FACTOR (IN COSN,SINN) ACCOUNTS FOR THE SUM IN KERNELM
!     ON ENTRY THE FIRST TIME, GRPMN IS SIN,COS * Kmn(analytic)
!
      allocate (g1(istore,nu2,0:nf,2), g2(istore,nu2,0:nf,2),
     1   kernels(istore), kernelc(istore), gcos(istore,2), 
     2   gsin(istore,2), stat = m)
      if (m .ne. 0) stop 'Allocation error in fourp'
      
      g1 = 0
      g2 = 0
      
      do 10 n = 0,nf
        do 10 kv = 1,nv
          cosn = 0.5_dp*onp*cosv(n,kv)
          sinn = 0.5_dp*onp*sinv(n,kv)
          iu = kv
          do ku = 1,nu2
            do ip = 1,istore
              kernels(ip) = grp(iu,ip) - grp(imirr(iu),ip)        !sin symmetry
              g1(ip,ku,n,1) = g1(ip,ku,n,1) + cosn*kernels(ip)
              g2(ip,ku,n,1) = g2(ip,ku,n,1) + sinn*kernels(ip)
              if (lasym) then
              kernelc(ip) = grp(iu,ip) + grp(imirr(iu),ip)        !cos symmetry
              g1(ip,ku,n,2) = g1(ip,ku,n,2) + cosn*kernelc(ip)
              g2(ip,ku,n,2) = g2(ip,ku,n,2) + sinn*kernelc(ip)
              end if
            end do
            iu = iu + nv
          end do   
 10   continue

!
!     PERFORM KU (POLOIDAL ANGLE) TRANFORM
!
      do 30 m = 0,mf
        do 30 ku = 1,nu2
          cosm = -cosui(m,ku)
          sinm =  sinui(m,ku)
          do 30 n= 0,nf
            do ip = 1,istore
              gcos(ip,1) = g1(ip,ku,n,1)*sinm
              gsin(ip,1) = g2(ip,ku,n,1)*cosm
              grpmn(ip+istart,m,n,1) = grpmn(ip+istart,m,n,1)
     1        + gcos(ip,1) + gsin(ip,1)
            end do

            if (n .ne. 0) then
              do ip = 1,istore
                grpmn(ip+istart,m,-n,1) = grpmn(ip+istart,m,-n,1)
     1        + gcos(ip,1) - gsin(ip,1)
              end do
            endif  

            if (.not.lasym) cycle

            do ip = 1, istore
              gcos(ip,2) =-g1(ip,ku,n,2)*cosm
              gsin(ip,2) = g2(ip,ku,n,2)*sinm
              grpmn(ip+istart,m,n,2) = grpmn(ip+istart,m,n,2)
     1           + gcos(ip,2) + gsin(ip,2)
            end do

            if (n .ne. 0) then
              do ip = 1,istore
                 grpmn(ip+istart,m,-n,2) = grpmn(ip+istart,m,-n,2)
     1           + gcos(ip,2) - gsin(ip,2)
              end do
            end if  
 30   continue

      istart = iend

      deallocate (g1, g2, kernels, kernelc, gcos, gsin, stat = m)

      end subroutine fourp

      
      subroutine greenf(delgr, delgrp, ip)
      use vacmod
      use vparams, only: one, zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ip
      real(rprec), dimension(nuv), intent(out) :: delgr, delgrp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(2) :: ilow, ihigh
      integer :: ivoff, iskip, iuoff, i, kp, nloop
      real(rprec), dimension(:), allocatable ::
     1    ftemp, gsave, htemp, ga1, ga2, dsave
      real(rprec):: z2p, rcosip, rsinip, cosp, sinp,
     1    sxsave, sysave
C-----------------------------------------------
!
!       ON EXIT, DELGR IS THE DIFFERENCE OF "GREEN'S FUNCTION"
!       AND ANALYTIC APPROXIMATION, SUMMED OVER FIELD PERIODS
!       DELGRP IS DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION"
!       AND ANALYTIC APPROXIMATION
!
!       COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS
!
      allocate (ftemp(nuv), gsave(nuv), htemp(nuv), ga1(nuv), ga2(nuv),
     1          dsave(nuv), stat=i)
      if (i .ne. 0) stop 'allocation error in greenf'
            
      ilow(1) = 1
      ilow(2) = ip + 1
      ihigh(1) = ip - 1
      ihigh(2) = nuv
      ivoff = nuv + 1 - ip
      iskip = (ip - 1)/nv
      iuoff = nuv - nv*iskip
      z2p = -2*z1b(ip)
      rcosip = -2*rcosuv(ip)
      rsinip = -2*rsinuv(ip)
      delgr(ip)  = zero
      delgrp(ip) = zero
!
!     INITIALIZE ANALYTIC APPROXIMATIONS AND COMPUTE FIELD-PERIOD
!     INVARIANT VECTORS
!
      do i = 1, nuv
         ga1(i) = tanu(i + iuoff)*(guu_b(ip)*tanu(i+iuoff)+guv_b(ip)*
     1      tanv(i+ivoff)) + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
         ga2(i) = tanu(i + iuoff)*(auu(ip)*tanu(i+iuoff)+auv(ip)*
     1      tanv(i+ivoff)) + avv(ip)*tanv(i+ivoff)*tanv(i+ivoff)
         gsave(i) = rb2(ip) + rb2(i) + z1b(i)*z2p
         dsave(i) = drv(ip) + z1b(i)*snz(ip)
      end do
!
!     SUM OVER FIELD-PERIODS
!
      do kp = 1, nfper
         cosp = rcosip*cosper(kp) + rsinip*sinper(kp)
         sinp = rsinip*cosper(kp) - rcosip*sinper(kp)
         sxsave = -0.5_dp*(snr(ip)*cosp-snv(ip)*sinp)/r1b(ip)
         sysave = -0.5_dp*(snr(ip)*sinp+snv(ip)*cosp)/r1b(ip)
         if (kp .le. 1) then
            do nloop = 1, 2
               do i = ilow(nloop), ihigh(nloop)
                 ga2(i) = ga2(i)/ga1(i)
                 ga1(i) = one/sqrt(ga1(i))
                 ftemp(i) = one/(gsave(i) + cosp*rcosuv(i) 
     1                    + sinp*rsinuv(i))
                 htemp(i) = sqrt(ftemp(i))
                 delgrp(i) = -ga2(i)*ga1(i) + ftemp(i)*htemp(i)*
     1              (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
                 delgr(i) = htemp(i) - ga1(i)
               end do
            end do
         else
            do i = 1,nuv        
              ftemp(i) =one/(gsave(i) + cosp*rcosuv(i) + sinp*rsinuv(i))
              htemp(i) = sqrt(ftemp(i))
              delgrp(i) = delgrp(i) + ftemp(i)*htemp(i)*
     1           (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
              delgr(i) = delgr(i) + htemp(i)
           end do  
         endif
      end do

      deallocate (ftemp, gsave, htemp, ga1, ga2, dsave, stat=i)

      end subroutine greenf

      
      subroutine scalpot(bvec, amatrix, wint, ns, istore_max, ivacskip)
      use vacmod
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ns, ivacskip, istore_max
      real(rprec), intent(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
      real(rprec), intent(in) :: wint(*)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ip, istore, istart
      real(rprec), allocatable :: grpmn(:), green(:), gstore(:)
      real(rprec), allocatable :: greenp(:,:)
C-----------------------------------------------
      if (.not.allocated(amatsav)) 
     1   stop 'AMATSAV: Allocation error in scalpot' 

      allocate (grpmn(nuv2*mnpd2), stat=ip)
      if (ip .ne. 0) stop 'GRPMN: Allocation error in scalpot'
!
!     INITIALIZE VECTORS
!
      if (ivacskip .eq. 0) then
         amatrix = zero
         grpmn   = zero
      endif

      bvec = zero
!
!       COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!       ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!       AND GRPMN CONTAINS SIN,COS * TRANSFORM OF NORMAL DERIVATIVE
!       OF THE "GREEN'S FUNCTION"
!
      call analyt (grpmn, bvec, ivacskip)

      if (ivacskip .ne. 0) then
         bvec = bvec + bvecsav
      else
         allocate (green(nuv), gstore(nuv), greenp(nuv,istore_max))
         bvecsav = bvec
         gstore  = zero
!
!       COMPUTE SURFACE INTEGRALS OF SOURCE, "GREEN'S FUNCTION" NEEDED
!       FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
!       NOTE: SOURCE IS THE RHS, KERNEL IS THE LHS OF EQN.
!
         istart = 0
         do ip = 1, nuv2
            istore = 1 + mod(ip-1,istore_max)
!
!       COMPUTE EXACT AND APPROXIMATE "GREEN'S" FUNCTION AND GRADIENT
!
            call greenf (green, greenp(1,istore), ip)
            gstore = gstore + bexni(ip)*green
!
!       COMPUTE FOURIER INTEGRAL OF GRADIENT KERNEL ON UNPRIMED MESH
!
            if (istore.eq.istore_max .or. ip.eq.nuv2) 
     1        call fourp (grpmn, greenp, istore, istart, ip)

         end do
 
!
!       COMPUTE FOURIER INTEGRAL OF GRADIENT (SOURCE) PRIMED (UNPRIMED)
!
         call fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, ns)
         deallocate (green, greenp, gstore)
      endif

      deallocate (grpmn)
      amatrix = amatsav

      if (ivacskip .ne. 0) return 
!
!     SAVE NON-SINGULAR CONTRIBUTION TO BVEC IN BVECSAV
!
      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

      end subroutine scalpot

      
      subroutine surface(rc, rs, zs, zc, xm, xn, mnmax)
      use vacmod
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mnmax
      real(rprec), dimension(mnmax) :: rc, rs, zs, zc, xm, xn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
      integer :: i, mn, m, n, n1
      real(rprec), allocatable, dimension(:) :: 
     1   ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1
C-----------------------------------------------
!
!       THIS ROUTINE COMPUTES THE SURFACE VALUES OF R,Z AND DERIVATIVES
!
!
!       Compute R & Z (and their derivatives) on surface
!
      allocate (ruu(nuv2), ruv(nuv2), rvv(nuv2), zuu(nuv2), zuv(nuv2),
     1          zvv(nuv2), cosmn1(nuv2), sinmn1(nuv2), stat = i)
      if (i .ne. 0) stop 'Allocation error in SURFACE'
       
      r1b = 0;   rub = 0;   rvb = 0;  ruu = 0; ruv = 0; rvv = 0
      z1b = 0;   zub = 0;   zvb = 0;  zuu = 0; zuv = 0; zvv = 0
      do mn = 1, mnmax
         m = nint(xm(mn))
         n = nint(xn(mn)/(nfper))
         n1 = abs(n)
         cosmn1(:) = cosu1(:,m)*cosv1(:,n1) + csign(n)*sinu1(:,m)*
     1               sinv1(:,n1)
         sinmn1(:) = sinu1(:,m)*cosv1(:,n1) - csign(n)*cosu1(:,m)*
     1               sinv1(:,n1)
         do i = 1, nuv2
            r1b(i) = r1b(i) + rc(mn) * cosmn1(i)
            rub(i) = rub(i) - xm(mn) * rc(mn) * sinmn1(i)
            rvb(i) = rvb(i) + xn(mn) * rc(mn) * sinmn1(i)
            z1b(i) = z1b(i) + zs(mn) * sinmn1(i)
            zub(i) = zub(i) + xm(mn) * zs(mn) * cosmn1(i)
            zvb(i) = zvb(i) - xn(mn) * zs(mn) * cosmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rc(mn) * cosmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rc(mn) * cosmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rc(mn) * cosmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zs(mn) * sinmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zs(mn) * sinmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zs(mn) * sinmn1(i)
         end do
 
         if (.not.lasym) cycle

         do i = 1, nuv2
            r1b(i) = r1b(i) + rs(mn) * sinmn1(i)
            rub(i) = rub(i) + xm(mn) * rs(mn) * cosmn1(i)
            rvb(i) = rvb(i) - xn(mn) * rs(mn) * cosmn1(i)
            z1b(i) = z1b(i) + zc(mn) * cosmn1(i)
            zub(i) = zub(i) - xm(mn) * zc(mn) * sinmn1(i)
            zvb(i) = zvb(i) + xn(mn) * zc(mn) * sinmn1(i)
            ruu(i) = ruu(i) - xm(mn)*xm(mn)*rs(mn) * sinmn1(i)
            ruv(i) = ruv(i) + xm(mn)*xn(mn)*rs(mn) * sinmn1(i)
            rvv(i) = rvv(i) - xn(mn)*xn(mn)*rs(mn) * sinmn1(i)
            zuu(i) = zuu(i) - xm(mn)*xm(mn)*zc(mn) * cosmn1(i)
            zuv(i) = zuv(i) + xm(mn)*xn(mn)*zc(mn) * cosmn1(i)
            zvv(i) = zvv(i) - xn(mn)*xn(mn)*zc(mn) * cosmn1(i)
         end do
      end do

!
!       COMPUTE METRIC COEFFICIENTS AND AREA ELEMENTS
!       NOTE: guv = .5*np GUV; gvv = np*np* GVV, where GUV, GVV are the
!            real metric elements
!
      do i = 1,nuv2
        guu_b(i) = rub(i)*rub(i) + zub(i)*zub(i)
        guv_b(i) = (rub(i)*rvb(i)+ zub(i)*zvb(i))*onp*2.0_dp
        gvv_b(i) = (rvb(i)*rvb(i)+ zvb(i)*zvb(i)+(r1b(i)*r1b(i)))*onp2
        rb2(i) = (r1b(i)*r1b(i)) + z1b(i)*z1b(i)
        snr(i) = -r1b(i)*zub(i)
        snv(i) = rvb(i)*zub(i) - rub(i)*zvb(i)
        snz(i) = r1b(i)*rub(i)
        drv(i) = (-r1b(i)*snr(i)) - z1b(i)*snz(i)
        auu(i) = (0.5_dp*r1b(i))*(zuu(i)*rub(i) - ruu(i)*zub(i))
        auv(i) = (snv(i)*rub(i)+(zuv(i)*rub(i) 
     1         - ruv(i)*zub(i))*r1b(i))*onp
        avv(i) = (snv(i)*rvb(i)+(0.5_dp*r1b(i))*(zvv(i)*rub(i)
     1         - zub(i)*rvv(i) - snr(i)))*onp2
      end do

      if (.not.lasym) then
         do i = 1 + nv, nuv2 - nv
            rb2(imirr(i)) = rb2(i)
            r1b(imirr(i)) = r1b(i)
            z1b(imirr(i)) =-z1b(i)
         end do
      end if

      do i = 1,nuv
        rcosuv(i) = r1b(i)*cosuv(i)
        rsinuv(i) = r1b(i)*sinuv(i)
      end do  

      deallocate (ruu, ruv, rvv, zuu, zuv, zvv, cosmn1, sinmn1, stat=i)

      end subroutine surface

      
      subroutine tolicu(rcurr, zcurr)
      use vacmod
      use vacwires
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(nv) :: rcurr, zcurr
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, kper, kv
C-----------------------------------------------
 
      if (nv*nfper .ne. nvp) stop ' nvp!=nv*nfper in tolicu'
      i = 0
      do kper = 1, nfper
         do kv = 1, nv
            i = i + 1
            xw(i) = rcurr(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*
     1         sinuv(kv))
            yw(i) = rcurr(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*
     1         sinuv(kv))
            zw(i) = zcurr(kv)
         end do
      end do
      xw(nvp+1) = xw(1)
      yw(nvp+1) = yw(1)
      zw(nvp+1) = zw(1)
      do i = 1,nvp
        dx(i) = xw(i+1) - xw(i)
        dy(i) = yw(i+1) - yw(i)
        dz(i) = zw(i+1) - zw(i)
        vx(i) = yw(i)*dz(i) - zw(i)*dy(i)
        vy(i) = zw(i)*dx(i) - xw(i)*dz(i)
        vz(i) = xw(i)*dy(i) - yw(i)*dx(i)
      end do  

      end subroutine tolicu

      
      subroutine precal
      use vparams, only: zero, one, epstan
      use vacmod
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: kp, ku, kuminus, kv, kvminus, i, m, n, mn, 
     1   jmn, kmn, l, istat1
      real(rprec), dimension(0:mf + nf,0:mf,0:nf) :: cmn
      real(rprec) :: argu, argv, dn1, smn, f1, f2, f3
C-----------------------------------------------
!
!     THIS ROUTINE COMPUTES INITIAL CONSTANTS AND ARRAYS
!
      pi2 = 8*atan(one)
      pi3 = 0.5_dp*pi2**3
      pi4 = 2*pi2
      alp = pi2/(nfper)
      alu = pi2/(nu)
      alv = pi2/(nv)
      onp = one/(nfper)
      onp2 = onp*onp
      alvp = onp*alv

!
!     ALLOCATE PERSISTENT ARRAYS. DEALLOCATED IN FILEOUT ROUTINE
!
      allocate( tanu(2*nuv), tanv(2*nuv), sin2v(2*nuv),
     1     sinper(nfper), cosper(nfper), sinuv(nuv), cosuv(nuv),
     2     sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),
     3     cosv(-nf:nf,nv), sinui(0:mf,nu), cosui(0:mf,nu),
     4     cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),
     5     sinu1(nuv2,0:mf), cosu1(nuv2,0:mf),
     6     sinv1(nuv2,0:nf), cosv1(nuv2,0:nf), imirr(nuv), 
     7     xmpot(mnpd), xnpot(mnpd), stat=istat1)
      if (istat1.ne.0) stop 'allocation error in precal'


!
!       IMIRR(I) GIVES THE INDEX OF THE POINT TWOPI-THETA(I),TWOPI-ZETA(I)
!
      do kp = 1, nfper
         cosper(kp) = cos(alp*(kp - 1))
         sinper(kp) = sin(alp*(kp - 1))
      end do
      do ku = 1, nu
         kuminus = mod(nu + 1 - ku,nu) + 1
         do kv = 1, nv
            kvminus = mod(nv + 1 - kv,nv) + 1
            i = kv + nv*(ku - 1)
            imirr(i) = kvminus + nv*(kuminus - 1)
            cosuv(i) = cos(alvp*(kv - 1))
            sinuv(i) = sin(alvp*(kv - 1))
         end do
      end do
!
!       NOTE: ACTUAL ANGLE DIFFERENCE IS (KUP-1) - (KU-1)
!
      i = 0
      do ku = 1, 2*nu
         argu = 0.5_dp*alu*(ku - 1)
         do kv = 1, nv
            i = i + 1
            argv = 0.5_dp*alv*(kv - 1)
            sin2v(i) = sin(argv)*sin(argv)
            if (abs(argu - 0.25_dp*pi2)<epstan .or.
     1      abs(argu - 0.75_dp*pi2) < epstan) then
               tanu(i) = 0.9e30_dp
            else
               tanu(i) = 2.0_dp*tan(argu)
            endif
            if (abs(argv - 0.25_dp*pi2) < epstan) then
               tanv(i) = 0.9e30_dp
            else
               tanv(i) = 2.0*tan(argv)
            endif
         end do
      end do
      do m = 0, mf
         l40: do ku = 1, nu
            cosu(m,ku) = cos(alu*(m*(ku - 1)))
            sinu(m,ku) = sin(alu*(m*(ku - 1)))
            cosui(m,ku) = cosu(m,ku)*alu*alv*2.0
            sinui(m,ku) = sinu(m,ku)*alu*alv*2.0
            if (ku.eq.1 .or. ku.eq.nu2) cosui(m,ku) = 0.5_dp*cosui(m,ku)
            do kv = 1, nv
               i = kv + nv*(ku - 1)
               if (i > nuv2) cycle  l40
               cosu1(i,m) = cosu(m,ku)
               sinu1(i,m) = sinu(m,ku)
            end do
         end do l40
      end do
      do n = -nf, nf
         dn1 = alvp*(n*nfper)
         csign(n) = sign(one,dn1)
         l50: do ku = 1, nu
            do kv = 1, nv
               i = kv + nv*(ku - 1)
               cosv(n,kv) = cos(dn1*(kv - 1))
               sinv(n,kv) = sin(dn1*(kv - 1))
               if (i.gt.nuv2 .or. n.lt.0) cycle  l50
               cosv1(i,n) = cosv(n,kv)
               sinv1(i,n) = sinv(n,kv)
            end do
         end do l50
      end do
      mn = 0
      do n = -nf, nf
         do m = 0, mf
            mn = mn + 1
            xmpot(mn) = (m)
            xnpot(mn) = (n*nfper)
         end do
      end do
!
!       COMPUTE "CMN'S" AND THEIR SUMS , EQ (A14 AND A13) IN J.COMP.PHYS PAPER
!
      do m = 0, mf
         do n = 0, nf
            jmn = m + n
            kmn = abs(m - n)
            smn = 0.5_dp*(jmn + kmn)
            f1 = 1.0_dp
            f2 = 1.0_dp
            f3 = 1.0_dp
            do i = 1, kmn
               f1 = f1*(smn + (1 - i))
               f2 = f2*(i)
            end do
            cmn(0:mf+nf,m,n) = zero
            do l = kmn, jmn, 2
               cmn(l,m,n) = f1/(f2*f3)*((-1)**((l - m + n)/2))
               f1 = f1*0.25_dp*((jmn + l + 2)*(jmn - l))
               f2 = f2*0.5_dp*(l + 2 + kmn)
               f3 = f3*0.5_dp*(l + 2 - kmn)
            end do
         end do
      end do
!
!       The ALP factor comes from integral over field periods
!
        do m = 1,mf
           do n = 1,nf
              cmns(0:mf+nf,m,n) = 0.5_dp*alp*(cmn(0:mf+nf,m,n) +
     1        cmn(0:mf+nf,m-1,n) + cmn(0:mf+nf,m,n-1) + 
     2        cmn(0:mf+nf,m-1,n-1))
           end do
        end do   
      cmns(0:mf+nf,1:mf,0) = (0.5_dp*alp)*(cmn(0:mf+nf,1:mf,0)
     1                       + cmn(0:mf+nf,:mf-1,0))
      cmns(0:mf+nf,0,1:nf) = (0.5_dp*alp)*(cmn(0:mf+nf,0,1:nf)
     1                    + cmn(0:mf+nf,0,:nf-1))
      cmns(0:mf+nf,0,0)    = (0.5_dp*alp)*(cmn(0:mf+nf,0,0)
     1                    + cmn(0:mf+nf,0,0)) 

      end subroutine precal
