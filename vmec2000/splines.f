      subroutine add_tension(amat, wten, hx, tens, tensv, fpoly, n, nb, 
     1   ioff, nmat)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, nb, ioff, nmat
      real(rprec) tens, tensv, fpoly
      real(rprec), dimension(nmat,*) :: amat
      real(rprec), dimension(nmat) :: wten
      real(rprec), dimension(n) :: hx
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k, i, koff
      real(rprec), dimension(n) :: wten0, tenshx, work
      real(rprec) :: delta_x, tension
C-----------------------------------------------
 
!
!       THIS SUBROUTINE ADDS THE TENSION TERM WTEN*(DEL(K) - DEL(K-1))
!       TO THE REST OF THE CHI-SQ AMAT, WHERE DEL(K) = (G(K+1)-G(K))/H(K)
!       AND G(K) IS THE SECOND DERIVATIVE AT THE K-TH KNOT
!
!       IOFF ALLOWS FOR ADDING TENSION SEPARATELY TO IOTA (IOFF=0)
!       AND PRESSURE (IOFF=NIOTA) SPLINES
!
!       tens:   spline tension (optionally, at the 1st pt only;see note)
!       tensv:  vbl spline tension for n-1th point (see note)
!       fpoly:  vbl spline tension form factor (note: if tens1<>tens2
!               then tension(i-th point) = tens+(tensv-tens)*(i/n-1))**fpoly)
 
!
!       BOUNDS CHECKING
!
      if (n + ioff > nmat) stop '(n+ioff>nmat)'
      if (fpoly < 0.) stop '(fpoly<0)'
      if (n < 1) stop '(n < 1)'
 
!
!       COMPUTE TENSION COEFFICIENTS
!
      
      delta_x = sum(hx(:n-1))
      tension = 0.5*(delta_x/(n))**3
      if (fpoly.eq.0. .or. tens.eq.tensv) then
         tenshx(:n-1) = tens
      else
         do i = 1, n - 1
            tenshx(i) = tens + (tensv - tens)*(real(i - 1,rprec)/
     1         (n - 1))**fpoly
         end do
      endif
 
      do i = 1,n-1
        tenshx(i) = tension * tenshx(i) / hx(i)
        work(i) = hx(i)*(wten(i+ioff) + wten(i+ioff+1))
      enddo
      do i = 2,n-1
        wten0(i) = 0.5 * ( work(i) + work(i-1) )/(hx(i) + hx(i-1))
      enddo
      wten0(1) = wten0(2)
      wten0(n) = wten(n+ioff)
!
!       COMPUTE, FOR K = 1,N, B(K,L)*JACOBIAN(L,I) = W(K,I),
!       WHERE JACOBIAN = D[G]/D[F] and B is TRIDIAGONAL
!       SEE EQN(27) IN PHYS.PLASMAS 1, p 2277.
!
      do k = 1, n
         koff = k + ioff
         work(:n-1) = 0
!       SET UP COEFFICIENTS IN [G(K+1)-G(K)]/h(k) - [G(K)-G(K-1)]/h(k-1)
         if (k .eq. 1) then
            work(2) = tenshx(1)*wten0(2)
            work(1) = -work(2)
         else if (k .eq. n) then
            work(n-1) = tenshx(n-1)*wten0(n-1)
         else
            work(k-1) = tenshx(k-1)*wten0(k-1)
            work(k) = -(tenshx(k)+tenshx(k-1))*wten0(k)
            work(k+1) = tenshx(k)*wten0(k+1)
         endif
         if (nb .eq. natur) work(1) = 0
         work(n) = 0
!
!       COMPUTE work(j) = work(i)*Jacobian(i,j) and add to amat(k,j)
!
         call jacprod (work, hx, n, nb)
         amat(koff,1+ioff:n+ioff) = amat(koff,1+ioff:n+ioff) + work(:n)
      end do
 
      end subroutine add_tension
      
      subroutine getspline(amat, splnots, hk, delse, hs, indexs, isort, 
     1   ndata0, nots)
      use vsvd0
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ndata0, nots
      real(rprec) hs
      integer, dimension(ndata0) :: indexs, isort
      real(rprec), dimension(nots,ndata0) :: amat
      real(rprec), dimension(nots) :: splnots, hk
      real(rprec), dimension(ndata0) :: delse
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(nots) :: nk
      integer :: i, js, ia, k, k1, ib, j, nb
      real(rprec), dimension(ndata0) :: w, w1, u, u1, snodes
      real(rprec), dimension(nots,ndata0) :: bmat
C-----------------------------------------------
 
!
!       ON EXIT, AMAT = AMAT + BMAT, WHERE BMAT IS THE 2nd
!       DERIVATIVE COEFFICIENT MATRIX ARRAY, MULTIPLIED BY JACOBIAN
!       AND AMAT (ON RHS) WAS FUNCTION COEFFICIENT MATRIX ARRAY
!
 
!
!       SORT KNOT POSITIONS IN ASCENDING ORDER IN S-SPACE
!       USE SQRT(S) KNOT POSITIONS FOR IMPROVED RESOLUTION
!       NOTE: SNODES(I) IS THE VALUE OF SQRT(S) CORRESPONDING TO
!       THE MESH VALUES CORRESPONDING TO DELSE, INDEXS (COMPUTED OUTSIDE
!       THIS PROGRAM)
!
 
      do i = 1, ndata0
         js = indexs(i)
         snodes(i) = sqrt(hs*((js - 1) + delse(i)))
!        if (snodes(i) .le. zero) snodes(i) = epstan
      end do

!     Avoid roundoff error in SPLININT
      if( snodes(ndata0) .gt. splnots(nots) )
     1    snodes(ndata0) = splnots(nots)
 
      call sort_data (snodes, isort, ndata0)
 
!
!       COMPUTE MATRIX COEFFICIENTS RELATING SPLINE AT SPLNOTS
!       TO REAL-SPACE FUNCTION AT SORTED MESH POINTS RMESH
!
      amat(:nots,:ndata0) = zero
      bmat(:nots,:ndata0) = zero
 
!
!       SETUP SPLINE PARAMETERS AT EACH TIME STEP, SINCE SNODES
!       MAY BE CHANGING DYNAMICALLY IN S-SPACE
!
      call setup_int(splnots,snodes,hk,w,w1,u,u1,nk,nots,ndata0)
 
      ia = 1
      do k = 1, nots - 1
         if (nk(k) .gt. 0) then
            k1 = k + 1
            ib = ia + nk(k) - 1
            amat(k,ia:ib) = amat(k,ia:ib) + w(ia:ib)
            bmat(k,ia:ib) = bmat(k,ia:ib) + u(ia:ib)
            amat(k1,ia:ib) = amat(k1,ia:ib) + w1(ia:ib)
            bmat(k1,ia:ib) = bmat(k1,ia:ib) + u1(ia:ib)
            ia = ib + 1
         endif
      end do
 
      if (ib .ne. ndata0) stop 'ib!=ndat'
      nb = ideriv
 
      do j = 1, ndata0
         bmat(nots,j) = 0.
         call jacprod (bmat(1,j), hk, nots, nb)
         amat(:nots,j) = amat(:nots,j) + bmat(:nots,j)
      end do
 
      end subroutine getspline
      
      subroutine gety2(y, y2, h, nots, nb)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb
      real(rprec), dimension(*) :: y, y2, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer jmax
      real(rprec), dimension(nots) :: aspline, bspline, dspline
C-----------------------------------------------
 
      jspmin(1) = 2
      if (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      y2(1) = 0.
      y2(nots) = 0.
      if (nb .eq. ideriv) y2(1) = 6.0*(y(2)-y(1))/h(1)
      aspline(2:jmax) = h(2:jmax)
      bspline(2:jmax) = h(:jmax-1)
      dspline(2:jmax) = 2.0*(h(2:jmax)+h(:jmax-1))
      y2(2:jmax) = 6.0*((y(3:jmax+1)-y(2:jmax))/h(2:jmax)
     1   -(y(2:jmax)-y(:jmax-1))/h(:jmax-1))
 
      call tridslv(aspline,dspline,bspline,y2,jspmin,jmax,0,nots,1)
 
      end subroutine gety2
      
      subroutine initspline(amat, splnot, h, weight, nots)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots
      real(rprec), dimension(nots,nots) :: amat
      real(rprec), dimension(*) :: splnot, h, weight
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer i
      real(rprec) :: eps
C-----------------------------------------------
 
      if (nots .lt. 3) stop 'nots<3'
      eps = 1.0/(splnot(nots)-splnot(1))
 
      amat = 0.
      do i = 1, nots
         amat(i,i) = weight(i)
      end do
 
      do i = 1, nots - 1
         h(i) = splnot(i+1) - splnot(i)
         if (eps*h(i) .le. 1.e-8_dp) stop 'h(i)<1.e-8'
      end do
 
      end subroutine initspline
      
      subroutine jacprod(c, h, nots, nb)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb, jmax
      real(rprec), dimension(*) :: c, h
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec), dimension(nots) :: 
     1   aspline, bspline, dspline, dum1
C-----------------------------------------------
 
!
!       THIS ROUTINE COMPUTES THE INNER PRODUCT COUT(I) = CIN(J)*JACOBIAN(J,I)
!       WHERE JACOBIAN(J,I) = D[G(J)]/D[F(I)]
!       HERE, G(J) ARE SECOND-DERIVATIVE KNOTS, F(I) FUNCTION KNOTS
!
!       COMPUTE COEFFICIENT ARRAY ELEMENTS A*X(I+1) + D*X(I) + B*X(I-1)
!       (TO BE SAFE, RECOMPUTE EACH TIME, SINCE IOTA, P SPLINES MAY
!        DIFFER FROM CALL TO CALL)
!
      aspline(1) = h(1)
      dspline(1) = 2.0*h(1)
      aspline(2:nots-1) = h(2:nots-1)
      bspline(2:nots-1) = h(:nots-2)
      dspline(2:nots-1) = 2.0*(h(2:nots-1)+h(:nots-2))
 
      jspmin(1) = 2
      if (nb .eq. ideriv) jspmin(1) = 1
      jmax = nots - 1
      call tridslv(aspline,dspline,bspline,c,jspmin,jmax,0,nots,1)
      dum1(1) = 6.0*(c(2)-c(1))/h(1)
      dum1(2:nots) = 6.0*(c(:nots-1)-c(2:nots))/h(:nots-1)
      c(2:nots-1) = dum1(2:nots-1) - dum1(3:nots)
      c(1) = dum1(1)
      c(nots) = dum1(nots)

      end subroutine jacprod
      
      subroutine set_dual(data, hi, yi, yi2, hp, yp, yp2, wten, alsq, 
     1   niota, npres, nots)
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer niota, npres, nots
      real(rprec), dimension(nots) :: data
      real(rprec), dimension(niota) :: hi, yi, yi2
      real(rprec), dimension(npres) :: hp, yp, yp2
      real(rprec), dimension(nots) :: wten
      real(rprec), dimension(nots,nots) :: alsq
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer nb, ioff
C-----------------------------------------------
!
!       ADD TENSION TO DIAGONAL BLOCKS
!
      nb = ideriv
      ioff = 0
      call add_tension (alsq, wten, hi, tensi, tensi2, fpolyi, niota, nb
     1   , ioff, nots)
      call add_tension (alsq, wten, hp, tensp, zero, zero, npres, nb, 
     1   niota, nots)
 
!
!       FREEZE EDGE PRESSURE IF NO PRESSURE SPECIED
!       COMPUTE SOLUTION FOR SPLINES
!
      call solver (alsq, data, nots)
      yi(:niota) = data(:niota)
      yp(:npres) = data(1+niota:npres+niota)
 
!
!       COMPUTE SECOND DERIVATIVES
!
      call gety2 (yi, yi2, hi, niota, nb)
      call gety2 (yp, yp2, hp, npres, nb)
 
      end subroutine set_dual
      
      subroutine setspline(x,weight,y,h,yfit,y2,wten,tens,nots,nb)
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nb
      real(rprec) tens
      real(rprec), dimension(*) :: x, weight
      real(rprec), dimension(nots) :: y
      real(rprec), dimension(*) :: h
      real(rprec), dimension(nots) :: yfit
      real(rprec), dimension(*) :: y2, wten
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ioff, istat
      real(rprec) :: alsq(nots,nots)
C-----------------------------------------------
!
!       x:        independent coordinate array
!       y:        dependent y(x) array
!       yfit:     fitted values (under tension) to y array
!       h:        x(i+1) - x(i) array
!       y2:       y'' array used for splines
!       wten:        weight array for tension (changed on exit)
!       alsq:        matrix elements for least squares fit (from s-integrations)
!       nots:        number of independent coordinates (knots)
!       nb:        = NATUR, use natural boundary condition at left knot
!                  = IDERIV, use derivative (dy/dx =0) boundary condition at left knot
 
!
!       IT IS ASSUMED THAT X,Y,WTEN ARE ALL SORTED (ON X(I) < X(I+1))
!
 
!
!       INITIALIZE ALSQ TO ZERO, COMPUTE H ELEMENTS
!
      call initspline (alsq, x, h, weight, nots)
 
!
!       SET UP SPLINE MATRIX ASPLINE AND NON-DIMENSIONLIZE TENSION
!
      ioff = 0
      call add_tension (alsq, wten, h, tens, zero, zero, nots, nb, 
     1   ioff, nots)
 
!
!       SOLVE FOR COEFFICIENTS
!
      yfit(:nots) = y(:nots)
      call solver (alsq, yfit, nots)
!
!       OBTAIN Y'' COEFFICIENTS AND STORE IN Y2
!
      call gety2 (yfit, y2, h, nots, nb)
 
      end subroutine setspline
      
      subroutine setup_int(xknots,smesh,hx,w,w1,u,u1,nk,nots,nmesh)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, nmesh
      integer, dimension(nots) :: nk
      real(rprec), dimension(nots) :: xknots
      real(rprec), dimension(nmesh) :: smesh
      real(rprec), dimension(nots) :: hx
      real(rprec), dimension(nmesh) :: w, w1, u, u1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: epstan = 1.d-10
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ksp1, k, i, k1
      real(rprec) :: smesh1, hk6
C-----------------------------------------------
!
!       FOR THE SPLINED FUNCTIONS CORRESPONDING TO XKNOTS (PRESSURE,
!       IOTA), THIS ROUTINE COMPUTES THE A AND B MATRIX ELEMENTS
!       (STORED IN W,U,W1,U1) WHICH ARE NEEDED TO EVALUATE THE FUNCTIONS
!       IN REAL-SPACE IN TERMS OF THEIR (VARIABLE) SPLINE KNOT VALUES.
!       THIS 'UNDOES' WHAT SPLINT ROUTINE DOES. LET Y(I) DENOTE THE
!       FUNCTION AT THE POINT SMESH(I) SUCH THAT
!
!                   XKNOTS(K) < SMESH(I) <= XKNOTS(K)
!
!       THEN,  Y(I) = W(I)*YK  + U(I)*GK  + W1(I)*YK1   + U1(I)*GK1
!
!       WHERE YK, GK ARE THE SPLINE AND 2ND DERIVATIVES AT KNOT K
!             YK1,GK1 ARE THE SAME AT KNOT K+1
!
      ksp1 = nots - 1
      smesh1 = smesh(1)
      if (smesh1 .le. xknots(1)) smesh(1) = xknots(1) + epstan
 
      nk = 0
 
      k = 1
      do i = 1, nmesh
  140    continue
         k1 = k + 1
!
!       XKNOTS = SQRT(HS*(JS-1)) DEFINED IN STARK,PRESSURE ROUTINE
!       (THIS CORRESPONDS TO APPROXIMATELY EQUAL SPACING ALONG MIDPLANE)
!
         if (smesh(i).gt.xknots(k) .and. smesh(i).le.xknots(k1)) then
            nk(k) = nk(k) + 1
            hk6 = hx(k)*hx(k)/6.0
            w1(i) = (smesh(i)-xknots(k))/hx(k)
            if (w1(i)<(-epstan) .or. w1(i)>1.0+epstan) stop 'w1(i)'
            w(i) = 1.0 - w1(i)
            u(i) = hk6*w(i)*(w(i)*w(i)-1.0)
            u1(i) = hk6*w1(i)*(w1(i)*w1(i)-1.0)
         else
            k = k + 1
            if (k .gt. ksp1) stop 'K>KSP1'
            go to 140
         endif
      end do
 
      smesh(1) = smesh1
 
      end subroutine setup_int
      
      subroutine sort_data (x, index_array, n)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
c-----------------------------------------------
      integer n
      integer, dimension(n) :: index_array
      real(rprec), dimension(n) :: x
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j
      integer, dimension(1) :: isamax
      real(rprec), dimension(n) :: dumx
C-----------------------------------------------
!
!       RETURNS INDEX(I) ARRAY, SO THAT X(INDEX(I)) IS SORTED, I=1,N
!       RETURNS Xin(INDEX(I)) = Xout(I)
!
 
      do i = n, 1, -1
         isamax = maxloc(abs(x))
         j = isamax(1)
         dumx(i) = x(j)
         x(j) = 0.
         index_array(i) = j
      end do
 
      x = dumx
 
      end subroutine sort_data
      
      subroutine splinint(grn, cm, jacob, h, u, u1, w, w1, nk, nots, 
     1   ifunc, nmesh)
      use vparams
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer nots, ifunc, nmesh
      integer, dimension(*) :: nk
      real(rprec), dimension(nmesh) :: grn, cm
      real(rprec), dimension(nots) :: jacob, h
      real(rprec), dimension(nmesh) :: u, u1, w, w1
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, k, nb, ia, ib, k1, ksp1, nmesh1
      real(rprec), dimension(nmesh) :: func
      real(rprec), dimension(nots) :: af, bs
C-----------------------------------------------
 
!
!       COMPUTES af,bs FACTORS IN
!       (ifunc=INTFUN)     Int[ GRN(X) * F(X) ds ]      or
!       (ifunc=INTDER)     Int[ GRN(X) * {d(CmF)/ds} ds ]
!                 =   af(k)*f(k) + bs(k)*f''(k)
!       WHERE f(k),f''(k) ARE SPLINE COEFFICIENTS OF F(X)
!
!       NOTE: FOR ifunc = INTDER, the OHS factor in CmF cancels
!             THE HS factor in the Integral.
!             FOR ifunc = INTFUN, GRN is assumed to be pre-multiplied
!             OUTSIDE this routine by HS factor
!       ALSO, COMPUTES af(k) + (SUM on i)bs(i)*J(i,k) = jacob(k),
!       WHERE J(i,k) = d[g(i)]/d[f(k)], g = f''
!
!       nk(k): Number of smesh-pts in k-th spline interval
!       xknots(k) < smesh <= xknots(k+1),   k = 1,nots-1
!
!       NOTE: The ifunc=INTDER case is done by integrating by parts,
!             so that the half-point integration (GRN at half mesh pts)
!             becomes a full-point integration in Cm*F.
!
      nb = ideriv           !Pressure, iota derivatives vanish at origin
      ksp1 = nots - 1
      nmesh1 = nmesh - 1
 
      if (ifunc .eq. intder) then
!
!       Integrate by parts (in s), func(1) and func(nmesh) are 'surface terms'
!
         func(1) = -cm(1)*grn(2)
         func(2:nmesh1) = cm(2:nmesh1)*(grn(2:nmesh1)-grn(3:nmesh1+1))
         func(nmesh) = cm(nmesh)*grn(nmesh)
      else
         func = grn
      endif
 
      af(:nots) = zero
      bs(:nots) = zero
 
      ia = 1
      do k = 1, ksp1
         if (nk(k) .ne. 0) then
            k1 = k + 1
            ib = ia + nk(k) - 1
            do j = ia,ib
              af(k)  = af(k)  + func(j)*w(j)
              bs(k)  = bs(k)  + func(j)*u(j)
              af(k1) = af(k1) + func(j)*w1(j)
              bs(k1) = bs(k1) + func(j)*u1(j)
            enddo
            ia = ib + 1
         endif
      end do
 
      if (ib .ne. nmesh) stop 'ib!=nmesh'
      if (nb .eq. natur) bs(1) = 0.           !Natural boundary conditions
      bs(nots) = 0.
      call jacprod (bs, h, nots, nb)         !Returns bs(i)=bs(j)*J(j,i)
      jacob(:nots) = af(:nots) + bs(:nots)
 
      end subroutine splinint
      
      subroutine splint(xa, ya, y2a, n, x, y, yp, ndim)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer n, ndim
      real(rprec), dimension(*) :: xa, ya, y2a, x, y, yp
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0, c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, i, k
      real(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, deriv
C-----------------------------------------------
!
!       SPLINE INTERPOLATION ROUTINE (Numerical Recipes, pg. 89)
!       XA: ordered array of length N of ordinates at which function YA=F(XA)
!           is tabulated
!       YA: array of length N , = F(XA)
!       Y2A: array of second derivatives at XA points
!       computed from call to SPLINE
!       X : value at which Y = F(X) is to be computed from splines
!       YP = dY/dX at X
!       NDIM: dimension of X, Y, YP arrays
 
 
      deriv = yp(1)
      klo = 1
      khi = n
      do i = 1, ndim
         do while(khi - klo .gt. 1)
            k = (khi + klo)/2
            if (xa(k) .gt. x(i)) then
               khi = k
            else
               klo = k
            endif
         end do
 
         h = xa(khi) - xa(klo)
         a = xa(khi) - x(i)
         b = x(i) - xa(klo)
         h2 = h*h
         a2 = a*a
         b2 = b*b
         y26lo = c1o6*y2a(klo)
         y26hi = c1o6*y2a(khi)
         y(i) = (a*(ya(klo)+(a2-h2)*y26lo)+b*(ya(khi)+(b2-h2)*y26hi))/h
         if (deriv .ne. zero) yp(i) = (ya(khi)-ya(klo)+y26hi*(3.0*b2-h2)
     1      -y26lo*(3.0*a2-h2))/h
         if (i.lt.ndim .and. x(i+1).gt.x(i)) then
            khi = n
         else
            klo = 1
         endif
      end do
 
      end subroutine splint
      
      function splints (x)
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, k
      real(rprec) :: h, a, a2, b, b2, h2, y26lo, y26hi, yx, splints
C-----------------------------------------------
 
      klo = 1
      khi = isnodes
 
    1 continue
      if (khi - klo .gt. 1) then
         k = (khi + klo)/2
         if (sknots(k) .gt. x) then
            khi = k
         else
            klo = k
         endif
         go to 1
      endif
 
      h = sknots(khi) - sknots(klo)
      a = sknots(khi) - x
      b = x - sknots(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2stark(klo)
      y26hi = c1o6*y2stark(khi)
      yx = (a*(ystark(klo)+(a2-h2)*y26lo)+b*(ystark(khi)+(b2-h2)*y26hi))
     1   /h
      splints = yx

      end function splints
       
      function splintx(x)
      use vparams
      use vsvd
      use csplinx
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) x
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1o6 = 1._dp/6._dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: klo, khi, k
      real(rprec) :: h, a, b, h2, a2, b2, y26lo, y26hi, qmidx0,
     1   splintx      
C-----------------------------------------------
 
 
       call setspline (rmidx, wmidx, qmidx, hmidx, ymidx, y2midx, 
     1    tenmidx, tenmidx(1), nptsx, natur)
 
      klo = 1
      khi = nptsx
 
    1 continue
      if (khi - klo .gt. 1) then
         k = (khi + klo)/2
         if (rmidx(k) .gt. x) then
            khi = k
         else
            klo = k
         endif
         go to 1
      endif
 
      h = rmidx(khi) - rmidx(klo)
      if( h.eq.zero )then
        splintx = zero
        return
      end if  
      a = rmidx(khi) - x
      b = x - rmidx(klo)
      h2 = h*h
      a2 = a*a
      b2 = b*b
      y26lo = c1o6*y2midx(klo)
      y26hi = c1o6*y2midx(khi)
      qmidx0 = (a*(ymidx(klo)+(a2-h2)*y26lo)+b*(ymidx(khi)+
     1   (b2-h2)*y26hi))/h
      splintx = qmidx0

      end function splintx
