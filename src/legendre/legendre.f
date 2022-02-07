       recursive function factorial(num) result(fact)
       use kind_spec
       implicit none
       integer, intent(IN) :: num
       real(rprec) :: fact
       integer:: nfact
       
       if(num == 1 .or. num == 0) then
         nfact = 1
       else
         nfact = num*factorial(num-1)
       endif

       fact = nfact

       end function factorial


       subroutine build_matrices_legendre(n, a, b, a_inv, b_inv)
       use kind_spec
       implicit none
       integer, intent(IN) :: n
       real(rprec), intent(OUT), dimension(0:n,0:n) :: a, b, 
     1  a_inv, b_inv
       integer :: idx, i, j
       real(rprec) :: factorial
       real(rprec), dimension(13), parameter :: norm_legend =
     1    ( / 1, 1, 3, 5, 35, 63, 231, 429, 6435, 
     2        12155, 46189, 88179, 676039 / )       
       real(rprec), dimension(13), parameter ::  norm_inv_legend =
     1    ( / 1, 1, 2, 2, 8, 8, 16, 16, 128, 
     2        128, 256, 256, 1024 / )
       real(rprec), dimension(13*13), parameter :: b_legend =
     1  ( / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     2  0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3  1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4  0, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     5  7, 0, 20, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0,
     6  0, 27, 0, 28, 0, 8, 0, 0, 0, 0, 0, 0, 0,
     7  33, 0, 110, 0, 72, 0, 16, 0, 0, 0, 0, 0, 0,
     8  0, 143, 0, 182, 0, 88, 0, 16, 0, 0, 0, 0, 0,
     9  715, 0, 2600, 0, 2160, 0, 832, 0, 128, 0, 0, 0, 0,
     A  0, 3315, 0, 4760, 0, 2992, 0, 960, 0, 128, 0, 0, 0,
     B  4199, 0, 16150, 0, 15504, 0, 7904, 0, 2176, 0,
     C  256, 0, 0,
     D  0, 20349, 0, 31654, 0, 23408, 0, 10080, 0, 2432,
     E  0, 256, 0,
     F  52003, 0, 208012, 0, 220248, 0, 133952, 0, 50048,
     G  0, 10752, 0, 1024 / )
       real(rprec), dimension(13*13), parameter :: b_legend_inv =
     1  ( / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     2  0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3  -1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4  0, -3, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     5  3, 0, -30, 0, 35, 0, 0, 0, 0, 0, 0, 0, 0,
     6  0, 15, 0, -70, 0, 63, 0, 0, 0, 0, 0, 0, 0,
     7  -5, 0, 105, 0, -315, 0, 231, 0, 0, 0, 0, 0, 0,
     8  0, -35, 0, 315, 0, -693, 0, 429, 0, 0, 0, 0, 0,
     9  35, 0, -1260, 0, 6930, 0, -12012, 0, 6435, 0, 0,
     A  0, 0,
     B  0, 315, 0, -4620, 0, 18018, 0, -25740, 0, 12155, 
     C  0, 0, 0,
     D  -63, 0, 3465, 0, -30030, 0, 90090, 0, -109395, 0,
     E  46189, 0, 0,
     F  0, -693, 0, 15015, 0, -90090, 0, 218790, 0, -230945,
     G  0, 88179, 0,
     H  231, 0, -18018, 0, 225225, 0, -1021020, 0, 2078505,
     I  0, -1939938, 0, 676039 / )
     
!--------------------------------------------------------------------------------------------------
!      P_n(y):: LEGENDRE polynomial in [-1,1]. Matrix coefficients from 
!      "Handbook of Mathematical Functions", Abramowitz and Stegun, 1973, p.795.
!      If higher degree polyn. needed, then data set must be extended.
!
!      DESCRIPTION OF MATRICES:
!
!      A, A_inv: convert Legendre Pol. to and from y-powers in [-1,1]
!
!              norm_legend(m)*y**m = SUM_n=0^m (b_legend{mn} L_n)
!              norm_legend_inv(m)*L_m = SUM_n=0^N (b_legend_inv{mn} y**n)
!
!           ==>   A_{mn} = b_legend_{mn}/norm_legend(m)                if  m<=n  or  0. otherwise
!           ==>   A_inv_{mn} = b_legend_inv_{mn}/norm_inv_legend(m)    if  m<=n  or  0. otherwise
!
!      B, B_inv: convert y-powers in [-1,1] to and from x-powers in [0,1]
! 
!             x = (y + 1)/2;   x**m =  SUM_n=0^m  2**(-m) m!/(n!(m-n)!)  y**n           
!             y = 2*x - 1;     y**m =  SUM_n=0^m  2**n (-1)**(m-n) m!/(n!(m-n)!)  x**n    
!
!           ==>   B_{mn} =  2**(-m) m!/(n!(m-n)!)                     if  m<=n  or  0. otherwise
!           ==>   B_inv_{mn} =  2**n (-1)**(m-n) m!/(n!(m-n)!)        if  m<=n  or  0. otherwise
!       
!
!      NOTICE that:   A_inv = (A)**(-1);  B_inv = (B)**(-1) 
!--------------------------------------------------------------------------------------------------

       if (n > 12) stop 'N(Legendre) CANNOT be larger than 12!!!'
   
       a = 0
       a_inv = 0
       b = 0
       b_inv = 0
   
       do i = 0, n
         idx = 13*i
         do j = 0, i
           a(i,j) = b_legend(idx+j+1)/norm_legend(i+1)
           a_inv(i,j) = b_legend_inv(idx+j+1)/norm_inv_legend(i+1)
           b(i,j) = 2._dp**(-i) * factorial(i)/
     1       (factorial(j)*factorial(i-j))
           b_inv(i,j) = 2._dp**j*(-1._dp)**(i-j)* factorial(i)/
     1       (factorial(j)*factorial(i-j))
         enddo
       enddo
   
       end subroutine build_matrices_legendre

   
       subroutine power_to_legendre(n, a, b, ac, tc)
       use kind_spec
       implicit none
       integer, intent(IN):: n
       real(rprec), dimension(0:n), intent(IN):: ac
       real(rprec), dimension(0:n), intent(OUT):: tc
       real(rprec), dimension(0:n,0:n), intent(IN):: a, b
       integer:: i, j, k
!------------------------------------------------------------------
!      Given the following notation:
!
!           TC == (tc(1), ...tc(n))==>  vector of coefficients
!                 for Legendre series in [-1,1]
!           AC == (ac(1), ...ac(n))==>  vector of coefficients
!                 for power series in [0,1]
!      then:
!                        TC = AC* B * A
!------------------------------------------------------------------     
       do i = 0, n
         tc(i) = 0
         do j= 0, n
           do k = 0, n
             tc(i) = tc(i) + ac(j) * b(j,k) * a(k,i)
           enddo
         enddo
       enddo
       
       end subroutine power_to_legendre

       
       subroutine legendre_to_power(n, a_inv, b_inv, tc, ac)
       use kind_spec
       implicit none
       integer, intent(IN):: n
       real(rprec), dimension(0:n), intent(IN):: tc
       real(rprec), dimension(0:n), intent(OUT):: ac
       real(rprec), dimension(0:n,0:n), intent(IN):: a_inv, b_inv
       integer:: i, j, k
!---------------------------------------------------------------------
!      Given the following notation:
!
!           AC == (ac(1), ...ac(n))==>  vector of coefficients for 
!                 power series in [0,1]
!           TC == (tc(1), ...tc(n))==>  vector of coefficients for 
!                 Legendre series in [-1,1]
!      then:
!                        AC = TC* A_INV * B_INV
!----------------------------------------------------------------------   
       do i = 0, n
         ac(i) = 0
         do j= 0, n
           do k = 0, n
             ac(i) = ac(i) + tc(j) * a_inv(j,k) * b_inv(k,i)
           enddo
         enddo
       enddo
       
       end subroutine legendre_to_power

          
       subroutine point_wise_to_power(npts, x, f, n, ac, i_full,
     1   nsp_pts)
       use kind_spec
       implicit none
       integer, intent(IN):: npts, i_full, nsp_pts
       integer, intent(INOUT):: n
       real(rprec), dimension(npts), intent(IN):: x, f
       real(rprec), dimension(0:n), intent(OUT):: ac
       integer :: i, j, k , n_count, n_int      
       character(LEN=2) :: cdum
       character(LEN=15) :: form
       real(rprec) :: res, ohs_sp
       real(rprec), dimension(:), allocatable :: tc, norm, y,
     1   fy, y_sp, fy_sp
       real(rprec), dimension(:,:), allocatable:: a_inv, b_inv,
     1   a, b, inorm, tleg
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
       REAL(rprec) , EXTERNAL :: integral
!--------------------------------------------------------------------------------
!      MEANING OF ARGUMENTS:
!      NPTS = number of points where the point-wise function is given
!      X = independent variable values for these points in [0,1]
!      F = dependent variable values on the X points
!      N = tentative number of Legendre polynomials used in the reconstruction
!      AC = final power series coefficients in x ([0,1])
!      I_FULL =0, if input data on full mesh; =1, if on half mesh
!      NSP_PTS = number of points for the spline fit of F on the half mesh.
!--------------------------------------------------------------------------------
!
!      ABSTRACT:
!      Conversion of a point-wise solution in [0,1] given on a discrete set
!      of NPTS radial points (either on the half or the full mesh, see below)
!      first to an expansion in Legendre polynomials in [-1,1] via:
!
!          F = sum_m=0^N  <F,P_m>/<P_m,P_m>   <f,g> == INT_[-1,1] f(x)g(x)dx
!
!      and defining the Legendre coefficient vector TC as:
!
!          TC = {<F,P_m>/<P_m,P_m>, m=1, N}
!
!      then calculates the equivalent power series coefficients vector in 
!      [0,1], AC:
!
!                AC = TC*A_inv*B_inv
!
!      The value of Legendre Polynomials used "N" depends CRITICALLY on the 
!      number of points NPTS in the discrete data set of the incoming point-wise
!      solution. Because if the integrals are evaluated using these points, the
!      wiggles in the polynomial may not be well sampled by the existing set of
!      points and may cause deviations of the ortogonality condition: 
!
!                        <P_m,P_n> = 2/(2*n+1)* delta_{mn}
!
!      This makes the reconstruction to deteriorate quickly for increasing "N". 
!      To avoid this, the incoming data values are SPLINED to a number of points
!      NSP_PTS on a half-like mesh using Hermite cubic splines. Then, the
!      value of "N" is chosen requiring:
!
!               | 1- .5*(2*n+1)*<P_N,P_n> | < RES
!
!      RES is prescribed to be <= .5%. Until this condition is achieved, "N" is
!      subsequently reduced. If "N" gets equal to 2 in this reduction, the 
!      calculation is stopped and an error message issued.
!
!      WARNING: The point-wise solution MAY BE given on the half mesh or full mesh
!      of a equally spaced radial grid, that is:
!
!      ***I_FULL = 1 beginning with x(3/2) and finishing with x(M-1/2), where  
!         x(1) = 0. and x(M) = 1. 
!      ***I_FULL = 0 beginning with x(1) = 0. and finishing with x(M) = 1.
!
!--------------------------------------------------------------------------------

       ac = 0
       res = 1
              
!...   Map [0,1] to [-1,1]
       
       allocate (y(npts), fy(npts), y_sp(nsp_pts), fy_sp(nsp_pts))

       y(:npts) = 2*x(:npts) - 1
       fy(:npts) = f(:npts)

!...   I_full = 0 if FY given on the full mesh, and =1 if on the half mesh.   
!      FY is splined on nsp_pts on the half_mesh
       
       ohs_sp = 2._dp/nsp_pts
       y_sp= (/( -1._dp+ (real(i,rprec)-0.5_dp)*ohs_sp,
     1    i = 1, nsp_pts)/)  

       call spline_it(npts, y, fy, nsp_pts, y_sp, fy_sp, i_full)
       
!...   Determine number of Legendre Polynomials to be used

       n_count = 0
       do
         if(res <= 0.005_dp) exit
           n_int = n - n_count    
           if(n_int < n) deallocate (a, b, a_inv, b_inv, tleg,
     1       inorm, norm)
           allocate(a(0:n_int,0:n_int), b(0:n_int,0:n_int),
     1        a_inv(0:n_int,0:n_int), b_inv(0:n_int,0:n_int),
     2        tleg(0:n_int, nsp_pts), inorm(0:n_int,0:n_int), 
     3        norm(0:n_int))
     
           norm(0:n_int) = (/(2._dp/(2*i+1),
     1       i=0, n_int)/)     

           call build_matrices_legendre (n_int, a, b, a_inv, b_inv) 
                
!...   Build Initial Legendre polynomials on grid
       
           do i = 0, n_int
             tleg(i,:nsp_pts) = 0
             do k = n_int, 0, -1
               tleg(i,:nsp_pts)= a_inv(i,k) + y_sp(:nsp_pts)
     1            *tleg(i,:nsp_pts)
             enddo
           enddo
              
!...   Compute orthogonality and norms 

           do j = 0, n_int
             do i = 0, n_int
               inorm(i,j) = integral(nsp_pts, y_sp, tleg(j,1:nsp_pts), 
     1           tleg(i,1:nsp_pts))
             enddo
           write(cdum,'(i2)')n_int+1
           cdum = adjustl(cdum)
           form = '('//trim(cdum)//'(x,f8.3))'
           form = adjustl(form)
!          write(20,trim(form)) (inorm(i,j), i=0, n)
           enddo

!...   Compute residual. If res > 0.02, reduce N and iterate
       
           res = abs((norm(n_int) - inorm(n_int,n_int))/(norm(n_int)))
           n_count = n_count + 1
           if(n_count == n-1) stop 'Precision bad for even 2 Polynomial'           

       enddo 
      
       allocate(tc(0:n_int))
       do j = 0, n_int
         tc(j) = integral(nsp_pts, y_sp, fy_sp, tleg(j,1:nsp_pts))
     1    /norm(j)
       enddo
              
       call legendre_to_power(n_int, a_inv, b_inv, tc, ac(0))
       
       deallocate(a_inv, b_inv, a, b, tc, tleg, inorm, norm, fy, 
     1   y, y_sp, fy_sp) 
       
       end subroutine point_wise_to_power


       subroutine spline_it(ndata, xdata, ydata, npts, x, y, i_full)
       use kind_spec
       implicit none
       integer, intent(IN) :: ndata, npts, i_full
       real(rprec), dimension(ndata), intent(IN):: xdata, ydata
       real(rprec), dimension(npts), intent(OUT):: x, y
       real(rprec), dimension(npts) :: dy      
       real(rprec), dimension(:), allocatable:: xfull, yfull, 
     1   dyfull, wk
       integer:: iwk, ierr
       logical:: spline      

       allocate(xfull(ndata+i_full), yfull(ndata+i_full), 
     1     dyfull(ndata+i_full), wk(2*(ndata+i_full)) )

       if (i_full .eq. 1) then                                            !! I_FULL= 1, data on HALF_mesh
         xfull(1) = -1
         xfull(ndata+1) = 1
         xfull(2:ndata) = 0.5_dp*(xdata(1:ndata-1)+xdata(2:ndata)) 
         yfull(1) = ydata(1) + (xdata(1)+1.)*
     1    (ydata(2)-ydata(1))/(xdata(2)-xdata(1))
         yfull(2:ndata) = 0.5_dp*(ydata(1:ndata-1)+ydata(2:ndata)) 
         yfull(ndata+1) = ydata(ndata) + (1 - xdata(ndata))*
     1    (ydata(ndata-1)-ydata(ndata))/(xdata(ndata-1)-xdata(ndata))
       else                                                              !!I_FULL =0, data on FULL mesh;
         xfull(1:ndata) = xdata(1:ndata)
         yfull(1:ndata) = ydata(1:ndata)
       endif
       
       spline = .false.
       wk = 0
       ierr = 0
       iwk = 2*ndata

       call PCHEZ(ndata+i_full, xfull, yfull, dyfull, spline,
     1   wk, iwk, ierr)
        if(ierr.lt.0) stop 'LEGENDRE: error in SPLINE'
       
       call PCHEV(ndata+i_full, xfull, yfull, dyfull,
     1   npts, x, y, dy, ierr)
        if(ierr.lt.0) stop 'LEGENDRE: error in EVAL_SPLINE'

       deallocate(xfull, yfull, dyfull, wk)        

       end subroutine spline_it


       function integral (n, x, y1, y2) 
       use kind_spec
       implicit none
       integer, intent(IN) :: n 
       real(rprec), intent(IN), dimension(n) :: y1, y2 , x
       real(rprec), parameter :: one = 1
       real(rprec), dimension(n) :: coef   
       real(rprec) :: h, integral
!------------------------------------------------------------------------------
!      WARNING: The integral formula assumes values on the half MESH,
!      that is, the independent variable is given from x(3/2) to x(M-1/2),
!      but the limits of the integral are x(1)=-1 and x(M)=1. N = M-1 is the
!      number of points on the half mesh.
!      In the current version, fulfilment of this condition is only checked
!      for the case in which the interval of integration is the interval
!      of orthogonality of Legendre polynomias, i.e., [-1,1]. 
!      Integration formula is a 2-nd order HALF-MESH formula from
!      "Numerical Recipes", W.Press et al, Chapter 4, page 110.
!------------------------------------------------------------------------------
       if (n < 10) stop 'Too few points to carry out integration.!'
       if (x(n) < x(1)) stop' B <<A in INTEGRAL!!!!!'
       if (x(n) == one .or. x(1) == -one) stop 'HALF MESH INTEGRAL!'
        
       coef = 1
       
       h = (x(n)-x(1))/(n-1)       
       integral = h* SUM(coef*y1*y2)
       
       end function integral
