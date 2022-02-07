      subroutine axisopt(fsq, r00, iresidue, ivac)
      use vsvd
      use vparams, only: zero, one, nthreed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iresidue, ivac
      real(rprec) :: fsq, r00
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: smax = 0.998_dp
      real(rprec), parameter :: smin = 0.985_dp
      character*(60), parameter :: optbegin =
     1   'Begin variation of Raxis to minimize total RMS error'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: delstep, dedrmax, factor, delerr,
     1   dedr, rstepx1
      real(rprec), save :: delstep_old, errmax, errold, rstepx,
     1   raxold, scale
C-----------------------------------------------
 
 
      if (iresidue.lt.1 .or. errsvd*1.e6_dp.lt.one .or.
     1    fsq.gt.fturnon_axis .or. ivac .le.2) return 
!
!     MOVE R-AXIS BASED ON dR/dt = (-dEsvd/dR)
!     LIMIT MAXIMUM RSTEPX TO RSTEPX0
!     TRY TO FIND ZERO-CROSSING IN dEsvd/dR (ESTIMATED NUMERICALLY)
!
 
      if (iresidue .eq. 1) then                    !First time through
         iresidue = 2
         raxold = r00
         errold = errsvd
         errmax = zero
         rstepx = rstepx0
         scale = smax
         if (iopt_raxis .gt. 0) then
            write (*, 115) optbegin
            write (nthreed, 115) optbegin
         endif
      else
         delerr = errsvd - errold                !delta E-svd
         delstep = r00 - raxold                  !delta R-axis
         if (delerr.ne.zero .and. abs(delstep).gt.1.e-3_dp*rstepx0) then
            dedr = delerr/delstep
            errmax = max(errmax,errsvd)
            dedrmax = 2.0*errmax/rwidth
            rstepx1 = min(one,abs(dedr)/dedrmax)*rsfac*rstepx0
            factor = sign(one,(-dedr))        !Move in -dE/dR direction
            rstepx = rstepx1*factor
            scale = smax
            if (delstep*delstep_old .le. zero) scale = smin
            delstep_old = delstep
            raxold = r00
            errold = errsvd
         endif
      endif
      rsfac = scale*rsfac
c-5/1/96 raxmse = raxmse + rstepx
      raxmse = raxold + rstepx
  115 format(2x,a)
 
      end subroutine axisopt


      subroutine chisq(amat_i, amat_p, data, idata, isize, itotal)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: itotal
      integer, dimension(*) :: idata, isize
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, n, i1, ispec, is(5), ip(5), j2
      real(rprec) :: delsq, delp, dels
      character*(130) label
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!
 
         chisqerr(:jchix) = zero
 
         do i = 1, itotal
            delsq = (sum(ystark(:isnodes)*amat_i(:isnodes,i)) +
     1               sum(ythom(:ipnodes)*amat_p(:ipnodes,i))-data(i))**2
            if (i.ge.idata(ithom0) .and. i<idata(ithom0)+isize(ithom0)) 
     1         then
               chisqerr(ithom0) = chisqerr(ithom0) + delsq
            else if (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) then
               i1 = i - idata(istark0) + 1
               if (i1 .eq. islope) then
                  chisqerr(islope0) = delsq
               else if (i1 .eq. icurrout) then
                  chisqerr(icurr0) = delsq
               else
                  chisqerr(istark0) = chisqerr(istark0) + delsq
               endif
            else if (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) then
               chisqerr(idiam0) = chisqerr(idiam0) + delsq
            else if (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) then
               chisqerr(iflxs0) = chisqerr(iflxs0) + delsq
            else if (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) then
               chisqerr(ibrzfld) = chisqerr(ibrzfld) + delsq
            endif
         end do
 
!
         errsvd = sum(chisqerr(:jchix))
         if (.not.lpprof) errsvd = errsvd - chisqerr(ithom0)
 
      if (iequi.ne.1 .or. .not.lrecon) return
         if (.not.lpprof) then
            write (nthreed, 15)
         else
            write (nthreed, 10)
         endif
         if (lpprof) then
            do n = 1, nchistp
               write (nthreed, 20) nchi2(n), chi2(ithom0,n), 
     1            chi2(istark0,n), chi2(icurr0,n), chi2(idiam0,n), 
     2            chi2(iflxs0,n), chi2(ibrzfld,n), chi2(jchix1,n)
            end do
         else
            do n = 1, nchistp
               write (nthreed, 20) nchi2(n), chi2(istark0,n), 
     1            chi2(icurr0,n), chi2(idiam0,n), chi2(iflxs0,n), 
     2            chi2(ibrzfld,n), chi2(jchix1,n)
            end do
         endif
 
!
!       PRINT OUT MATRIX ELEMENTS (5 EACH FOR PRESSURE, IOTA)
!
         write (nthreed, 200)
         delp = (ipnodes - 1)/4.
         dels = (isnodes - 1)/4.
         ip(1) = 1
         is(1) = 1
         ip(2:4) = ip(1) + int((((/(j2,j2=2,4)/)) - 1)*delp)
         is(2:4) = is(1) + int((((/(j2,j2=2,4)/)) - 1)*dels)
         ip(5) = ipnodes
         is(5) = isnodes
         write (label, 210) ip(1), ip(2), ip(3), ip(4), ip(5), is(1), 
     1      is(2), is(3), is(4), is(5)
         write (nthreed, 220) label
         ispec = 0
         do i = 1, itotal
            if (i.ge.idata(ithom0) .and. 
     1         i.lt.idata(ithom0)+isize(ithom0)) then
               i1 = i - idata(ithom0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, '   PRES  (')
            else if (i.ge.idata(istark0) .and. i<idata(istark0)
     1            +isize(istark0)) then
               i1 = i - idata(istark0) + 1
               if (i1 .eq. islope) then
                  ispec = ispec - 1
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, 1, ip, is, '  IOTA0  (')
               else if (i1 .eq. icurrout) then
                  ispec = ispec - 1
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, 1, ip, is, ' CURRENT (')
               else
                  call printmatrix (amat_p(1,i), amat_i(1,i), data(i), 
     1               i, i1 + ispec, ip, is, '   MSE   (')
               endif
            else if (i.ge.idata(idiam0) .and. i<idata(idiam0)
     1            +isize(idiam0)) then
               i1 = i - idata(idiam0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, ' DIAMAG  (')
            else if (i.ge.idata(iflxs0) .and. i<idata(iflxs0)
     1            +isize(iflxs0)) then
               i1 = i - idata(iflxs0) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, ' FLUXES  (')
            else if (i.ge.idata(ibrzfld) .and. i<idata(ibrzfld)
     1            +isize(ibrzfld)) then
               i1 = i - idata(ibrzfld) + 1
               call printmatrix (amat_p(1,i), amat_i(1,i), data(i), i, 
     1            i1, ip, is, '  BR-BZ  (')
            endif
         end do
   20 format(i6,1p8e12.4)
   10 format(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER   Thomscat      Stark',5x,'Current',5x,'Diamag.',5x,
     4   'Saddle',6x,' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
   15 format(/,30x,'ABSOLUTE CHI-SQ ERROR BY DATA TYPE'/,30x,
     1   '(NOT NORMED BY NUMBER DATA POINTS)'/,20x,
     2   'NOTE: STARK CHISQ MAY BE EVALUATED AT REDISTRIBUTED KNOTS'/,
     3   '  ITER      Stark',5x,'Current',5x,'Diamag.',5x,'Saddle',6x,
     4   ' B-Loops',6x,'TOTAL',/,1x,5('-'),7(2x,10('-')))
  200 format(//,38x,'SPLINE MATRIX ELEMENTS BY DATA TYPE'/,30x,
     1   ' AI(i,j)*iota(j) + AP(i,j)*[mu0*pres(j)] = DATA(i)'/,30x,
     2   ' NOTE: DATA(I) IS THE RAW DATA NORMED TO SIGMA(I)'//)
  210 format('   I   TYPE         DATA(I)','  AP(I,',i2,')  AP(I,',i2,
     1   ')  AP(I,',i2,')  AP(I,',i2,')','  AP(I,',i2,')  AI(I,',i2,
     2   ')  AI(I,',i2,')  AI(I,',i2,')','  AI(I,',i2,')  AI(I,',i2,')')
  220 format(a,/,3x,'-',3x,4('-'),9x,7('-'),10(2x,8('-')))
 
      end subroutine chisq

      
      subroutine printmatrix(amatp, amati, data, i, i1, ip, is, type)
      use vparams, only: rprec, nthreed
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer i, i1
      real(rprec) data
      character*(10) type
      integer, dimension(*) :: ip, is
      real(rprec), dimension(*) :: amatp, amati
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: k
C-----------------------------------------------
 
      write (nthreed, 10) i, type, i1, data, (amatp(ip(k)),k=1,5), (
     1   amati(is(k)),k=1,5)
 
   10 format(1x,i3,a10,i2,')',1p11e10.2)
 
      end subroutine printmatrix

       
      subroutine store_chisq
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
!
!       COMPUTES CHI**2 FROM DIFFERENT SUBROUTINES AT VARIOUS TIME-STEPS
!       WRITTEN BY D.K. LEE (3/93)
!
 
      if (mod(iter2,nstep).ne.10 .and. iequi.eq.0) return
         chisqerr(jchix1) = sum(chisqerr(:jchix))
         if (.not.lpprof) chisqerr(jchix1) = chisqerr(jchix1) 
     1      - chisqerr(ithom0)
         nchistp = nchistp + 1
         if (nchistp .gt. mstp) return
         chi2(:,nchistp) = chisqerr
         nchi2(nchistp) = iter2 - 10
         if (iequi .eq. 1) nchi2(nchistp) = iter2
         if (iter2 .eq. 10) nchi2(nchistp) = 1

      end subroutine store_chisq

      
      subroutine findphi(reven, rodd, rmeas, dse, dso, rmid, ismeas, 
     1   iumeas, indexr, npts)
      use vmec_main
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer npts
      integer, dimension(npts) :: ismeas, iumeas
      integer, dimension(2*ns) :: indexr
      real(rprec), dimension(ns,nzeta,*) :: reven, rodd
      real(rprec), dimension(npts) :: rmeas, dse, dso
      real(rprec), dimension(2*ns) :: rmid
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, k, itemp, jtemp, i1, j1, ns2
C-----------------------------------------------
 
 
!
!       THIS ROUTINE FINDS THE LOWER FLUX INDEX [=INDEXS] CORRESPONDING
!       TO THE MEASURED VALUE [=RMEAS] OF (R,Z=0) ALONG THE MIDPLANE
!       (THETA=0 OR PI).
!       THE QUANTITIES RMEAS(K=1,NPTS) ARE INTERPOLATED AS FOLLOWS:
!
!       RMEAS(K) = R[ISMEAS(K),IUMEAS(K)]*[1-DSO(K)] +
!                  R[ISMEAS(K)+1,IUMEAS(K)]*DSO(K)
!
!       BECAUSE OF THE SQRT(S) BEHAVIOUR OF R IN THE FIRST RADIAL ZONE,
!       THE ACTUAL S-INTERPOLAND IN THE FIRST ZONE IS DSO(K)**2 = DSE(K).
!       IN ALL OTHER ZONES, DSE(K) = DSO(K).
!
 
      if (npts .le. 0) return 
      ns2 = 2*ns
!
!     COMPUTE THE GRID VALUES (S-COORDINATE) OF R ALONG THE MIDPLANE,
!     STARTING AT THETA=PI (I=NTHETHA2) AND ENDING AT THETA=0 (I=1)
!
 
      rmid(:ns) = reven(indexr(:ns),1,ntheta2) + sqrts(indexr(:ns))
     1   *rodd(indexr(:ns),1,ntheta2)
      rmid(ns+1:ns2) = reven(indexr(ns+1:ns2),1,1) + 
     1   sqrts(indexr(ns+1:ns2))*rodd(indexr(ns+1:ns2),1,1)
 
!
!     FIND THE RADIAL ZONE INDEX [=ITEMP], WHICH BRACKETS THE MEASURED R-VALUE
!
!     RMID(ITEMP-1) .le. RMEAS .le. RMID(ITEMP)
!
 
      do k = 1, npts
         itemp = 0
         do i = 1, ns2 - 1
            if (rmeas(k) .lt. rmid(i)) then
               itemp = i
               go to 100
            endif
         end do
         itemp = ns2
!
!         FIND FLUX-COORDINATE S-INDEX [=ISMEAS], POLOIDAL ANGLE
!         INDEX [=IUMEAS], AND INTERPOLAND [=DSO]
!
 
  100    continue
         if (itemp.gt.1 .and. itemp.lt.ns2) then
            i1 = itemp - 1
            jtemp = indexr(itemp)
            j1 = indexr(i1)
            dso(k) = (rmeas(k)-rmid(i1))/(rmid(itemp)-rmid(i1))
            if (j1 .lt. jtemp) then                 !THETA = 0
               ismeas(k) = j1
               iumeas(k) = 1
            else
               ismeas(k) = jtemp
               dso(k) = 1.0 - dso(k)
               iumeas(k) = ntheta2
            endif
         else
            dso(k) = 1.0
            ismeas(k) = indexr(itemp) - 1
!           IGNORE MEASURED POINTS OUTSIDE GRID
            if (itemp.eq.1 .or. rmeas(k).gt.rmid(ns2-1)) dso(k) = -1.0
            if (itemp .eq. ns2) iumeas(k) = 1
         endif
!
!        ACCOUNT FOR SQRT(S) SINGULARITY IN 1st ZONE
!        DSE IS THE S-FRACTIONAL INTERPOLAND
!
         dse(k) = dso(k)
         if (ismeas(k) .eq. 1) dse(k) = dso(k)*dso(k)
      end do
 
      end subroutine findphi

      
      subroutine fixrecon(ier)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ier
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: angle_variance  = 0.3_dp
      real(rprec), parameter :: radial_variance = 0.3_dp
      real(rprec), parameter :: p_threshold = 1.e-3_dp
      real(rprec), parameter :: c1p5 = 1.5_dp
      integer, parameter :: inode_max = 15
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(itse) :: isortp
      integer, dimension(1) :: isamax
      integer :: istat1 = 0, istat2 = 0, istat3 = 0, ii, i, ineg
     1   , ipos, ipmax, ioff, ileft, n, n1, icount, index, ind, m, m1
      real(rprec) :: delstark, datapos, dataneg
      real(rprec), dimension(imse+1) :: datalsq_s, stark_temp
      real(rprec) :: datamax, datamin, t1, rneg,
     1   rpos, datanorm, presmax, presin, presout, presmin, tsign
      real(rprec), dimension(itse) ::
     1   datalsq_p, ythom0, y2thom0, qtemp
      logical :: l1v(imse), l4v(itse)
C-----------------------------------------------
!
!
!                INDEX OF LOCAL VARIABLES
!
!         needflx    =NEEDIT, loop required for flux match
!                    >=ISYMCOIL, loop required but flux computed by
!                               invoking symmetry in Z
!                    =IDONTNEED, loop not required for flux match
!        needbfld    =NEEDIT, loop required for B-field match
!                    =ISAMECOIL, loop at same position as previous loop
!                    >=ISYMCOIL, loop required but B-field computed
!                               by invoking symmetry in Z
!                    =IDONTNEED, loop not required in B-field match
!          dsiext    connected flux loop signals due to external coils
!          plflux    array of measured (inferred) plasma contrib. to flux loops
!          plbfld    array of measured (inferred) plasma contrib. to B-loops
!

      call free_mem_recon

      ier = 0
!
!       ONLY QUANTITIES INDEPENDENT OF RADIAL MESH CAN GO HERE
!
!       STARK-DATA CONSISTENCY CHECK
!
      delstark = tan(dcon*angle_variance)

!
!       SORT STARK DATA IN ASCENDING ORDER IN R-SPACE
!       AND RE-INDEX DATASTARK, SIGMA_STARK ARRAYS ON ASCENDING RSTARK ARRAY
!
!       SCALE MOTIONAL STARK DATA TO RADIANS
!       RSTARK = R position along Z=0 plane of measurement
!       DATASTARK = ARCTAN ( Bpol/Btor ) at RSTARK, in degrees
!       AND IS CONVERTED TO Bpol/Btor
!

      ii = 0
      l1v(:imse) = (rstark(:imse) .gt. zero)
      do i = 1, imse
         if (l1v(i)) ii = ii + 1
      end do
      if (ii .ne. imse) then
         print *, 'There is a zero in the RSTARK array ?!'
         ier = 1
         return
      endif

      allocate (indexs1(imse+3), indexu1(imse+3), isortr(imse+3),
     1   isorts(imse+3), delse1(imse+3), delso1(imse+3),
     2   starkcal(imse+3), qmeas(imse+3), qcalc(imse+3),
     3   fpsical(imse+3), stark_weight(imse+3),
     4   rsort(imse+3), rsort0(imse+3), stat=istat1)
      if (istat1 .ne. 0) then
         print *, ' ISTAT1 = ', istat1, ' ALLOCATING INDEXS1'
         stop
      endif

      datalsq_s(:imse) = datastark(:imse)
      stark_temp(:imse) = sigma_stark(:imse)
      call sort_data (rstark, isorts, imse)
      do i = 1, imse
         datastark(i) = datalsq_s(isorts(i))
         sigma_stark(i) = stark_temp(isorts(i))
         if (sigma_stark(i) .ge. cbig) then
            print *, 'SIGMA_STARK missing'
            ier = 1
            return
         endif
         if (sigma_stark(i) .lt. zero) sigma_stark(i) =
     1     abs(sigma_stark(i)*datastark(i))
          !CONVERT TO Bpol/Btor ... applying profile offsets along the way!
         datastark(i) = tan(dcon*(datastark(i)+mseangle_offset+
     1      mseangle_offsetm*mseprof(i)))
         sigma_stark(i) = tan(dcon*sigma_stark(i))
         if (sigma_stark(i) .eq. zero) sigma_stark(i) = 0.1*delstark
      end do
      rstarkmin = rstark(1)
      rstarkmax = rstark(imse)

!
!     NEED FOR SCALING EDGE IOTA
!     HERE, SINCE RSTARK IS ORDERED, DATAMAX -> RIGHT OF AXIS
!     AND DATAMIN -> LEFT OF AXIS
!
!       GET RC0MSE (ESTIMATE FOR MAGNETIC AXIS)
!
      rwidth = radial_variance*(rbc(0,1)+abs(rbs(0,1)))

      if (imse .gt. 0) then
         datamax = datastark(imse)
         datamin = datastark(1)
         rc0mse = 0.
         ineg = 1
c                                     !NO ZERO CROSSING: FIND MIN ANYHOW
          if (datamax * datamin.gt.zero) then      !NO ZERO CROSSING: FIND MIN ANYHOW
            datamin = abs(datamin)
            do i = 2,imse
              t1 = abs(datastark(i))
              if (t1.lt.datamin) then
                datamin = t1
                ineg = i
              endif
            end do
             if (ineg.eq.imse) ineg = ineg-1
            ipos = ineg + 1
            goto 310
          else if (( datamax*signiota.lt.zero ) .or.
     >             (datamin*signiota.gt.zero)) then
            datastark = -datastark
          endif

!
!       ALLOW FOR POSSIBLE MULTIPLE ZERO CROSSINGS (WIGGLES) IN DATA
!

         do i = 1, imse
            if (datastark(i)*signiota .le. zero) then
               ineg = i                          !LEFT OF MAGNETIC AXIS
            else
               exit                              !RIGHT OF MAGNETIC AXIS
            endif
         end do
         do i = imse, ineg + 1, -1
            if (datastark(i)*signiota .le. zero) then
               exit                              !LEFT OF MAGNETIC AXIS
            else
               ipos = i                          !RIGHT OF MAGNETIC AXIS
            endif
         end do

  310    continue

         rneg = rstark(ineg)
         rpos = rstark(ipos)
         dataneg = datastark(ineg)
         datapos = datastark(ipos)
      endif                                      !End of if(imse>0)
      if (datapos .ne. dataneg) then
         rc0mse = (datapos*rneg - dataneg*rpos)/(datapos - dataneg)
         rwidth=delstark*abs((rneg-rpos)/(datapos-dataneg))+rwidth
      endif
      if (ipos .gt. ineg + 1) rc0mse = 0.5_dp*(rpos + rneg)

!
!       ESTIMATE MAGNETIC AXIS FROM RAXIS
!
      raxmse = rc0mse
      if (rc0mse.eq.0.0_dp .or. iopt_raxis.ne.1) raxmse = raxis(0,1)
      rstepx0 = 0.005_dp*rwidth

!
!       COMPUTE SPLINES IN R-SPACE FOR MATCHING IOTA(0)
!

      datanorm = zero
      delse1(:imse) = one
      qcalc(:imse) = one/sigma_stark(:imse)**2
      datalsq_s(:imse) = datastark(:imse)*qcalc(:imse)
      scstark = sum(abs(datastark(:imse)/sigma_stark(:imse)))
      datanorm = sum(abs(datastark(:imse)))
      scstark = datanorm/scstark
c04-96        call setspline(rstark,qcalc,datalsq_s,stark_temp,ystark0,
c04-96     >  y2stark0,delse1,0.1*tensi/scstark**2,imse,NATUR)

!
!       DETERMINE NUMBER OF IOTA KNOTS IN SQRT(S)-SPACE
!
      if (isnodes .le. 0) then
         isnodes = min(inode_max,imse + 1)
         isnodes = max(5,isnodes)               !At least 5 spline knots
      endif
      if (isnodes .lt. 5) stop 'MUST PICK ISNODES > 4'
      write (nthreed, *)
     1   'Number of iota-spline knots (in s-space):     ', isnodes

      allocate (nk_ia(isnodes), nk_ib(isnodes), hstark(isnodes),
     1  y2stark(isnodes), ystark(isnodes), sknots(isnodes), stat=istat1)
      if (istat1 .ne. 0) then
         print *, ' ISTAT1 = ', istat1, ' ALLOCATING NK_IA'
         stop
      endif

!
!       COMPUTES NODES IN SQRT(S) SPACE FOR SPLINES
!       THIS ASSUMES A FIXED NO - ISNODES - OF KNOTS
!       THIS MAY NOT PRESERVE ISNODES.
!       ALSO, IT IS NOT NECESSARY TO TAKE EQUALLY SPACED KNOTS IN
!       SQRT(S) SPACE. INDEED, THE FOLLOWING CHOICES ARE POSSIBLE:
!
!       SKNOTS(I) = HNODES*(I-1)   .eq.>   EQUAL-SPACED IN SQRT(S)
!
!       SKNOTS(I) = SQRT(HNODES*(I-1)) .eq.>  EQUAL-SPACED IN S
!
!       DO NOT - UNDER ANY CIRCUMSTANCES - CHANGE THE ARGUMENTS TO
!       THE SPLINT, GETSPLINE, SETUP_INT ROUTINES FROM SQRTS,SHALF
!       TO SQRTS**2, SHALF**2 TO DO S-INTERPOLATION. RATHER, CHANGE
!       SKNOTS (AND PKNOTS) ACCORDING TO THE ABOVE FORMULA. THIS IS
!       ABSOLUTELY CRUCIAL, SINCE ONLY IN SQRT(S) SPACE DO THE
!       FIRST DERIVATIVE BOUNDARY CONDITIONS, d IOTA/d SQRT(S) = 0
!       (SIMILAR FOR P) APPLY AT THE AXIS, S=0.
!
!
      do i = 1, isnodes
         sknots(i) = real(i - 1,rprec)/(isnodes - 1)
      end do

      hstark(:isnodes-1) = sknots(2:isnodes) - sknots(:isnodes-1)


!
!       SET UP DATA ARRAY SCALE FACTORS
!       ACTUAL PRESSURE = PRESPEAK * PFAC * P(INTERNAL)
!       IF (LPOFR) THEN DATA ARE INPUT vs R (REAL SPACE)
!       IF (.NOT.LPOFR),DATA ARE INPUT vs S (FLUX SPACE)
!
      if (itse .eq. 0) call getpresprofile         !!Simulate 'data'
      if (itse .gt. 0) then

         allocate (sthom(itse), delse2(itse), delso2(itse), pcalc(itse),
     1      indexs2(itse), indexu2(itse), stat=istat1)
         if (istat1 .ne. 0) then
            print *, ' ISTAT1 = ', istat1, ' ALLOCATING STHOM'
            stop
         endif


         datathom(:itse) = datathom(:itse)*presfac
         presmax = maxval(datathom(:itse))

!
!       SORT DATA IN ASCENDING ORDER IN R-SPACE (LPOFR) OR S-SPACE(.NOT.LPOFR)
!       AND RE-INDEX DATATHOM, SIGMA_THOM ARRAYS ON ASCENDING RTHOM ARRAY

         datalsq_p(:itse) = datathom(:itse)
         qtemp(:itse) = sigma_thom(:itse)
         call sort_data (rthom, isortp, itse)
         if (lpofr) then
            do i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               rthom(i) = rthom(i) + pres_offset
               if (rthom(i) .le. zero) then
                  print *, 'Check units of PRES_OFFSET: rthom < 0!'
                  ier = 1
                  return
               endif
               if (datathom(i) .eq. presmax) ipmax = i
               sigma_thom(i) = qtemp(isortp(i))
               if (sigma_thom(i) .ge. cbig) then
                  print *, 'SIGMA_THOM missing'
                  ier = 1
                  return
               endif
               if (sigma_thom(i) .lt. zero) then
                  sigma_thom(i) = abs(sigma_thom(i)*datathom(i))
               else
                  if (sigma_thom(i) .gt. zero) then
                     sigma_thom(i) = presfac*sigma_thom(i)
                  else
                     sigma_thom(i) = p_threshold*presmax
                  endif
               endif
            end do
         else
            do i = 1, itse
               datathom(i) = datalsq_p(isortp(i))
               sthom(i) = rthom(i)
               if (datathom(i) .eq. presmax) ipmax = i
               sigma_thom(i) = qtemp(isortp(i))
               if (sigma_thom(i) .ge. cbig) then
                  print *, 'SIGMA_THOM missing'
                  ier = 1
                  return
               endif
               if (sigma_thom(i) .lt. zero) then
                  sigma_thom(i) = abs(sigma_thom(i)*datathom(i))
               else
                  if (sigma_thom(i) .gt. zero) then
                     sigma_thom(i) = presfac*sigma_thom(i)
                  else
                     sigma_thom(i) = p_threshold*presmax
                  endif
               endif
            end do
         endif

!
!       THROW AWAY NOISY (SMALL) PRESSURE DATA BELOW P_THRESHOLD
!       STARTING FROM PEAK WORKING TO LARGER, SMALLER R
!
         ineg = ipmax
         ipos = ipmax
         do while(ineg.gt.1 .and.
     1            datathom(ineg-1).ge.p_threshold*presmax)
            ineg = ineg - 1
         end do
         do while(ipos.lt.itse .and.
     1            datathom(ipos+1).ge.p_threshold*presmax)
            ipos = ipos + 1
         end do
         itse = ipos - ineg + 1
         ioff = ineg - 1
         do i = 1, itse
            datathom(i) = datathom(ioff+i)
         end do
         do i = 1, itse
            rthom(i) = rthom(ioff+i)
         end do
         do i = 1, itse
            sigma_thom(i) = sigma_thom(ioff+i)
         end do
!
!       COMPUTE PRESSURE AND 1/SIGMA SPLINES IN R-SPACE (OR S-SPACE)
!       a. PRESSURE SPLINE
!
         datanorm = zero
         delse2(:itse) = one
         pcalc(:itse) = one/sigma_thom(:itse)**2
         datalsq_p(:itse) = datathom(:itse)*pcalc(:itse)
         scthom = sum(abs(datathom(:itse)/sigma_thom(:itse)))
         datanorm = sum(abs(datathom(:itse)))
         scthom = datanorm/scthom
         call setspline (rthom, pcalc, datalsq_p, qtemp, ythom0,
     1      y2thom0, delse2, 0.1*tensp/scthom**2, itse, natur)


!
!       FIND PRESSURE PEAK USING SMOOTHED DATA
!
         isamax = maxloc(ythom0(:itse))
         i = isamax(1)
         pthommax = ythom0(i)
         rthompeak = rthom(i)

         ileft = 0                    !Count data points to left of peak
         l4v(:itse) = rthom(:itse) < rthompeak
         do i = 1, itse
            if (l4v(i)) ileft = ileft + 1
         end do

         if (ipnodes .le. 0) then
            ipnodes = max(ileft + 1,itse - ileft)
            ipnodes = min(inode_max,ipnodes)
            ipnodes = max(5,ipnodes)            !At least 5 spline knots
            if (.not.lpprof) ipnodes = 7
         endif
         if (ipnodes < 5) stop 'MUST PICK IPNODES > 4'
         write (nthreed, *)
     1      'Number of pressure-spline knots (in s-space): ', ipnodes

         allocate( nk_pa(ipnodes), nk_pb(ipnodes), ythom(ipnodes),
     1      y2thom(ipnodes), hthom(ipnodes), pknots(ipnodes) )
         if (istat1 .ne. 0) then
            print *, ' ISTAT1 = ', istat1, ' ALLOCATION NK_PA'
            stop
         endif

!
!       COMPUTE NODES IN SQRT(S) SPACE FOR SPLINES
!       (SEE COMMENTS ABOVE PRECEDING SKNOTS(I) CALCULATION)
!
         do i = 1, ipnodes
            pknots(i) = real(i - 1,rprec)/(ipnodes - 1)
         end do
         hthom(:ipnodes-1) = pknots(2:ipnodes) - pknots(:ipnodes-1)
!
!       COMPUTE MINOR RADII FOR DETERMINING PHIEDGE
!
         if (lpofr) then
            rthommax = rthom(itse)
            rthommin = rthom(1)
            presin = datathom(1)
            presout = datathom(itse)
            presmin = min(presin,presout)
            ipresin = 0
            ipresout = 0
            if (presin .eq. presmin) then
               ipresin = 1
               if (presout.le.c1p5*presmin .or. presout<=0.1_dp*presmax)
     1            ipresout = 1
            else
               ipresout = 1
               if(presin.le.c1p5*presmin .or. presin<=0.1_dp*presmax)
     1            ipresin=1
            endif
         else
            ipresin = 0                        !Only use theta=0 in pofs
            ipresout = 1
         endif
      endif                                   !End of if(itse.gt.0) test

!
!       COMPUTE INDICES OF FLUX LOOP MEASUREMENTS NEEDED FOR
!       LOOP SIGNAL MATCHING
!       ALSO COMPUTE CONNECTED EXTERNAL POLOIDAL FLUX ARRAY

      if (.not.lfreeb) then
         nflxs = 0
         nobser = 0
         nobd  = 0
         nbsets = 0
      end if

      nmeasurements = imse + itse + 2 + nflxs   !!1 for diamag, 1 for edge MSE
      do n = 1, nobser
         needflx(n) = idontneed
         iconnect(1,nobd+n) = n            !For outputting absolute flux

         if (lasym) cycle

!        Save Index of up-down symmetric spatial observation points
         do n1 = 1,n-1
           if ((xobser(n1).eq. xobser(n)) .and.
     >         (zobser(n1).eq.-zobser(n))) needflx(n) = n1
         enddo
      end do

      do n = 1, nobd
         dsiext(n) = zero
         if (sigma_flux(n) .lt. zero) sigma_flux(n) =
     1     abs(sigma_flux(n)*dsiobt(n))
         if (sigma_flux(n) .eq. zero) sigma_flux(n) = 0.0001
         do icount = 1, 4
            index = iconnect(icount,n)
            tsign = sign(1,index)
            if(index.ne.0)dsiext(n) = dsiext(n)+psiext(abs(index))*tsign
         end do
      end do

      do n = 1,nflxs
        index = indxflx(n)
        if (index.gt.0) then
          plflux(index) = dsiobt(index) - dsiext(index)        !n-th connected PLASMA flux
          do icount = 1,4
            ind = abs(iconnect(icount,index))
            if ((ind.gt.0).and.(needflx(ind).eq.IDONTNEED))
     1      needflx(ind) = NEEDIT
          enddo
        endif
      enddo

!
!       COMPUTE INDICES OF EXTERNAL BFIELD MEASUREMENTS NEEDED
!       FOR SIGNAL MATCHING
!       FOR MULTIPLE ANGLES AT THE SAME R,Z,PHI LOCATION, IT IS
!       ASSUMED THE LOOP DATA ARE CONSECUTIVELY ORDERED
!
      do n = 1, nbsets
         nmeasurements = nmeasurements + nbfld(n)
         do m = 1,nbcoils(n)
            needbfld(m,n) = IDONTNEED
            if( (m.gt.1).and.(rbcoil(m,n).eq.rbcoil(m-1,n)).and.
     1        (zbcoil(m,n).eq.zbcoil(m-1,n)) )needbfld(m,n)=ISAMECOIL
            if( sigma_b(m,n).lt.zero )sigma_b(m,n) =
     1      abs(sigma_b(m,n) * bbc(m,n))
            if( sigma_b(m,n).eq.zero )sigma_b(m,n) = 0.0001

            if (lasym) cycle
!       CHECK FOR ANTISYMMETRIC SITUATED COIL FOR M1 < M, N <= NSETS
            do n1 = 1, nbsets
               do m1 = 1, m-1
               if( (rbcoil(m1,n1).eq.rbcoil(m,n)) .and.
     1            (zbcoil(m1,n1).eq.-zbcoil(m,n)).and.
     2            (abcoil(m1,n1).eq.abcoil(m,n)) )
     3         needbfld(m,n) = n1 + nbsets*(m1-1)
               enddo
            enddo

         enddo
      enddo

      do n1 = 1, nbsets
        do m1 = 1,nbfld(n1)
          index = indxbfld(m1,n1)
          if( index.gt.0 )then
!m-th PLASMA B-field
            plbfld(index,n1) = bbc(index,n1) - bcoil(index,n1)
            if( needbfld(index,n1).eq.IDONTNEED )
     1      needbfld(index,n1) = NEEDIT
          endif
        enddo
      enddo

      end subroutine fixrecon

      
      subroutine getpresprofile
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      real(rprec) :: sigmin
C-----------------------------------------------
 
!
!       WRITE OVER THOMPSON DATA
!
      lpofr = .false.                   !!these data are at s-half nodes
      lpprof = .false.
      itse = 10
      pthommax = datathom(1)                  !!compute in final version
      sigmin = 0.03*pthommax
 
      do i = 1, itse
         rthom(i) = real(i - 1,rprec)/(itse - 1)
         datathom(i) = (pthommax - sigmin)*(1. - rthom(i))**2 + sigmin
         sigma_thom(i) = 0.2*pthommax
      end do
 
      end subroutine getpresprofile

      
      subroutine getgreen
      use vsvd
      use vparams, only: twopi
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: jnm, i
      real(rprec), dimension(10) :: ak, bk
      real(rprec), dimension(9) :: ae, be
      real(rprec), dimension(jngrn) :: ye, yk, sqrt1u
      real(rprec):: dqk2, qk2, eta, alg, sum1, suma, sum2,
     1     sumb, sum3, sumc, sum4, sumd
C-----------------------------------------------
 
      data ak/3.0072519903686507E-04_dp, 3.9684709020989806E-03_dp,
     1   1.0795990490591656E-02_dp, 1.0589953620989356E-02_dp,
     2   7.5193867218083799E-03_dp, 8.9266462945564728E-03_dp,
     3   1.4942029142282098E-02_dp, 3.0885173001899746E-02_dp,
     4   9.6573590301742396E-02_dp, 1.3862943611198872e+0_dp/
      data bk/6.6631752464607272E-05_dp, 1.7216147097986537E-03_dp,
     1   9.2811603829686118E-03_dp, 2.0690240005100891E-02_dp,
     2   2.9503729348688723E-02_dp, 3.7335546682286003E-02_dp,
     3   4.8827155048118076E-02_dp, 7.0312495459546653E-02_dp,
     4   1.2499999999764055e-1_dp, 5.0000000000000000e-1_dp/
      data ae/3.2519201550638976E-04_dp, 4.3025377747931137E-03_dp,
     1   1.1785841008733922E-02_dp, 1.1841925995501268E-02_dp,
     2   9.0355277375409049E-03_dp, 1.1716766944657730E-02_dp,
     3   2.1836131405486903E-02_dp, 5.6805223329308374E-02_dp,
     4   4.4314718058336844E-1_dp/
      data be/7.2031696345715643E-05_dp, 1.8645379184063365E-03_dp,
     1   1.0087958494375104E-02_dp, 2.2660309891604169E-02_dp,
     2   3.2811069172721030E-02_dp, 4.2672510126591678E-02_dp,
     3   5.8592707184265347E-02_dp, 9.3749995116366946E-02_dp,
     4   2.4999999999746159E-1_dp/
 
!
!       Compute "Green's Functions" for Poloidal Flux, 2*pi*R*A-sub-phi,
!       BR, and BZ at point (XT,ZT) due to unit current (mu0*I = 1) at (XS,ZS) ...
!       modified to interpolate on k**2 - 3-34-92 - sph
!
      jnm = jngrn - 1
      odqk2 = (jnm)
      dqk2 = 1.0_dp/odqk2
      do i = 2, jnm
         qk2 = dqk2*(i - 1)
         qsq(i) = qk2
         eta = 1 - qk2
         alg = log(eta)
         sum1 = ((((ak(1)*eta+ak(2))*eta+ak(3))*eta+ak(4))*eta+ak(5))*
     1      eta + ak(6)
         suma = (((sum1*eta + ak(7))*eta+ak(8))*eta+ak(9))*eta + ak(10)
         sum2 = ((((bk(1)*eta+bk(2))*eta+bk(3))*eta+bk(4))*eta+bk(5))*
     1      eta + bk(6)
         sumb = (((sum2*eta + bk(7))*eta+bk(8))*eta+bk(9))*eta + bk(10)
         yk(i) = suma - alg*sumb
         sum3 = (((ae(1)*eta+ae(2))*eta+ae(3))*eta+ae(4))*eta
         sumc = (((((sum3 + ae(5))*eta+ae(6))*eta+ae(7))*eta+ae(8))*eta+
     1      ae(9))*eta
         sum4 = (((be(1)*eta+be(2))*eta+be(3))*eta+be(4))*eta
         sumd = (((((sum4 + be(5))*eta+be(6))*eta+be(7))*eta+be(8))*eta+
     1      be(9))*eta
         ye(i) = sumc - alg*sumd + 1
         yf(i) = ((1 + eta)*yk(i)-2*ye(i))/qk2
      end do
      ye(1) = 0.25_dp*twopi
      ye(jngrn) = 1
      yk(1) = ye(1)
      yk(jngrn) = 2*yk(jnm) - yk(jngrn-2)
      yf(1) = 0.
      yf(jngrn) = 2*yf(jnm) - yf(jngrn-2)
      qsq(1) = 0
      qsq(jngrn) = 1
 
      sqrt1u = sqrt(qsq(:jngrn))/twopi
c                                      !Factor of 1/2 from sqrt(4*xs*xt)
      yek(:jngrn) = 0.5_dp*sqrt1u*(ye(:jngrn)-yk(:jngrn))
      yeq(:jngrn) = 0.25_dp*qsq(:jngrn)*sqrt1u*ye(:jngrn)
c                                 !Factor of 2 absorbed by sqrt(4 xt xs)
      yf(:jngrn) = twopi*sqrt1u*yf(:jngrn)
      dyek(:jnm) = (yek(2:jnm+1)-yek(:jnm))*odqk2
      dyeq(:jnm) = (yeq(2:jnm+1)-yeq(:jnm))*odqk2
      dyf(:jnm) = (yf(2:jnm+1)-yf(:jnm))*odqk2
 
      end subroutine getgreen

      
      subroutine getlim
      use vmec_main
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: resup = 3.0_dp
      real(rprec), parameter :: resdn = 0.5_dp
      real(rprec), parameter :: eps = 0.005_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ntheta_2pi, nphi_plane, i,  
     1   nthtot, iexpand, ishrink, ionlim, n,
     2   limpts, nonlim, nexpand, nshrink, ilim0, nlim0
      real(rprec), dimension(2*ntheta1) ::
     1   rbdy, zbdy, rubdy, zubdy
      real(rprec) :: fshrink, distmax, fexpand
C-----------------------------------------------
 
c
c     DETERMINES WHEN PLASMA TOUCHES LIMITER
c     USE DOUBLE THE NO. OF THETA POINTS FOR INCREASED RESOLUTION
c
      ntheta_2pi = ntheta1
      nphi_plane = 1                 !Pick a phi plane (phi = 0 for now)
 
      rbdy(:ntheta3*2-1:2) = r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + r1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zbdy(:ntheta3*2-1:2) = z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + z1(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      rubdy(:ntheta3*2-1:2) = ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + ru(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
      zubdy(:ntheta3*2-1:2) = zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     1   nphi_plane,0) + zu(ns*nphi_plane:ns*ntheta3*nphi_plane:ns*
     2   nphi_plane,1)
 
 
      if (.not.lasym) then
!     FOR NOW, THIS ONLY WORKS FOR NZETA=1 (PHI=0 PLANE)
!     TO EXTEND TO OTHER PHI PLANES, MUST USE IREFLECT(NZETA)
      do i = 1, ntheta_2pi - ntheta2
         rbdy(2*(ntheta2+i)-1) = rbdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         zbdy(2*(ntheta2+i)-1) = -zbdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         rubdy(2*(ntheta2+i)-1) = -rubdy(2*(ntheta1-ntheta2-i)+3)
      end do
      do i = 1, ntheta_2pi - ntheta2
         zubdy(2*(ntheta2+i)-1) = zubdy(2*(ntheta1-ntheta2-i)+3)
      end do
      end if
 
!
!     FIND EVEN INDEXED POINTS BY INTERPOLATION
!
      nthtot = 2*ntheta_2pi
      rbdy(nthtot) = .5_dp*(rbdy(1)+rbdy(nthtot-1))
      zbdy(nthtot) = .5_dp*(zbdy(1)+zbdy(nthtot-1))
      rubdy(nthtot) = .5_dp*(rubdy(1)+rubdy(nthtot-1))
      zubdy(nthtot) = .5_dp*(zubdy(1)+zubdy(nthtot-1))
      rbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rbdy(3:ntheta_2pi*2-1:2) +
     1      rbdy(:ntheta_2pi*2-3:2))
      zbdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zbdy(3:ntheta_2pi*2-1:2) +
     1      zbdy(:ntheta_2pi*2-3:2))
      rubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(rubdy(3:ntheta_2pi*2-1:2)+
     1   rubdy(:ntheta_2pi*2-3:2))
      zubdy(2:(ntheta_2pi-1)*2:2) = .5_dp*(zubdy(3:ntheta_2pi*2-1:2)+
     1   zubdy(:ntheta_2pi*2-3:2))
 
      fshrink = 0.0_dp
      distmax = sum(rbdy(:nthtot)**2) + sum(zbdy(:nthtot)**2)
      fexpand = distmax
      iexpand = 0
      ishrink = 0
      ionlim = 0
 
      do n = 1, nlim
         limpts = limitr(n)
         call cauchy (rbdy, zbdy, rubdy, zubdy, rlim(:,n), zlim(:,n),
     1      reslim(:,n), seplim(:,n), distmax, nthtot, limpts)

           do i = 1,limpts
c       LIMITER POINT ON PLASMA
            if( (abs(reslim(i,n)-resdn).lt.eps) )then
!    .gt.      .or. (abs(reslim(i,n)).gt.resup) )then
              ionlim = i
              nonlim = n
c       LIMITER POINT OUTSIDE PLASMA
            else if( reslim(i,n).lt.RESDN )then
              if( seplim(i,n).le.fexpand )then
                fexpand = seplim(i,n)
                iexpand = i
                nexpand = n
              endif
c       LIMITER POINT INSIDE PLASMA
            else if( reslim(i,n).ge.RESDN )then
              if( seplim(i,n).gt.fshrink )then
                fshrink =  seplim(i,n) 
                ishrink = i
                nshrink = n
              endif
            endif
          enddo
      end do
 
c
c       LOGIC: IF THERE IS A LIMITER POINT INSIDE PLASMA, THEN MUST
c       SHRINK CONTOUR. OTHERWISE, IF THERE IS AT LEAST ONE LIMITER
c       POINT ON CONTOUR, AND ALL THE REST OUTSIDE, DO NOT CHANGE ANYTHING.
c       FINALLY, IF ALL LIMITER POINTS ARE OUTSIDE PLASMA, EXPAND PLASMA
c       TO OSCULATE WITH LIMITER
c
      if (ishrink .gt. 0) then
         gphifac = -sqrt(fshrink)
         ilim0 = ishrink
         nlim0 = nshrink
      else if (ionlim .gt. 0) then
         gphifac = 0.
         ilim0 = ionlim
         nlim0 = nonlim
      else if (iexpand .gt. 0) then
         gphifac = sqrt(fexpand)
         ilim0 = iexpand
         nlim0 = nexpand
      endif
 
      dlim_min = gphifac
      rlim_min = rlim(ilim0,nlim0)
      zlim_min = zlim(ilim0,nlim0)
 
c     overall damping in time, rsfac/sqrt(rsfac) = sqrt(rsfac)
      gphifac = gphifac/r01
      if (abs(gphifac) .gt. 0.04*one)
     1   gphifac = 0.04*gphifac/abs(gphifac)
      gphifac = 0.20*gphifac/sqrt(rsfac)
 
      end subroutine getlim


      subroutine cauchy(rbdy, zbdy, rubdy, zubdy, rlim, zlim, residue, 
     1   sep, distmax, ntheta, nlim)
      use vparams, only: twopi, dp, rprec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer, intent(in) :: ntheta, nlim
      real(rprec), intent(in) :: distmax
      real(rprec), dimension(ntheta), intent(in) :: 
     1     rbdy, zbdy, rubdy, zubdy
      real(rprec), dimension(nlim) :: rlim, zlim, residue, sep
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero=0, p5=0.5_dp, two=2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: n, i, imin, imax
      real(rprec), dimension(ntheta) :: dsq, dsepdu
      real(rprec), dimension(ntheta) :: x1u, y1u
      real(rprec) :: delu, dmin, delta_d, alpha, gam0
C-----------------------------------------------
c       Check that the points (rlim(i),zlim(i)) are inside boundary surface
c       using Cauchys theorem in "complex"-plane (for a fixed
c       toroidal plane, nphi=const. It is assumed that rbdy, zbdy are
c       extended around the full interval, 0-2pi, in theta.
c       with rbdy(1) = rbdy(ntheta+1) (i.e., ntheta intervals)
c
c       Because of numerical inaccuracies, instead of testing on
c       res = 0 (outside), res = 1 (inside), we use the test:
c       res >= .5, inside;  res < .5, outside
c**********************************************************************
c
c       LOCAL VARIABLE ARRAYS
c
c       dsq:    Distance squared between limiter point and plasma boundary
c       sep:    Minimum dsq for each limiter point
c    dsepdu:    .5* d(dsq)/d(theta)
c   residue:    Contour integral of 1/(X-rlim)in complex X=(R,Z) plane
c
      delu = twopi/ntheta
      dmin = 1.E-20_DP*distmax
 
      do n = 1, nlim
         residue(n) = zero
         x1u = rbdy(:ntheta) - rlim(n)
         y1u = zbdy(:ntheta) - zlim(n)
         dsq(:ntheta) = x1u*x1u + y1u*y1u
         dsepdu(:ntheta) = x1u*rubdy(:ntheta) + y1u*zubdy(:ntheta)
         residue(n) = residue(n) + sum((x1u*zubdy(:ntheta) - 
     1      y1u*rubdy(:ntheta))/(dsq(:ntheta)+dmin))
 
         residue(n) = residue(n)/ntheta
!
!        Find actual minimum distance from nth limiter point to boundary
!
         sep(n) = distmax
         do i = 1,ntheta
            if( dsq(i).le.sep(n) )then
               imin = i
               sep(n) = dsq(i)
            endif
         enddo
 
!        gamu = two*abs(dsepdu(imin))*delu
 
         if (dsepdu(imin) .le. zero) then
            imax = 1 + mod(imin,ntheta)
         else
            imax = imin
            imin = imax - 1
            if (imin .eq. 0) imin = ntheta
         endif
 
         delta_d = two*(dsepdu(imax)-dsepdu(imin))
         alpha = delta_d/delu
!        gamu = gamu/delta_d
!        sep(n) = sep(n) - p5*alpha*gamu**2
         if (alpha .ne. zero) 
     1      gam0 = 0.5_dp - (dsq(imax)-dsq(imin))/(alpha*delu**2)
         sep(n) = dsq(imin) - p5*alpha*(gam0*delu)**2
         if (sep(n) .lt. zero) sep(n) = zero
 
      end do
 
      end subroutine cauchy

           
      subroutine newprofil(phipog)
      use vmec_main
      use vacmod
      use realspace
      use vforces, lu => czmn, lv => crmn
      use vsvd
      use vspline
      use xstuff
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: icount, inodes, js
      integer, save :: idata(jchix), isize(jchix)       !INCREMENTAL UPDATES
      real(rprec), dimension(ipnodes + isnodes) :: datalsq, wten
      real(rprec), dimension(1+nrzt) :: r12sqr
      real(rprec) :: amat_lsq(isnodes+ipnodes,isnodes+ipnodes),
     1   djac_p(ipnodes,nmeasurements), djac_i(isnodes,nmeasurements)
      real(rprec), dimension(nmeasurements) :: datainput
      real(rprec) :: treconon, delt1, pfac0, aminout,
     1   aminin, ymin, pfactor, treconoff

C----------------------------------------------- 
!     IDATA: GIVES STARTING INDEX IN DATAINPUT ARRAY
!            FOR STORING EACH OF THE DATA TYPES
!     ISIZE: GIVES NUMBER OF DATA STARTING AT IDATA
!     INDEXING OF DJAC ARRAY (ASSUMES K STARTS AT 0)
 
      inodes = isnodes + ipnodes
 
      call second0 (treconon)
 
!     Unfreeze magnetic axis 
      if (iresidue.eq.0 .and. fsq*1.e6_dp.lt.one) iresidue = 1
      delt1 = one/real(ipedsvd,rprec)
      if (iresidue .eq. 0) delt1 = one
 
!
!       COMPUTE AVERAGE RADIAL FORCE BALANCE CONSTRAINT
!
      call radfor (pfac0)
      pfac = pfac + delt1*(pfac0 - pfac)
!
!       UPDATE PHI SCALE FACTOR (PHIFAC)
!       SCALE TOROIDAL FLUX TO MATCH PRESSURE WIDTH OR LIMITER
!
      aminout = max(rthommax,rstarkmax) - r00
      aminin = r00 - min(rthommin,rstarkmin)
      apres = (aminin*ipresin + aminout*ipresout)/(ipresin +ipresout)
      aminor = ((r00 - rinner)*ipresin + (router - r00)*ipresout)/
     1       (ipresin + ipresout)
      if (imatch_phiedge.ne.1 .and. ivac.gt.1 .or. imatch_phiedge.eq.3 
     1   .and. (.not.lfreeb)) then
         call newphi (phipog)
         call gettflux
      endif
 
      icount = 0
 
      if (.not.(mod(iter2 - iter1,ipedsvd).ne.0 .and. iequi.eq.0
     1     .and. iresidue.gt.0)) then
 
!
!       SETUP COMMON BLOCKS FOR FLUX-MATCHING ROUTINES
!

         if (iphidiam + nflxs + nbfldn.gt.0 .or. iequi.gt.0) then
            r12sqr(2:nrzt) = sqrt(armn_o(2:nrzt))
            call flux_init (phipog)
         endif
 
!
!       COMPUTE MATRIX ELEMENTS FOR THOMPSON SCATTERING DATA
!
         idata(ithom0) = icount + 1
         call getthom(djac_i(1,idata(ITHOM0)), djac_p(1,idata(ITHOM0)),
     1     datainput(idata(ITHOM0)), r1(1:,0), r1(1:,1), isize(ITHOM0))
         icount = icount + isize(ithom0)
 
!
!       COMPUTE MOTIONAL STARK EFFECT. THIS CALL ALSO INITIALIZES
!       THE ALSQ, DATALSQ ARRAYS AND SETS UP THE SPLINE NODES.
!
         idata(istark0) = icount + 1
         call getmse(djac_i(1,idata(ISTARK0)), djac_p(1,idata(ISTARK0)),
     1     datainput(idata(ISTARK0)), r1(1:,0), r1(1:,1),lu,
     2     lu(1+nrzt), zu(1:,0), zu(1:,1), phipog, isize(ISTARK0))
         icount = icount + isize(istark0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR DIAMAGNETIC FLUX LOOP
!
         idata(idiam0) = icount + 1
         call getdiam (djac_i(1,idata(idiam0)), djac_p(1,idata(idiam0))
     1      , datainput(idata(idiam0)), isize(idiam0))
         icount = icount + isize(idiam0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL POLOIDAL FLUXES
!
         idata(iflxs0) = icount + 1
         call getflux (djac_i(1,idata(iflxs0)), djac_p(1,idata(iflxs0))
     1      , datainput(idata(iflxs0)), r12sqr, clmn_e(1), clmn_o(1), 
     2      blmn_o(1), armn_o(1), blmn_e(1), azmn_o(1), isize(iflxs0))
         icount = icount + isize(iflxs0)
 
!
!       COMPUTE MATRIX ELEMENTS FOR EXTERNAL MAGNETIC FIELD MATCHING
!
         idata(ibrzfld) = icount + 1
         call getbfld (djac_i(1,idata(ibrzfld)), djac_p(1,idata(ibrzfld)
     1      ), datainput(idata(ibrzfld)), r12sqr, azmn_o(1), blmn_o(1), 
     2      clmn_e(1), clmn_o(1), armn_o(1), blmn_e(1), isize(ibrzfld))
         icount = icount + isize(ibrzfld)
 
!
!       SQUARE DATA MATRIX ELEMENTS AND STORE IN ALSQ
!
 
         if (icount .gt. nmeasurements) stop 'icount>nmeasurements'
         if (iequi .eq. 0) then
 
            call sgemvmm (djac_i, djac_p, amat_lsq, datainput, datalsq, 
     1         wten, icount, isnodes, ipnodes, inodes)
 
!
!       COMPUTE IOTA, PRESSURE SPLINE COEFFICIENTS
!
            call set_dual (datalsq, hstark, ystark, y2stark, hthom, 
     1         ythom, y2thom, wten, amat_lsq, isnodes, ipnodes, inodes)
 
            if (.not.lpprof) then
               ymin = minval(ythom(1:ipnodes))
               ythom(:ipnodes) = ythom(:ipnodes) - ymin
            endif
 
!
!       COMPUTE IOTA, PRESSURE AT R(js) FROM SPLINED INPUT
!       DATA ALONG THE MIDPLANE
!
            call splint (sknots, ystark, y2stark, isnodes, sqrts, 
     1         isplinef, zero, ns)
            call splint (sknots, ystark, y2stark, isnodes, shalf(2), 
     1         isplineh(2), zero, ns1)
            call splint (pknots, ythom, y2thom, ipnodes, sqrts, psplinef
     1         , zero, ns)
            call splint (pknots, ythom, y2thom, ipnodes, shalf(2), 
     1         psplineh(2), zero, ns1)
 
            pfactor = dmu0*pthommax             !!*pfac moved to getthom
            do js = 1,ns
              isplinef(js) = isplinef(js) - iotaf(js)
              isplineh(js) = isplineh(js) - iotas(js)
              psplinef(js) = pfactor*psplinef(js) - presf(js)
              psplineh(js) = pfactor*psplineh(js) - mass(js)
            end do  
 
         endif                     ! iequi>0
!
!       COMPUTE CHISQ
!
         call chisq (djac_i, djac_p, datainput, idata, isize, icount)

      endif                        ! mod(iter2-iter1,ipedsvd) == 0

      if (iequi .eq. 0) then
!
!       UPDATE PRESSURE SPLINE AND ESTABLISH INTERNAL
!       (CODE) UNITS OF PRESSURE. P(real) = dmu0 * pthommax * P(splined)
!       WHERE P(code-units) = mu0 * P(real)
!       SMOOTH TIME VARIATION OF PROFILES
!
         do js = 1,ns
           iotaf(js) = iotaf(js) + delt1*isplinef(js)
           iotas(js) = iotas(js) + delt1*isplineh(js)
           presf(js) = presf(js) + delt1*psplinef(js)
           mass(js)  = mass(js)  + delt1*psplineh(js)
         end do
      endif
 
!
!     STORE CHISQ
!
      call store_chisq
 
!
!       OPTIMIZE MAGNETIC AXIS POSITION BY MINIMIZING RMS ERROR
!       RETURNS RAXMSE AS UPDATED GUESS FOR NEW AXIS POSITION
!       TO BE USED IN SUBROUTINE RESIDUE
!
      call axisopt (fsq, r00, iresidue, ivac)
 
!       Compute force to fix axis at RAXMSE
      grmse = -0.05*(r00 - raxmse)
 
      call second0 (treconoff)
      timer(6) = timer(6) + (treconoff - treconon)
 
      end subroutine newprofil
      

      subroutine flux_init(phipog)
      use vmec_main
      use vmec_params, only: signgs
      use vforces, only : r12=>armn_o, gsqrt=>azmn_o, orsq=>blmn_o
      use vsvd
      use realspace
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: l, js
      real(rprec), external :: dot_g
C-----------------------------------------------
 
!
!     COMPUTE OBSERVATION POINT - INVARIANT FUNCTIONS OF RADIUS
!     CURRENT = PHIP * INTEGRAL(0,2pi)[ guu / gsqrt]   (on full mesh)
!     RM2     = < R**(-2) >
!     VRM2    = V` * RM2    (V` = 2pi * VP)
!     ORSQ    = SQRT(G) * R**(-2)   (on half mesh)
!
!     MUST HAVE GONE THROUGH NEWPROFILE DETERMINATION OF IOTA AT
!     LEAST ONCE, OTHERWISE IOTAS IS UNDEFINED!
!
      if (iresidue .le. 0) return 
      current(1) = zero
      presint(1) = one
 
      do l = 2,nrzt-1
        orsq(l) = p5*( phipog(l) + phipog(l+1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      enddo
      do l = ns,nrzt,ns
        orsq(l) = ( c1p5*phipog(l) - p5*phipog(l-1) ) *
     1  (ru0(l)*ru0(l) + zu0(l)*zu0(l))
      enddo

      do js = 2, ns
         current(js) = twopi*DOT_G(nznt,orsq(js),ns,wint(js),ns)
         presint(js) = one
      end do
 
      do l = 2, nrzt
         orsq(l) = gsqrt(l)/r12(l)**2
      end do
      do js = 2, ns
         vrm2(js) = twopi*DOT_G(nznt,orsq(js),ns,wint(js),ns)
         rm2(js) = vrm2(js)/(twopi*signgs*vp(js))
         ovrm2(js) = one/vrm2(js)
         ochip(js) = one/(phip(js)*iotas(js))
         presph(js) = presf(js) - presf(js - 1)
      end do
 
      end subroutine flux_init

      
      subroutine getbfld(amat_i, amat_p, data_array, r12sqr, 
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12, kcbfld)
      use vmec_main
      use vmec_params, only: signgs
      use vsvd
      use realspace, only: wint
      use vspline, only: hthom, hstark
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcbfld
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(nrzt) :: r12sqr,
     1  gsqrt, orsq, gobsr1, gobsz1, r12, z12
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, n1, m1, iobs, isym, iobsmax, iloop, 
     1   msym, nsym, indexbfld, l, lk
      real(rprec), dimension(:), allocatable :: gobsz2, gobsr2
      real(rprec):: tpisg, sumir, sumiz, wscaleb, sumpr,
     1    sumpz, deltab, plasbfld, coscoil, sincoil, t2, t1
C----------------------------------------------
       
!     IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET
      kcbfld = 0
      if (iresidue.le.0 .or. (nbfldn.eq.0 .and. iequi.eq.0))return
      
      allocate( gobsz2(nrzt), gobsr2(nrzt), stat=l)
      if (l .ne. 0) stop 'allocation problem in getbfld'
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!       GOBSR1 = BR, GOBSZ1 = BZ 
!       IF A (BR,BZ) OR (BRHO,BTHETA) PAIR, JUST CALL GRNBFLD ONCE
!
      if (any(rm2(2:ns) .eq. zero)) stop 'rm2 = 0'
      tpisg = twopi * signgs                !Positive volume integral
      do n1 = 1, nbsets
         do m1 = 1, nbcoils(n1)
            isym = needbfld(m1,n1)
            if (isym.eq.needit .or. iequi.eq.1) then
               call grnbfld(r12sqr,r12,z12,gobsr1,gobsz1,nrzt,m1,n1)
               do l = 2,nrzt
                 gobsr2(l) = gobsr1(l)*orsq(l)
                 gobsr1(l) = gobsr1(l)*gsqrt(l)
                 gobsz2(l) = gobsz1(l)*orsq(l)
                 gobsz1(l) = gobsz1(l)*gsqrt(l)
               end do
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!
               do js = 2, ns
                  sumir = zero
                  sumiz = zero
                  sumpr = zero
                  sumpz = zero
                  do lk = js,nrzt,ns
                    sumir = sumir + gobsr2(lk)*wint(lk)
                    sumpr = sumpr + gobsr1(lk)*wint(lk)
                    sumiz = sumiz + gobsz2(lk)*wint(lk)
                    sumpz = sumpz + gobsz1(lk)*wint(lk)
                  enddo
                  imb(js,m1,n1,1) = tpisg*sumir
                  pmb(js,m1,n1,1) = (-tpisg*sumpr) + imb(js,m1,n1,1)
     1               /rm2(js)
                  imb(js,m1,n1,2) = tpisg*sumiz
                  pmb(js,m1,n1,2) = (-tpisg*sumpz) + imb(js,m1,n1,2)
     1               /rm2(js)
               end do

            else if (isym .eq. ISAMECOIL) then       !Same coil position as previous coil
              do js = 2,ns
                imb(js,m1,n1,1) = imb(js,m1-1,n1,1)
                pmb(js,m1,n1,1) = pmb(js,m1-1,n1,1)
                imb(js,m1,n1,2) = imb(js,m1-1,n1,2)
                pmb(js,m1,n1,2) = pmb(js,m1-1,n1,2)
              enddo
            endif
          enddo   !m1
        enddo   !n1

!
!       CHECK FOR SYMMETRIC COIL (MAY BE IN DIFFERENT COIL SET,
!       SO HAD TO MOVE OUT OF M1,N1 LOOP ABOVE)
!
        do n1 = 1,nbsets
          do m1 = 1,nbcoils(n1)
            isym = needbfld(m1,n1)
            if (isym .ge. ISYMCOIL) then
              msym = 1 + (isym-1)/nbsets
              nsym = 1 + mod(isym-1,nbsets)
              do js = 2,ns                     !BR(-Z) = -BR(Z), BZ(-Z) = BZ(Z)
                imb(js,m1,n1,1) =-imb(js,msym,nsym,1)
                pmb(js,m1,n1,1) =-pmb(js,msym,nsym,1)
                imb(js,m1,n1,2) = imb(js,msym,nsym,2)
                pmb(js,m1,n1,2) = pmb(js,msym,nsym,2)
              enddo
            endif
          enddo
        enddo

!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        do 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at end
          do n1 = 1, nbsets
            iobsmax = nbfld(n1)
            if (iloop .eq. 1) iobsmax = nbcoils(n1)
            if (iobsmax .gt. 0) then
              do 1000 iobs = 1, iobsmax
                indexbfld = indxbfld(iobs,n1)
                if (iloop .eq. 1) indexbfld = iobs
                if (indexbfld .le. 0) goto 1000
                coscoil = cos( abcoil(indexbfld,n1) )
                sincoil = sin( abcoil(indexbfld,n1) )
                do js = 2,ns
                  pmb(js,0,n1,1) = ochip(js) *
     >            (pmb(js,indexbfld,n1,1)*coscoil +
     >             pmb(js,indexbfld,n1,2)*sincoil)
                  imb(js,0,n1,1) = ovrm2(js) *
     >            (imb(js,indexbfld,n1,1)*coscoil +
     >             imb(js,indexbfld,n1,2)*sincoil) 
                end do

                     if (iloop .eq. 0) then
                        deltab = plbfld(indexbfld,n1)
                        kcbfld = kcbfld + 1
 
                        call splinint (imb(1,0,n1,1), current, 
     1                     amat_i(1,kcbfld), hstark, u_ib, u1_ib, 
     2                     w_ib, w1_ib, nk_ib, isnodes, intder, ns)
 
                        call splinint (pmb(1,0,n1,1), presint, 
     1                     amat_p(1,kcbfld), hthom, u_pb, u1_pb, w_pb, 
     2                     w1_pb, nk_pb, ipnodes, intder, ns)
 
                        wscaleb = one/sigma_b(indexbfld,n1)
                        data_array(kcbfld) = wscaleb*deltab
                        t2 = dmu0*pthommax      !!*pfac moved to getthom
 
                        amat_i(:,kcbfld) = wscaleb*amat_i(:,kcbfld)
                        wscaleb = wscaleb*t2
                        amat_p(:,kcbfld) = wscaleb*amat_p(:,kcbfld)
 
                     else      !Store plasma fluxes in EXTFLX for output
                        plasbfld = zero
                        do js = 2, ns
                           t1 = current(js)*iotaf(js) - current(js-1)*
     1                        iotaf(js - 1)
                           plasbfld = plasbfld + pmb(js,0,n1,1)*
     1                        presph(js) + imb(js,0,n1,1)*t1
                        end do
                        plbfld(iobs,n1) = plasbfld
                     endif
 1000         continue
            endif
          enddo       !n1
 2000   continue

      deallocate( gobsz2, gobsr2, stat=l)

      end subroutine getbfld

      
      subroutine getdiam(amat_i, amat_p, data_array, kcdiam)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vforces, only : r12=>armn_o, ru12=>azmn_e
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcdiam
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ilimit = 2
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, lk, l
      real(rprec), dimension(ns) :: gp, gi, gip
      real(rprec), dimension(isnodes) :: amat2_i
      real(rprec) :: wdiam, z12, tv, ti, t2, sum1
C-----------------------------------------------
 
 
      kcdiam = 0
      if (iphidiam.eq.0 .or. iresidue.lt.ilimit) return 
!
!       COMPUTE FIT TO DIAMAGNETIC SIGNAL, USING EQUILIBRIUM RELATION
!       (modified 7/96 by SPH)
!
!       PHI-DIAMAG = 2*pi*INT[ Gp dp/ds + Gi d(<Bu>)/ds ]
!
!       where
!
!       Gp = Ru * Z * <sqrt(g)> /(R * phip)
!       Gi = Ru * Z * iota / R
!
 
      kcdiam = kcdiam + 1
c-7/96  dNewPhiedge = signgs*twopi*hs*SSUM_1(ns1,phip(2),1)
c-7/96  VacPhiedge  = signgs*bsubvvac*hs*SSUM_1(ns1,vrm2(2),1)
c-7/96  delphid0    = VacPhiedge - dNewPhiedge
 
      wdiam = one/sigma_delphid
      gp(1) = zero
      gi(1) = zero
      do js = 2, ns
         gp(js) = zero
         do lk = 1, nznt
            l = js + ns*(lk - 1)
            z12 = .5_dp*(z1(l,0)+z1(l-1,0)+shalf(l)*(z1(l,1)+z1(l-1,1)))
            gp(js) = gp(js) + ru12(l)*z12/r12(l)*wint(l)
         end do
      end do
!
!       NOTE: gip terms comes from linearizing the iota*d/ds[current*iota]
!             terms
!
      do js = 2, ns
         tv = twopi*vp(js)/phip(js)
         ti = -gp(js)*signgs*wdiam
         gi(js) = ti*iotas(js)
         gp(js) = -gp(js)*tv*wdiam
         gip(js) = ti*(current(js)*iotaf(js)-current(js-1)*iotaf(js-1))
      end do
 
      call splinint (gi, current, amat_i(1,kcdiam), hstark, u_ib, u1_ib
     1   , w_ib, w1_ib, nk_ib, isnodes, intder, ns)
 
      call splinint (gip(2), current(2), amat2_i, hstark, u_ia, u1_ia, 
     1   w_ia, w1_ia, nk_ia, isnodes, intfun, ns1)
 
      call splinint (gp, presint, amat_p(1,kcdiam), hthom, u_pb, u1_pb, 
     1   w_pb, w1_pb, nk_pb, ipnodes, intder, ns)
 
      amat_i(:isnodes,kcdiam) = amat_i(:isnodes,kcdiam) + amat2_i(:
     1   isnodes)
 
      t2 = dmu0*pthommax                        !!*pfac moved to getthom
      amat_p(:ipnodes,kcdiam) = t2*amat_p(:ipnodes,kcdiam)
      sum1 = sum(iotas(2:ns)*gip(2:ns))
      data_array(kcdiam) = wdiam*phidiam + sum1
      if (iequi .eq. 0) then
!
!       Eliminate p variation until well-converged
!
!@        do i = 1,ipnodes
!@          data_array(kcdiam) = data_array(kcdiam) -
!@     >    amat_p(i,kcdiam)*ythom(i)
!@          amat_p(i,kcdiam) = 0.
!@        end do
 
!
!       FINAL OUTPUT (ALSO USE FOR DEBUGGING)
!
      else
!
!       Integrate by parts
!
         delphid = gp(ns)*presf(ns) + gi(ns)*iotaf(ns)*current(ns) - 
     1      gp(2)*presf(1) - gi(2)*current(1)*iotaf(1)
         do js = 2, ns1
            delphid = delphid - presf(js)*(gp(js+1)-gp(js)) - iotaf(js)*
     1         current(js)*(gi(js+1)-gi(js))
         end do
         delphid = delphid/wdiam
      endif
!@        do js = 2,ns
!@        end do
!@
!@        sumi = sum(amat_i(:isnodes,kcdiam)*ystark(:isnodes))
!@     >       - sum(amat2_i(isnodes)*ystark(:isnodes))
!@        sump = DOT_G(ipnodes,amat_p(1,kcdiam),1,ythom,1)
!@
!@        write(*,1212)delphid,(sumi+sump)/wdiam,delphid0
!@ 1212   format(' DelPhid = ',1pe10.3,' PhiD = ',1pe10.3,
!@     >   ' DelPhid0 = ',1pe10.3)
 
      end subroutine getdiam

      
      subroutine getflux(amat_i, amat_p, data_array, r12sqr, 
     1  gobser1, gobser2, orsq, r12, z12, gsqrt, kcflux)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vsvd
      use vspline
      use vparams, only: zero
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcflux
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(nrzt) :: r12sqr,
     1  gobser1, gobser2, orsq, r12, z12, gsqrt     
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js, l, iloop, iobsmax, iobs, index, indexflx, 
     1   isym, n1, lk
      real(rprec) :: t1, t2, tpisg, sign0, sumi, sump, delta
C-----------------------------------------------
 
!       IRESIDUE > 0, OTHERWISE FLUX_INIT NOT CALLED YET!
        kcflux = 0
        if (iresidue.le.0 .or. (nflxs.eq.0 .and. iequi.eq.0)) return

!
!       COMPUTES MATRIX ELEMENTS NEEDED TO RELATE OBSERVED FLUX
!       MEASUREMENTS TO CURRENT AND PRESSURE EXPANSION COEFFICIENTS
!       R12,Z12 ARE THE PLASMA R,Z COORDINATES AT THE HALF
!       RADIAL NODE POINTS
!


!
!       COMPUTE SYMMETRIZED PSI(R,Z)+PSI(R,-Z) FLUX "GREEN'S FUNCTION"
!
!
!       COMPUTE "GREEN'S FUNCTION" KERNEL ONLY FOR NEEDED OBSERVATION POINTS
!       (OR FOR FINAL OUTPUT IF IEQUI=1)
!
      tpisg = twopi*signgs                     !Positive volume integral
      do n1 = 1, nobser
         isym = needflx(n1)
         if (isym.eq.needit .or. iequi.eq.1) then
            call grnflx (r12sqr, r12, z12, gobser1, nrzt, n1)
            do l = 2,nrzt
              gobser2(l) = gobser1(l)*orsq(l)
              gobser1(l) = gobser1(l)*gsqrt(l)
            end do  
!
!       DO INTEGRAL OVER ANGLES (ALL INTEGRALS ARE FROM THETA=0,TWOPI)
!       IM = <G/R**2>, PM = <G(1/R**2/<R**-2> - 1)>
!
            do js = 2, ns
              sumi = zero
              sump = zero
              do lk = js ,nrzt, ns
                sumi = sumi + gobser2(lk)*wint(lk)
                sump = sump + gobser1(lk)*wint(lk)
              enddo
              im(js,n1) = tpisg*sumi
              pm(js,n1) = (-tpisg*sump) + im(js,n1)/rm2(js)
            end do

          else if( isym.ge.ISYMCOIL )then    !only for up-down symmetric plasmas
            do js = 2,ns
              im(js,n1) = im(js,isym)
              pm(js,n1) = pm(js,isym)
            enddo
          endif
        enddo    !n1 loop
!
!       COMPUTE SPLINE MATRIX ELEMENTS BY INTEGRATING OVER RADIUS
!
        do 2000 iloop = 0,iequi                !iequi = 0 normally, = 1 at end
          iobsmax = nflxs
          if( iloop.eq.1 )iobsmax = nobd + nobser
          do 1000 iobs = 1,iobsmax
            indexflx = indxflx(iobs)
            if( iloop.eq.1 )indexflx = iobs
            if( indexflx.le.0 )go to 1000
            do js = 2,ns
              pm(js,0) = zero
              im(js,0) = zero
            enddo

            do l = 1,4                !This could be halved by using symmetry
              index = iconnect(l,indexflx)
              if( index.ne.0 )then
                sign0 = 1.0
                if( index.lt.0 )then
                  sign0 = -sign0
                  index = -index
                endif
                do js = 2,ns
                  pm(js,0) = pm(js,0) + sign0*pm(js,index)
                  im(js,0) = im(js,0) + sign0*im(js,index)
                enddo
              endif
            enddo
          
            do js = 2,ns
              pm(js,0) = pm(js,0)*ochip(js)
              im(js,0) = im(js,0)*ovrm2(js)
            enddo

               if (iloop .eq. 0) then
                  kcflux = kcflux + 1
                  delta = plflux(indexflx)
 
                  call splinint (im, current, amat_i(1,kcflux), hstark, 
     1               u_ib, u1_ib, w_ib, w1_ib, nk_ib, isnodes, intder, 
     2               ns)
 
                  call splinint (pm, presint, amat_p(1,kcflux), hthom, 
     1               u_pb, u1_pb, w_pb, w1_pb, nk_pb, ipnodes, intder, 
     2               ns)
 
                  t1 = one/sigma_flux(indexflx)
                  data_array(kcflux) = t1*delta
                  amat_i(:,kcflux) = t1*amat_i(:,kcflux)
                  t2 = t1*dmu0*pthommax         !!*pfac moved to getthom
                  amat_p(:,kcflux) = t2*amat_p(:,kcflux)
 
 
               else            !Store plasma fluxes in PLFLUX for output
                  plflux(indexflx) = zero
                  do js = 2, ns
                     plflux(indexflx) = plflux(indexflx) + pm(js,0)*
     1                  presph(js) + im(js,0)*(current(js)*iotaf(js)-
     2                  current(js-1)*iotaf(js-1))
                  end do
               endif
 1000     continue
 2000   continue

      end subroutine getflux

      
      subroutine getmse(amat_i, amat_p, data_array, re, ro, lue, luo, 
     1   zue, zuo, phipog, kcstark)
      use vmec_main
      use vmec_params, only: signgs
      use realspace
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcstark
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(ns,nzeta,*) :: 
     1   re, ro, lue, luo, zue, zuo
      real(rprec), dimension(*) :: phipog
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: lt, i, js, ks, j, irnodes
      real(rprec), dimension(ns,ntheta2) :: luef, luof
      real(rprec) :: dlu, dzu, guu1, edgeiota, guu2, edgefactor, 
     1   wedge, t1
      real(rprec), save :: facedge 
C-----------------------------------------------
 
!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE MSE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (ISNODES, INCLUDING S=0 AND S=1)
!
!       THE RANGE OF MSE DATA IS EXTENDED TO INCLUDE S=1 BY USING THE
!       CURRENT MATCHING CONDITION TO CONSTRAIN IOTA(S=1)
!
!       COMING INTO THIS ROUTINE, THE RSTARK, DATASTARK HAVE BEEN
!       PREVIOUSLY ORDERED SO RSTARK(1) < RSTARK(2) < ....
!
 
!
!       COMPUTE FULL MESH VALUES OF LAMBDA
!
      call lamfull (luef, luof, lue, luo)
!
!       COMPUTE OUTER EDGE PITCH (IOTA) TO MATCH TOTAL CURRENT
!       IOTA(EDGE) = MU0 * IPLASMA/ 2 * PI *< Guu/SQRTg > PHIP(s=1)
!       NEED THIS TO SPLINE IOTA OVER FULL S-RANGE ( 0 .le. S .le.1 )
!
      guu1 = dot_product(c1p5*wint(ns:nrzt:ns)*guu(ns:nrzt:ns),
     1   phipog(ns:nrzt:ns))
      guu2 = dot_product(cp5*wint(ns-1:nrzt-1:ns)*guu(ns-1:nrzt-1:ns),
     1   phipog(ns-1:nrzt-1:ns))
 
      if (iresidue.eq.0 .or. iotaf(ns).eq.zero) then
         facedge = one
      else if (mod(iter2 - iter1,ipedsvd) .eq. 0) then
         facedge = (guu1*iotas(ns) - guu2*iotas(ns1))/(iotaf(ns)*(guu1
     1       - guu2))
      endif
      edgefactor = facedge*(guu1 - guu2)*signgs*twopi
      edgeiota = currv/edgefactor
 
      irnodes = max(0,imse) + 1
      lt = 1                                     !Outer R edge
      dlu = luef(ns,lt) + luof(ns,lt)
      wedge = (zue(ns,1,lt)+zuo(ns,1,lt))/(dlu*router)
      rstark(irnodes) = router
      datastark(irnodes) = wedge*edgeiota        !Edge pitch
      sigma_stark(irnodes) = abs(sigma_current*wedge/edgefactor)
 
!
!       THROW AWAY POINTS OUTSIDE GRID
!       NOTE: IF ONLY OUTER POINT KEPT, THE CALL TO SORT IS UNNECESSARY
!
      rsort0(:irnodes) = rstark(:irnodes)
      call sort_data (rsort0,isortr,irnodes)
      kcstark = 0
      do i = 1,irnodes
        j = isortr(i)
        if( ((rsort0(i).gt.rinner) .and.
     >       (rsort0(i).le.router)) .or. (iequi.ne.0) )then
           kcstark = kcstark+1
           rsort(kcstark) = rsort0(i)                        !kcstark <= i
           starkcal(kcstark) = datastark(j)                 !sorted data array
           qcalc(kcstark) = one/sigma_stark(j)                !qcalc = sorted weight array
        endif
      enddo
 
!
!       COMPUTE IOTA(0) FROM SLOPE AT RSTARK=R00
!
c04-96        kcstark = kcstark+1
c04-96        rsort(kcstark) = r00                !Magnetic axis (s=0)
c04-96        qcalc(kcstark) = 1.0/scstark
c04-96
c04-96        if( imse.gt.0 )then
c04-96        slope0 = 1.0
c04-96        call splint(rstark,ystark0,y2stark0,
c04-96     >  imse,r00,dum,slope0,1)
c04-96        starkcal(kcstark) = r00*slope0*luef(1,1)/dkappa
c04-96        else
c04-96c       EXTEND BOUNDARY POINTS TO INCLUDE FULL RANGE IN THETA
c04-96        starkcal(kcstark) = ai(0)
c04-96        endif
 
!
!       FIND S,THETA INDICES FOR RSORT ARRAY
!
      call findphi (re, ro, rsort, delse1, delso1, rmid, indexs1, 
     1   indexu1, indexr, kcstark)
 
!
!       COMPUTE MATRIX ELEMENTS FOR IOTA SPLINE NODES CORRESPONDING
!       TO ORDERED RSORT ARRAY ( = RSORT S )
!
      if (kcstark .gt. nmse) stop 'kcstark>nmse'
      call getspline (amat_i, sknots, hstark, delse1, hs, indexs1, 
     1   isorts, kcstark, isnodes)
 
!
!       MATCH TO STARK MSE DATA ACCORDING TO THE FORMULA:
!
!       Bz/Btor = IOTA*Zu/[ R*(1+LAMu) ]
!
!       NOTE: QCALC, DATA = STARKCAL CORRESPOND TO RSORT_S(I)
!       WITH INDEX KS = ISORTS(I) (INDEXED ON RSORT BEFORE IT WAS SORTED)
!       SAME IS TRUE FOR DELSE,O1, INDEXS1, INDEXU1
!
 
      islope = 0
      do i = 1, kcstark
c                     !Index BEFORE sorting on sknots (indexed on rsort)
         ks = isorts(i)
         js = indexs1(ks)
         lt = indexu1(ks)
!
!       COMPUTE WEIGHT(J) = Zu / (R * [1 + LAMu]), WHICH IS THE FACTOR
!       RELATING MSE PITCH = WEIGHT(J) * IOTA(J) TO ROTATIONAL TRANSFORM.
!
!       ON ENTRY INTO THIS LOOP,
!       QCALC = 1/SIGMA_STARK
!
         dlu = (one - delse1(ks))*luef(js,lt) + delse1(ks)*luef(js+1,lt
     1      ) + (one - delso1(ks))*sqrts(js)*luof(js,lt) + delso1(ks)*
     2      sqrts(js+1)*luof(js+1,lt)
         dzu = (one - delse1(ks))*zue(js,1,lt) + delse1(ks)*zue(js+1,1,
     1      lt) + (one - delso1(ks))*sqrts(js)*zuo(js,1,lt) + delso1(ks
     2      )*sqrts(js+1)*zuo(js+1,1,lt)
         stark_weight(ks) = dzu/(rsort(ks)*dlu)
 
         if (rsort(ks) .eq. router) icurrout = i
 
c04-96        if( rsort(ks).eq.r00 )then                        !IOTA(0)
c04-96          islope = i
c04-96          stark_weight(ks) = abs(wedge)                  !Need in
c04-96          starkcal(ks) = weight(ks)*starkcal(ks)        !Need for
c04-96          data_array(i) = starkcal(ks) * qcalc(ks)
c04-96          amat_i(1,i) = amat_i(1,i) * stark_weight(ks) * qcalc(ks)
c04-96         else
         data_array(i) = starkcal(ks)*qcalc(ks)
         t1 = stark_weight(ks)*qcalc(ks)
         amat_i(:isnodes,i) = t1*amat_i(:isnodes,i)
c04-96        endif
         if (i.eq.icurrout .and. qcalc(ks).eq.zero) stop 'CURR ERR'
      end do
 
      imse2 = kcstark
      amat_p(:,:kcstark) = zero
 
      end subroutine getmse

      
      subroutine gettflux
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: p01=1.e-2_dp, zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real(rprec) :: tfac, phifac0
C-----------------------------------------------
 
!
!       UPDATE PHIEDGE SCALE FACTOR BY THE FORMULA
!       FDOT/F = -2*(1 - apres/aminor), WHERE F = PHISCALE
!
      tfac = p01*rsfac
 
      if (imatch_phiedge .eq. 0) then
         if (aminor .eq. zero) then
            phifac0 = phifac
         else
           phifac0 = phifac*(apres/aminor)
         end if  
         gphifac = tfac*(phifac0 - phifac)
 
      else if (imatch_phiedge .eq. 2) then
         call getlim
         gphifac = tfac*phifac*gphifac
      endif
 
      end subroutine gettflux

      
      subroutine getthom(amat_i, amat_p, data_array, re, ro, kcthom)
      use vmec_main
      use vsvd
      use vspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer kcthom
      real(rprec), dimension(isnodes,*) :: amat_i
      real(rprec), dimension(ipnodes,*) :: amat_p
      real(rprec), dimension(*) :: data_array
      real(rprec), dimension(ns,nzeta,*) :: re, ro
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, ks
      integer, dimension(itse + 2) :: isortp
      real(rprec), dimension(itse) :: datalsq_p
      real(rprec), dimension(itse + 2) :: rgrid
      real(rprec):: sig_fac, t1
      real(rprec) :: rmat1u(itse)
      logical :: l1v(itse)
C-----------------------------------------------
 
!
!       THIS SUBROUTINE COMPUTES THE LEAST-SQUARES AMATRIX AND DATA
!       ARRAYS FOR MATCHING TO THE PRESSURE DATA AT EQUALLY-SPACED KNOTS
!       IN SQRT(PHI-FLUX) SPACE (IPNODES, INCLUDING S=0 AND S=1)
!
!       COMING INTO THIS ROUTINE, THE RTHOM, DATATHOM HAVE BEEN
!       PREVIOUSLY ORDERED SO RTHOM(1) < RTHOM(2) < ....
!
 
!       IF (LPOFR), user has input P(R,Z)
!       If (.NOT.LPOFR),  then user has input P(s), not P(R)
!
 
      if (.not.lpofr) then                       !p(R) or p(s) ?
 
!       CONSTRUCT RTHOM BASED ON P(s)
         rthom(:itse) = sthom(:itse)
 
         call pofs (re, ro, ns, rthom, itse)
         rthommax = rthom(itse)
         rthommin = rthom(1)
 
      endif
 
!
!       IF NO PRESSURE DATA, MAKE SURE CHISQ-THOM <= CHI_TARGET
!
      sig_fac = one
 
!
!       CONSTRUCT EVERYTHING BASED ON P(R)
!       FOR POINTS OUTSIDE GRID, SET R = either rmin,rmax
!
      kcthom = 0
      if (itse .gt. 0) then
         l1v(:itse) = .false.
         datalsq_p(:itse) = datathom(:itse)*pfac   !sorted data array
         sigma_thom(:itse) = sigma_thom(:itse)/sig_fac
         pcalc(:itse) = 1.0/sigma_thom(:itse)      !pcalc = sorted sigma array
         rmat1u(:itse) = rthom(:itse)
         where (rmat1u(:itse) .lt. rinner) 
            rmat1u(:itse) = rinner
         elsewhere
            l1v(:itse) = rmat1u(:itse) .gt. router
         end where
         where (l1v(:itse)) rmat1u(:itse) = router
         rgrid(kcthom+1:itse+kcthom) = rmat1u(:itse)
         kcthom = itse + kcthom
      endif
 
 
!
!       FIND S,THETA INDICES FOR GRIDDED R-THOM ARRAY (RGRID)
!
      call findphi (re, ro, rgrid, delse2, delso2, rmid, indexs2, 
     1   indexu2, indexr, kcthom)
 
!
!       COMPUTE MATRIX ELEMENTS FOR PRESSURE SPLINE NODES CORRESPONDING
!       TO ORDERED RGRID ARRAY
!
      call getspline (amat_p, pknots, hthom, delse2, hs, indexs2,
     1   isortp, kcthom, ipnodes)
 
!
!       MATCH PRESSURE SPLINE KNOTS TO THOMSON SCATTERING DATA
!
!       ON ENTRY INTO THIS LOOP, PCALC = 1/SIGMA_THOM
!
 
      do i = 1, kcthom
         ks = isortp(i)         !Index BEFORE sorting on pknots (indexed
 
         data_array(i) = datalsq_p(ks)*pcalc(ks)
         t1 = pthommax*pcalc(ks)
         amat_p(:ipnodes,i) = t1*amat_p(:ipnodes,i)
      end do
      if (.not.lpofr) rthompeak = rgrid(isortp(1))
 
      itse2 = kcthom
      amat_i(:,:itse2) = zero
 
      end subroutine getthom

      
      subroutine grnbfld(xsqr, xs, zs, br, bz, idim, nobs1, nobs2)
      use vsvd
      use vparams, only: one
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: four=4.0_dp, p5=0.5_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer idim, nobs1, nobs2
      real(rprec), dimension(idim) :: xsqr, xs, zs, br, bz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer j, i1, i2
      real(rprec) :: xt, zt, xtsqr, xt4, oxt, zrp, zrm, xvv, 
     1 fxu, sqrxu, qqp, qqm, delqp, delqm, yeqp, yeqm, brp, brm
C-----------------------------------------------
!
!       COMPUTE BR = (ZT-ZS)/RT/SQRT(4*RT*RS) * F1(k)
!               BZ = 1/SQRT(4*RT*RS)*[RS/RT * F2(k) - F1(k)]
!       WHERE   F1 = k/2pi[ (E(k) - K(k)) + q1(k)*E(k) ]
!               F2 = k/(2pi)  [ q1(k)*E(k) ]
!               q1 = .5*k**2/(1. - k**2)   [Most singular piece near k=1]
!               k**2 = 4*RT*RS/[(RT+RS)**2 + (ZT-ZS)**2]
!
      xt = rbcoil(nobs1,nobs2)
      zt = zbcoil(nobs1,nobs2)
      xtsqr = p5/rbcoilsqr(nobs1,nobs2)            !1/2 from symmetrizing
      xt4 = four*xt
      oxt = one/xt

      do j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xtsqr/xsqr(j)
        qqp = fxu/(xvv + zrp*zrp)
        qqm = fxu/(xvv + zrm*zrm)
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = int(qqp*odqk2) + 1
        i2 = int(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        delqp = qqp - qsq(i1)
        delqm = qqm - qsq(i2)
        yeqp = (yeq(i1) + delqp*dyeq(i1))/(one - qqp)
        yeqm = (yeq(i2) + delqm*dyeq(i2))/(one - qqm)
        brp = yek(i1) + delqp*dyek(i1)
        brm = yek(i2) + delqm*dyek(i2)
        br(j) = sqrxu*oxt*(zrp*(brp+yeqp) + zrm*(brm+yeqm))
        bz(j) = sqrxu*((xs(j)*oxt-1.0)*(yeqp + yeqm) - (brp + brm))
      enddo

      end subroutine grnbfld

      
      subroutine grnflx(xsqr, xs, zs, ansp, idim, nobs)
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: p5 =0.5_dp, four=4.0_dp
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer idim, nobs
      real(rprec), dimension(idim) :: xsqr, xs, zs, ansp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i1, i2, j
      real(rprec) :: xt, zt, xtsqr, xt4, 
     4   zrp, zrm, xvv, fxu, sqrxu, qqp, qqm
C-----------------------------------------------
!
!       EVALUATES "GREEN'S" FUNCTION FOR POLOIDAL FLUX AT INTERIOR
!       POINTS XS,ZS AND OBSERVATION POINT XT,ZT (ANSP)
!       (RECALL THETA INTEGRATION ONLY FROM ZERO TO PI, SO NEED
!        TO REFLECT ZS TO -ZS, AT LEAST IN UP-DOWN SYMMETRIC CASE)
!
      xt = xobser(nobs)
      zt = zobser(nobs)
      xtsqr = p5*xobsqr(nobs)        !1/2 factor from averaging up,down
      xt4 = four*xt
      do j = 2,idim
        zrp = zt - zs(j)
        zrm = zt + zs(j)
        xvv =(xt + xs(j))**2
        fxu = xs(j)*xt4
        sqrxu = xsqr(j)*xtsqr
        qqp = fxu/(xvv + zrp*zrp)                !k**2 for zplasma > 0
        qqm = fxu/(xvv + zrm*zrm)                !k**2 for zplasma < 0
!
!       WHICH INDEX LIES BELOW ?
!
        i1 = int(qqp*odqk2) + 1
        i2 = int(qqm*odqk2) + 1
!
!       LINEAR INTERPOLATION
!
        ansp(j)  =  sqrxu *( ( yf(i1)+(qqp-qsq(i1))*dyf(i1) )
     >   + ( yf(i2)+(qqm-qsq(i2))*dyf(i2) ) )
      enddo
 
      end subroutine grnflx

      
      subroutine pofs(re, ro, ns, rthom, itse)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer ns, itse
      real(rprec), dimension(ns) :: re, ro
      real(rprec), dimension(*) :: rthom
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: j, js, i
      real(rprec), dimension(ns) :: s
      real(rprec) :: sqrjs, rt2, ds
C-----------------------------------------------
!
!       Interpolate Rmid = Re + Ro to get R(s)
!
!       THIS ROUTINE INTERPOLATES THE RTHOM "S-SPACE" ARRAY
!       ONTO THE INSTANTANEOUS RMID ARRAY
!       ON INPUT, RTHOM IS THE S-VALUE ARRAY FOR THOMPSON DATA
!       ON OUTPUT,RTHOM IS THE CORRESPONDING (THETA=0) RMID ARRAY
!
      do j = 1, ns
         s(j) = real(j - 1,rprec)/(ns - 1)
      end do
 
      js = 1
      do i = 1, itse
         rt2 = rthom(i)
  100    continue
         if (rt2.ge.s(js) .and. rt2.le.s(js+1)) then
            ds = (rt2 - s(js))/(s(js+1)-s(js))
            sqrjs = sqrt(rt2)
            rthom(i) = re(js) + (re(js+1)-re(js))*ds + sqrjs*
     1         (ro(js)+(ro(js+1)-ro(js))*ds)
         else
            js = js + 1
            if (js < ns) go to 100
         endif
      end do
 
      end subroutine pofs

      
      subroutine radfor(pfac0)
      use vmec_main
      use vmec_params, only: signgs
      use vforces, only : r12=>armn_o, gsqrt=>azmn_o, gor=>clmn_o
      use realspace
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) pfac0
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec) :: p05 = 0.05, p5 = 0.5_dp, c1p5 = 1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: js
      real(rprec), dimension(ns) :: vpres
      real(rprec) :: delpres, pedge, t1
      real(rprec), external :: dot_g
C-----------------------------------------------
 
 
!
!       COMPUTE VPRES, NEEDED FOR F00 PRESSURE BALANCE
!
      gor(2:nrzt) = gsqrt(2:nrzt) / r12(2:nrzt)
      do js = 2, ns
         vpres(js) =signgs*DOT_G(nznt,gor(js),ns,wint(js),ns)
      end do
 
      pedge = c1p5*pres(ns) - p5*pres(ns1)
      pressum0 = dot_product(wint(ns:nrzt:ns)*zu0(ns:nrzt:ns),
     1   r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      pressum0 = signgs*pedge*pressum0
      pressum0 = pressum0 + hs*dot_product(vpres(2:ns),pres(2:ns))
 
      if (pressum0 .eq. zero) pressum0 = one
 
      pfac0 = pfac
      if (iresidue .ge. 3) return                    !!AXIS MOVED BY FSQR IN RESIDUE
!
!       COMPUTE AVERAGE FORCE BALANCE CONSTRAINT FOR FIXING R(0)
!       (INTEGRAL OF RADIAL FORCE BALANCE,M=0,N=0, OVER S)
!
 
      if (1.e6_dp*fsq .le. one) then
         delpres = 0.
         delpres = -fsqsum0/pressum0
         t1 = abs(delpres)
         if (t1 .gt. p05) delpres = p05*delpres/t1   !!Wait til close
         pfac0 = pfac*(one + delpres)
      endif
 
      end subroutine radfor

      
      subroutine readrecon
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      character*(50), dimension(2), save :: raxis_message
      character*(50), dimension(4), save :: phiedge_message
C-----------------------------------------------
      data raxis_message/'Magnetic axis position fixed', 
     1   'Magnetic axis position optimized'/
      data phiedge_message/
     1   'Phiedge determined to match pressure minor radius', 
     2   'Phiedge matched to input value', 
     3   'Phiedge determined by limiter position', 
     4   'Phiedge determined by Ip'/
 
      iphidiam = 0
      if (imse > nmse) stop 'IMSE>NMSE'
      if (itse > ntse) stop 'ITSE>NTSE'
      if ((imse>0 .or. nflxs>0 .or. nbfldn>0) .and. itse.ge.0) then
         iresidue = 0
      else
         lrecon = .false.
      end if   
 
      if (.not.lrecon) return                   !No reconstruction matching
      ncurr = 0                                  !Just to be safe
      if (sigma_current .ge. cbig) stop 'SIGMA_CURRENT missing'
      if (sigma_delphid .ge. cbig) print *, ' SIGMA_DELPHID missing'
      if (sigma_current < zero) sigma_current = abs(sigma_current*
     1   curtor)
      if (sigma_delphid < zero) sigma_delphid = abs(sigma_delphid*
     1   phidiam)
      write (nthreed, 150)
  150 format(/' DATA MATCHING PARAMETERS: ',/,1x,35('-'))
      write (nthreed, 155) imse, itse, nflxs, nobser, nobd, nbfldn, 
     1   nbcoilsn, sigma_current, 1.e3_dp*sigma_delphid, tensp, tensi,
     2   tensi2, fpolyi, mseangle_offset, presfac, pres_offset, lpofr
      write (nthreed, 152) mseangle_offsetm
  152 format('mse-angleM offset',/,f13.3)
  155 format('   imse       itse      nflxs     nobser       nobd',
     1'     nbfldn    nbcoils  sigma_current(A)   sigma_delphid(mWb)',/,
     2   i7,6i11,3x,1pe15.3,4x,0pf17.3,/,
     3   '    tension(p)   tension(i)  tension2(i)  fpolyi  ',
     4   'mse-angle offset  pres scale factor pressure offset  lpofr',/,
     5   3f13.3,f9.3,f18.3,f19.3,f16.3,6x,l1)
      write (nthreed, 200)
  200 format(/,' LEGEND',/,1x,6('-'))
 
      if (curtor < cbig) then
         write (nthreed, 210) 1.e-6_dp*curtor
      else
         write (nthreed, *) 'Need toroidal plasma current'
         stop 15
      endif
  210 format(' Matching to toroidal current = ',f10.3,' [MA]')
      sigma_current = dmu0*sigma_current
      if (nflxs > 0) then
         write (nthreed, *) 'Fitting ', nflxs, 
     1      ' external flux loop measurements'
      else
         write (nthreed, *) 
     1      'Not fitting external flux loop measurements.'
      endif
      if (phidiam<cbig .and. sigma_delphid<cbig) then
         iphidiam = 1
         write (nthreed, 220) 1.e3_dp*phidiam
      else
         write (nthreed, *) 'No fit to diamagnetic flux'
      endif
  220 format(' Fitting diamagnetic flux     = ',f10.3,' [mWb]')
      write (nthreed, *) raxis_message(iopt_raxis+1)
      write (nthreed, *) phiedge_message(imatch_phiedge+1)

      end subroutine readrecon

      
      subroutine sgemvmm(amat_i, amat_p, amatsq, b, bsq, wten,
     1   mdata, niota, npres, nots)
      use kind_spec
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer mdata, niota, npres, nots
      real(rprec), dimension(niota,mdata) :: amat_i
      real(rprec), dimension(npres,mdata) :: amat_p
      real(rprec), dimension(nots,nots) :: amatsq
      real(rprec), dimension(mdata) :: b
      real(rprec), dimension(nots) :: bsq, wten
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      real(rprec), parameter :: zero = 0.0_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, j, ioff, joff
C-----------------------------------------------

      amatsq = zero
      bsq    = zero
 
!
!       INITIALIZE IOTA, PRESSURE DIAGONAL ELEMENTS (ALREADY 'SQUARED')
!
!
!       COMPUTE LOWER DIAGONAL HALF OF SQUARE OF MATRIX
!       A(trans)*A and A(trans)*B
!       BY ADDING CONTRIBUTIONS FROM EXTERNAL MAGNETICS SIGNALS
!
 
!
!       FIRST UPPER LEFT NIOTA X NIOTA BLOCK
!
      do i = 1, niota
         bsq(i) = bsq(i) + sum(b*amat_i(i,:))
         do j = 1, i
            amatsq(i,j) = amatsq(i,j) + sum(amat_i(i,:)*amat_i(j,:))
         end do
      end do
 
!
!       LOWER NPRES X NIOTA BLOCK, NPRES X NPRES BLOCK
!
      do i = 1, npres
         ioff = i + niota
         bsq(ioff) = bsq(ioff) + sum(b*amat_p(i,:))
         do j = 1, niota
            amatsq(ioff,j) = amatsq(ioff,j) + 
     1                       sum(amat_p(i,:)*amat_i(j,:))
         end do
         do j = 1, i
            joff = j + niota
            amatsq(ioff,joff) = amatsq(ioff,joff) +
     1                          sum(amat_p(i,:)*amat_p(j,:))
         end do
      end do
 
      do i = 1, nots
         wten(i) = amatsq(i,i)
         amatsq(1:i-1,i) = amatsq(i,1:i-1)
      end do
 
 
      end subroutine sgemvmm
      

      subroutine smoothdata(nwout)
      use vmec_main
      use vsvd
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: nwout
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      integer, parameter :: ndata_elems = 11
c-----------------------------------------------
c   L o c a l   V a r i a b l e s
c-----------------------------------------------
      integer :: npts, i, ndata1
      real(rprec), dimension(2*ns-1) :: hmid,ymid,y2mid,
     1  wmid,tenmid,requal
      real(rprec), dimension(:), pointer :: data

c-----------------------------------------------
!
!       spline output data onto equally-spaced r-mesh for plotting
!
      npts = 2*ns - 1
      if (npts .le. 1) return
      
      do i = 1, npts
         if (i .le. ns) curint(i) = twopi/dmu0*buco(ns + 1 - i)
         if (i .gt. ns) curint(i) = twopi/dmu0*buco(i - ns + 1)
      end do
 
      do i = 1, npts
         wmid(i) = 1.0
         tenmid(i) = 0.1
         requal(i) = rmid(1) + ((i - 1)*(rmid(npts)-rmid(1)))/
     1      (npts - 1)
      end do
 
      do ndata1 = 1, ndata_elems
         select case (ndata1) 
         case (1) 
            data => datamse
         case (2) 
            data => qmid
         case (3) 
            data => shear
         case (4) 
            data => presmid
         case (5) 
            data => alfa
         case (6) 
            data => curmid
         case (7) 
            data => curint
         case (8) 
            data => psimid
         case (9) 
            data => ageo
         case (10) 
            data => volpsi
         case (11) 
            data => phimid
         end select
         call setspline (rmid, wmid, data, hmid, ymid, y2mid, tenmid, 
     1      tenmid(1), npts, natur)
         call splint (rmid, ymid, y2mid, npts, requal, data(1), 
     1      zero, npts)
      end do
!
!     write out splined data
!
      write (nwout, 703) (datamse(i),requal(i),qmid(i),shear(i),presmid
     1   (i),alfa(i),curmid(i),i=1,npts)
      write (nwout, 703) (rsort(i),atan(datastark(isortr(i)))/dcon,abs(
     1   qmeas(i)),i=1,imse2 - 1)
      write (nwout, 703) (rthom(i),datathom(i),i=1,itse)
  703 format(5e20.13)
 
      if (lmac) then
        write (nmac, 705)
  705 FORMAT(//,' FOLLOWING DATA EQUALLY SPACED IN R-MIDPLANE'//,
     1   '        RMID       J-PHI       SHEAR        QMID',
     2   '   MSE-PITCH     PRESMID         PSI        AMID',
     3   '      VOLUME         PHI',/,
     4   '         [M]    [A/M**2]                        ',
     5   '       [Deg]        [Pa]        [Wb]         [M]',
     6   '      [M**3]        [Wb]',/)
        write (nmac, 707) (requal(i),curmid(i),shear(i),qmid(i),
     1    datamse(i),presmid(i),psimid(i),ageo(i),volpsi(i),
     2    phimid(i),i=1,npts)
        write (nmac, 709) phimid(2*ns-1), psimid(2*ns-1)
      end if
  707 format(1p10e12.3)
  709 format(/,' phi-edge =',t16,1pe12.3,t40,'psi-edge =',t56,1pe12.3)
       
      end subroutine smoothdata
