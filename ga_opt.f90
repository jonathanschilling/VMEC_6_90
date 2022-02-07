#!/bin/sh
#---------------------------------------------------------------------
cat > temp.c << "EOC"
cat > ga_optim.f << "EOF"
      subroutine GA_preset
      use kind_spec
      use ga_mod
      implicit none

c
      kountmx=5
      irestrt=0
      itourny=0
      ielite=0
      iunifrm=0
      iniche=0
      iskip=0
      iend=0
      nchild=1
      unique_ind=0
      ibound=0
      nichflg(1:nparmax)=1
      parmax = 0
      parmin = 0
      nposibl = 0
      microga=0
      save_space = .false.

      end subroutine GA_preset


       subroutine GA_driver(fcn, n_opt, n_var, x, fvec, tol, eps,
     1     num_iter_opt, max_processors, filename, info, lwa, lrestart )

       use kind_spec
       use ga_mod
       use gade_mod, ONLY: ga_de
c      use optim_params, ONLY : ldiag_opt
       use system_mod
c      use read_namelist_mod
       use safe_open_mod
       implicit none
#ifdef MPI_OPT
       include 'mpif.h'                                       !mpi stuff
#endif
       integer :: n_opt, n_var, info, lwa, num_iter_opt, max_processors 
       real(rprec), dimension(n_opt), target :: fvec
       real(rprec), dimension(n_var) :: x
       real(rprec), dimension(n_var) :: partemp
       real(rprec) :: tol, eps, chi_sq, tmp                   !ga_evaluate
       external fcn
       character*(*) ::  filename
       logical :: lrestart

       character*(100) ::  temp
       integer, parameter :: master=0
       integer :: num_iter_max, iunit, myid=master
       integer :: i, j, istat, iflag, nfev, ierr

c ******************************************************************
c  entries for the 'ga' namelist
c
c npopsiz    -  population size
c idum       -  if < 0, then |idum| is used as seed for random-number gen.
c pmutate    -  probability for random jump mutation
c pcross     -  crossover probability
c ielite     -  /=0  make sure best parent is preserved into decendent populations
c icreep     -  creep mutation flag:  only do creep mutations if .ne. 0
c pcreep     -  probability for random creep mutation
c iunifrm    -  =0 single point crossover at random chromosome point
c              /=0 uniform crossover
c iniche     - /=0 turn on niching
c nichflg    - array of flags for the free-parameters,
c              each non-zero entry enables niching for that free-parameter
c iskip
c iend
c nchild     - default=1; if =2, then each crossover creates 2 children, 
c              the second child having the second parents genes
c parmin     - array specifying minimum value for each free-parameter, 
c parmax     - array specifying maximum value for each free-parameter, 
c ibound     - =1 then interpret parmin and parmax as scale-factors to be
c              multiplied by the initial guess values for each parameter
c nposibl
c nowrite    - =0 then write output during optimization
c microga    - =0 perform random mutations
c             /=0 perform micro-GA
c unique_ind
c itourny
c
c ******************************************************************     
#ifdef MPI_OPT
!     Get mpi parameters
      call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)       !mpi stuff
#endif
      info = 0
      num_iter_max = num_iter_opt
c
c      temp="cp ../ga_params.in" // " ./ga_params." //
c     >            trim(filename)
c      call system(temp)
c      iunit = 72
c      call safe_open(iunit, istat, 'ga_params.'//trim(filename),
c     1               'old', 'formatted')
c      if(istat .ne. 0) 
c     1   write(6,*) "Error opening ga_params.in: istat=", istat
c      READ (iunit, NML = ga, iostat=istat) 
c      close(iunit)

      itourny=1
      maxgen=ngen
      kountmx=maxgen
      nparam=n_var
      num_obj = n_opt

      if( ibound .eq. 1 ) then
         par_max(:n_var) = x(:n_var)*parmax(:n_var)
         par_min(:n_var) = x(:n_var)*parmin(:n_var)

         where (par_max(:n_var) < par_min(:n_var) )
            partemp(:n_var) = par_max(:n_var)
            par_max(:n_var) = par_min(:n_var)
            par_min(:n_var) = partemp(:n_var)
         end where

      else
         par_max = parmax
         par_min = parmin
      endif

      if (all(nposibl .eq. 0 )) nposibl=15

      nfev=1
      f_obj => fvec

!     if (myid .eq. master) write(6, nml = ga_de)

      irestrt = 0
      if( lrestart ) irestrt = 1
c      if(irestrt .ne. 0) then
c        temp="cp ../ga_restart." // trim(filename) // " ."
c        call system(temp)
c      endif
      
c
c      store initial parameter values as a unique individual
       parent = 0
       child = 0
       iparent = 0
       ichild = 0
       if(unique_ind .gt. 0 ) then
         unique_ind=min(unique_ind, npopsiz)
         parent(1:nparam,unique_ind) = x(1:nparam)
       endif


       call ga_sp(fcn, n_opt, fvec, chi_sq, filename, nfev, iflag, 
     >            max_processors, myid)

!       iflag=1
!       chi_sq = ga_evaluate(fcn, n_opt, fvec, nparam, parent(1,jbest),
!     >           iflag, nfev)

       if (myid .eq. master) then
          write(6,*) "final solution: "
          write(6,*) "best individual : ", jbest
          write(6,*) "x ", (parent(i,jbest),i=1,nparam)
c         write(6,*) "fvec ",(fvec(i),i=1,n_opt)
          write(6,*) "y ", chi_sq
       end if

       x(1:n_var) = parent(1:n_var,jbest)

       iflag=-100
       call fcn(n_opt, npopsiz, parent(1,jbest), fvec, iflag, nfev)
     
       end subroutine GA_driver


      subroutine ga_sp(fcn, nopt, fvec, best, filename, nfev, iflag,
     >                 max_num_processors, myid)
      use kind_spec
      use ga_mod
      use safe_open_mod
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif
      external fcn

      integer :: nopt
      real(rprec), dimension(nopt) :: fvec

      integer, parameter :: master = 0
      integer :: kount, npossum,ig2sum, istart, istore
      integer :: ncross, ipick, mate1, mate2, ierr, istat
      integer :: i, j, nfev, iflag, max_num_processors, myid
      real(rprec) :: fbar, best, evals
      character*(100) filename
      save
c
c     call input
c
c  Perform necessary initialization and read the ga.restart file.
      call ga_initial(istart,npossum,ig2sum,filename,myid)
c
c  $$$$$ Main generational processing loop. $$$$$
      kount=0
      nfit_eval=nfev
      istore=0
      iunit_ga_out = 24
      if (myid .eq. master)
     1  call safe_open(iunit_ga_out, istat, 'ga_out.'//trim(filename),
     2                'unknown', 'formatted')


      do 20 i=istart,maxgen+istart-1
         iflag=-1
         if (myid .eq. master) then
            write (6,1111) i
            write (iunit_ga_out,1111) i
c           write (iunit_ga_out,1050)
c
c  Evaluate the population, assign fitness, establish the best
c  individual, and write output information.
            write(6,*) 'pre ga_evalout', max_num_processors
            write(6,*) fbar,best,nopt,nfev,max_num_processors,iflag
         end if
         call ga_evalout(fbar, best, fcn, nopt, fvec, nfev,
     >        max_num_processors, iflag, myid)
         istore=istore+1
         geni(istore)=float(i)
         genavg(istore)=fbar
         genmax(istore)=best
         if (npopsiz.eq.1 .or. iskip.ne.0) then
            if (myid .eq. master) close(iunit_ga_out)
            call ga_restart(i,istart,kount,filename, myid)
            return
         endif
c
c  niching
         if (iniche.ne.0) call ga_niche(myid)
c
c  selection, crossover and mutation
         ncross=0
         ipick=npopsiz
         do 45 j=1,npopsiz,nchild
c
c  Perform selection.
            call ga_selectn(ipick,j,mate1,mate2)
c
c  Now perform crossover between the randomly selected pair.
            call crosovr(ncross,j,mate1,mate2)
 45      continue

         if (myid .eq. master) then
            write(6,1225) ncross
            write(iunit_ga_out,1225) ncross
         end if
c
c  Now perform random mutations.  If running micro-GA, skip mutation.
         if (microga.eq.0) call ga_mutate (myid)
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated.
         call ga_newgen(npossum,ig2sum,myid)
c
c  Implement micro-GA if enabled.
         if (microga.ne.0) call ga_micro(i,npossum,ig2sum,myid)
c
c  Write to restart file.
         call ga_restart(i,istart,kount,filename,myid)
 20   continue

c  $$$$$ End of main generational processing loop. $$$$$

      if (myid .eq. master) then
         write(iunit_ga_out,3000)
         do 100 i=1,maxgen
            evals=float(npopsiz)*geni(i)
            write(iunit_ga_out,3100) geni(i),evals,genavg(i),genmax(i)
 100     continue
         CLOSE (iunit_ga_out)
      end if
      
 1050 format(1x,'    Binary Code',16x,'Parameter Values and  Fitness')
 1111 format(//'#################  Generation',i5,'  #################')
 1225 format(/'  Number of Crossovers      =',i5)
 3000 format(2x//'Summary of Output'/
     +       2x,'Generation   Evaluations   Avg.Fitness   Best Fitness')
 3100 format(2x,3(e10.4,4x),e11.5)
 
      end subroutine ga_sp


      subroutine ga_initial(istart,npossum,ig2sum,filename,myid)
c#######################################################################
c
c  This subroutine sets up the program by generating the g0, g1 and
c  ig2 arrays, and counting the number of chromosomes required for the
c  specified input.  The subroutine also initializes the random number
c  generator, parent and iparent arrays (reads the ga.restart file).
      use kind_spec
      use ga_mod
      use safe_open_mod
      implicit none
#ifdef MPI_OPT
       include 'mpif.h'                                       !mpi stuff
#endif
      integer, parameter :: master = 0
      integer :: istart, npossum, ig2sum, myid
      integer :: i, j, k, l, itemp, istat, ierr
      character*(100) :: filename
      real(rprec) :: rand
      save
c
c
      do  i=1,nparam
         g0(i)=par_min(i)
         pardel(i)=par_max(i)-par_min(i)
         itemp=2**nposibl(i)
         g1(i)=pardel(i)/(itemp-1)
      enddo

      do  i=1,nparam
         ig2(i)=nposibl(i)
      enddo
c
c  Count the total number of chromosomes (bits) required
      nchrome=0
      npossum=0
      ig2sum=0
      do 9 i=1,nparam
         nchrome=nchrome+ig2(i)
         npossum=npossum+2**nposibl(i)
         ig2sum=ig2sum+(2**ig2(i))
 9    continue
      if (nchrome.gt.nchrmax) then
         if (myid .eq. master) then
            write(6,1800) nchrome
            write(iunit_ga_out,1800) nchrome
            close(iunit_ga_out)
         end if
         stop
      endif
c
      if (npossum.lt.ig2sum .and. microga.ne.0 
     1   .and. myid.eq.master) then
         write(6,2100)
         write(iunit_ga_out,2100)
      endif
c
c  Initialize random number generator
      call ran3(idum,rand)
c
      if(irestrt.eq.0) then
c  Initialize the random distribution of parameters in the individual
c  parents when irestrt=0.
         istart=1
         do 10 i=1,npopsiz
            do 15 j=1,nchrome
               call ran3(1,rand)
               iparent(j,i)=1
               if(rand.lt.0.5d0) iparent(j,i)=0
 15         continue
 10      continue
         if (npossum.lt.ig2sum) call ga_possibl(parent,iparent,myid)
c  insert unique individual
         if (unique_ind .gt. 0) then
            do i=1, nparam
               call ga_code(unique_ind, i, parent, iparent)
            enddo
         endif

      else
c  If irestrt.ne.0, read from restart file.
         if (myid .eq. master) then
            iunit_ga_restart = 25
            call safe_open(iunit_ga_restart, istat, 
     1               '../ga_restart.'//trim(filename),
     1               'unknown', 'formatted')
            read (iunit_ga_restart,*) istart,npopsiz
            do j=1,npopsiz
               read(iunit_ga_restart,*) k,(iparent(l,j),l=1,nchrome)
            enddo
            CLOSE (iunit_ga_restart)
         end if
#ifdef MPI_OPT
         call MPI_BCAST(istart, 1, MPI_INTEGER, master, 
     1     MPI_COMM_WORLD, ierr)
         call MPI_BCAST(npopsiz, 1, MPI_INTEGER, master, 
     1     MPI_COMM_WORLD, ierr)
      do l = 1, nchrome 
         if (myid .eq. master) fitness = iparent(l,:)
         call MPI_BCAST(fitness, indmax, MPI_REAL8, master, 
     1        MPI_COMM_WORLD, ierr)
         if (myid .ne. master) iparent(l,:) = fitness
      end do
#endif       
         
      endif
c
      if(irestrt.ne.0) call ran3(idum-istart,rand)
c
 1800 format(1x,'ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
 2000 format(1x,'ERROR: You have a parameter with a number of '/
     +       1x,'   possibilities > 2**30!  If you really desire this,'/
     +       1x,'   change the DO loop 7 statement and recompile.'//
     +       1x,'   You may also need to alter the code to work with'/
     +       1x,'   REAL numbers rather than INTEGER numbers; Fortran'/
     +       1x,'   does not like to compute 2**j when j>30.')
 2100 format(1x,'WARNING: for some cases, a considerable performance'/
     +       1x,'   reduction has been observed when running a non-'/
     +       1x,'   optimal number of bits with the micro-GA.'/
     +       1x,'   If possible, use values for nposibl of 2**n,'/
     +       1x,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
c
      return
      end subroutine ga_initial


      subroutine ga_evalout(fbar, best, fcn, nopt, fvec, nfev, 
     >                   num_processors, iflag, myid)
c#######################################################################
c
c  This subroutine evaluates the population, assigns fitness,
c  establishes the best individual, and outputs information.
      use kind_spec
      use ga_mod
c      use optim, ONLY : ldiag_opt
      implicit none
      external fcn
#ifdef MPI_OPT
#else
      external ga_fitness_parallel
#endif      
      integer, parameter :: master=0
      integer :: nopt, n, j, k, kk, icross, ncreep, iflag, myid
      real(rprec), dimension(nopt) :: fvec
      real(rprec), dimension(nparmax) :: paramsm,paramav
      integer :: nfev, num_processors 
      real(rprec) :: fitsum, funcval, fbar, best, rand
      integer :: jstart, jend, istat, jstat
      logical :: ldiag_opt
      
      save
c
c

      fitsum = 0 
      best=-1.0e30_dp

      ldiag_opt = .false.

      if (myid .eq. master) write(6,*) 'in ga_evalout',num_processors
         write(6,*) fbar,best,nopt,nfev,iflag
c  ,iflag
c                   fbar,best,nopt,nfev,num_processors,iflag

      do 29 n=1,nparam
         paramsm(n)=0
 29   continue
      jstart=1
      jend=npopsiz
      if(iskip.ne.0) jstart=iskip
      if(iend.ne.0) jend=iend
c
      do  j=jstart,jend

         call ga_decode(j,parent,iparent)

c        if(iskip.ne.0 .and. iend.ne.0 .and. iskip.eq.iend) then
c        if(lscreen) then
c        if(nchrome .le. 120) then
c        write(6,1075) j,(iparent(k,j),k=1,nchrome)
c        else
c        write(6,1075) j,(iparent(k,j),k=1,120)
c        write(6,1077) (iparent(k,j),k=121,nchrome)
c        endif
c        write(6,1076)   (parent(kk,j),kk=1,nparam),0.0
c        endif
c        endif
         if(LDIAG_OPT .and. myid.eq.master) then
            if(nchrome .le. 120) then
               write(iunit_ga_out,1075) j,(iparent(k,j),k=1,nchrome)
            else
               write(iunit_ga_out,1075) j,(iparent(k,j),k=1,120)
               write(iunit_ga_out,1077) (iparent(k,j),k=121,nchrome)
            endif
            write(iunit_ga_out,1076) (parent(kk,j),kk=1,nparam)
         endif

      enddo
#if defined(MPI_OPT)
         call ga_fitness_mpi (jend-jstart+1, f_obj, num_obj, 
     1        fcn, nfev, fitness)
#else
         if (myid .eq. master) write(6,'(1x,i4,a,i4,a)') jend-jstart+1,
     1         ' processes started on ',num_processors, ' processors'
c
c        flush out buffer before multiprocessing
#if defined(RISC)
         call flush_(6)
         call flush_(iunit_ga_out)
#else
         call flush(6)
         call flush(iunit_ga_out)
#endif 

         call multiprocess(jend-jstart+1, num_processors, 
     >                  ga_fitness_parallel, fcn )
#endif 
        nfev=nfev+jend-jstart+1
        nfit_eval=nfev

!       Clean up...
        iflag=-100
        call fcn(nopt, npopsiz, parent(1,jbest), fvec, iflag, nfev)
#ifndef MPI_OPT
        do j=jstart, jend
c
c  Call function evaluator, write out individual and fitness, and add
c  to the summation for later averaging.
c        iflag=j
c        funcval = ga_evaluate(fcn, nopt, fvec, nparam, parent(1,j),
c    1                      iflag, nfev)

         read(j+1000, iostat=istat) jstat, iflag
         if( jstat .ne. j ) then
         write(6,*) "wrong index read in evalout"
         iflag=-14
         exit
         endif

         read(j+1000, iostat=istat) funcval
         fitness(j)=funcval
         close(j+1000, status='delete')
        end do
#endif        
        do 30 j = jstart, jend
         fitsum=fitsum+fitness(j)
         do 22 n=1,nparam
            paramsm(n)=paramsm(n)+parent(n,j)
 22      continue
c
c  Check to see if fitness of individual j is the best fitness.
         if (fitness(j).gt.best) then
            best=fitness(j)
            jbest=j
            do 24 k=1,nchrome
               ibest(k)=iparent(k,j)
 24         continue
         endif
 30   continue
     
c  Compute parameter and fitness averages.
      fbar=fitsum/npopsiz
      do 23 n=1,nparam
         paramav(n)=paramsm(n)/npopsiz
 23   continue

c
c  Write output information
      if (myid.eq.master) then
      if (ldiag_opt) then
         if (npopsiz.eq.1) then
            if(nchrome .le. 120) then
               write(iunit_ga_out,1075) 1,(iparent(k,1),k=1,nchrome)
            else
            write(iunit_ga_out,1075) 1,(iparent(k,1),k=1,120)
            write(iunit_ga_out,1077) (iparent(k,j),k=121,nchrome)
            end if
            write(iunit_ga_out,1076)   (parent(k,1),k=1,nparam)
            write(iunit_ga_out,1078)  fitness(1)
            write(iunit_ga_out,*) ' Average Values:'
            write(iunit_ga_out,1275) (parent(k,1),k=1,nparam)
            write(iunit_ga_out,1276) fbar
         else
            write(iunit_ga_out,1275) (paramav(k),k=1,nparam)
            write(iunit_ga_out,1276) (fitness(j),j=1,npopsiz)
         end if
      end if
      write(6,1100) fbar
      write(iunit_ga_out,1100) fbar
      write(6,1200) best
      write(iunit_ga_out,1200) best
      end if
      
 1075 format(i3,1x,(120i1))
 1077 format(3x,1x,(120i1))
 1076 format(3x,1x,10(1x,e10.4))
 1078 format(10x,e12.5)
 1100 format(1x,'Average Function Value of Generation=',e12.5)
 1200 format(1x,'Maximum Function Value              =',e12.5/)
 1275 format(/' Average Values:',18x,10(1x,e10.4))
 1276 format(10x,10e12.5)

      end subroutine ga_evalout


      subroutine ga_niche(myid)
c#######################################################################
c
c  Implement "niching" through Goldberg''s multidimensional phenotypic
c  sharing scheme with a triangular sharing function.  To find the
c  multidimensional distance from the best individual, normalize all
c  parameter differences.
c
      use kind_spec
      use ga_mod
      implicit none
      integer, parameter :: master=0
      real(rprec) :: alpha, del, del2, sigshar, sumshar, share
      integer :: nniche, jj, ii, j, k, myid
      save
c
c   Variable definitions:
c
c  alpha   = power law exponent for sharing function; typically = 1.0
c  del     = normalized multidimensional distance between ii and all
c            other members of the population
c            (equals the square root of del2)
c  del2    = sum of the squares of the normalized multidimensional
c            distance between member ii and all other members of
c            the population
c  nniche  = number of niched parameters
c  sigshar = normalized distance to be compared with del; in some sense,
c            1/sigshar can be viewed as the number of regions over which
c            the sharing function should focus, e.g. with sigshar=0.1,
c            the sharing function will try to clump in ten distinct
c            regions of the phase space.  A value of sigshar on the
c            order of 0.1 seems to work best.
c  share   = sharing function between individual ii and j
c  sumshar = sum of the sharing functions for individual ii
c
c      alpha=1.0
      sigshar=0.1_dp
      nniche=0
      do 33 jj=1,nparam
         nniche=nniche+nichflg(jj)
 33   continue
      if (nniche.eq.0) then
         if (myid .eq. master) then
            write(6,1900)
            write(iunit_ga_out,1900)
            close(iunit_ga_out)
         end if
         stop
      endif
      do 34 ii=1,npopsiz
         sumshar=0
         do 35 j=1,npopsiz
            del2=0
            do 36 k=1,nparam
               if (nichflg(k).ne.0) then
                  del2=del2+((parent(k,j)-parent(k,ii))/pardel(k))**2
               endif
 36         continue
            del=sqrt(del2)/nniche
            if (del.lt.sigshar) then
c               share=1.0-((del/sigshar)**alpha)
               share=1-(del/sigshar)
            else
               share=0
            endif
            sumshar=sumshar+share/npopsiz
 35      continue
         if (sumshar.ne.0.0_dp) fitness(ii)=fitness(ii)/sumshar
 34   continue
 
 1900 format(1x,'ERROR: iniche=1 and all values in nichflg array = 0'/
     +       1x,'       Do you want to niche or not?')
 
      end subroutine ga_niche


      subroutine ga_selectn(ipick,j,mate1,mate2)
c#######################################################################
c
c  Subroutine for selection operator.  Presently, tournament selection
c  is the only option available.
c
      use ga_mod
      implicit none
      
      integer :: ipick, j, mate1, mate2
      integer :: n
      save
c
      if(itourny.eq.1) then
         call ga_select(mate1,ipick)
         call ga_select(mate2,ipick)
c        write(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
         do 46 n=1,nchrome
            ichild(n,j)=iparent(n,mate1)
            if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate2)
 46      continue
      endif
c
      end subroutine ga_selectn


      subroutine crosovr(ncross,j,mate1,mate2)
c#######################################################################
c
c  Subroutine for crossover between the randomly selected pair.
      use kind_spec
      use ga_mod
      implicit none 
      integer :: ncross, j, mate1, mate2
      integer :: n, icross
      real(rprec) :: rand
      save
c
      if (iunifrm.eq.0) then
c  Single-point crossover at a random chromosome point.
         call ran3(1,rand)
         if(rand.gt.pcross) goto 69
         ncross=ncross+1
         call ran3(1,rand)
         icross=2+int((nchrome-1)*rand)
         do 50 n=icross,nchrome
            ichild(n,j)=iparent(n,mate2)
            if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
 50      continue
      else
c  Perform uniform crossover between the randomly selected pair.
         do 60 n=1,nchrome
            call ran3(1,rand)
            if(rand.le.pcross) then
               ncross=ncross+1
               ichild(n,j)=iparent(n,mate2)
               if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
            endif
 60      continue
      endif
 69   continue
 
      end subroutine crosovr

 
      subroutine ga_mutate (myid)
c#######################################################################
c
      use kind_spec
      use ga_mod
      implicit none 
      integer, parameter :: master = 0
      integer :: nmutate, ncreep, i, j, k, myid
      real(rprec) :: rand, creep
      save
c
c  This subroutine performs mutations on the children generation.
c  Perform random jump mutation if a random number is less than pmutate.
c  Perform random creep mutation if a different random number is less
c  than pcreep.
      nmutate=0
      ncreep=0
      do 70 j=1,npopsiz
         do 75 k=1,nchrome
c  Jump mutation
            call ran3(1,rand)
            if (rand.le.pmutate) then
               nmutate=nmutate+1
               if(ichild(k,j).eq.0) then
                  ichild(k,j)=1
               else
                  ichild(k,j)=0
               endif
               if (nowrite.eq.0 .and. myid.eq.master) then
                  write(6,1300) j,k
                  write(iunit_ga_out,1300) j,k
               end if
            endif
 75      continue
c  Creep mutation (one discrete position away).
         if (icreep.ne.0) then
            do 76 k=1,nparam
               call ran3(1,rand)
               if(rand.le.pcreep) then
                  call ga_decode(j,child,ichild)
                  ncreep=ncreep+1
                  creep=1
                  call ran3(1,rand)
                  if (rand.lt.0.5_dp) creep=-1
                  child(k,j)=child(k,j)+g1(k)*creep
                  if (child(k,j).gt.par_max(k)) then
                     child(k,j)=par_max(k)-1.0d0*g1(k)
                  elseif (child(k,j).lt.par_min(k)) then
                     child(k,j)=par_min(k)+g1(k)
                  endif
                  call ga_code(j,k,child,ichild)
                  if (nowrite.eq.0 .and. myid.eq.master) then
                     write(6,1350) j,k
                     write(iunit_ga_out,1350) j,k
                  end if   
               endif
 76         continue
         endif
 70   continue
      if (myid .eq. master) then
         write(6,1250) nmutate,ncreep
         write(iunit_ga_out,1250) nmutate,ncreep
      end if
 
 1250 format(/'  Number of Jump Mutations  =',i5/
     +        '  Number of Creep Mutations =',i5)
 1300 format('*** Jump mutation performed on individual  ',i4,
     +       ', chromosome ',i3,' ***')
 1350 format('*** Creep mutation performed on individual ',i4,
     +       ', parameter  ',i3,' ***')
c
      end subroutine ga_mutate


      subroutine ga_newgen(npossum,ig2sum,myid)
c#######################################################################
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated; if not, and if ielite=1,
c  then reproduce the best parent into a random slot.
c
      use kind_spec
      use ga_mod
      implicit none 
      integer, parameter :: master = 0
      integer :: npossum, ig2sum, kelite, jelite
      integer :: irand, i, j, n, myid
      real(rprec) :: rand
      save
c
      if (npossum.lt.ig2sum) call ga_possibl(child,ichild,myid)
      kelite=0
      do 94 j=1,npopsiz
         jelite=0
         do 95 n=1,nchrome
            iparent(n,j)=ichild(n,j)
            if (iparent(n,j).eq.ibest(n)) jelite=jelite+1
            if (jelite.eq.nchrome) kelite=1
 95      continue
 94   continue
      if (ielite.ne.0 .and. kelite.eq.0) then
         call ran3(1,rand)
         irand=1+int(npopsiz*rand)
         iparent(1:nchrome,irand)=ibest(1:nchrome)
         if (myid .eq. master) write(iunit_ga_out,1260) irand
      endif
c
 1260 format('  Elitist Reproduction on Individual ',i4)
c
      end subroutine ga_newgen


      subroutine ga_micro(i,npossum,ig2sum, myid)
c#######################################################################
c
c  Micro-GA implementation subroutine
c
      use kind_spec
      use ga_mod
      implicit none 
      integer, parameter :: master=0
      integer :: i, npossum, ig2sum, myid
      integer :: icount, j, n
      real(rprec) :: diffrac, rand
      save
c
c
c  First, check for convergence of micro population.
c  If converged, start a new generation with best individual and fill
c  the remainder of the population with new randomly generated parents.
c
c  Count number of different bits from best member in micro-population
      icount=0
      do 81 j=1,npopsiz
         do 82 n=1,nchrome
            if(iparent(n,j).ne.ibest(n)) icount=icount+1
 82      continue
 81   continue
c
c  If icount less than 5% of number of bits, then consider population
c  to be converged.  Restart with best individual and random others.
      diffrac=real(icount,rprec)/((npopsiz-1)*nchrome)
      if (diffrac.lt.0.05_dp) then
         do 87 n=1,nchrome
            iparent(n,1)=ibest(n)
 87      continue
         do 88 j=2,npopsiz
            do 89 n=1,nchrome
               call ran3(1,rand)
               iparent(n,j)=1
               if(rand.lt.0.5_dp) iparent(n,j)=0
 89         continue
 88      continue
         if (npossum.lt.ig2sum) call ga_possibl(parent,iparent,myid)
         if (myid .eq. master) then
           write(6,1375) i
            write(iunit_ga_out,1375) i
         end if   
      endif
 
 1375 format(//'%%%%%%%  Restart micro-population at generation',
     +       i5,'  %%%%%%%')
c
      return
      end subroutine ga_micro


      subroutine ga_select(mate,ipick)
c#######################################################################
c
c  This routine selects the better of two possible parents for mating.
c
      use kind_spec
      use ga_mod
      implicit none 
      integer :: mate, ipick, ifirst, isecond
      save
c
      if(ipick+1.gt.npopsiz) call ga_shuffle(ipick)
      ifirst=ipick
      isecond=ipick+1
      ipick=ipick+2
      if(fitness(ifirst).gt.fitness(isecond)) then
         mate=ifirst
      else
         mate=isecond
      endif
c     write(3,*)'select',ifirst,isecond,fitness(ifirst),fitness(isecond)
c
      end subroutine ga_select


      subroutine ga_shuffle(ipick)
c#######################################################################
c
c  This routine shuffles the parent array and its corresponding fitness
c
      use kind_spec
      use ga_mod
      implicit none 
      integer :: ipick, j, n, itemp, iother
      real(rprec) :: rand, temp
      save
c
      ipick=1
      do 10 j=1,npopsiz-1
         call ran3(1,rand)
         iother=j+1+int((npopsiz-j)*rand)
         do 20 n=1,nchrome
            itemp=iparent(n,iother)
            iparent(n,iother)=iparent(n,j)
            iparent(n,j)=itemp
 20      continue
         temp=fitness(iother)
         fitness(iother)=fitness(j)
         fitness(j)=temp
 10   continue
c
      end subroutine ga_shuffle


      subroutine ga_decode(i,array,iarray)
c#######################################################################
c
c  This routine decodes a binary string to a real number.
c
      use kind_spec
      use ga_mod
      implicit none 
      integer :: i, j, m, l, k, iparam
      real(rprec), dimension(nparmax,indmax) :: array
      integer, dimension(nchrmax,indmax) :: iarray
      save
 
      l=1
      do 10 k=1,nparam
         iparam=0
         m=l
         do 20 j=m,m+ig2(k)-1
            l=l+1
            iparam=iparam+iarray(j,i)*(2**(m+ig2(k)-1-j))
 20      continue
         array(k,i)=g0(k)+g1(k)*iparam
 10   continue

      end subroutine ga_decode


      subroutine ga_code(j,k,array,iarray)
c#######################################################################
c
c
c  This routine codes a parameter into a binary string.
c
      use kind_spec
      use ga_mod
      implicit none 
      integer :: j, k, i, istart, iparam, m
      real(rprec), dimension(nparmax,indmax) :: array
      integer, dimension(nchrmax,indmax) :: iarray

      save
c
c  First, establish the beginning location of the parameter string of
c  interest.
      istart=1
      do 10 i=1,k-1
         istart=istart+ig2(i)
 10   continue
c
c  Find the equivalent coded parameter value, and back out the binary
c  string by factors of two.
      m=ig2(k)-1
      if (g1(k).eq.0.0_dp) return
      iparam=nint((array(k,j)-g0(k))/g1(k))
      do 20 i=istart,istart+ig2(k)-1
         iarray(i,j)=0
         if ((iparam+1).gt.(2**m)) then
            iarray(i,j)=1
            iparam=iparam-2**m
         endif
         m=m-1
 20   continue
c     write(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
c
      end subroutine ga_code


      subroutine ga_possibl(array, iarray, myid)
c#######################################################################
c
c  This subroutine determines whether or not all parameters are within
c  the specified range of possibility.  If not, the parameter is
c  randomly reassigned within the range.  This subroutine is only
c  necessary when the number of possibilities per parameter is not
c  optimized to be 2**n, i.e. if npossum < ig2sum.
c
      use kind_spec
      use ga_mod
      implicit none

      integer, parameter :: master = 0
      real(rprec), dimension(nparmax,indmax) :: array
      integer, dimension(nchrmax,indmax) :: iarray
      integer :: i, j, n2ig2j, irand, myid
      real(rprec) :: rand

      save
c
      do 10 i=1,npopsiz
         call ga_decode(i,array,iarray)
         do 20 j=1,nparam
            n2ig2j=ig2(j)
            if(nposibl(j).ne.n2ig2j .and. array(j,i).gt.par_max(j)) then
               call ran3(1,rand)
               irand=int((2**nposibl(j))*rand)
               array(j,i)=g0(j)+irand*g1(j)
               call ga_code(i,j,array,iarray)
               if (nowrite.eq.0 .and. myid.eq.master) then
                  write(6,1000) i,j
                  write(iunit_ga_out,1000) i,j
               end if
            endif
 20      continue
 10   continue

 1000 format('*** Parameter adjustment to individual     ',i4,
     1       ', parameter  ',i3,' ***')
 
      end subroutine ga_possibl


      subroutine ga_restart(i,istart,kount,filename, myid)
c#######################################################################
c
c  This subroutine writes restart information to the ga.restart file.
c
      use kind_spec
      use ga_mod
      use safe_open_mod
      implicit none 
      integer, parameter :: master = 0
      integer :: kount, i, j, l, istart, istat, myid
       character*(100) :: filename
      save

      if (myid .ne. master) return
      
      kount=kount+1
      if(i.eq.maxgen+istart-1 .or. kount.eq.kountmx) then
         iunit_ga_restart = 25
         call safe_open(iunit_ga_restart, istat, 
     1                'ga_restart.'//trim(filename),
     2               'unknown', 'formatted')
         rewind iunit_ga_restart
         write(iunit_ga_restart,*) i+1,npopsiz
         do 80 j=1,npopsiz
            write(iunit_ga_restart,*) j,(iparent(l,j),l=1,nchrome)
c        if(nchrome .le. 60) then
c           write(iunit_ga_restart,1500) j,(iparent(l,j),l=1,nchrome)
c        else
c           write(iunit_ga_restart,1500) j,(iparent(l,j),l=1,60)
c           write(iunit_ga_restart,1501) (iparent(l,j),l=61,nchrome)
c        endif
 80      continue
         CLOSE (iunit_ga_restart)
         kount=0
      endif
c
 1500 format(i5,3x,60i2)
 1501 format(5x,3x,60i2)

      end subroutine ga_restart


      subroutine ran3(idum,rand)
c#######################################################################
c
c  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
c  any negative value to initialize or reinitialize the sequence.
c  This function is taken from W.H. Press, "Numerical Recipes" p. 199.
c
      use kind_spec
      implicit none 
      integer :: idum, iff, inext, inextp, i, ii, k
      real(rprec) :: rand, ma, mj, mk
      save
      real(rprec), parameter :: mbig=4000000, mseed=1618033,
     1           mz=0, fac=1._dp/mbig
c     parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
c
c  According to Knuth, any large mbig, and any smaller (but still large)
c  mseed can be substituted for the above values.
      dimension ma(55)
      data iff /0/
      if (idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=mseed-abs(idum)
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      rand=mj*fac

      end subroutine ran3

      function ga_evaluate(fcn, nopt, fvec, n, x, iflag, nfev) 
c#######################################################################
        use kind_spec
        integer :: n, nfev
        real(rprec) :: ga_evaluate
        real(rprec), dimension(n) :: x
        integer :: nopt
        real(rprec), dimension(nopt) :: fvec
        external  fcn

      integer   i, j, iflag 

      call fcn(nopt, n, x, fvec, iflag, nfev)
      ga_evaluate = -sum(fvec(:nopt)**2) 

      end function ga_evaluate


#if defined(MPI_OPT)
      subroutine ga_fitness_mpi (np, fvec, nopt, fcn, nfev, funcval)
      use ga_mod
      implicit none
      include 'mpif.h'                                       !mpi stuff

      integer :: nopt
      real(rprec), dimension(nparam) :: x
      real(rprec), dimension(nopt) :: fvec

      integer, parameter :: master = 0                       !mpi stuff
      integer :: np, nfev, j, iflag, istat
      external fcn
      real(rprec) :: funcval(np)

      integer :: status(MPI_STATUS_SIZE)                     !mpi stuff
      integer :: myid, numprocs, i                           !mpi stuff
      integer :: numsent, sender, ierr                       !mpi stuff
      integer :: anstype, column                             !mpi stuff

!******************************************
!
!  mpi setup calls; set barrier so all processors get here before starting
!
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )       !mpi stuff
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )   !mpi stuff
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                 !mpi stuff

!******************************************
!
!     ****Master portion of the code****
!
      if (myid .eq. master) then
         numsent = 0    !numsent is a counter used to track how many
                        !jobs have been sent to workers
c
c     Send forward difference displacements from master to each
c           worker process. Tag with these with the column number.
c
         do j = 1,min(numprocs-1,np)
            x(:) = parent(:,j)
            call MPI_SEND(x, nparam, MPI_REAL8, j, 
     1                  j, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(1) in ga_fitness_mpi'
            numsent = numsent+1
         end do          !j = 1,min(numprocs-1,n)
c
c      Looping through the columns, collect answers from the workers.
c      As answers are received, new uncalculated columns are sent
c      out to these same workers.
c
         do j = 1,np
            call MPI_RECV(fvec, nopt, MPI_REAL8, 
     1           MPI_ANY_SOURCE, MPI_ANY_TAG, 
     2           MPI_COMM_WORLD, status, ierr)
            if (ierr .ne. 0) stop 'MPI_RECV error(1) in ga_fitness_mpi'
            sender     = status(MPI_SOURCE)    
            anstype    = status(MPI_TAG)       ! column is tag value
            if (anstype .gt. np) stop 'ANSTYPE > NP IN ga_fitness_mpi'
            
            funcval(anstype) = -sum(fvec(:nopt)**2)
            write(6,'(a,1pe10.3,a,i3,a,i3)')' FUNCVAL = ', 
     1       -funcval(anstype),
     2      ' for iteration ', anstype+nfev,' processor = ', sender

c
c           If more columns are left, then send another column to the worker(sender)
c           that just sent in an answer
c
            if (numsent .lt. np) then
               numsent = numsent+1
               x(:) = parent(:,numsent)
               
               call MPI_SEND(x, nparam, MPI_REAL8, 
     1                       sender, numsent, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0) 
     1            stop 'MPI_SEND error(2) in ga_fitness_mpi'

            else                ! Tell worker that there is no more work to do
               
               call MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0)stop 'MPI_END error(3) in ga_fitness_mpi'
            endif      ! if( myid .eq. master ) then
         end do     ! do j = 1,n
c
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
c
      else if (myid .le. np) then        ! i.e., if( myid .ne. master )
c
c        Otherwise accept the next available column, check the tag,
c        and if the tag is non-zero call subroutine fcn.
c        If the tag is zero, there are no more columns
c        and worker skips to the end.
c
 90      call MPI_RECV(x, nparam, MPI_REAL8, master, 
     1                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         if (ierr .ne. 0) stop 'MPI_RECV error(2) in ga_fitness_mpi'
        
         column = status(MPI_TAG)                !!ID of pseudo-processor issuing this message
         if (column .eq. 0) then
            go to 200
         else
            iflag = column
c           Call the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA stores the local fvec_min array

            call fcn(nopt, nparam, x, fvec, iflag, nfev)
            if (iflag.ne.0) go to 300
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            call MPI_SEND(fvec, nopt, MPI_REAL8, master, 
     1                    column, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(4) in ga_fitness_mpi'
            go to 90    !Return to 90 and check if master process has sent any more jobs
         endif
 200     continue
      endif       ! if( myid .ne. master )

!
!     Broadcast the funcval array to all processors FROM master
!
      call MPI_BCAST(funcval, np, MPI_REAL8, master,
     1     MPI_COMM_WORLD, ierr)
      if (ierr .ne. 0) go to 100 

      return

 100  continue
      print *,' MPI_BCAST error in ga_fitness_mpi: IERR=', ierr
      
      return

 300  continue
      print *,' IFLAG = ', iflag, ' in ga_fitness_mpi call to fcn'
      stop

      end subroutine ga_fitness_mpi
#else
      subroutine ga_fitness_parallel(j,fcn)
c#######################################################################
      use ga_mod
      implicit none

      integer :: j, iflag, istat
      external fcn
      real(rprec) :: funcval, ga_evaluate

      iflag=j
      funcval = ga_evaluate(fcn, num_obj, f_obj, nparam, parent(1,j),
     >                   iflag, nfit_eval)
      write (j+1000) j, iflag
      write (j+1000) funcval
      close (j+1000)

      end subroutine ga_fitness_parallel

#endif
      subroutine write_gade_nml(iunit)
      use kind_spec
      use gade_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: iunit
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i
      external write_array
C-----------------------------------------------
      write (iunit,'(a)')'&GA_DE'
      write (iunit,200) 'NPOPSIZ = ', npopsiz
      write (iunit,200) 'NGEN    = ', ngen
      write (iunit,200) 'IDUM    = ', idum
      write (iunit,200) 'IBOUND  = ', ibound
      write (iunit,200) 'NOWRITE = ', nowrite
      write (iunit,200) 'MICROGA = ', microga
      write (iunit,200) 'ISKIP   = ', iskip
      write (iunit,200) 'IEND    = ', iend
      write (iunit,200) 'NCHILD  = ', nchild
      write (iunit,200) 'ITOURNY = ', itourny
      write (iunit,200) 'IELITE  = ', ielite
      write (iunit,200) 'ICREEP  = ', icreep
      write (iunit,200) 'IUNIFRM = ', iunifrm
      write (iunit,200) 'INICHE  = ', iniche
      write (iunit,200) 'STRATEGY= ', strategy
      write (iunit,200) 'CR_STRATEGY = ', cr_strategy
      write (iunit,200) 'OUT_ITER= ', out_iter
      write (iunit,200) 'UNIQUE_IND = ', unique_ind
      write (iunit,210) 'PCROSS  = ', pcross
      write (iunit,210) 'F_CROSS = ', f_cross
      write (iunit,210) 'PMUTATE = ', pmutate
      write (iunit,210) 'PCREEP  = ', pcreep
      write (iunit,101) 'SAVE_SPACE = ', save_space
      write (iunit,100) 'NICHFLG = '
      write (iunit,110) (nichflg(i), i=1,nparmax)
      write (iunit,100) 'NPOSIBL = '
      write (iunit,110) (nposibl(i), i=1,nparmax)
      call write_array (iunit, 'PARMIN', parmin, nparmax)
      call write_array (iunit, 'PARMAX', parmax, nparmax)
      write (iunit,'(a)') '/'

 100  format (2x,a)
 101  format (2x,a,l2)
 110  format (2x, 16i5)
 200  format (2x,a,i5)
 210  format (2x,a,1pe20.14)

      end subroutine write_gade_nml
EOF
EOC
