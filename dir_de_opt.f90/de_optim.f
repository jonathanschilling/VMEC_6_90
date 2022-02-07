      subroutine DE_preset
      use kind_spec
      use de_mod
      implicit none

      npopsiz = 0
      ngen = 0
      idum = 0
      strategy = 2
      CR_strategy = 0
      f_cross = 0.5_dp
      pcross = 0.3_dp
      parmin = 0.5_dp
      parmax = 2
      ibound = 1
      out_iter = 1
      save_space = .false.

      n_pop = 0
      n_free = 0
      nopt = 0
      if( allocated(ui_xc) ) deallocate( ui_xc)

      end subroutine DE_preset
 

      subroutine DE_driver(fcn, n_opt, n_var, x, fvec, tol, eps,  
     1     num_iter_opt, max_processors, filename, info, lwa, lrestart )
c ******************************************************************
      use kind_spec
      use de_mod
      use gade_mod, ONLY: ga_de
      use system_mod
c     use read_namelist_mod
      use safe_open_mod
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif
      integer, parameter :: iflag_cleanup = -100, master = 0
      integer :: n_opt, n_var, info, lwa, num_iter_opt, max_processors
      integer :: myid = master 
      real(rprec), dimension(n_opt), target :: fvec
      real(rprec), dimension(n_var) :: x
      real(rprec) :: tol, eps, chi_sq, ga_evaluate, tmp
      external fcn
      character*(*) ::  filename
      logical :: lrestart
c
c     local variables
c

      character*(100) ::  temp
      integer :: iunit, ierr
      integer :: i, j, istat, iflag, irestart = 25
      real(rprec), dimension(n_var) :: par_max, par_min, partemp
      real(rprec) ::  tgt

c ******************************************************************
c  entries for the 'de' namelist
c
c npopsiz    -  population size
c ngen       -  number of generations
c idum       -  if < 0, then |idum| is used as seed for random-number gen.
c strategy   - The strategy of the mutation operations is used in HDE,
c              except that cross-over strategy is separated into separate
c              specifier CR_strategy below
c CR_strategy - cross-over strategy 0=exponential, 1=binomial
c f_cross    - crossover scaling factor
c pcross     -  crossover probability (CR)
c parmin     - array specifying minimum value for each free-parameter, 
c parmax     - array specifying maximum value for each free-parameter, 
c ibound     - =1 then interpret parmin and parmax as scale-factors to be
c              multiplied by the initial guess values for each parameter
c              =0 interpret parmin and parmax as absolute values
c out_iter   - The intermediate output will be produced after "out_iter"
c              iterations. No intermediate output will be produced if
c              "out_iter < 1".
c
c ******************************************************************     
#ifdef MPI_OPT
!     Get mpi parameters
      call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)       !mpi stuff
#endif
      info = 0

      n_free = n_var
      nopt = n_opt
      n_pop = npopsiz
      if( n_pop <= 0) n_pop = 10*n_free
      if( (pcross < 0._dp) .or. (pcross > 1.0_dp) ) then
         if (myid .eq. master)
     1       write(6,*) 'Illegal value of pcross for de:', pcross
         stop
      endif
!     out_iter = max(1, out_iter)


      if( ibound .eq. 1 ) then
         par_max(:n_var) = x(:n_var)*parmax(:n_var)
         par_min(:n_var) = x(:n_var)*parmin(:n_var)

         where (par_max(:n_var) < par_min(:n_var) )
            partemp(:n_var) = par_max(:n_var)
            par_max(:n_var) = par_min(:n_var)
            par_min(:n_var) = partemp(:n_var)
         end where

      else
         par_max(:n_var) = parmax(:n_var)
         par_min(:n_var) = parmin(:n_var)
      endif

      nfev=1

      if (myid .eq. master) then
!         write(6, nml = ga_de)

! ------setup restart file -------------

         if( lrestart ) 
     1      call system("cp ../de_restart."//trim(filename)//" .")

         call safe_open(irestart, istat, 'de_restart.'//trim(filename),
     1               'unknown', 'formatted')
    
      end if
      
      iflag=1
      tgt = 0
      call DE_Evolve(fcn,n_free,par_min,par_max, tgt, n_pop, ngen,
     1        f_cross,pcross,strategy,CR_strategy,out_iter,iunit,
     2        irestart, max_processors,lrestart, x, chi_sq,nfev)
 
      if (myid .eq. master) then
         write(6,*) "final solution: "
c        write(6,*) "x ", x(:n_var)
c        write(6,*) "fvec ",(fvec(i),i=1,n_opt)
         write(6,*) "x " 
         write(6,990) (x(i), i=1,n_var)
990   format(5(1pe22.14))
         write(6,*) "y ", chi_sq
      end if
      
      iflag = iflag_cleanup
      call fcn(n_opt, npopsiz, x, fvec, iflag, nfev)
     
      end subroutine DE_driver


      subroutine DE_Evolve(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, 
     1     F_XC, 
     2     CR_XC, strategy, CR_strategy, refresh, iwrite, irestart,
     3     num_proc, lrestart,
     4     best_XC, bestval, nfeval)
!.......................................................................
!    
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB 
!  version of differential evolution (DE). This FORTRAN 90 code 
!  has been tested on Compaq Visual Fortran v6.1. 
!  Any users new to the DE are encouraged to read the article of Storn and Price. 
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the 
!    ICEC''96 contest by differential evolution. IEEE conf. on Evolutionary 
!    Computation, 842-844.
!
!  This Fortran 90 program is written by Dr. Feng-Sheng Wang 
!  Department of Chemical Engineering, National Chung Cheng University, 
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!
!  Modified by M. Zarnstorff, PPPL, to implement all strategies and make random
!  selection process be compatible with c versions and restart capability.
!  Also, added multi-processor interface to STELLOPT suite.
!  Modified by S. Hirshman, ORNL to implement MPI interface to STELLOPT.
!.........................................................................
!                obj : The user provided file for evlauting the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!        CR_strategy : cross-over strategy 0=exponential, 1=binomial
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
!           irestart : unit number for writing (and reading) the restart file
!           num_proc : number of processors to use in evaluating population
!           lrestart : =T then read in restart file, =F ignore any restart files

!    best_XC(Dim_XC) : The best real decision parameters.  If non-zero, this
!                      array also contains an initial guess that is included
!                      in the population
!            bestval : The best objective function.
!             nfeval : The number of function call.
!

      use kind_spec
      use de_mod, ONLY: ui_XC, nopt 
      use gade_mod, ONLY : save_space
      implicit none
#ifdef MPI_OPT
      include 'mpif.h'                                       !mpi stuff
#endif

      integer, intent(in) :: NP, Dim_XC, itermax, strategy,  
     1                       CR_strategy, iwrite, irestart, refresh
      real(rprec), intent(in) :: VTR, CR_XC, F_XC
      real(rprec), dimension(Dim_XC), intent(in) :: XCmin, XCmax
      real(rprec), dimension(Dim_XC), intent(inout) :: best_XC
      real(rprec), dimension(Dim_XC) :: temp_var
      real(rprec), intent(out) :: bestval
      integer, intent(out) :: nfeval, num_proc  
      logical :: lrestart
      external  obj

      logical, dimension(NP,Dim_XC) :: lcross
      integer, parameter :: master = 0
      integer :: myid = master, ierr
      integer :: i, j, iter, k, n,istat, r_free, seed_size
      integer, dimension(1) :: inumber, ibest
      integer, dimension(NP) :: a1, a2, a3, a4, a5
      INTEGER, ALLOCATABLE :: seed (:) 

      real(rprec) :: tempval, s, diversity, xxval, val_mean, val_max
      real(rprec), dimension(NP) :: val, rand, new_val
      real(rprec), dimension(Dim_XC) :: bestmemit_XC
      real(rprec), dimension(Dim_XC) :: rand_C1
      real(rprec), allocatable, dimension(:,:) :: pop_XC, rand_XC
      real(rprec), dimension(nopt) :: fvec

      intrinsic max, min, random_number, mod, abs, any, all, maxloc, 
     1          minloc, maxval, minval

#ifdef MPI_OPT
!-----Get mpi parameters 
      call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)        !mpi stuff
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)                 !mpi stuff
#endif

      allocate (ui_XC(NP, Dim_XC), pop_XC(NP, Dim_XC),
     1          rand_XC(NP, Dim_XC), stat=ierr)

!-----Initialize the random-number seed stuff

      CALL RANDOM_SEED ( )  ! Processor reinitializes the seed
      if (myid .eq. master) then               
         if (ierr .ne. 0) stop 'Allocation error 1 in DE_Evolve'
         CALL RANDOM_SEED (SIZE = seed_size)  ! size of seed array
         ALLOCATE (seed(seed_size), stat=ierr)
         if (ierr .ne. 0) stop 'Allocation error 1 in DE_Evolve'
         CALL RANDOM_SEED (GET=seed(1:seed_size)) ! Gets the current seed


!!-----Initialize a population --------------------------------------------!!

         pop_XC=0
         do i=1,NP
            call random_number(rand_C1)
            pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
         end do
      
         val = 0

         if( any(best_XC /= 0._dp)) pop_XC(1,:) = best_XC   ! init with given case

!!----- Restart if requested ----

         if (lrestart) then
            write(6,*) 'Restarting from saved file'

            read (irestart,*) r_free
            if (r_free /= Dim_XC) then
               write(6,*) '***DE Restart Error***'
               write(6,*)'***Mismatch on number of free parameters! ***'
               stop
            end if
          
            read (irestart, *) seed(:seed_size)
            CALL RANDOM_SEED (PUT=seed(1:seed_size)) ! Sets seed from array

            do i=1,NP
               read (irestart,*,end=100) 
     1              k,val(i),(pop_XC(i,j),j=1, Dim_XC)
            enddo

100         REWIND (unit=irestart)
         end if
      end if

#ifdef MPI_OPT
!!all processors need same random variables and initial values
      call MPI_BCAST(seed_size, 1, MPI_INTEGER, master, 
     1     MPI_COMM_WORLD, ierr)
      if (myid .ne. master) ALLOCATE (seed(seed_size), stat=ierr)
      call MPI_BCAST(seed, seed_size, MPI_REAL8, master, 
     1     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(val, np, MPI_REAL8, master, 
     1     MPI_COMM_WORLD, ierr)
      do i = 1, NP 
         if (myid .eq. master) temp_var = pop_XC(i,:)
         call MPI_BCAST(temp_var, Dim_XC, MPI_REAL8, master, 
     1        MPI_COMM_WORLD, ierr)
         if (myid .ne. master) pop_XC(i,:) = temp_var
      end do
#endif       

!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
      nfeval=0

      ui_XC = pop_XC    ! prepare for initial evaluation

      if( any(val(:NP) == 0._dp) ) then
         call DE_Evaluate( num_proc, obj, val, NP, Dim_XC, nfeval)
      else if (myid .eq. master) then
         write(6,*) 'Using evaluation from restart file'
      endif

c     do i=1,NP
c        call obj(pop_XC(i,:), val(i))
c        nfeval=nfeval+1
c     end do      

      ibest = minloc(val)
      bestval = val(ibest(1))
      bestmemit_XC=pop_XC(ibest(1),:)
      best_XC=bestmemit_XC

      if (myid .eq. master) 
     1    write(6,*)'initial best=',bestval,' at loc ',ibest

!!--------------------------------------------------------------------------!!
!!------Perform evolutionary computation------------------------------------!! 
!!
      do iter=1, itermax

        if (myid .ne. master) goto 1000
!---- prep for crossover  -----------------------------------------!
        if( CR_XC == 1.0_dp ) then
           lcross = .true.

        else if( CR_strategy == 1) then
           call random_number(rand_XC)
           lcross = rand_XC < CR_XC   

           do i=1, NP
              if( .not. any(lcross(i,:))) then
                 n = 1 + rand_XC(i,1)*(Dim_XC-1)
                 lcross(i,n) = .true.
              endif
           enddo

        else   ! if cr_strategy = 0 or anything else
           lcross = .false.
           call random_number(rand)
           do i=1, NP
              n = 1 + rand(i)*(Dim_XC-1)
              do
                 lcross(i,n) = .true.
                 n = mod(n,Dim_XC) + 1
                 call random_number(s)
                 if(s >= CR_XC .or. lcross(i,n) ) exit
              enddo
           enddo
        endif


!!------Setup Mutation parents----------------------------------------------!!
        if( NP < 2) stop 'population too small for Diff. Evolution'

        call random_number(rand)
        a1 = 1 + rand*(NP-1)

        a2 = a1        
        do while (any(a2 == a1) )
           call random_number(rand)
           where (a2 == a1) a2 = 1 + rand*(NP-1)
        enddo

        if( strategy == 2 .or. strategy == 4 .or. strategy == 5) then
           if( NP < 3) stop 'population too small for DE strategy'
           a3 = a1        
           do while (any(a3 == a1 .or. a3 == a2))
              call random_number(rand)
              where (a3 == a1 .or. a3 == a2)
     1            a3 = 1 + rand*(NP-1)
           enddo
        endif

        if( strategy == 4 .or. strategy == 5) then
           if( NP < 4) stop 'population too small for DE strategy'
           a4 = a1        
           do while (any(a4 == a1 .or. a4 == a2 .or. a4 == a3))
              call random_number(rand)
              where (a4 == a1 .or. a4 == a2 .or. a4 == a3)
     1            a4 = 1 + rand*(NP-1)
           enddo
        endif

        if( strategy == 5) then
           if( NP < 5) stop 'population too small for DE strategy'
           a5 = a1        
           do while (any(a5 == a1 .or. a5 == a2 .or. a5 == a3 
     1             .or. a5 == a4))
              call random_number(rand)
              where (a5 == a1 .or. a5 == a2 .or. a5 == a3 .or. 
     1             a5 == a4)
     2            a5 = 1 + rand*(NP-1)
           enddo
        endif


!---- preload new descendents   -----------------------------------------!

        ui_XC = pop_XC       ! start with new pop. = old pop.


!---- for strategies using best previous case, preload best case into new pop.

        if( strategy == 1 .or. strategy == 3 .or. strategy == 4) then
           do i=1,NP
              where(lcross(i,:)) ui_XC(i,:) = bestmemit_XC(:)
           enddo
        endif

!---- select and do a  mutation strategy-----------------------------------!

!=======Choice of strategy=================================================
!  We have tried to come up with a sensible naming-convention: DE/x/y/z
!    DE :  stands for Differential Evolution
!    x  :  a string which denotes the vector to be perturbed
!    y  :  number of difference vectors taken for perturbation of x
!    z  :  crossover method (exp = exponential, bin = binomial)
!
!  There are some simple rules which are worth following:
!  1)  F is usually between 0.5 and 1 (in rare cases > 1)
!  2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first
!  3)  To start off NP = 10*D is a reasonable choice. Increase NP if 
!      misconvergence happens.
!  4)  If you increase NP, F usually has to be decreased
!  5)  When the DE/best... schemes fail DE/rand... usually works and vice versa
!
!
! Here are Storn''s comments on the different strategies:
!
! (1) DE/best/1/'z'
!     Our oldest strategy but still not bad. However, we have found several
!     optimization problems where misconvergence occurs.
!
! (2) DE/rand/1/'z'
!     This is one of my favourite strategies. It works especially well when the
!     "bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5
!     as a first guess.
!
! (3) DE/rand-to-best/1/'z'
!     This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.
!     If you get misconvergence try to increase NP. If this does not help you
!     should play around with all three control variables.
!
! (4) DE/best/2/'z' is another powerful strategy worth trying
!
! (5) DE/rand/2/'z' seems to be a robust optimizer for many functions
!
!===========================================================================

        select case (strategy)

        case (1)
           where(lcross)     ! note: in these locations ui_XC = best case
     1     ui_XC=ui_XC+F_XC*(pop_XC(a1,:)-pop_XC(a2,:))

        case default
           where(lcross)
     1     ui_XC=pop_XC(a3,:)+F_XC*(pop_XC(a1,:)-pop_XC(a2,:))

        case (3)
           where(lcross)     ! note: in these locations ui_XC = best case 
     1     ui_XC=pop_XC+F_XC*(ui_XC-pop_XC+pop_XC(a1,:)-pop_XC(a2,:))

        case (4)
           where(lcross)     ! note: in these locations ui_XC = best case 
     1     ui_XC=ui_XC+F_XC*(pop_XC(a1,:)-pop_XC(a2,:)+pop_XC(a3,:)
     2                     -pop_XC(a4,:))

        case (5)
           where(lcross) 
     1     ui_XC=pop_XC(a5,:)+F_XC*(pop_XC(a1,:)-pop_XC(a2,:)+
     2                              pop_XC(a3,:)-pop_XC(a4,:))

        end select


!------Confine each of feasible individuals in the lower-upper bound-------!!

        do i=1,NP
           ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
        enddo

 1000   continue
#ifdef MPI_OPT
!!all processors need same ui_XC values going into DE_Evaluate
      do i = 1, NP 
         if (myid .eq. master) temp_var = ui_XC(i,:)
         call MPI_BCAST(temp_var, Dim_XC, MPI_REAL8, master, 
     1        MPI_COMM_WORLD, ierr)
         if (myid .ne. master) ui_XC(i,:) = temp_var
      end do
#endif    


!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!

        call DE_Evaluate( num_proc, obj, new_val, NP, Dim_XC, nfeval)

        j = 0
        do i=1,NP
           if (new_val(i) < val(i)) then
              j = j + 1
              pop_XC(i,:)=ui_XC(i,:)
              val(i)=new_val(i)
              if (val(i) < bestval) then
                  bestval=val(i)
                  best_XC=ui_XC(i,:)
              end if
           end if
        end do

        if (myid .eq. master) write(6,*) 'selecting,',j,' improvements'
!
        if (save_space) then
        if (myid .eq. master) then
        j = 1
        nfeval = nfeval + 1
        call obj(nopt, Dim_XC, best_XC, fvec, j, nfeval)
        end if
        end if


!------Write restart file ----------------------------------------

        if (myid .eq. master) then
           write(irestart,*) Dim_XC

           CALL RANDOM_SEED (GET=seed(1:seed_size)) ! Sets seed from array

           write(irestart,*) seed(:seed_size)

           do i=1,NP
              write(irestart,*) i,val(i),(pop_XC(i,j),j=1, Dim_XC)
           enddo

           REWIND (unit=irestart)
        end if

!------Write output summary if neede------------------------------

        if( refresh > 0 ) then
           if (mod(iter,refresh) == 0) then
              val_mean = sum(val(:NP))/NP
              val_max = maxval(val(:NP))

              if (myid .eq. master) then
c          write(unit=iwrite,FMT=203) iter
                 write(6, FMT=203) iter        
                 do i=1,Dim_XC
c                write(unit=iwrite, FMT=202) i, best_XC(i)
                    write(6,FMT=202) i,best_XC(i)
                 end do
c          write(unit=iwrite, FMT=201) bestval
                 write(6, FMT=201) bestval,val_mean,val_max 
              end if
           end if
        end if

        if ( bestval <= VTR ) then
c          write(iwrite, *) ' The best fitness is smaller than VTR'
           if (myid .eq. master)
     1        write(6, *) 'The best fitness is smaller than VTR' 
           exit
        endif

      end do

      deallocate (seed, ui_XC, pop_XC, rand_XC)

!!------end the evolutionary computation------------------------------!!
201   format(2x, 'bestval=', ES14.7, ' mean=',ES14.7,' max=',ES14.7, /)
202   format(5x, 'best_XC(', I3, ')=', ES12.5)
203   format(2x, 'No. of iteration =', I8)

      end subroutine DE_Evolve


      subroutine DE_Evaluate( num_proc, fcn, val, NP, Dim_XC, nfeval)
      use kind_spec
      use DE_mod
      implicit none

      integer, parameter :: iflag_cleanup = -100
      integer :: num_proc, NP, Dim_XC, nfeval
      real(rprec), dimension(NP) :: val

      real(rprec), dimension(nopt) :: fvec
      real(rprec) :: funcval

      integer :: iflag, i, j, istat, jstat
      external fcn
#ifdef MPI_OPT
#else
      external de_parallel
#endif      

      n_pop = NP
      n_free = Dim_XC
      nfev = nfeval

#ifdef MPI_OPT
      call de_mpi(np, fcn, val)
#else
      call multiprocess(NP, num_proc, de_parallel, fcn)

!
!  gather results here from all processors
!
      do  j=1, NP

c
c  read in the results from the individual evaluations
c
         read(j+1000, iostat=istat) jstat, iflag, funcval
         if( istat .ne. 0) write(6,*) 'Iostat =',istat,' for case ',j

         if( jstat .ne. j ) then
            write(6,*) "wrong index read in de_evaluate"
            iflag=-14
            exit
         endif

         val(j)=funcval

         close(j+1000, status='delete')

      enddo
#endif
      nfeval = nfeval + NP
      iflag = iflag_cleanup
      call fcn (nopt, n_free, ui_XC(1,:), fvec, iflag, nfeval)

      end subroutine DE_Evaluate

#if defined(MPI_OPT)
      subroutine de_mpi(np, fcn, funcval)
      use de_mod
      implicit none
      include 'mpif.h'                                       !mpi stuff

      real(rprec), dimension(n_free) :: x
      real(rprec), dimension(nopt) :: fvec

      integer, parameter :: master = 0                       !mpi stuff
      integer :: np, j, iflag, istat
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
            x(:) = ui_XC(j,:)
            call MPI_SEND(x, n_free, MPI_REAL8, j, 
     1                  j, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(1) in de_mpi' 
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
            if (ierr .ne. 0) stop 'MPI_RECV error(1) in de_mpi'
            sender     = status(MPI_SOURCE)    
            anstype    = status(MPI_TAG)       ! column is tag value
            if (anstype .gt. np) stop 'ANSTYPE > NP IN de_mpi'
            
            funcval(anstype) = sum(fvec(:nopt)**2)
c           write(6,'(a,1pe10.3,a,i3,a,i3)')' FUNCVAL = ', 
c    1       funcval(anstype),
c    2      ' for iteration ', anstype+nfev,' processor = ', sender

c
c           If more columns are left, then send another column to the worker(sender)
c           that just sent in an answer
c
            if (numsent .lt. np) then
               numsent = numsent+1
               x(:) = ui_XC(numsent,:)
               
               call MPI_SEND(x, n_free, MPI_REAL8, 
     1                       sender, numsent, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0) stop 'MPI_SEND error(2) in de_mpi'

            else                ! Tell worker that there is no more work to do
               
               call MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_WORLD, ierr)
               if (ierr .ne. 0) stop 'MPI_END error(3) in de_mpi'
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
 90      call MPI_RECV(x, n_free, MPI_REAL8, master, 
     1                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         if (ierr .ne. 0) stop 'MPI_RECV error(2) in de_mpi'
        
         column = status(MPI_TAG)                !!ID of pseudo-processor issuing this message
         if (column .eq. 0) then
            go to 200
         else
            iflag = column
c           Call the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA stores the local fvec_min array

            call fcn(nopt, n_free, x, fvec, iflag, nfev)
            if (iflag.ne.0) go to 300
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            call MPI_SEND(fvec, nopt, MPI_REAL8, master, 
     1                    column, MPI_COMM_WORLD, ierr)
            if (ierr .ne. 0) stop 'MPI_SEND error(4) in de_mpi'
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
      print *,' MPI_BCAST error in de_mpi: IERR=', ierr
      
      return

 300  continue
      print *,' IFLAG = ', iflag, ' in de_mpi call to fcn'
      stop

      end subroutine de_mpi
#else
      subroutine de_parallel(j,fcn)
      use de_mod
      implicit none

      real(rprec), dimension(n_free) :: x
      real(rprec), dimension(nopt) :: fvec

      integer :: j, iflag, istat
      external fcn
      real(rprec) :: funcval

      iflag=j
      x(:) = ui_XC(j,:)

      call fcn(nopt, n_free, x, fvec, iflag, nfev)
      funcval = sum(fvec(:nopt)**2) 

      write (j+1000) j, iflag, funcval
      close (j+1000)

      write(6,'(a,f12.5,a,i3)')' FUNCVAL = ', funcval,' for iteration ',
     1   j+nfev

      end subroutine de_parallel
#endif
