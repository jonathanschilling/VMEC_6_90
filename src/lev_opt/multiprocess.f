      subroutine multiprocess(numprocess, maxprocess, wrapperfcn, fcn)
      use system_mod
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer :: numprocess, maxprocess
      real :: fcn
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: i, status, iretpid, ierror
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL wrapperfcn, fcn, myfork
C-----------------------------------------------
 
 
!   This Program was written by S. P. Hirshman (7/1/99), under contract
!   the U.S. DOE, and should not be used without explicit permission.
!
!   Input Variable Names
!   MaxProcess:   user defined constant, maximum number of processes that this routine will
!                 request. It should be about max[NCPU/(1+NSPAWN), 1], where NCPU is the number of
!                 machine cpus, and NSPAWN is the maximum number of processes spawned by the
!                 call to loc_function (for stellopt code, NSPAWN = 1)
!   NumProcess:   the TOTAL number of processes to be launched in parallel (if possible). If
!                 this exceeds the number of available processors (max_process, then the next
!                 available processor scheduled by the operating system
!                 returns only after ALL processes are completed.
!   WrapperFcn:   Fortran subroutine that performs the desired task.
!                 It takes for arguments
!                    (1) the index of the process to execute, where 1 <= index <= NumProcess
!                    (2) the subroutine (Fcn) which is called ito evaluate specific information
!                        (such as functional minima, etc.).
!
!   Calling convention from a Fortran program:
!
!   call MultiProcess(nprocess, maxProcess, Wrapper_Subroutine, Worker_Subroutine)
 
      iretpid = 0; ierror = 0

      if (maxprocess .gt. 1) then
         write(6,*)
         write(6,*) ' Begin multi-processing: request ', numprocess, 
     1      ' processes distributed among ', maxprocess, ' processors'

         call vmec_flush(6)

         do i = 1, numprocess
            call myfork (i, maxprocess, wrapperfcn, fcn)
         end do
 
!        Wait for all processes to finish...
         do i = 1, numprocess
            call pxfwait (status, iretpid, ierror)
         end do
 
      else
 
         do i = 1, numprocess
            call wrapperfcn (i, fcn)
         end do
      endif
 
      end subroutine multiprocess


      subroutine myfork(i, maxprocess, wrapper, fcn)
      use system_mod
      implicit none      
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      integer i, maxprocess
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: pid, status, iretpid, ierror, werror
      integer, save :: nprocess = 0
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL wrapper, fcn
C-----------------------------------------------
 
      if (i .eq. 1) nprocess = 0
      ierror = -1
 
!     Child process: limit number to max_process to avoid potential system hang-up
 
      do while(ierror .ne. 0)
 
         if (nprocess .lt. maxprocess) call pxffork (pid, ierror)
 
         if (ierror .ne. 0) then
!           wait for next available processor
            call pxfwait (status, iretpid, werror)
!           if (status.gt.0 .and. nprocess.ge.1) then
            if (nprocess .ge. 1) nprocess = nprocess - 1
!           else
!              nprocess = 0
!           endif
         endif
      end do
 
      if (pid .eq. 0) then
         call wrapper (i, fcn)
#if defined(CRAY)
         call exit(1)
#else         
         stop
#endif         
      endif
 
      nprocess = nprocess + 1
 
      end subroutine myfork
