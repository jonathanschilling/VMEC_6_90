      module gade_mod
      use kind_spec
      implicit none

!    shared
      integer, parameter :: indmax=500, nchrmax=5000, nparmax=500
      integer, parameter :: max_gen=1000

      integer :: npopsiz, ngen, idum, ibound
      real(rprec)  :: pcross
      real(rprec), dimension(nparmax) :: parmax, parmin
      logical :: save_space
!
!   GA Specific
      integer :: nowrite,microga,unique_ind
      real(rprec)  :: pmutate,pcreep
      integer :: iskip,iend,nchild,itourny,ielite,icreep,iunifrm,
     +           iniche
      integer, dimension(nparmax) :: nposibl, nichflg
!
!   DE specific
      integer :: strategy, CR_strategy, out_iter
      real(rprec) :: f_cross

       namelist /ga_de/ npopsiz,ngen,idum,ibound,pcross,
     +          strategy, CR_strategy, out_iter, f_cross,
     +          nowrite,microga,unique_ind,pmutate,pcreep,
     +          iskip,iend,nchild,itourny,ielite,icreep,iunifrm,
     +          iniche,save_space,
     +          parmin,parmax,nposibl,nichflg

      contains

      subroutine read_gade_namelist (iunit, istat)
      integer :: iunit, istat

      read (iunit, nml=ga_de, iostat=istat)

      end subroutine read_gade_namelist

      end module gade_mod
