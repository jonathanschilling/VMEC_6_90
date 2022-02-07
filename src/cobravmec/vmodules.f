       module normalize_data
       use kind_spec
       implicit none
       integer :: nfp_v
       real(rprec), parameter :: mu_0=1.2566368e-6_dp
       real(rprec) :: r0,b0_v,amin,beta0, twopi
       end module normalize_data
 
       module general_dimensions
       integer :: ns_cob, mnmax_v, mndim, nlist 
       end module general_dimensions

       module readin_data
       use kind_spec
       use general_dimensions
       integer,dimension(:),allocatable :: list   
       real(rprec), dimension(:), allocatable :: hiota, hphip,
     1   hpres
       real(rprec), dimension(:), allocatable :: xn_v, xm_v
       real(rprec), dimension(:), allocatable :: lmnsh,bmnch,
     1   bsupvmnh, bsupumnh
       logical :: lscreen
       end module readin_data
       
       module ballooning_data
       use kind_spec
       implicit none
       integer :: tsymm = 1
       real(rprec), allocatable, dimension(:) :: 
     1    init_zeta_v, init_theta_v
       real(rprec) :: init_zeta, init_theta
       real(rprec), parameter :: tole = 1.e-3_dp, x0 = 0
       integer, parameter :: np0_in = 101, kth = 1, k_w = 7
       integer,parameter :: lmax=8, krich=3, km=krich-1
       integer, parameter :: jnewton_max=200
       integer :: np0
       real(rprec), parameter :: newton_tol = tole
       
!      lmax=#max iterations in Richardson's; krich= min.#samples to interpolate
!      x0=radial wave number (default=0.)
!      k_w: number of potential wells included in integration domain!
!      tole= tolerance for eigenvalue fit error
!      np0_in=coarsest grid used in Richardson's; kth-largest eigenvalue obtained
!      jnewton_max= maximum number of iterations in Newton-Raphson to get theta
!      newton_tol= tolerance in Newton_Raphson

       end module ballooning_data
       
       module fmesh_quantities
       use kind_spec
       use general_dimensions
       use readin_data, ONLY: xn_v, xm_v
       implicit none
       real(rprec), dimension(:), allocatable:: iotaf, phipf,
     1   presf, mercierf
       real(rprec), dimension(:), allocatable:: iotapf, prespf
       real(rprec), dimension(:), allocatable:: rmncf, zmnsf,
     1   lmnsf, bmncf, bsupvmnf, bsupumnf
       real(rprec), dimension(:), allocatable:: rmncpf, zmnspf,
     1   lmnspf, bmncpf
       end module fmesh_quantities

       module summod
       use kind_spec
       implicit none 
              
       real(rprec), dimension(:), allocatable :: bfields, bfieldze,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, aux, 
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, lam2, cks_vmec, ckth_vmec, bfieldi, bfield2i,
     7  lam3

       contains

       subroutine alloc_summod(npt)
       integer :: npt, istat

       allocate(bfields(npt), bfieldze(npt), bfieldth(npt), rboo(npt),
     1  rs(npt), rze(npt), rth(npt), zs(npt), zze(npt), zth(npt), 
     2  ccosi(npt), ssine(npt), lam1(npt), zetang(npt), thetang(npt), 
     3  arg(npt), lambdaze(npt), lambdath(npt), bsubze(npt), bsubs(npt),
     4  bsubth(npt), lambdas(npt), rboo2(npt),jacob2(npt),bfield2i(npt), 
     5  rjac2i(npt), bfield(npt), bsupth(npt), bsupze(npt), gtssub(npt),
     6  gstsub(npt), gzssub(npt), gszsub(npt), gtzsub(npt), gztsub(npt),
     7  gttsub(npt), gzzsub(npt), gsssub(npt), gttsup(npt), gzzsup(npt),
     8  gtssup(npt), gstsup(npt), gzssup(npt), gszsup(npt), gtzsup(npt),
     9  gztsup(npt), gsssup(npt), lam2(npt),  bfieldi(npt), aux(npt),
     A  lam3(npt), cks_vmec(npt), ckth_vmec(npt), stat= istat)

       if (istat .ne. 0) stop 'Allocation error in COBRA summod...'

       end subroutine alloc_summod

       subroutine free_summod
       integer :: istat
       
       deallocate( bfields, bfieldze, bfieldi, bfield2i, aux, lam3,
     1  bfieldth, rboo, rs, rze, rth, zs, zze, zth, ccosi, ssine, lam1,
     2  zetang, thetang, arg, lambdaze, lambdath, bsubze, bsubs, bsubth,
     3  lambdas, rboo2, jacob2, rjac2i, bfield, bsupth, bsupze, lam2,
     4  gtssub, gstsub, gzssub, gszsub, gtzsub, gztsub, gttsub, gzzsub,
     5  gsssub, gttsup, gzzsup, gtssup, gstsup, gzssup, gszsup, gtzsup,
     6  gztsup, gsssup, cks_vmec, ckth_vmec, stat=istat)

       end subroutine free_summod

       end module summod
