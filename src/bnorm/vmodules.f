      module meshes
      integer ::  nu,nv,np,mb,nb,mf,nf,md,nd,nuv,nuvp,nvp,nuvh1
      end module meshes

      module  bnvariables
      use kind_spec
      real(rprec) ::  pi,pi2,alp,alu,alv,alvp,fnuv,curpol
      real(rprec), allocatable :: sinu(:,:), conu(:,:), sinv(:,:),
     1    conv(:,:)
      real(rprec), allocatable :: tanu(:), tanv(:)
      integer, allocatable :: indu(:), indv(:)
      real(rprec), allocatable ::   x(:), y(:), z(:), djx(:), 
     1    djy(:), djz(:)
      real(rprec), allocatable ::  xu(:), yu(:), zu(:), xv(:), 
     1    yv(:), zv(:)
      real(rprec), allocatable :: guu(:), guv(:), gvv(:), dju(:),
     1    djv(:), sqf(:)
      real(rprec), allocatable ::  au(:), av(:)
      real(rprec), allocatable :: bsubu(:,:), bsubv(:,:), bu(:), 
     1    bv(:)
      real(rprec), allocatable :: cr(:,:), cz(:,:), cl(:, :)

      end module bnvariables
      
      module neswrite
      use kind_spec
      integer :: ntheta, nzeta, nfp, mnmax, mnmax_ws
      integer, dimension(:), allocatable :: ixm, ixn, ixm_ws, ixn_ws
      real(rprec), allocatable, dimension(:) :: raxis, zaxis, rmnc, zmns
      real(rprec), allocatable, dimension(:) :: rmnc_ws, zmns_ws
      real(rprec) :: coil_separation               
      real(rprec) :: iota_edge, phip_edge
      end module neswrite      

      module normal_info
      use kind_spec
      real(rprec) :: u0, v0, rb_ws, zb_ws, vb_ws, surf_norm(3)
      end module normal_info
