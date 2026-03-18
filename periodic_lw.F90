
module mlbm_stream_lw
use precision, only: wp
implicit none
private

public :: lw2_fem_weights
public :: stream_lw2_fem

real(wp), parameter :: C(-1:1,-1:1,1:8) = reshape([&
    [1,0,-1,4,0,-4,1,0,-1], &
    [-1,-4,-1,0,0,0,1,4,1], &
    [-1,0,1,-4,0,4,-1,0,1], &
    [1,4,1,0,0,0,-1,-4,-1], &
    [2,4,0,4,0,-4,0,-4,-2], &
    [0,-4,-2,4,0,-4,2,4,0], &
    [-2,-4,0,-4,0,4,0,4,2], &
    [0,4,2,-4,0,4,-2,-4,0]] / 12.0_wp, [3,3,8])

real(wp), parameter :: D(-1:1,-1:1,1:8) = reshape([&
    [-1,2,-1,-4,8,-4,-1,2,-1], &
    [-1,-4,-1,2,8,2,-1,-4,-1], &
    [-1,2,-1,-4,8,-4,-1,2,-1], &
    [-1,-4,-1,2,8,2,-1,-4,-1], &
    [-5,-2,1,-2,16,-2,1,-2,-5], &
    [1,-2,-5,-2,16,-2,-5,-2,1], &
    [-5,-2,1,-2,16,-2,1,-2,-5], &
    [1,-2,-5,-2,16,-2,-5,-2,1]] / 6.0_wp, [3,3,8])

integer :: nhalo = 1

contains

function lw2_fem_weights(dt,imeth) result(w)
   real(wp), intent(in) :: dt
   integer, intent(in) :: imeth

   real(wp) :: w(-1:1,-1:1,1:8)

   real(wp), parameter :: cx(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]

   ! FDM weights (mainly for testing purposes)
   if (imeth == 1) then
      central_fdm: block
         real(wp), dimension(-1:1,-1:1) :: fx, fy, fxx, fxy, fyy
         integer :: k

         fx  = 0
         fy  = 0
         fxx = 0
         fxy = 0
         fyy = 0

         !fy( :,0) = [1,0,-1]/2.0_wp
         !fx(0, :) = [1,0,-1]/2.0_wp

         fy(:,-1) = [1,0,-1]/12.0_wp
         fy(:, 0) = [1,0,-1]/3.0_wp
         fy(:, 1) = [1,0,-1]/12.0_wp

         fx(-1,:) = [1,0,-1]/12.0_wp
         fx( 0,:) = [1,0,-1]/3.0_wp
         fx( 1,:) = [1,0,-1]/12.0_wp


         fyy(:,0) = [1,-2,1]
         fxx(0,:) = [1,-2,1]

         fxy( 1, 1) =  1
         fxy(-1, 1) = -1
         fxy(-1,-1) =  1
         fxy( 1,-1) = -1

         fxy = fxy/4.0_wp

         do k = 1, 8
            w(:,:,k) = -dt*(cx(k)*fx + cy(k)*fy) +&
               0.5_wp*dt**2*(cx(k)**2*fxx + 2.0_wp*cx(k)*cy(k)*fxy + cy(k)**2*fyy)
         end do
      end block central_fdm
      return
   end if

   ! Default weights are the FEM ones
   fem: block
      real(wp), dimension(-1:1,-1:1,1:8) :: C, D
      integer :: k
      do k = 1, 8

         C( 0, 0, k) = 0
         C( 1, 0, k) =  cx(k)/3.0_wp
         C( 0, 1, k) =  cy(k)/3.0_wp
         C(-1, 0, k) = -cx(k)/3.0_wp
         C( 0,-1, k) = -cy(k)/3.0_wp
         C( 1, 1, k) = ( cx(k) + cy(k))/12.0_wp
         C(-1, 1, k) = (-cx(k) + cy(k))/12.0_wp
         C(-1,-1, k) = (-cx(k) - cy(k))/12.0_wp
         C( 1,-1, k) = ( cx(k) - cy(k))/12.0_wp

         D( 0, 0, k) = ( 4*cx(k)**2 + 4*cy(k)**2) / 3.0_wp
         D( 1, 0, k) = (-2*cx(k)**2 +   cy(k)**2) / 3.0_wp
         D( 0, 1, k) = (   cx(k)**2 - 2*cy(k)**2) / 3.0_wp
         D(-1, 0, k) = (-2*cx(k)**2 +   cy(k)**2) / 3.0_wp
         D( 0,-1, k) = (   cx(k)**2 - 2*cy(k)**2) / 3.0_wp
         D( 1, 1, k) = (-cx(k)**2 - 3*cx(k)*cy(k) - cy(k)**2) / 6.0_wp
         D(-1, 1, k) = (-cx(k)**2 + 3*cx(k)*cy(k) - cy(k)**2) / 6.0_wp
         D(-1,-1, k) = (-cx(k)**2 - 3*cx(k)*cy(k) - cy(k)**2) / 6.0_wp
         D( 1,-1, k) = (-cx(k)**2 + 3*cx(k)*cy(k) - cy(k)**2) / 6.0_wp

         ! x- and y-indexes are reversed below
         C(:,:,k) = transpose(C(:,:,k))
         D(:,:,k) = transpose(D(:,:,k))
      end do

      w = -dt*C - 0.5_wp*dt**2*D

   end block fem

end function

subroutine stream_lw2_fem(nx,ny,fsrc,fdst,w)
   integer, intent(in), value :: nx, ny
   real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
   real(wp), intent(in) :: w(-1:1,-1:1,1:8)
   real(wp), intent(inout) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)

   integer :: i, j, k

   ! Rest particles
   fdst(1:ny,1:nx,0) = fsrc(1:ny,1:nx,0)

   ! We implicitly assume that h = 1

   !$omp parallel default(private) shared(nx,ny,fsrc,fdst) firstprivate(w)

   !$omp do collapse(2)
   do k = 1, 8
      !$omp tile sizes(1,128)
      do j = 1, nx
      do i = 1, ny
         block
            real(wp) :: delta, tmp(-1:1,-1:1)

            tmp = w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k)
!            delta = sum(w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k))

            delta = (tmp(1,1) + tmp(-1,-1)) + (tmp(-1,1) + tmp(1,-1))
            delta = delta + ((tmp(1,0) + tmp(-1,0)) + (tmp(0,1) + tmp(0,-1)))
            delta = delta + tmp(0,0)

            fdst(i,j,k) = fsrc(i,j,k) + delta
         end block
      end do
      end do
   end do

   !$omp end parallel

end subroutine

end module






#ifndef NCOLLAPSE
#define NCOLLAPSE 2
#endif

#ifndef TILE_1
#define TILE_1 1
#endif

#ifndef TILE_2
#define TILE_2 256
#endif

module periodic_lw

   use precision, only: wp
   use lattice, only: lattice_grid, &
      cx, cy, csqr, w0, ws, wd, w, rho0, shift


   implicit none
   private

   public :: lw_stream
   public :: lw_collide
   public :: lw_bc

   integer, parameter :: nhalo = 1

   abstract interface
      subroutine lw_stream_kernel(nx,ny,fsrc,fdst,dt)
         import wp, nhalo
         integer, intent(in), value :: nx, ny
         real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
         real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
         real(wp), intent(in) :: dt
      end subroutine
   end interface

   integer, parameter :: nlimit = 100**2

contains

   subroutine lw_collide(grid)
      class(lattice_grid), intent(inout) :: grid

#if 0
      call lw_collision_kernel(grid%nx,grid%ny,&
         grid%f(:,:,:,grid%iold), &
         grid%omega)
#else
      ! (nx,ny,f,rho,ux,uy,omega,indp)
      call d2q9_collision(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%rho, grid%ux, grid%uy, grid%omega, &
         grid%wrk)
#endif

   end subroutine

   subroutine lw_stream(grid)
      use mlbm_stream_lw
      class(lattice_grid), intent(inout) :: grid

      procedure(lw_stream_kernel), pointer :: kernel => null()

      real(wp), allocatable :: wfem(:,:,:)

#if 1
      if (.not. associated(kernel)) then
         ! TODO: allow runtime selection
         kernel => lw_stream_kernel_fdm_v2
      end if

      call kernel(grid%nx,grid%ny,&
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         grid%dt)
#else
      if (.not. allocated(wfem)) then
         wfem = lw2_fem_weights(grid%dt, IMETH = 0)
      end if

      call stream_lw2_fem(grid%nx,grid%ny,&
         grid%f(:,:,:,grid%iold), &
         grid%f(:,:,:,grid%inew), &
         wfem)
#endif

   end subroutine


   subroutine lw_bc(grid)
      class(lattice_grid), intent(inout) :: grid

      call lw_pbc_kernel(grid%nx,grid%ny,grid%f(:,:,:,grid%iold))

   end subroutine

   subroutine lw_collision_kernel(nx,ny,pdf,omega)
      implicit none

      integer, intent(in), value :: nx, ny
      real(wp), intent(inout) :: pdf(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: omega

      ! Collision variables
      real(wp) :: drho, rho, irho, ux, uy, feq(0:8), f(0:8)
      real(wp) :: indp, uxx, uyy, uxpy, uxmy
      real(wp) :: omegabar
      real(wp) :: omega_w0, omega_ws, omega_wd
      integer :: y, x

      omegabar = 1.0_wp - omega
      !
      ! Collision
      !
      omega_w0 = 3.0_wp * omega * w0
      omega_ws = 3.0_wp * omega * ws
      omega_wd = 3.0_wp * omega * wd

      !$omp parallel do collapse(2) default(private) shared(pdf,nx,ny) &
      !$omp                         firstprivate(omega,omegabar, omega_w0, omega_ws, omega_wd)
      do x = 1, nx
         do y = 1, ny

            ! Gather PDFs
            f(0) = pdf(y,x,0)
            f(1) = pdf(y,x,1)
            f(2) = pdf(y,x,2)
            f(3) = pdf(y,x,3)
            f(4) = pdf(y,x,4)
            f(5) = pdf(y,x,5)
            f(6) = pdf(y,x,6)
            f(7) = pdf(y,x,7)
            f(8) = pdf(y,x,8)

            ! density
            drho = f(0) + (((f(5) + f(7)) + (f(6) + f(8))) + &
                   ((f(1) + f(3)) + (f(2) + f(4))))

            if (shift) then
               rho = drho + rho0
            else
               rho = drho
            end if

            !irho = 1.0_wp/rho

            ! velocity
            ux = (((f(5) - f(7)) + (f(8) - f(6))) + (f(1) - f(3))) / rho
            uy = (((f(5) - f(7)) + (f(6) - f(8))) + (f(2) - f(4))) / rho

            uxx = ux*ux
            uyy = uy*uy

            if (shift) then
               indp = -0.5_wp * (uxx + uyy)
!               indp = -1.5_wp * (uxx + uyy)
            else
               indp = 1.0_wp/3.0_wp - 0.5_wp * (uxx + uyy)
!               indp = 1.0_wp - 1.5_wp * (uxx + uyy)
            endif

            feq(0) = rho*(indp)
            feq(1) = rho*(indp + ux + 1.5_wp*uxx)
            feq(2) = rho*(indp + uy + 1.5_wp*uyy)
            feq(3) = rho*(indp - ux + 1.5_wp*uxx)
            feq(4) = rho*(indp - uy + 1.5_wp*uyy)

            uxpy = ux + uy
            feq(5) = rho*(indp + uxpy + 1.5_wp*uxpy*uxpy)
            feq(7) = rho*(indp - uxpy + 1.5_wp*uxpy*uxpy)

            uxmy = ux - uy
            feq(6) = rho*(indp - uxmy + 1.5_wp*uxmy*uxmy)
            feq(8) = rho*(indp + uxmy + 1.5_wp*uxmy*uxmy)

            if (shift) feq = drho/3.0_wp + feq

            pdf(y,x,0) = omegabar*f(0) + omega_w0*feq(0)
            pdf(y,x,1) = omegabar*f(1) + omega_ws*feq(1)
            pdf(y,x,2) = omegabar*f(2) + omega_ws*feq(2)
            pdf(y,x,3) = omegabar*f(3) + omega_ws*feq(3)
            pdf(y,x,4) = omegabar*f(4) + omega_ws*feq(4)
            pdf(y,x,5) = omegabar*f(5) + omega_wd*feq(5)
            pdf(y,x,6) = omegabar*f(6) + omega_wd*feq(6)
            pdf(y,x,7) = omegabar*f(7) + omega_wd*feq(7)
            pdf(y,x,8) = omegabar*f(8) + omega_wd*feq(8)

         end do
      end do

   end subroutine

   subroutine lw_stream_kernel_fvm(nx,ny,fsrc,fdst,dt)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      fdst(1:ny,1:nx,0) = fsrc(1:ny,1:nx,0)

      !$omp parallel do collapse(2)
      do k = 1, 8
         do j = 1, nx
            do i = 1, ny

            block
               real(wp), parameter :: p2 = 0.5_wp, p8 = 0.125_wp
               real(wp) :: fc, fw, fe, fn, fs, fnw, fne, fsw, fse
               real(wp) :: cfw, cfe, cfn, cfs
               real(wp) :: vx, vy
               real(wp) :: delta

               vx = dt*cx(k)
               vy = dt*cy(k)

               fc  = fsrc(i,j,k)
               fn  = fsrc(i+1,j,k)
               fs  = fsrc(i-1,j,k)
               fe  = fsrc(i,j+1,k)
               fw  = fsrc(i,j-1,k)
               fne = fsrc(i+1,j+1,k)
               fnw = fsrc(i+1,j-1,k)
               fse = fsrc(i-1,j+1,k)
               fsw = fsrc(i-1,j-1,k)

               ! WEST cell face
               cfw = p2*(fc + fw) - p2*vx*(fc - fw) - &
                          p8*vy*(fnw + fn - fsw - fs)

               ! NORTH cell face
               cfn = p2*(fc + fn) - p2*vy*(fn - fc) - &
                          p8*vx*(fne + fe - fnw - fw)

               ! EAST cell face
               cfe = p2*(fc + fe) - p2*vx*(fe - fc) - &
                          p8*vy*(fne + fn - fse - fs)

               ! SOUTH cell face
               cfs = p2*(fc + fs) - p2*vy*(fc - fs) - &
                          p8*vx*(fse + fe - fsw - fw)

               delta = -vx*(cfe - cfw) - vy*(cfn - cfs)

               fdst(i,j,k) = fsrc(i,j,k) + delta

            end block

            end do
         end do
      end do

   end subroutine

   subroutine lw_stream_kernel_fdm_v1(nx,ny,fsrc,fdst,dt)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      real(wp), parameter, dimension(0:8) :: vx = cx, vy = cy
      real(wp), parameter, dimension(0:8) :: vxx = cx**2, vyy = cy**2
      real(wp), parameter, dimension(0:8) :: vxy = 2*cx*cy

      !
      ! Streaming
      !
      fdst(1:ny,1:nx,0) = fsrc(1:ny,1:nx,0)

      !$omp parallel do collapse(NCOLLAPSE) schedule(static)
      do k = 1, 8
         !$omp tile sizes(TILE_1,TILE_2)
         do j = 1, nx
            do i = 1, ny

            block
               real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
               real(wp) :: first, second, delta

               dfy = 0.5_wp*(fsrc(i+1,j,k) - fsrc(i-1,j,k))
               dfx = 0.5_wp*(fsrc(i,j+1,k) - fsrc(i,j-1,k))

               dfyy = fsrc(i+1,j,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i-1,j,k)
               dfxx = fsrc(i,j+1,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i,j-1,k)
               dfxy = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + &
                               fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k))

               first = -dt*(vx(k)*dfx + vy(k)*dfy)
               second = 0.5_wp*dt**2*(vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)
               delta = first + second

               fdst(i,j,k) = fsrc(i,j,k) + delta
            end block

            end do
         end do
      end do

   end subroutine


   subroutine lw_stream_kernel_fdm_v2(nx,ny,fsrc,fdst,dt)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      real(wp), dimension(0:8) :: vx, vy, vxx, vxy, vyy

      vx = -dt*cx
      vy = -dt*cy
      vxx = 0.5_wp*vx*vx
      vxy = vx*vy
      vyy = 0.5_wp*vy*vy

      !
      ! Streaming
      !
      fdst(1:ny,1:nx,0) = fsrc(1:ny,1:nx,0)

      ! From the difference stencil (!)
      vx = vx * 0.5_wp
      vy = vy * 0.5_wp
      vxy = vxy * 0.25_wp

      !$omp parallel if(nx*ny > nlimit) default(private) shared(nx,ny,fsrc,fdst) &
      !$omp     firstprivate(vx,vy,vxx,vxy,vyy)

      !$omp do collapse(NCOLLAPSE) schedule(static)
      do k = 1, 8
         !$omp tile sizes(TILE_1,TILE_2)
         do j = 1, nx
            do i = 1, ny

            block
               real(wp) :: dE, dW, dN, dS, dNE, dNW, delta

               dE = fsrc(i,j+1,k) - fsrc(i,j,k)
               dW = fsrc(i,j,k) - fsrc(i,j-1,k)
               dN = fsrc(i+1,j,k) - fsrc(i,j,k)
               dS = fsrc(i,j,k) - fsrc(i-1,j,k)

               ! Mixed differences
               dNE = (fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k))
               dNW = (fsrc(i+1,j-1,k) - fsrc(i-1,j-1,k))

               fdst(i,j,k) = fsrc(i,j,k) + vx(k) *(dE + dW) + vy(k) * (dN + dS) + &
                     + (vxx(k) * (dE - dW) + vyy(k) * (dN - dS) + vxy(k) * (dNE - dNW))
            end block

            end do
         end do
      end do

      !$omp end parallel

   end subroutine


   subroutine lw_pbc_kernel(nx,ny,fdst)
      implicit none
      integer, intent(in) :: nx, ny
      real(wp), intent(inout) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)

      integer :: k
      ! In the Lax-Wendroff method, we need to couple all of the
      ! PDFs and not just the incoming ones!

      do k = 0, 8
         ! SOUTH HALO
         fdst(1:ny,   0,k) = fdst(1:ny,nx,k)

         ! NORTH HALO
         fdst(1:ny,nx+1,k) = fdst(1:ny, 1,k)

         ! WEST HALO
         fdst(   0,1:nx,k) = fdst(ny,1:nx,k)

         ! EAST HALO
         fdst(ny+1,1:nx,k) = fdst( 1,1:nx,k)

         !
         ! CORNERS
         !
         fdst(0,0,k) = fdst(ny, nx, k)
         fdst(ny+1,nx+1,k) = fdst(1, 1,k)

         fdst(0, nx+1, k) = fdst(ny,1, k)
         fdst(ny+1,0, k) = fdst(1,nx, k)
      end do

   end subroutine


   subroutine d2q9_collision(nx,ny,f,rho,ux,uy,omega,indp)

   !$ use omp_lib, only: omp_get_num_threads, omp_get_thread_num

      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: omega
      real(wp), intent(inout) :: f(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: rho(ny,nx), ux(ny,ny), uy(ny,nx), indp(ny,nx)

      real(wp) :: omegabar, omega_w0, omega_ws, omega_wd
      real(wp) :: t0,t1,t2,t3,t4,t5,t6,t7,t8
      real(wp) :: cu, usqr, invrho, drho
      real(wp) :: vel_trm_13, vel_trm_24
      real(wp) :: vel_trm_57, vel_trm_68
      real(wp) :: velxpy, velxmy

      integer :: i, j, k

      real(wp), parameter :: one_third = 1.0_wp / 3.0_wp
      real(wp), parameter :: one_half = 1.0_wp / 2.0_wp
      real(wp), parameter :: th = 3.0_wp / 2.0_wp

      ! TODO: THIS KERNEL DOESN'T SUPPORT SHIFTING...
      logical, parameter :: shift = .false.

      omegabar = 1.0_wp - omega
      omega_w0 = 3.0_wp * omega * w0
      omega_ws = 3.0_wp * omega * ws
      omega_wd = 3.0_wp * omega * wd

      !$omp parallel if (nx*ny > nlimit) default(private) shared(nx,ny,f,rho,ux,uy,indp) &
      !$omp    firstprivate(omega,omegabar,omega_w0,omega_ws,omega_wd)

      !$omp do schedule(static)
      xloop: do j = 1, nx

      !$omp simd private(drho, t0,t1,t2,t3,t4,t5,t6,t7,t8)
      do i = 1, ny

         ! pull pdfs travelling in different directions
         t0 = f(i,j,0)
         t1 = f(i,j,1)
         t2 = f(i,j,2)
         t3 = f(i,j,3)
         t4 = f(i,j,4)
         t5 = f(i,j,5)
         t6 = f(i,j,6)
         t7 = f(i,j,7)
         t8 = f(i,j,8)

         ! density
         drho = (((t5 + t7) + (t6 + t8)) + &
                ((t1 + t3) + (t2 + t4))) + t0

         if (shift) then
            rho(i,j) = drho + rho0
         else
            rho(i,j) = drho
         end if

         ! velocity
         ux(i,j) = (((t5 - t7) + (t8 - t6)) + (t1 - t3)) / rho(i,j)
         uy(i,j) = (((t5 - t7) + (t6 - t8)) + (t2 - t4)) / rho(i,j)

         if (shift) then
            indp(i,j) = -0.5_wp * (ux(i,j)**2 + uy(i,j)**2)
         !   indp(i,j) = -1.5_wp * (uxx + uyy)
         else
            indp(i,j) = 1.0_wp/3.0_wp - 0.5_wp * (ux(i,j)**2 + uy(i,j)**2)
         !   indp(i,j) = 1.0_wp - 1.5_wp * (uxx + uyy)
         endif

         ! direction independent part
         !indp(i,j) = one_third - 0.5_wp * (ux(i,j)**2 + uy(i,j)**2)

         ! pdf zero
         f(i,j,0) = omegabar*f(i,j,0) + omega_w0*rho(i,j)*indp(i,j)

      end do

      !$omp simd private(vel_trm_13,vel_trm_24)
      do i = 1, ny
         vel_trm_13 = indp(i,j) + th * ux(i,j) * ux(i,j)

         f(i,j,1) = omegabar*f(i,j,1) + omega_ws * rho(i,j) * (vel_trm_13 + ux(i,j))
         f(i,j,3) = omegabar*f(i,j,3) + omega_ws * rho(i,j) * (vel_trm_13 - ux(i,j))

         vel_trm_24 = indp(i,j) + th * uy(i,j) * uy(i,j)

         f(i,j,2) = omegabar*f(i,j,2) + omega_ws * rho(i,j) * (vel_trm_24 + uy(i,j))
         f(i,j,4) = omegabar*f(i,j,4) + omega_ws * rho(i,j) * (vel_trm_24 - uy(i,j))
      end do

      !$omp simd private(velxpy,vel_trm_57,velxmy,vel_trm_68)
      do i = 1, ny
         velxpy = ux(i,j) + uy(i,j)
         vel_trm_57 = indp(i,j) + th * velxpy * velxpy
         f(i,j,5) = omegabar*f(i,j,5) + omega_wd * rho(i,j) * (vel_trm_57 + velxpy)
         f(i,j,7) = omegabar*f(i,j,7) + omega_wd * rho(i,j) * (vel_trm_57 - velxpy)

         velxmy = ux(i,j) - uy(i,j)
         vel_trm_68 = indp(i,j) + th * velxmy * velxmy

         f(i,j,6) = omegabar*f(i,j,6) + omega_wd * rho(i,j) * (vel_trm_68 - velxmy)
         f(i,j,8) = omegabar*f(i,j,8) + omega_wd * rho(i,j) * (vel_trm_68 + velxmy)
      end do

      end do xloop
      !$omp end do

      !$omp end parallel

   end subroutine

end module

