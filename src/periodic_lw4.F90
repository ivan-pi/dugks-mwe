module periodic_lw4

   use precision, only: wp

   ! The abstract parent class
   use lattice, only: lattice_grid, &
      cx, cy, csqr, w0, ws, wd, w, rho0, shift

   implicit none
   private

   public :: lattice_grid_lw4

   ! Lax-Wendroff based discretization
   type, extends(lattice_grid) :: lattice_grid_lw4
      private
   contains
      procedure, pass(grid) :: alloc => alloc_grid_lw4
      procedure, pass(grid) :: collision => lw_collide
      procedure, pass(grid) :: streaming => lw_stream
      procedure, pass(grid) :: bc => lw_bc
   end type

   integer, parameter :: nhalo = 2

   ! Limit for OpenMP multi-threading
   integer, parameter :: nlimit = 100**2

contains

   !
   ! Overriden methods
   !

   subroutine alloc_grid_lw4(grid,nx,ny,nu,dt,log)
      class(lattice_grid_lw4), intent(inout) :: grid
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: nu, dt
      logical, intent(in), optional :: log

      call grid%alloc_mem(nx,ny,nhalo=2)
      call grid%set_properties(nu,dt)

   end subroutine

   subroutine lw_collide(grid)
      class(lattice_grid_lw4), intent(inout) :: grid

      call d2q9_collision(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%rho, grid%ux, grid%uy, grid%omega, &
         grid%wrk)

   end subroutine

#include "d2q9_collision.fi"

   subroutine lw_stream(grid)
      class(lattice_grid_lw4), intent(inout) :: grid

      call stream_fdm_4th_order(grid%nx,grid%ny,&
         grid%f(:,:,:,grid%iold),&
         grid%f(:,:,:,grid%inew),&
         grid%dt)

   end subroutine

   subroutine lw_bc(grid)
      class(lattice_grid_lw4), intent(inout) :: grid

      call lw_pbc_kernel(grid%nx,grid%ny,grid%f(:,:,:,grid%iold))

   contains

      subroutine lw_pbc_kernel(nx,ny,fdst)
         implicit none
         integer, intent(in) :: nx, ny
         real(wp), intent(inout) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)

         integer :: k

         if (nhalo /= 2) error stop
         ! In the Lax-Wendroff method, we need to couple all of the
         ! PDFs and not just the incoming ones!

         ! TODO: does multi-threading help here?
         !       we could use nested parallelism (do + section)

         do k = 0, 8

            ! SOUTH HALO
            fdst(1:ny,-1:0,k) = fdst(1:ny,nx-1:nx,k)

            ! NORTH HALO
            fdst(1:ny,nx+1:nx+2,k) = fdst(1:ny,1:2,k)

            ! WEST HALO
            fdst(-1:0,1:nx,k) = fdst(ny-1:ny,1:nx,k)

            ! EAST HALO
            fdst(ny+1:ny+2,1:nx,k) = fdst(1:2,1:nx,k)

            !
            ! CORNERS (2-by-2 blocks)
            !
            fdst(-1:0,-1:0,k) = fdst(ny-1:ny,nx-1:nx,k)
            fdst(ny+1:ny+2,nx+1:nx+2,k) = fdst(1:2,1:2,k)

            fdst(-1:0,nx+1:nx+2,k) = fdst(ny-1:ny,1:2,k)
            fdst(ny+1:ny+2,-1:0,k) = fdst(1:2,nx-1:nx,k)

         end do

      end subroutine

   end subroutine


   subroutine stream_fdm_4th_order(nx,ny,fsrc,fdst,dt)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in)  :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      real(wp), dimension(0:8) :: vx, vy, vxx, vxy, vyy

      ! Streaming variables
      real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
      real(wp) :: fc, fn, fs, fw, fe, fne, fnw, fse, fsw
      real(wp) :: fnn, fss, fww, fee, fne2, fnw2, fse2, fsw2


      vx = dt*cx
      vy = dt*cy
      vxx = 0.5_wp*vx*vx
      vxy = vx*vy
      vyy = 0.5_wp*vy*vy


      !$omp parallel if (nx*ny>nlimit) default(private) &
      !$omp    shared(nx,ny,fsrc,fdst) firstprivate(vx,vy,vxx,vxy,vyy)

      !$omp do nowait
      do j = 1, nx
         do i = 1, ny
           fdst(i,j,0) = fsrc(i,j,0)
         end do
      end do

      ! TODO: consider collapsing loops and loop blocking here

      !$omp do collapse(2)
      do k = 1, 8
         do j = 1, nx
            !$omp simd
            do i = 1, ny


               fc  = fsrc(i  , j, k)

               fe  = fsrc(i+1, j, k)
               fw  = fsrc(i-1, j, k)

               fee  = fsrc(i+2, j, k)
               fww  = fsrc(i-2, j, k)

               fn  = fsrc(i, j+1, k)
               fs  = fsrc(i, j-1, k)

               fnn  = fsrc(i, j+2, k)
               fss  = fsrc(i, j-2, k)

               fne = fsrc(i+1, j+1, k)
               fnw = fsrc(i-1, j+1, k)
               fsw = fsrc(i-1, j-1, k)
               fse = fsrc(i+1, j-1, k)

               fne2 = fsrc(i+2, j+2, k)
               fnw2 = fsrc(i-2, j+2, k)
               fsw2 = fsrc(i-2, j-2, k)
               fse2 = fsrc(i+2, j-2, k)

! ATTENTION; x and y may be swapped here...

               dfy = (1.0_wp/12.0_wp)*(fww - fee) + (2.0_wp/3.0_wp)*(fe - fw)
               dfx = (1.0_wp/12.0_wp)*(fss - fnn) + (2.0_wp/3.0_wp)*(fn - fs)
#if 1
               dfyy = -(1.0_wp/12.0_wp)*fww + (4.0_wp/3.0_wp)*fw - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fe - (1.0_wp/12.0_wp)*fee
               dfxx = -(1.0_wp/12.0_wp)*fss + (4.0_wp/3.0_wp)*fs - (5.0_wp/2.0_wp)*fc  + (4.0_wp/3.0_wp)*fn - (1.0_wp/12.0_wp)*fnn
               dfxy = (1.0_wp/3.0_wp)*(fne - fnw + fsw - fse) - (1.0_wp/48.0_wp)*(fne2 - fnw2 + fsw2 - fse2)
#else
               dfyy = fe - 2.0_wp*fc + fw
               dfxx = fn - 2.0_wp*fc + fs
               dfxy = 0.25_wp*(fne - fse - fnw + fsw)
#endif
               fdst(i,j,k) = fc - vx(k)*dfx - vy(k)*dfy + (vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy)

            end do
         end do

      end do

      !$omp end parallel

   end subroutine

end module

