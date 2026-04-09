module periodic_lw6

   use precision, only: wp

   ! The abstract parent class
   use lattice, only: lattice_grid, &
      cx, cy, csqr, w0, ws, wd, w, rho0, shift

   implicit none
   private

   public :: lattice_grid_lw6

   ! Lax-Wendroff based discretization
   type, extends(lattice_grid) :: lattice_grid_lw6
      private
   contains  
      procedure, pass(grid) :: alloc => alloc_grid_lw6
      procedure, pass(grid) :: collision => lw_collide
      procedure, pass(grid) :: streaming => lw_stream
      procedure, pass(grid) :: bc => lw_bc
   end type

   integer, parameter :: NHALO = 3

   ! Limit for OpenMP multi-threading
   integer, parameter :: nlimit = 100**2

contains

   !
   ! Overriden methods
   !

   subroutine alloc_grid_lw6(grid,nx,ny,nu,dt,log)
      class(lattice_grid_lw6), intent(inout) :: grid
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: nu, dt
      logical, intent(in), optional :: log

      call grid%alloc_mem(nx,ny,nhalo=nhalo)
      call grid%set_properties(nu,dt)

   end subroutine

   subroutine lw_collide(grid)
      class(lattice_grid_lw6), intent(inout) :: grid

      call d2q9_collision(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%rho, grid%ux, grid%uy, grid%omega, &
         grid%wrk)

   end subroutine

#include "d2q9_collision.fi"

   subroutine lw_stream(grid)
      class(lattice_grid_lw6), intent(inout) :: grid

      call stream_fdm_6th_order(grid%nx,grid%ny,&
         grid%f(:,:,:,grid%iold),&
         grid%f(:,:,:,grid%inew),&
         grid%dt)

   end subroutine

   subroutine stream_fdm_6th_order(nx,ny,fsrc,fdst,dt)
      integer, intent(in), value :: nx, ny
      real(wp), intent(in)  :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(out) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
      real(wp), intent(in) :: dt

      integer :: i, j, k

      real(wp), dimension(0:8) :: vx, vy, vxx, vxy, vyy

      ! Streaming variables
      real(wp) :: dfx, dfy, dfxx, dfxy, dfyy
      real(wp) :: dx1, dx2, dx3
      real(wp) :: dy1, dy2, dy3

      real(wp), parameter :: mb =  3.0_wp/8.0_wp
      real(wp), parameter :: mc = -3.0_wp/80.0_wp
      real(wp), parameter :: md =  1.0_wp/360.0_wp


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
            do i = 1, ny

               dy1 = fsrc(i+1,j,k) - fsrc(i-1,j,k)
               dy2 = fsrc(i+2,j,k) - fsrc(i-2,j,k)
               dy3 = fsrC(i+3,j,k) - fsrc(i-3,j,k)

               dx1 = fsrc(i,j+1,k) - fsrc(i,j-1,k)
               dx2 = fsrc(i,j+2,k) - fsrc(i,j-2,k)
               dx3 = fsrc(i,j+3,k) - fsrc(i,j-3,k)

! First order difference
! −1/60 3/20  −3/4  0  3/4   −3/20 1/60

               dfy = (3.0_wp/4.0_wp) * dy1 - (3.0_wp/20.0_wp) * dy2 + &
                    (1.0_wp/60.0_wp) * dy3

               dfx = (3.0_wp/4.0_wp) * dx1 - (3.0_wp/20.0_wp) * dx2 + &
                    (1.0_wp/60.0_wp) * dx3

#if 1
! Second order difference
! 1/90   −3/20 3/2   −49/18   3/2   −3/20 1/90

               dfyy =           1.5_wp * (fsrc(i+1,j,k) + fsrc(i-1,j,k)) - &
                      (3.0_wp/20.0_wp) * (fsrc(i+2,j,k) + fsrc(i-2,j,k)) + &
                      (1.0_wp/90.0_wp) * (fsrc(i+3,j,k) + fsrc(i-3,j,k)) - &
                      (49.0_wp/18.0_wp) * fsrc(i,j,k)

               dfxx =           1.5_wp * (fsrc(i,j+1,k) + fsrc(i,j-1,k)) - &
                      (3.0_wp/20.0_wp) * (fsrc(i,j+2,k) + fsrc(i,j-2,k)) + &
                      (1.0_wp/90.0_wp) * (fsrc(i,j+3,k) + fsrc(i,j-3,k)) - &
                      (49.0_wp/18.0_wp) * fsrc(i,j,k)

! Mixed second derivative
               dfxy = mb*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k)) &
                    + mc*(fsrc(i+2,j+2,k) - fsrc(i-2,j+2,k) + fsrc(i-2,j-2,k) - fsrc(i+2,j-2,k)) &
                    + md*(fsrc(i+3,j+3,k) - fsrc(i-3,j+3,k) + fsrc(i-3,j-3,k) - fsrc(i+3,j-3,k))
!               dfxy = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k))
#else
!               dfy = 0.5_wp*(fsrc(i+1,j,k) - fsrc(i-1,j,k))
!               dfx = 0.5_wp*(fsrc(i,j+1,k) - fsrc(i,j-1,k))

!               dfyy = fsrc(i+1,j,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i-1,j,k)
!               dfxx = fsrc(i,j+1,k) - 2.0_wp*fsrc(i,j,k) + fsrc(i,j-1,k)
               dfxy = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k) + &
                               fsrc(i-1,j-1,k) - fsrc(i+1,j-1,k))
#endif
               fdst(i,j,k) = fsrc(i,j,k) &
                           - vx(k)*dfx - vy(k)*dfy &
                           + vxx(k)*dfxx + vxy(k)*dfxy + vyy(k)*dfyy

            end do
         end do

      end do

      !$omp end parallel

   end subroutine


   subroutine lw_bc(grid)
      class(lattice_grid_lw6), intent(inout) :: grid

      call lw_pbc_kernel(grid%nx,grid%ny,grid%f(:,:,:,grid%iold))

   contains

      subroutine lw_pbc_kernel(nx,ny,fdst)
         implicit none
         integer, intent(in) :: nx, ny
         real(wp), intent(inout) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)

         integer :: k

         ! In the Lax-Wendroff method, we need to couple all of the
         ! PDFs and not just the incoming ones!

         ! TODO: does multi-threading help here?
         !       we could use nested parallelism (do + section)

         do k = 0, 8

            ! SOUTH HALO
            fdst(1:ny,-2:0,k) = fdst(1:ny,nx-2:nx,k)

            ! NORTH HALO
            fdst(1:ny,nx+1:nx+3,k) = fdst(1:ny,1:3,k)

            ! WEST HALO
            fdst(-2:0,1:nx,k) = fdst(ny-2:ny,1:nx,k)

            ! EAST HALO
            fdst(ny+1:ny+3,1:nx,k) = fdst(1:3,1:nx,k)

            !
            ! CORNERS (2-by-2 blocks)
            !
            fdst(-2:0,-2:0,k) = fdst(ny-2:ny,nx-2:nx,k)
            fdst(ny+1:ny+3,nx+1:nx+3,k) = fdst(1:3,1:3,k)

            fdst(-2:0,nx+1:nx+3,k) = fdst(ny-2:ny,1:3,k)
            fdst(ny+1:ny+3,-2:0,k) = fdst(1:3,nx-2:nx,k)

         end do

      end subroutine

   end subroutine

end module

