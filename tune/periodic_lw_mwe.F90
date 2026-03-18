
#ifndef NCOLLAPSE
#define NCOLLAPSE 2
#endif

#ifndef TILE_1
#define TILE_1 1
#endif

#ifndef TILE_2
#define TILE_2 256
#endif

#ifndef KERN
#define KERN 2
#endif

module periodic_lw

   implicit none
   private

   public :: lw_stream_kernel_fdm_v1, time

   integer, parameter :: wp = kind(1.0d0)


   integer, parameter :: nhalo = 1

   real(wp), parameter :: cx(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]

   real(wp), parameter :: w0 = 4.0_wp / 9.0_wp, &
                          ws = 1.0_wp / 9.0_wp, &
                          wd = 1.0_wp / 36.0_wp
   real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

   real(wp), parameter :: csqr = 1.0_wp/3.0_wp

   integer, parameter :: gx = 512, gy = 512

contains

   function time(fsrc,fdst,dt) result(elapsed) bind(c,name="time")
      use omp_lib, only: omp_get_wtime
      real(wp), intent(in) :: fsrc(gy,gx,0:8)
      real(wp), intent(out) :: fdst(gy,gx,0:8)
      real(wp), intent(in), value :: dt

      integer, parameter :: dp = kind(1.0d0), sp = kind(1.0e0)
      real(dp) :: tbegin, tend
      real(sp) :: elapsed
      integer :: k

      tbegin = omp_get_wtime()
      do k = 1, 10
         select case(KERN)
         case(1)
            call lw_stream_kernel_fdm_v1(gx-2,gy-2,fsrc,fdst,dt)
         case(2)
            call lw_stream_kernel_fdm_v2(gx-2,gy-2,fsrc,fdst,dt)
         end select
      end do
      tend = omp_get_wtime()
      elapsed = (tend - tbegin)/10*1.0E3_sp

      print '(A, es15.8e2)', "elapsed (ms) = ", elapsed

   end function


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

      !$omp parallel do collapse(NCOLLAPSE)
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

      !$omp parallel do collapse(NCOLLAPSE) firstprivate(vx,vy,vxx,vxy,vyy)
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
               dNE = 0.25_wp*(fsrc(i+1,j+1,k) - fsrc(i-1,j+1,k))
               dNW = 0.25_wp*(fsrc(i+1,j-1,k) - fsrc(i-1,j-1,k))

               delta = (vx(k) * 0.5_wp *(dE + dW) + vy(k) * 0.5_wp * (dN + dS)) &
                     + (vxx(k) * (dE - dW) + vyy(k) * (dN - dS) &
                     +  vxy(k) * (dNE - dNW))

               fdst(i,j,k) = fsrc(i,j,k) + delta
            end block

            end do
         end do
      end do

   end subroutine





end module

