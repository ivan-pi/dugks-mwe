module periodic_lw_fem
use precision, only: wp
implicit none
private

public :: lw2_fem_weights
public :: stream_lw2_fem
public :: stream_lw2_fem_rhs

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
   real(wp) :: delta, tmp(-1:1,-1:1)

   ! Rest particles
   fdst(1:ny,1:nx,0) = fsrc(1:ny,1:nx,0)

   ! We implicitly assume that h = 1

   !$omp parallel if(nx*ny > 100**2) default(private) shared(nx,ny,fsrc,fdst) firstprivate(w)

#ifdef __flang__
   !$omp do collapse(2) schedule(static)
   do k = 1, 8
      do j = 1, nx
      do i = 1, ny
#else
   !$omp do collapse(2) schedule(static)
   do k = 1, 8
      !$omp tile sizes(1,256)
      do j = 1, nx
      do i = 1, ny
#endif
            tmp = w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k)
!            delta = sum(w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k))

            delta = (tmp(1,1) + tmp(-1,-1)) + (tmp(-1,1) + tmp(1,-1))
            delta = delta + ((tmp(1,0) + tmp(-1,0)) + (tmp(0,1) + tmp(0,-1)))
            delta = delta + tmp(0,0)

            fdst(i,j,k) = fsrc(i,j,k) + delta
      end do
      end do
   end do

   !$omp end parallel

end subroutine



subroutine stream_lw2_fem_rhs(nx,ny,fsrc,fdst,w)
   integer, intent(in), value :: nx, ny
   real(wp), intent(in) :: fsrc(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)
   real(wp), intent(in) :: w(-1:1,-1:1,1:8)
   real(wp), intent(inout) :: fdst(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8)

   integer :: i, j, k
   real(wp) :: delta, tmp(-1:1,-1:1)

   ! Rest particles
   fdst(1:ny,1:nx,0) = 0.0_wp

   ! We implicitly assume that h = 1

   !$omp parallel if(nx*ny > 100**2) default(private) shared(nx,ny,fsrc,fdst) firstprivate(w)

#ifdef __flang__
   !$omp do collapse(2) schedule(static)
   do k = 1, 8
      do j = 1, nx
      do i = 1, ny
#else
   !$omp do collapse(2) schedule(static)
   do k = 1, 8
      !$omp tile sizes(1,256)
      do j = 1, nx
      do i = 1, ny
#endif
            tmp = w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k)
!            delta = sum(w(-1:1,-1:1,k)*fsrc(i-1:i+1,j-1:j+1,k))

            delta = (tmp(1,1) + tmp(-1,-1)) + (tmp(-1,1) + tmp(1,-1))
            delta = delta + ((tmp(1,0) + tmp(-1,0)) + (tmp(0,1) + tmp(0,-1)))
            delta = delta + tmp(0,0)

            fdst(i,j,k) = delta
      end do
      end do
   end do

   !$omp end parallel

end subroutine

end module
