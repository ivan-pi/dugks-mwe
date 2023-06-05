module lattice

   use precision
   use gnuplot_io, only: output_gnuplot_grid

   implicit none
   private

   public :: wp
   public :: lattice_grid

   public :: alloc_grid, dealloc_grid

   public :: set_properties

   public :: update_macros
   
   public :: output_gnuplot

   public :: set_pdf_to_equilibrium
   public :: equilibrium

   public :: cx, cy, csqr

   character(*), parameter :: FMT_REAL_SP = '(es15.8e2)'

   type :: lattice_grid

      integer :: nx, ny

      real(wp), allocatable :: f(:,:,:,:)

      real(wp), allocatable :: rho(:,:)
      real(wp), allocatable ::  ux(:,:)
      real(wp), allocatable ::  uy(:,:)
      
      real(wp) :: nu, dt, tau
      real(wp) :: omega, trt_magic
      real(wp) :: csqr
      
      integer :: iold, inew, imid
      
      procedure(collision_interface), pointer, pass(grid) :: collision => null()
      procedure(streaming_interface), pointer, pass(grid) :: streaming => null()

      character(len=:), allocatable :: filename, foldername

      character(len=:), allocatable :: logfile
      procedure(gridlog_interface), pointer, pass(grid) :: logger => null()
      integer :: logunit
   contains
      procedure :: set_output_folder
   end type

   abstract interface
      subroutine collision_interface(grid)
         import lattice_grid
         class(lattice_grid), intent(inout) :: grid
      end subroutine
      subroutine streaming_interface(grid)
         import lattice_grid
         class(lattice_grid), intent(inout) :: grid
      end subroutine
      subroutine gridlog_interface(grid, step)
         import lattice_grid
         class(lattice_grid), intent(in) :: grid
         integer, intent(in) :: step
      end subroutine
   end interface

   real(wp), parameter :: cx(0:8) = [real(wp) :: 0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [real(wp) :: 0, 0, 1, 0, -1, 1, 1, -1, -1]

   real(wp), parameter :: w0 = 4._wp / 9._wp, &
                          ws = 1._wp / 9._wp, &
                          wd = 1._wp / 36._wp

   real(wp), parameter :: csqr = 1._wp/3._wp
   real(wp), parameter :: invcsqr = 1._wp/csqr

contains

   pure function equilibrium(rho,ux,uy) result(feq)
      real(wp), intent(in) :: rho, ux, uy
      real(wp) :: feq(0:8)

      real(wp) :: uxx, uyy, uxy, uxpy, uxmy
      real(wp) :: indp

      uxx = ux*ux
      uyy = uy*uy
      uxy = ux*uy

      indp = 1.0_wp - 1.5_wp * (uxx + uyy)

      feq(0) = w0*rho*(indp)
      feq(1) = ws*rho*(indp + 3.0_wp*ux + 4.5_wp*uxx)
      feq(2) = ws*rho*(indp + 3.0_wp*uy + 4.5_wp*uyy)
      feq(3) = ws*rho*(indp - 3.0_wp*ux + 4.5_wp*uxx)
      feq(4) = ws*rho*(indp - 3.0_wp*uy + 4.5_wp*uyy)

      uxpy = ux + uy
      feq(5) = wd*rho*(indp + 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      feq(7) = wd*rho*(indp - 3.0_wp*uxpy + 4.5_wp*uxpy*uxpy)
      
      uxmy = ux - uy
      feq(6) = wd*rho*(indp - 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)
      feq(8) = wd*rho*(indp + 3.0_wp*uxmy + 4.5_wp*uxmy*uxmy)

   end function


   subroutine alloc_grid(grid,nx,ny,nf,log)
      type(lattice_grid), intent(out) :: grid
      
      integer, intent(in) :: nx, ny
      integer, intent(in), optional :: nf
      logical, intent(in), optional :: log

      integer :: nf_
      logical :: log_
      character(len=:), allocatable :: logfile_
      integer :: ny_, rmd

      grid%nx = nx
      grid%ny = ny

      ! Round to closest dimension of 16
      ny_ = ny
      rmd = mod(ny_, 16)      
      if (rmd > 0) ny_ = ny_ + (16 - rmd)

      ! Number of pdf fields
      nf_ = 2
      if (present(nf)) nf_ = nf

      ! PDF memory
      allocate(grid%f(ny_,nx,0:8,nf_))

      ! Macroscopic fields
      allocate(grid%rho(ny,nx))
      allocate(grid%ux(ny,nx))
      allocate(grid%uy(ny,nx))

      !
      ! Initialize field pointers
      !
      grid%inew = 1
      grid%iold = 2
      if (nf_ > 2) then
         grid%imid = 3
      else
         grid%imid = -1
      end if

      !
      ! Initialize logging file
      !
      log_ = .true.
      if (present(log)) log_ = log

      if (log_) then
         logfile_ = "lattice_grid_log.txt"
         if (allocated(grid%logfile)) then
            logfile_ = grid%logfile
         end if
         open(newunit=grid%logunit,file=logfile_,status='unknown')
      end if

   end subroutine

   subroutine dealloc_grid(grid)
      type(lattice_grid), intent(inout) :: grid
      
      logical :: isopen

      !
      ! close log file
      !
      inquire(grid%logunit,opened=isopen)
      if (isopen) then
         close(grid%logunit)
      end if

      !
      ! deallocate allocatable storage
      !
      call deallocate_all_(grid)

   contains
      
      !> Deallocate all allocatable objects by applying intent(out)
      subroutine deallocate_all_(grid)
         type(lattice_grid), intent(out) :: grid
         ! Use an associate to prevent spurious warnings
         associate(nx => grid%nx)
            return
         end associate
      end subroutine
   
   end subroutine

   subroutine set_properties(grid, nu, dt, magic)
      type(lattice_grid), intent(inout) :: grid
      real(wp), intent(in) :: nu, dt
      real(wp), optional :: magic

      real(wp) :: tau

      grid%nu = nu
      grid%dt = dt

      ! TODO: lattice dependent logic 
      grid%csqr = csqr
      tau = invcsqr*nu

      grid%tau = tau
      grid%omega = dt/(tau + 0.5_wp*dt)

      if (present(magic)) then
         grid%trt_magic = magic
      else
         ! bgk by default
         !grid%trt_magic = (2.0_wp - grid%omega)**2 / (4.0_wp * (grid%omega)**2)
         grid%trt_magic = (tau/dt)**2
      end if

      print *, "trt magic = ", grid%trt_magic

   end subroutine


   subroutine set_pdf_to_equilibrium(grid)
      type(lattice_grid), intent(inout) :: grid

      real(wp) :: rho_, ux_, uy_
      integer :: x, y

      associate( nx => grid%nx, &
                 ny => grid%ny, &
                 rho => grid%rho, &
                 ux => grid%ux, &
                 uy => grid%uy, &
                 pdf => grid%f(:,:,:,grid%iold))

      do x = 1, nx
         do y = 1, ny

            rho_ = rho(y,x)
             ux_ =  ux(y,x)
             uy_ =  uy(y,x)

            pdf(y,x,:) = equilibrium(rho_, ux_, uy_)

         end do
      end do

      end associate

!      if (size(grid%f,4) > 2) then
!         grid%f(:,:,:,grid%imid) = grid%f(:,:,:,grid%iold)
!         grid%f(:,:,:,grid%inew) = grid%f(:,:,:,grid%iold)
!         call grid%collision()
!      end if

   end subroutine

   subroutine perform_step(grid)
      type(lattice_grid), intent(inout) :: grid

      call grid%streaming()  ! write from iold to inew
      call grid%collision()  ! update inew in place

      swap: block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp         
      end block swap

   end subroutine

   subroutine update_macros(grid)
      type(lattice_grid), intent(inout) :: grid

      integer :: ld

      ld = size(grid%f,1)

      call update_macros_kernel(grid%nx, grid%ny, ld, &
         grid%f(:,:,:,grid%inew), &
         grid%rho, grid%ux, grid%uy)

   contains

      subroutine update_macros_kernel(nx,ny,ld,f,grho,gux,guy)
         integer, intent(in) :: nx, ny, ld
         real(wp), intent(in) :: f(ld, nx, 0:8)
         real(wp), intent(inout), dimension(ny,nx) :: grho, gux, guy

         real(wp) :: rho, invrho, ux, uy, fs(0:8)
         integer :: x, y

         !$omp parallel do collapse(2) default(private) shared(nx,ny,f,grho,gux,guy)
         do x = 1, nx
            do y = 1, ny

               fs = f(y,x,:)

               ! density
               rho = fs(0) + (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                      ((fs(1) + fs(3)) + (fs(2) + fs(4)))) 

               grho(y,x) = rho

               ! velocity
               invrho = 1.0_wp/rho
               ux = invrho * (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3)))
               uy = invrho * (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4)))

               gux(y,x) = ux
               guy(y,x) = uy

            end do
         end do
         !$omp end parallel do

      end subroutine

   end subroutine


   subroutine output_gnuplot(grid, step)
      type(lattice_grid), intent(in) :: grid
      integer, intent(in), optional :: step

      character(len=64) :: istr
      character(len=:), allocatable :: fullname
      integer :: istat

      istr = ''
      if (present(step)) then
         write(istr,'(I0.9)') step
      end if

      fullname = ''
      if (allocated(grid%foldername)) then
         call execute_command_line("mkdir -p "//trim(grid%foldername), &
            exitstat=istat, wait=.true.)

         if (istat /= 0) then
            write(*,'(A)') "[output_grid] error making directory "//grid%foldername
            error stop
         end if
         fullname = trim(grid%foldername)
      end if

      fullname = fullname//'/'//grid%filename//trim(istr)//'.txt'
      
      call output_gnuplot_grid(fullname, &
         grid%nx,grid%ny, &
         grid%rho,grid%ux,grid%uy)

   end subroutine

   subroutine set_output_folder(grid,foldername,verbose)
      class(lattice_grid), intent(inout) :: grid
      character(len=*), intent(in) :: foldername
      logical, intent(in), optional :: verbose

      integer :: istat
      character(len=:), allocatable :: mkdir_opts

      mkdir_opts = '-p'
      if (present(verbose)) then
         if (verbose) then
            mkdir_opts = mkdir_opts // 'v'
         end if
      end if

      call execute_command_line( &
         'mkdir ' // mkdir_opts // ' ' // trim(foldername), &
         exitstat=istat, wait=.true.)

      if (istat /= 0) then
         write(*,'(A)') "[set_output_folder] error making directory "//grid%foldername
         error stop
      end if

      grid%foldername = foldername

   end subroutine

end module
