module lattice

   use precision, only: wp, dp, walltime
   use gnuplot_io, only: output_gnuplot_grid

   implicit none
   private

   public :: wp
   public :: lattice_grid
   public :: perform_step

   public :: alloc_grid, dealloc_grid

   public :: set_properties
   public :: report_bw
   public :: total_bw

   public :: update_macros
   
   public :: output_gnuplot

   public :: set_pdf_to_equilibrium

   public :: cx, cy, csqr, w0, ws, wd, w, rho0, shift

   character(*), parameter :: FMT_REAL_SP = '(es15.8e2)'

   type, abstract :: lattice_grid
      integer :: nx, ny

      real(wp), allocatable :: f(:,:,:,:)

      real(wp), allocatable :: rho(:,:), wrk(:,:)
      real(wp), allocatable ::  ux(:,:)
      real(wp), allocatable ::  uy(:,:)
      
      real(wp) :: nu, dt, tau
      real(wp) :: omega, trt_magic
      real(wp) :: csqr
      
      integer :: iold, inew, imid
      

      character(len=:), allocatable :: filename, foldername

      character(len=:), allocatable :: logfile
      procedure(gridlog_interface), pointer, pass(grid) :: logger => null()
      integer :: logunit

      real(dp) :: ts = 0, tc = 0, tb = 0
      integer :: total_steps = 0
   contains
      procedure(collision_interface), pass(grid), deferred :: collision
      procedure(streaming_interface), pass(grid), deferred :: streaming
      procedure(bc_interface),        pass(grid), deferred :: bc
      procedure :: set_output_folder
      procedure :: alloc => alloc_grid
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
      subroutine bc_interface(grid)
         import lattice_grid
         class(lattice_grid), intent(inout) :: grid
      end subroutine
      subroutine gridlog_interface(grid, step)
         import lattice_grid
         class(lattice_grid), intent(in) :: grid
         integer, intent(in) :: step
      end subroutine
   end interface

   real(wp), parameter :: cx(0:8) = [0, 1, 0, -1, 0, 1, -1, -1, 1]
   real(wp), parameter :: cy(0:8) = [0, 0, 1, 0, -1, 1, 1, -1, -1]

   real(wp), parameter :: w0 = 4.0_wp / 9.0_wp, &
                          ws = 1.0_wp / 9.0_wp, &
                          wd = 1.0_wp / 36.0_wp
   real(wp), parameter :: w(0:8) = [w0,ws,ws,ws,ws,wd,wd,wd,wd]

   real(wp), parameter :: csqr = 1.0_wp/3.0_wp
   real(wp), parameter :: invcsqr = 3.0_wp

#ifdef DUGKS
   ! DUGKS implementation does not support distribution shifting
   logical, parameter :: shift = .false.
#else
   logical, parameter :: shift = .false.
#endif
   real(wp), parameter :: rho0 = 1.0_wp

   integer, parameter :: nhalo = 1
contains

   subroutine alloc_grid(grid,nx,ny,nu,dt,log)
      class(lattice_grid), intent(out) :: grid
      
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: nu, dt
      logical, intent(in), optional :: log

      integer, parameter :: nf = 2
      logical :: log_
      character(len=:), allocatable :: logfile_

      grid%nx = nx
      grid%ny = ny

      ! PDF memory
#ifdef DUGKS
      allocate(grid%f(ny,nx,0:8,nf))
#else
      allocate(grid%f(1-nhalo:ny+nhalo,1-nhalo:nx+nhalo,0:8,2))
#endif

      ! Macroscopic fields
      allocate(grid%rho(ny,nx))
      allocate(grid%wrk(ny,nx))
      allocate(grid%ux(ny,nx))
      allocate(grid%uy(ny,nx))

      block
         real(wp) :: macro_mem, pdf_mem

         macro_mem = 4 * product(shape(grid%rho,kind=wp)) * wp
         pdf_mem = product(shape(grid%f,kind=wp)) * wp

         write(*,'(A,G0)') "Mem. macros (MB) = ", macro_mem / 1024.0_wp**2
         write(*,'(A,G0)') "Mem. pdf (MB)    = ", pdf_mem / 1024.0_wp**2
      end block

      !
      ! Initialize field pointers
      !
      grid%inew = 1
      grid%iold = 2

      !
      ! Initialize logging file
      !
      if (present(log)) then
         if (log) then
            logfile_ = "lattice_grid_log.txt"
            if (allocated(grid%logfile)) then
               logfile_ = grid%logfile
            end if
            open(newunit=grid%logunit,file=logfile_,status='unknown')
         end if
      end if

      !
      ! Initialize fluid settings
      !
      call set_properties(grid,nu,dt,magic=0.25_wp)

   end subroutine

   subroutine dealloc_grid(grid)
      class(lattice_grid), intent(inout) :: grid
      
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
         class(lattice_grid), intent(out) :: grid
         ! Use an associate to prevent spurious warnings
         associate(nx => grid%nx)
            return
         end associate
      end subroutine
   
   end subroutine

   subroutine set_properties(grid, nu, dt, magic)
      class(lattice_grid), intent(inout) :: grid
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
      class(lattice_grid), intent(inout) :: grid

      real(wp) :: rho_, ux_, uy_
      integer :: x, y

      do x = 1, grid%nx
         do y = 1, grid%ny

            rho_ = grid%rho(y,x)
             ux_ = grid%ux(y,x)
             uy_ = grid%uy(y,x)

            grid%f(y,x,:,grid%iold) = &
               equilibrium(grid%rho(y,x), grid%ux(y,x), grid%uy(y,x))

         end do
      end do

   contains

      pure function equilibrium(rho,ux,uy) result(feq)
         real(wp), intent(in) :: rho, ux, uy
         real(wp) :: feq(0:8)

         real(wp) :: uxx, uyy, uxpy, uxmy
         real(wp) :: indp

         uxx = ux*ux
         uyy = uy*uy

         if (shift) then
            indp = -1.5_wp * (uxx + uyy)
         else
            indp = 1.0_wp - 1.5_wp * (uxx + uyy)
         endif

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

         if (shift) feq = feq + w*(rho - rho0)

      end function

   end subroutine


#if 1
   subroutine perform_step(grid)
      class(lattice_grid), intent(inout) :: grid

      real(dp) :: t(4)

      t(1) = walltime()
      call grid%collision()  ! update iold in place
      t(2) = walltime()
      call grid%bc()
      t(3) = walltime()
      call grid%streaming()
      t(4) = walltime()

      grid%tc = grid%tc + (t(2) - t(1))
      grid%tb = grid%tb + (t(3) - t(2))
      grid%ts = grid%ts + (t(4) - t(3))
      grid%total_steps = grid%total_steps + 1

      block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp
      end block

   end subroutine
#else
   subroutine perform_step(grid)
      class(lattice_grid), intent(inout) :: grid

      call grid%streaming()  ! write from iold to inew
      call grid%collision()  ! update inew in place

      swap: block
         integer :: itmp
         itmp = grid%iold
         grid%iold = grid%inew
         grid%inew = itmp         
      end block swap

   end subroutine
#endif

   subroutine update_macros(grid)
      class(lattice_grid), intent(inout) :: grid

      call update_macros_kernel(grid%nx, grid%ny, &
         grid%f(:,:,:,grid%iold), &
         grid%rho, grid%ux, grid%uy)

   contains

      subroutine update_macros_kernel(nx,ny,f,rho,ux,uy)
         integer, intent(in) :: nx, ny
#ifdef DUGKS
         real(wp), intent(in), contiguous :: f(:,:,0:)
#else
         real(wp), intent(in) :: f(0:ny+1,0:nx+1,0:8)
#endif
         real(wp), intent(out), dimension(ny,nx) :: rho, ux, uy

         real(wp) :: invrho, fs(0:8)
         integer :: x, y

         !$omp parallel do default(private) shared(nx,ny,f,rho,ux,uy)
         do x = 1, nx
            do y = 1, ny

               fs = f(y,x,:)

               ! density
               rho(y,x) = fs(0) + (((fs(5) + fs(7)) + (fs(6) + fs(8))) + &
                      ((fs(1) + fs(3)) + (fs(2) + fs(4))))

               if (shift) rho(y,x) = rho(y,x) + rho0

               ! velocity
               ux(y,x) = (((fs(5) - fs(7)) + (fs(8) - fs(6))) + (fs(1) - fs(3))) / rho(y,x)
               uy(y,x) = (((fs(5) - fs(7)) + (fs(6) - fs(8))) + (fs(2) - fs(4))) / rho(y,x)

            end do
         end do

      end subroutine

   end subroutine


   subroutine output_gnuplot(grid, step)
      class(lattice_grid), intent(in) :: grid
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


! TODO: generalize bandwidth computations for DUGKS scheme

   real(wp) function total_bw(nx,ny,elapsed_per_step) ! in GB/s
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: elapsed_per_step

      real(wp) :: stream_bw, collision_bw, macro_bw

      stream_bw = 2 * product(real([nx,ny,9,wp],wp))
      collision_bw = 2 * product(real([nx,ny,9,wp],wp))
      macro_bw = product(real([nx,ny,3],wp))

      total_bw = (stream_bw + collision_bw + macro_bw) / &
         elapsed_per_step / 1024.0_wp**3

   end function

   subroutine report_bw(grid)
      class(lattice_grid), intent(in) :: grid

      real(wp) :: stream_bw, collision_bw
      integer :: sz

      sz = wp
      stream_bw = 2 * product(real([grid%nx,grid%ny,9,sz],wp))

      ! PDFs
      collision_bw = 2 * product(real([grid%nx,grid%ny,9,sz],wp))
      collision_bw = collision_bw + product(real([grid%nx,grid%ny,4,sz],wp))

      write(*,'("Perfomance metrics:")')
      write(*,'("  Eff. BW (streaming) = ",G0.4," GB/s")') &
         stream_bw / (grid%ts/grid%total_steps) / 1024.0_wp**3
      write(*,'("  Eff. BW (collision) = ",G0.4," GB/s")') &
         collision_bw / (grid%tc/grid%total_steps) / 1024.0_wp**3

      associate(tt => grid%ts + grid%tc + grid%tb)
      write(*,'("  Streaming Time = ",G0.3," s (",F4.1," %)")') grid%ts, 100*grid%ts/tt
      write(*,'("  Collision Time = ",G0.3," s (",F4.1," %)")') grid%tc, 100*grid%tc/tt
      write(*,'("  Boundary Time  = ",G0.3," s (",F4.1," %)")') grid%tb, 100*grid%tb/tt
      end associate

   end subroutine

end module
