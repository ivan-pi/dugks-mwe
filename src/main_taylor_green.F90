program main_taylor_green

   use precision, only: wp, dp, walltime, walltime_remaining
   use lattice, only: lattice_grid, &
      set_pdf_to_equilibrium, &
      alloc_grid, &
      update_macros, &
      output_gnuplot, &
      set_properties, &
      dealloc_grid, &
      perform_step, &
      report_bw

   use taylor_green, only: taylor_green_t, pi
!   use collision_bgk, only: collide_bgk
!   use periodic_lbm, only: perform_lbm_step, lbm_stream
   use periodic_dugks, only: dugks_collide, dugks_stream, dugks_bc
   use periodic_lw, only: lw_collide, lw_stream, lw_bc

!$ use omp_lib

! amdflang: F90-S-0034-Syntax error
!   implicit none (type,external)

   implicit none

   integer, parameter :: nprint = 10000
   integer :: nx, ny

   integer :: step, nsteps

   type(taylor_green_t) :: tg

   type(lattice_grid) :: grid
   real(wp) :: cfl, dt, nu, tau
   real(wp) :: kx, ky, umax, nrm
   real(wp), parameter :: cs = 1.0_wp/sqrt(3.0_wp)

   real(wp) :: t, tmax

   real(wp) :: dt_over_tau
   character(len=64) :: arg

   real(dp) :: sbegin, send
   integer :: cases(6), ic

   logical :: with_omp
   integer :: maxthr

   write(*,'(A)') "=== TAYLOR-GREEN-VORTEX ==="

   ! Print configuraiton
   write(*,'("Work precision: ",I0)') wp

#ifdef DUGKS
   write(*,'(A)') "Algorithm: DUGKS"
#else
   write(*,'(A)') "Algorithm: Lax-Wendroff streaming"
#endif

   block
      use lattice, only: shift
      write(*,'("Distribution shifting: ", L1)') shift
   end block

   with_omp = .false.
   maxthr = 1
!$ with_omp = .true.
!$ maxthr = omp_get_max_threads()

   write(*,'("OpenMP: ", L1)') with_omp
   write(*,'("Max. threads: ",I0)') maxthr

!cases = [32,64,96,128,160,192]
!do ic = 1, size(cases)

!   nx = cases(ic)
!   ny = nx
   call read_env(nx, ny)
   print *, "nx, ny = ", nx, ny

   call alloc_grid(grid, nx, ny)

   grid%filename = "results"
   call grid%set_output_folder(foldername="taylor_green")

#ifdef DUGKS
   grid%collision => dugks_collide
   grid%streaming => dugks_stream
   grid%bc => dugks_bc
#else
   grid%collision => lw_collide
   grid%streaming => lw_stream
   grid%bc => lw_bc
#endif

   grid%logger => my_logger

   ! umax = Mach * cs
   ! nu = (umax * L) / Re

   select case("kraemer")
   case("guo")
      ! Guo (2013)
      ! WARNING: some parameters are uncertain in this case
      umax = 0.1_wp / sqrt(3._wp)
      nu = (umax * real(nx,wp)) / 1000._wp
   case("wu")
      ! Wu et al. (2018)
      umax = 0.015811388300841896 / sqrt(3.0_wp)
      nu = (umax * real(nx,wp)) / 100._wp
   case("kraemer")
      !  Kraemer et al. (2017)
      umax = 0.008 / sqrt(3._wp)
      nu = (umax * real(nx,wp)) / 10.0_wp
   case default
      ! Zhu (2017) settings; these are REPRODUCIBLE!
      umax = 0.01_wp / sqrt(3._wp)
      nu = (umax * real(nx,wp)) / 100._wp
   end select

   ! tau = nu / cs**2, cs**2 = 1/3
   tau = 3.0_wp * nu 

   ! read value for dt/tau or cfl
   call get_command_argument(1, arg)

   read(arg,*) cfl
!   cfl = 0.5_Wp
!   dt = 0.00554256
!   cfl = 0.1_wp
   dt = cfl / sqrt(2.0_wp)
   dt_over_tau = dt / tau

!   read(arg,*) dt_over_tau
!   dt_over_tau = 4.0_wp / 3.0_wp
!   dt = dt_over_tau*tau
!   cfl = sqrt(2._wp)*dt
   
   print *, "Ma     = ", umax / cs
   print *, "Re     = ", umax * max(nx, ny) / nu
   print *, "tau    = ", tau
   print *, "dt     = ", dt
   print *, "cfl    = ", cfl
   print *, "dt/tau = ", dt / tau

   if (cfl > 1.0) then
      write(*,'(A)') "Error: CFL > 1, exiting."
      stop 1
   end if

   call set_properties(grid, nu, dt, magic=1._wp/4._wp)

   print *, "omega = ", grid%omega

   ! ---- prepare flow case ----

   kx = 2*pi/real(nx,wp)
   ky = 2*pi/real(ny,wp)

   tg = taylor_green_t(nx,ny,kx,ky,umax,nu)
   print *, "umax = ", umax
   print *, "tc   = ", tg%td
   call write_gnuplot_include()

   tmax = log(10._wp)*tg%decay_time()
!   tmax = 0.7895683520871487 * tg%decay_time() ! Wu et al. 2018, or 10^5 timesteps
   nsteps = ceiling(tmax / dt)

   nsteps = 1000000
!   tmax = nsteps *dt

!   tmax =
!   nsteps = ceiling(tmax/dt)

   print *, "tmax = ", tmax
   print *, "nsteps = ", nsteps

   t = 0._wp
   call apply_initial_condition(tg, grid)

   call output_gnuplot(grid,step=0)
!   call grid%logger(step=0)

   sbegin = walltime()

   write(*,'(/,A)') "Starting time-stepping ... "

   time: do step = 1, nsteps

      call perform_step(grid)
      t = t + dt

      if (mod(step,nprint) == 0) then
         ! --- Output ---
         call update_macros(grid)
         write(*,*) "step = ", step, ", max(|u|) = ", &
            maxval(hypot(grid%ux,grid%uy)), calc_L2_norm(tg,grid,t), &
            ", remaining: ", walltime_remaining(step,nsteps,sbegin)

!         call output_gnuplot(grid,step)
!         call grid%logger(step)
      end if

      if (t >= tmax) then
         ! --- Exit Timeloop ---
         call update_macros(grid)
         call output_gnuplot(grid,step)
!         call grid%logger(step)
         exit time
      end if

   end do time

   send = walltime()

   print *, "MLUPS ", (real(nx,wp) * real(ny,wp) * real(step,wp) * 1.e-6_wp) / (send - sbegin)

   block
      real(wp) :: elapsed_per_step, stream_bw, collision_bw
      elapsed_per_step = (send - sbegin) / step

      stream_bw = 2 * product(real([nx,ny,9,wp],wp))
      collision_bw = 2 * product(real([nx,ny,9,wp],wp)) + product(real([nx,ny,3],wp))

      print *, "Eff. BW (global) = ", (stream_bw + collision_bw) / elapsed_per_step / 1024.0_wp**3

   end block
   call report_bw(grid)


   ! calculate average L2-norm
   call update_macros(grid)
   nrm = calc_L2_norm(tg, grid, t)

   print *, "L2-norm = ", nrm
   print *, "Final time = ", t
   print *, "Final u/umax = ", maxval(hypot(grid%ux,grid%uy)) / umax

   call dealloc_grid(grid)

!end do


contains

   subroutine read_env(nx,ny)
      integer, intent(out) :: nx, ny

      character(len=32) :: value
      integer :: status, n

      CALL GET_ENVIRONMENT_VARIABLE('N', VALUE, status=STATUS)

      if (status == 1) then
         ! Variable does not exist, use default values
         nx = 100
         ny = 100
      else if (status == 2) then
         write(*,*) "Error: processor doesn't support environment variables"
         stop
      else if (status == -1) then
         write(*,*) "Error: value is not long enough to hold variable"
         stop
      else
         ! Everything should be okay
         read(value,*) n
         nx = n
         ny = n
      end if

   end subroutine

   subroutine apply_initial_condition(case, grid)
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(inout) :: grid

      real(wp), parameter :: rho0 = 1.0_wp

      call case%eval(t=0.0_wp, &
                   p=grid%rho, &
                   ux=grid%ux, &
                   uy=grid%uy)

      ! convert pressure to lattice density
      grid%rho = grid%rho/grid%csqr + rho0

      call set_pdf_to_equilibrium(grid)
      call grid%bc() ! periodic bc

   end subroutine

   subroutine my_logger(grid,step)
      class(lattice_grid), intent(in) :: grid
      integer, intent(in) :: step

      write(grid%logunit, *) step, step*dt, maxval(hypot(grid%ux,grid%uy))
      flush(grid%logunit)

   end subroutine

   subroutine write_gnuplot_include()

      integer :: unit

      open(newunit=unit,file="lattice_grid_log.incl",status='unknown')

      write(unit,*) "dt = ", dt
      write(unit,*) "umax = ", umax
      write(unit,*) "tc = ", tg%td

      close(unit)
   end subroutine


   function calc_L2_norm(case, grid, t) result(nrm)
      ! TODO: move this function to the taylor_green module
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(in) :: grid
      real(wp), intent(in) :: t
      real(wp) :: nrm

      real(wp), allocatable :: pa(:,:), uxa(:,:), uya(:,:)
      real(wp) :: above, below

      allocate(pa,  mold=grid%rho) ! not actually used
      allocate(uxa, mold=grid%ux)
      allocate(uya, mold=grid%uy)

      call case%eval(t=t, &
                   p=pa, &
                   ux=uxa, &
                   uy=uya)

!      above = norm2(hypot(grid%ux - uxa, grid%uy - uya))
!      below = norm2(hypot(uxa, uya))
!      nrm = above/below

      ! For Kraemer test
!      nrm = above / umax

      ! # Code used in the Python version
      ! def L2_error(u,v,ua,va):
      !     return np.sqrt(np.sum((u-ua)**2 + (v-va)**2)/np.sum(ua**2 + va**2))

      above = sum((grid%ux - uxa)**2 + (grid%uy - uya)**2)
      !below = sum(uxa**2 + uya**2)
      !nrm = sqrt(above/below)
      nrm = sqrt(above) / umax

   end function


   real(wp) function nrm2(n,a)
      integer, intent(in) :: n
      real(wp), intent(in) :: a(n)

      integer :: tid, nthr, chunk, lo, hi

      tid = 0
      nthr = 1

      nrm = 0.0_wp

      !$omp parallel default(private) shared(n,a) reduction(+:nrm)
      !$ tid = omp_get_thread_num()
      !$ nthr = omp_get_num_threads()
      chunk = (n + nthr - 1)/nthr
      lo = tid*chunk + 1
      hi = min((tid+1)*chunk,n)
      nrm = dot_product(a(lo:hi),a(lo:hi))
      !$omp end parallel

      nrm = sqrt(nrm)
   end function


end program