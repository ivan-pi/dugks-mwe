program main_taylor_green

   use precision, only: wp
   use lattice, only: lattice_grid, &
      set_pdf_to_equilibrium, &
      alloc_grid, &
      update_macros, &
      output_gnuplot, &
      set_properties, &
      dealloc_grid

   use taylor_green, only: taylor_green_t, pi
!   use collision_bgk, only: collide_bgk
!   use periodic_lbm, only: perform_lbm_step, lbm_stream
   use periodic_dugks, only: perform_step, dugks_collide, dugks_stream

!$ use omp_lib

! amdflang: F90-S-0034-Syntax error
!   implicit none (type,external)

   implicit none

   integer, parameter :: nprint = 20000
   integer :: nx, ny

   integer :: step, nsteps

   type(taylor_green_t) :: tg

   type(lattice_grid) :: grid
   real(wp) :: cfl, dt, nu, tau
   real(wp) :: kx, ky, umax, nrm

   real(wp) :: t, tmax

   real(wp) :: dt_over_tau
   character(len=64) :: arg

   integer, parameter :: dp = kind(1.0d0)
!$ real(dp) :: sbegin, send

!$ print *, "--- In OpenMP mode ---"
!$ print *, "Maximimum number of threads: ", omp_get_max_threads()

   
   call read_env(nx,ny)
   print *, "nx, ny = ", nx, ny

   call alloc_grid(grid, nx, ny, nf = 2)

   grid%filename = "results"
   call grid%set_output_folder(foldername="taylor_green")

   grid%collision => dugks_collide
   grid%streaming => dugks_stream

   grid%logger => my_logger

   ! umax = Mach * cs
   umax = 0.01_wp / sqrt(3._wp)

   ! nu = (umax * L) / Re
   nu = (umax * real(nx,wp)) / 100._wp
   
   ! tau = nu / cs**2, cs**2 = 1/3
   tau = 3.0_wp * nu 

   ! read value for dt/tau or cfl
   call get_command_argument(1,arg)
   
   read(arg,*) cfl
   dt = cfl
   dt_over_tau = dt/tau

   !read(arg,*) dt_over_tau
   !dt = dt_over_tau*tau
   !cfl = sqrt(2._wp)*dt
   
   print *, "tau = ", tau
   print *, "dt/tau = ", dt/tau
   print *, "cfl = ", cfl

   call set_properties(grid, nu, dt, magic=1._wp/4._wp)

   print *, "omega = ", grid%omega

   ! ---- prepare flow case ----

   kx = 2*pi/real(nx,wp)
   ky = 2*pi/real(ny,wp)

   tg = taylor_green_t(nx,ny,kx,ky,umax,nu)
   print *, "umax = ", umax
   print *, "tc   = ", tg%td
   call write_gnuplot_include()

   tmax = log(2._wp)*tg%decay_time()
   nsteps = int(1.1_wp*tmax/dt)
   !nsteps = 5000

   print*, "nsteps = ", nsteps

   t = 0._wp
   call apply_initial_condition(tg, grid)

   call output_gnuplot(grid,step=0)
!   call grid%logger(step=0)

!$ sbegin = omp_get_wtime()

   time: do step = 1, nsteps

      call perform_step(grid)
      t = t + dt

      if (mod(step,nprint) == 0) then
         ! --- Output ---
         call update_macros(grid)
         
         print '(A,G0,2X,A,G0)', "step = ", step, ", max(|u|) = ", &
            maxval(hypot(grid%ux,grid%uy))

         call output_gnuplot(grid,step)
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

!$ send = omp_get_wtime()
!$ print *, "MLUPS ", (real(nx,wp) * real(ny,wp) * real(step,wp) * 1.e-6_wp) / (send - sbegin)

   ! calculate average L2-norm
   nrm = calc_L2_norm(tg, grid, t)

   print *, "L2-norm = ", nrm
   print *, "Final time = ", t

   call dealloc_grid(grid)

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
      type(taylor_green_t), intent(in) :: case
      type(lattice_grid), intent(in) :: grid
      real(wp), intent(in) :: t
      real(wp) :: nrm

      real(wp), allocatable :: pa(:,:), uxa(:,:), uya(:,:)
      real(wp) :: above, below

      allocate(pa,  mold=grid%rho) ! not needed
      allocate(uxa, mold=grid%ux)
      allocate(uya, mold=grid%uy)

      call case%eval(t=t, &
                   p=pa, &
                   ux=uxa, &
                   uy=uya)

      associate(nx => grid%nx, ny => grid%ny)


      above = norm2(hypot(&
         grid%ux(1:ny,1:nx) - uxa(1:ny,1:nx), &
         grid%uy(1:ny,1:nx) - uya(1:ny,1:nx)))
      below = norm2(hypot(uxa(1:ny,1:nx), uya(1:ny,1:nx)))

      end associate

      nrm = above/below

      ! # Code used in the Python version
      ! def L2_error(u,v,ua,va):
      !     return np.sqrt(np.sum((u-ua)**2 + (v-va)**2)/np.sum(ua**2 + va**2))

      !above = sum((grid%ux - uxa)**2 + (grid%uy - uya)**2)
      !below = sum(uxa**2 + uya**2)
      !nrm = sqrt(above/below)

   end function

end program