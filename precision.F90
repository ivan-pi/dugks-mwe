module precision
   implicit none
   private
   public :: sp, dp, wp

   public :: walltime
   public :: walltime_remaining

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
#if WITH_SP
   integer, parameter :: wp = sp
#else
   integer, parameter :: wp = dp
#endif

contains

   real(dp) function walltime()
   !$ use omp_lib, only: omp_get_wtime

   !$ if (.true.) then
   !$    walltime = omp_get_wtime()
   !$ else
         call cpu_time(walltime)
   !$ endif

   end function


   function walltime_remaining(current_step, total_steps, begin) result(t_rem)
        integer, intent(in) :: current_step, total_steps
        real(dp), intent(in) :: begin

        real(dp) :: t_rem, elapsed

        elapsed   = walltime() - begin

        if (current_step > 0) then
            ! (Time per step) * (Remaining steps)
            t_rem = (elapsed / real(current_step,kind=dp)) * &
               real(total_steps - current_step,kind=dp)
        else
            t_rem = 0.0_dp
        end if
    end function

end module
