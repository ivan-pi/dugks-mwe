module fftw3
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
end module

program test_fftw
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none

    integer :: npad, seed_val, i, n_args
    character(len=32) :: arg, val_str
    logical :: custom_seed

    ! --- Default Values ---
    npad = 0
    custom_seed = .false.

! --- CLI Parsing ---
    n_args = command_argument_count()
    i = 1
    do while (i <= n_args)
        call get_command_argument(i, arg)
        select case (trim(arg))
        case ("--help", "-h")
            call print_help()
            stop
        case ("--npad", "-p")
            if (i < n_args) then
                call get_command_argument(i+1, val_str)
                read(val_str, *) npad
                i = i + 1
            end if
        case ("--seed", "-s")
            if (i < n_args) then
                call get_command_argument(i+1, val_str)
                read(val_str, *) seed_val
                custom_seed = .true.
                i = i + 1
            end if
        case default
            print *, "Unknown argument: ", trim(arg)
            print *, "Use -h or --help for usage instructions."
            stop
        end select
        i = i + 1
    end do

    ! Seed random number generator
    if (custom_seed) call set_reproducible_seed(seed_val)

    print *, "Starting FFTW R2C Interface Tests..."
    print "(A, I2)", " Padding (npad): ", npad
    print "(A, I12)", " Random Seed:    ", seed_val
    print *, "------------------------------------"

    call test_1d_differentiation(npad)
    print *, "------------------------------------"
    call test_1d_periodic_fe_solver(npad)
    print *, "------------------------------------"
    call test_1d_cyclic_shift(npad)
    print *, "------------------------------------"
    call test_1d_parseval_energy(npad)
    print *, "------------------------------------"

    block
        integer, parameter :: im = 64, jm = 64, lm = 10
        call test_fftw_2d_single(im, jm, lm, npad)
        print "(' ------------------------------------')"
        call test_fftw_2d_many(im, jm, lm, npad)
        print "(' ------------------------------------')"
    end block

    print *, "Tests Completed."

contains

    subroutine print_help()
        print *, "FFTW R2C many_dft Test Suite"
        print *, "Usage: ./test_fftw [OPTIONS]"
        print *, ""
        print *, "Options:"
        print *, "  -h, --help      Show this help message and exit"
        print *, "  -p, --npad INT  Set number of ghost/padding cells (default: 0)"
        print *, "                  This tests FFTW's 'nembed' and stride logic."
        print *, "  -s, --seed INT  Set fixed random seed for reproducible tests."
        print *, "                  If omitted, system clock is used (non-deterministic)."
        print *, ""
        print *, "Note: All tests perform a Plan -> Initialize -> Execute cycle"
        print *, "      to ensure data integrity during FFTW planning."
    end subroutine print_help

    subroutine set_reproducible_seed(s_val)
        integer, intent(in) :: s_val
        integer :: n
        integer, allocatable :: seed_arr(:)
        call random_seed(size = n)
        allocate(seed_arr(n))
        seed_arr = s_val
        call random_seed(put = seed_arr)
    end subroutine set_reproducible_seed

    subroutine test_1d_differentiation(pad)
        integer, intent(in) :: pad
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: phys_in(:), phys_out(:), phys_ref(:)
        complex(DCP), allocatable :: spec_data(:)
        type(C_PTR) :: plan_fwd, plan_bwd
        integer :: i, n_size
        real(DP) :: L, dx, k_val, err
        real(DP), parameter :: pi = 4.0_DP * atan(1.0_DP)

        n_size = 16
        n(1) = n_size
        inembed(1) = n_size + 2*pad
        allocate(phys_in(1-pad:n_size+pad), phys_out(1-pad:n_size+pad))
        allocate(phys_ref(n_size), spec_data(n_size/2 + 1))

        ! Planning first with dummy/empty arrays
        plan_fwd = fftw_plan_many_dft_r2c(1, n, 1, phys_in(1), inembed, 1, n(1), &
                    spec_data(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        plan_bwd = fftw_plan_many_dft_c2r(1, n, 1, spec_data(1), [n(1)/2+1], 1, n(1)/2+1, &
                    phys_out(1), inembed, 1, n(1), FFTW_ESTIMATE)

        ! Init
        L = 2.0_DP * pi; dx = L / n_size
        phys_in = -999.0_DP
        do i = 1, n_size
            phys_in(i)  = sin((i-1) * dx)
            phys_ref(i) = cos((i-1) * dx)
        end do

        call fftw_execute_dft_r2c(plan_fwd, phys_in(1), spec_data(1))
        do i = 1, n_size/2 + 1
            k_val = (2.0_DP * pi / L) * real(i-1, DP)
            spec_data(i) = spec_data(i) * (0.0_DP, 1.0_DP) * k_val
        end do
        spec_data(n_size/2 + 1) = cmplx(real(spec_data(n_size/2 + 1)), 0.0_DP, kind=DCP)
        call fftw_execute_dft_c2r(plan_bwd, spec_data(1), phys_out(1))

        phys_out(1:n_size) = phys_out(1:n_size) / real(n_size, DP)
        err = maxval(abs(phys_out(1:n_size) - phys_ref))
        print *, "Test 1: 1D R2C Differentiation"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        print *, "   Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)

        call fftw_destroy_plan(plan_fwd)
        call fftw_destroy_plan(plan_bwd)
    end subroutine test_1d_differentiation

    subroutine test_1d_periodic_fe_solver(pad)
        integer, intent(in) :: pad
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x_orig(:), y_phys(:), x_final(:), row_phys(:)
        complex(DCP), allocatable :: y_spec(:), a_spec(:)
        type(C_PTR) :: fwd_p, bwd_p
        integer :: i, n_size
        real(DP) :: err

        n_size = 16
        n(1) = n_size
        inembed(1) = n_size + 2*pad
        allocate(x_orig(1-pad:n_size+pad), y_phys(1-pad:n_size+pad), &
                 x_final(1-pad:n_size+pad), row_phys(1-pad:n_size+pad))
        allocate(y_spec(n_size/2+1), a_spec(n_size/2+1))

        fwd_p = fftw_plan_many_dft_r2c(1, n, 1, y_phys(1), inembed, 1, n(1), &
                                      y_spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        bwd_p = fftw_plan_many_dft_c2r(1, n, 1, y_spec(1), [n(1)/2+1], 1, n(1)/2+1, &
                                      x_final(1), inembed, 1, n(1), FFTW_ESTIMATE)

        x_orig(1:n_size) = [(real(i, DP), i=1, n_size)]
        row_phys = 0.0_DP
        row_phys(1) = 4.0_DP/6.0_DP; row_phys(2) = 1.0_DP/6.0_DP; row_phys(n_size) = 1.0_DP/6.0_DP
        do i = 1, n_size
            y_phys(i) = (4.0_DP/6.0_DP) * x_orig(i) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i-2+n_size, n_size)+1) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i, n_size)+1)
        end do

        call fftw_execute_dft_r2c(fwd_p, row_phys(1), a_spec(1))
        call fftw_execute_dft_r2c(fwd_p, y_phys(1), y_spec(1))
        y_spec = y_spec / a_spec
        call fftw_execute_dft_c2r(bwd_p, y_spec(1), x_final(1))

        x_final(1:n_size) = x_final(1:n_size) / real(n_size, DP)
        err = maxval(abs(x_orig(1:n_size) - x_final(1:n_size)))
        print *, "Test 3: 1D Periodic FE Solver"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        print *, "   Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)

        call fftw_destroy_plan(fwd_p)
        call fftw_destroy_plan(bwd_p)
    end subroutine test_1d_periodic_fe_solver

    subroutine test_1d_cyclic_shift(pad)
        integer, intent(in) :: pad
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x_o(:), x_s(:), x_r(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: pf, pb
        integer :: i, n_size, shift = 5
        real(DP) :: arg, err, pi = 4.0_DP * atan(1.0_DP)

        n_size = 32
        n(1) = n_size
        inembed(1) = n_size + 2*pad
        allocate(x_o(1-pad:n_size+pad), x_s(1-pad:n_size+pad), &
                 x_r(1-pad:n_size+pad), spec(n_size/2+1))

        pf = fftw_plan_many_dft_r2c(1, n, 1, x_o(1), inembed, 1, n(1), spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        pb = fftw_plan_many_dft_c2r(1, n, 1, spec(1), [n(1)/2+1], 1, n(1)/2+1, x_r(1), inembed, 1, n(1), FFTW_ESTIMATE)

        do i = 1, n_size
            x_o(i) = exp(-((real(i-n_size/2, DP)/4.0_DP)**2))
            x_s(i) = exp(-((real(mod(i-1-shift+n_size, n_size)-n_size/2+1, DP)/4.0_DP)**2))
        end do

        call fftw_execute_dft_r2c(pf, x_o(1), spec(1))
        do i = 1, n_size/2 + 1
            arg = -2.0_DP * pi * real(i-1, DP) * real(shift, DP) / real(n_size, DP)
            spec(i) = spec(i) * cmplx(cos(arg), sin(arg), kind=DCP)
        end do
        call fftw_execute_dft_c2r(pb, spec(1), x_r(1))

        x_r(1:n_size) = x_r(1:n_size) / real(n_size, DP)
        err = maxval(abs(x_r(1:n_size) - x_s(1:n_size)))
        print *, "Test 4: 1D Cyclic Shift"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        print *, "   Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)

        call fftw_destroy_plan(pf)
        call fftw_destroy_plan(pb)
    end subroutine test_1d_cyclic_shift

    subroutine test_1d_parseval_energy(pad)
        integer, intent(in) :: pad
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: plan
        real(DP) :: e_phys, e_spec, err
        integer :: i, n_size

        n_size = 64
        n(1) = n_size
        inembed(1) = n_size + 2*pad
        allocate(x(1-pad : n_size+pad), spec(n_size/2 + 1))

        plan = fftw_plan_many_dft_r2c(1, n, 1, x(1), inembed, 1, n(1), spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)

        call random_number(x(1 : n_size))

        call fftw_execute_dft_r2c(plan, x(1), spec(1))
        e_phys = sum(x(1:n_size)**2)
        e_spec = abs(spec(1))**2
        do i = 2, n_size/2
            e_spec = e_spec + 2.0_DP * abs(spec(i))**2
        end do
        e_spec = (e_spec + abs(spec(n_size/2 + 1))**2) / real(n_size, DP)

        err = abs(e_phys - e_spec)
        print *, "Test 5: Parseval's Theorem"
        print '(A, E12.4)', "  Energy Difference: ", err
        print *, "   Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)

        call fftw_destroy_plan(plan)
    end subroutine test_1d_parseval_energy

    subroutine test_fftw_2d_single(im, jm, lm, pad)
        integer, intent(in) :: im, jm, lm, pad
        real(c_double), allocatable, target :: u(:,:,:), recov(:,:,:)
        complex(c_double_complex), allocatable, target :: cu(:,:)
        type(c_ptr) :: planf, planb
        integer :: k
        integer(c_int) :: n(2), nr(2), nc(2)
        real(c_double) :: norm, err, tolerance = 1.0e-12_c_double

        allocate(u(1-pad:im+pad, 1-pad:jm+pad, lm))
        allocate(recov(1-pad:im+pad, 1-pad:jm+pad, lm))
        allocate(cu(im/2 + 1, jm))

        ! Logic dimensions (Reversed for C-API)
        n = int([jm, im], c_int)

        ! Real Embedding: Physical dimensions of the padded plane
        ! nr(1) is the 'slow' dim, nr(2) is the 'fast' dim (with padding)
        nr = int([jm, im + 2*pad], c_int)

        ! Complex Embedding: Contiguous (no padding)
        nc = int([jm, im/2 + 1], c_int)

        ! Create plans using 'many' interface for a single (howmany=1) 2D transform
        ! This tells FFTW: "The row length is im+2*pad, but only transform im elements."
        planf = fftw_plan_many_dft_r2c(2, n, 1, u(1,1,1), nr, 1, 0, &
                                      cu(1,1), nc, 1, 0, FFTW_MEASURE)

        planb = fftw_plan_many_dft_c2r(2, n, 1, cu(1,1), nc, 1, 0, &
                                      recov(1,1,1), nr, 1, 0, FFTW_MEASURE)

        ! --- INIT ---
        u = -999.0_c_double
        call random_number(u(1:im, 1:jm, :))
        recov = 0.0_c_double

        ! --- EXECUTE ---
        do k = 1, lm
            ! We execute using the pointers to the specific k-slice
            call fftw_execute_dft_r2c(planf, u(1,1,k), cu(1,1))
            call fftw_execute_dft_c2r(planb, cu(1,1), recov(1,1,k))
        end do

        norm = 1.0_c_double / (real(im, c_double) * real(jm, c_double))
        recov(1:im, 1:jm, :) = recov(1:im, 1:jm, :) * norm
        err = maxval(abs(u(1:im, 1:jm, :) - recov(1:im, 1:jm, :)))

        print "(' Test 6: 2D R2C Identity')"
        print "('  Max Absolute Error:  ', E12.4)", err
        print "('   Result: ', A)", merge("SUCCESS", "FAILURE", err < tolerance)

        call fftw_destroy_plan(planf)
        call fftw_destroy_plan(planb)
    end subroutine test_fftw_2d_single

    subroutine test_fftw_2d_many(im, jm, lm, pad)
        integer, intent(in) :: im, jm, lm, pad
        real(c_double), allocatable, target :: u(:,:,:), recov(:,:,:)
        complex(c_double_complex), allocatable, target :: cu(:,:,:)
        type(c_ptr) :: planf, planb
        integer(c_int) :: n(2), nr(2), nc(2)
        integer(c_int) :: idist, odist
        real(c_double) :: err, tolerance = 1.0e-12_c_double

        allocate(u(1-pad:im+pad, 1-pad:jm+pad, lm))
        allocate(recov(1-pad:im+pad, 1-pad:jm+pad, lm))
        allocate(cu(im/2 + 1, jm, lm))

        ! Logic dimensions (Reversed for C-API)
        n = int([jm, im], c_int)

        ! Embedding: The full physical width of the leading dimensions
        ! The distance between rows in Real is im + 2*pad
        nr = int([jm, im + 2*pad], c_int)
        nc = int([jm, im/2 + 1], c_int)

        ! Distance between 2D planes
        idist = (im + 2*pad) * (jm + 2*pad)
        odist = (im/2 + 1) * jm

        planf = fftw_plan_many_dft_r2c(2, n, int(lm, c_int), &
                u(1,1,1), nr, 1, idist, &
                cu(1,1,1), nc, 1, odist, FFTW_ESTIMATE)

        planb = fftw_plan_many_dft_c2r(2, n, int(lm, c_int), &
                cu(1,1,1), nc, 1, odist, &
                recov(1,1,1), nr, 1, idist, FFTW_ESTIMATE)

        u = -777.0_c_double
        call random_number(u(1:im, 1:jm, :))

        call fftw_execute_dft_r2c(planf, u(1,1,1), cu(1,1,1))
        call fftw_execute_dft_c2r(planb, cu(1,1,1), recov(1,1,1))

        recov(1:im, 1:jm, :) = recov(1:im, 1:jm, :) / (real(im, c_double) * real(jm, c_double))
        err = maxval(abs(u(1:im, 1:jm, :) - recov(1:im, 1:jm, :)))

        print "(' Test 7: 2D Many R2C (Full Volume)')"
        print "('  Max Absolute Error:  ', E12.4)", err
        print "('   Result: ', A)", merge("SUCCESS", "FAILURE", err < tolerance)

        call fftw_destroy_plan(planf)
        call fftw_destroy_plan(planb)
    end subroutine test_fftw_2d_many

end program test_fftw