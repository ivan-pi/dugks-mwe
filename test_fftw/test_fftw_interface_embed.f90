program test_fftw_interface
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    logical :: use_embedding
    character(len=32) :: arg

    ! --- Simple CLI check ---
    use_embedding = .false.
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        if (trim(arg) == "--embedded" .or. trim(arg) == "1") use_embedding = .true.
    end if

    print *, "Starting FFTW R2C Interface Tests..."
    print *, "Mode: ", merge("EMBEDDED (Custom Bounds)", "STANDARD (Contiguous)   ", use_embedding)
    print *, "------------------------------------"
    
    call test_1d_differentiation(use_embedding)
    print *, "------------------------------------"
    call test_1d_periodic_fe_solver(use_embedding)
    print *, "------------------------------------"
    call test_1d_cyclic_shift(use_embedding)
    print *, "------------------------------------"
    call test_1d_parseval_energy(use_embedding)
    print *, "------------------------------------"
    print *, "Tests Completed."

contains

    subroutine test_1d_differentiation(embedded)
        logical, intent(in) :: embedded
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        ! Use deferred-shape arrays with custom bounds
        real(DP), allocatable :: phys_in(:), phys_out(:), phys_ref(:)
        complex(DCP), allocatable :: spec_data(:)
        type(C_PTR) :: plan_fwd, plan_bwd
        integer :: i, n_size, pad
        real(DP) :: L, dx, k_val, err

        real(DP), parameter :: pi = 4.0_DP * atan(1.0_DP)

        n_size = 16
        n(1) = n_size
        pad = merge(2, 0, embedded)
        
        ! The "Active" data is 1:n_size. The "Ghost" data is 1-pad:0 and n_size+1:n_size+pad.
        allocate(phys_in(1-pad : n_size+pad))
        allocate(phys_out(1-pad : n_size+pad))
        allocate(phys_ref(1 : n_size))
        allocate(spec_data(1 : n_size/2 + 1))

        ! inembed must reflect the TOTAL length of the allocated dimension
        inembed(1) = n_size + 2*pad
        L = 2.0_DP * pi
        dx = L / n_size

        ! Fill ghosts with garbage to ensure they are ignored
        phys_in = -999.0_DP
        do i = 1, n_size
            phys_in(i)  = sin((i-1) * dx)
            phys_ref(i) = cos((i-1) * dx)
        end do

        ! We pass phys_in(1) to FFTW. Because of custom bounds, 
        ! this is exactly where the transform should start.
        plan_fwd = fftw_plan_many_dft_r2c(1, n, 1, phys_in(1), inembed, 1, n(1), &
                                          spec_data(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        plan_bwd = fftw_plan_many_dft_c2r(1, n, 1, spec_data(1), [n(1)/2+1], 1, n(1)/2+1, &
                                          phys_out(1), inembed, 1, n(1), FFTW_ESTIMATE)

        call fftw_execute_dft_r2c(plan_fwd, phys_in(1), spec_data(1))

        ! Differentiation
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
        print *, "  Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)
        
        call fftw_destroy_plan(plan_fwd); call fftw_destroy_plan(plan_bwd)
    end subroutine test_1d_differentiation

    subroutine test_1d_periodic_fe_solver(embedded)
        logical, intent(in) :: embedded
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x_orig(:), y_phys(:), x_final(:), row_phys(:)
        complex(DCP), allocatable :: y_spec(:), a_spec(:)
        type(C_PTR) :: fwd_p, bwd_p
        integer :: i, n_size, pad
        real(DP) :: err

        n_size = 16
        n(1) = n_size
        pad = merge(3, 0, embedded)
        inembed(1) = n_size + 2*pad
        
        allocate(x_orig(1-pad:n_size+pad), y_phys(1-pad:n_size+pad), &
                 x_final(1-pad:n_size+pad), row_phys(1-pad:n_size+pad))
        allocate(y_spec(n_size/2+1), a_spec(n_size/2+1))

        x_orig(1:n_size) = [(real(i, DP), i=1, n_size)]
        
        row_phys = 0.0_DP
        row_phys(1) = 4.0_DP/6.0_DP
        row_phys(2) = 1.0_DP/6.0_DP
        row_phys(n_size) = 1.0_DP/6.0_DP

        do i = 1, n_size
            y_phys(i) = (4.0_DP/6.0_DP) * x_orig(i) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i-2+n_size, n_size)+1) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i, n_size)+1)
        end do

        fwd_p = fftw_plan_many_dft_r2c(1, n, 1, y_phys(1), inembed, 1, n(1), &
                                      y_spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        
        call fftw_execute_dft_r2c(fwd_p, row_phys(1), a_spec(1))
        call fftw_execute_dft_r2c(fwd_p, y_phys(1), y_spec(1))

        y_spec = y_spec / a_spec

        bwd_p = fftw_plan_many_dft_c2r(1, n, 1, y_spec(1), [n(1)/2+1], 1, n(1)/2+1, &
                                      x_final(1), inembed, 1, n(1), FFTW_ESTIMATE)
        call fftw_execute_dft_c2r(bwd_p, y_spec(1), x_final(1))

        x_final(1:n_size) = x_final(1:n_size) / real(n_size, DP)
        err = maxval(abs(x_orig(1:n_size) - x_final(1:n_size)))
        
        print *, "Test 3: 1D Periodic FE Solver"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        print *, "  Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)
        call fftw_destroy_plan(fwd_p); call fftw_destroy_plan(bwd_p)
    end subroutine test_1d_periodic_fe_solver

    subroutine test_1d_cyclic_shift(embedded)
        logical, intent(in) :: embedded
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x_o(:), x_s(:), x_r(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: pf, pb
        integer :: i, n_size, pad, shift = 5
        real(DP) :: arg, err

        real(DP), parameter :: pi = 4.0_DP * atan(1.0_DP)

        n_size = 32
        n(1) = n_size
        pad = merge(5, 0, embedded)
        inembed(1) = n_size + 2*pad
        
        allocate(x_o(1-pad:n_size+pad), x_s(1-pad:n_size+pad), &
                 x_r(1-pad:n_size+pad), spec(1:n_size/2+1))

        do i = 1, n_size
            x_o(i) = exp(-((real(i-n_size/2, DP)/4.0_DP)**2))
            x_s(i) = exp(-((real(mod(i-1-shift+n_size, n_size)-n_size/2+1, DP)/4.0_DP)**2))
        end do

        pf = fftw_plan_many_dft_r2c(1, n, 1, x_o(1), inembed, 1, n(1), spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        pb = fftw_plan_many_dft_c2r(1, n, 1, spec(1), [n(1)/2+1], 1, n(1)/2+1, x_r(1), inembed, 1, n(1), FFTW_ESTIMATE)

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
        print *, "  Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)
        call fftw_destroy_plan(pf); call fftw_destroy_plan(pb)
    end subroutine test_1d_cyclic_shift

    subroutine test_1d_parseval_energy(embedded)
        logical, intent(in) :: embedded
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), inembed(1)
        real(DP), allocatable :: x(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: plan
        real(DP) :: e_phys, e_spec, err
        integer :: i, n_size, pad

        n_size = 64
        n(1) = n_size
        pad = merge(10, 0, embedded)
        inembed(1) = n_size + 2*pad
        
        allocate(x(1-pad : n_size+pad), spec(1 : n_size/2 + 1))
        x = -1.0_DP
        call random_number(x(1 : n_size))

        plan = fftw_plan_many_dft_r2c(1, n, 1, x(1), inembed, 1, n(1), spec(1), [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
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
        print *, "  Result: ", merge("SUCCESS", "FAILURE", err < 1.0e-12_DP)
        call fftw_destroy_plan(plan)
    end subroutine test_1d_parseval_energy

end program test_fftw_interface