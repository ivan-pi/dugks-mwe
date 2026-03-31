program test_fftw_interface
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    ! --- Main Driver ---
    print *, "Starting FFTW R2C Interface Tests..."
    print *, "------------------------------------"
    call test_1d_differentiation()
    print *, "------------------------------------"
    call test_1d_periodic_fe_solver()
    print *, "------------------------------------"
    call test_1d_cyclic_shift()
    print *, "------------------------------------"
    call test_1d_parseval_energy()
    print *, "------------------------------------"
    print *, "Tests Completed."

contains

    subroutine test_1d_periodic_fe_solver()
        integer, parameter :: DP = C_DOUBLE
        integer, parameter :: DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), rank, howmany
        real(DP), allocatable :: x_orig(:), y_phys(:), x_final(:), row_phys(:)
        complex(DCP), allocatable :: y_spec(:), a_spec(:)
        type(C_PTR) :: forward_plan, backward_plan
        integer :: i, N_size
        real(DP) :: err

        rank = 1
        N_size = 16
        n(1) = N_size
        howmany = 1
        
        allocate(x_orig(N_size), y_phys(N_size), x_final(N_size), row_phys(N_size))
        allocate(y_spec(N_size/2 + 1), a_spec(N_size/2 + 1))

        ! 1. Define the Finite Element Operator (Circulant Row)
        ! Row entries: 4/6 (diag), 1/6 (off-diag). Periodic boundary wraps 1/6 to the end.
        row_phys = 0.0_DP
        row_phys(1) = 4.0_DP/6.0_DP
        row_phys(2) = 1.0_DP/6.0_DP
        row_phys(N_size) = 1.0_DP/6.0_DP

        ! 2. Define a known input x (a simple pulse or ramp)
        do i = 1, N_size
            x_orig(i) = real(i, DP)
        end do

        ! 3. Perform Matrix Multiplication (Direct Convolution in physical space)
        ! y_i = 1/6*x_{i-1} + 4/6*x_i + 1/6*x_{i+1} (with periodic wrapping)
        do i = 1, N_size
            y_phys(i) = (4.0_DP/6.0_DP) * x_orig(i) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i-2+N_size, N_size)+1) + &
                        (1.0_DP/6.0_DP) * x_orig(mod(i, N_size)+1)
        end do

        ! 4. Setup FFTW plans to solve the inverse
        forward_plan = fftw_plan_many_dft_r2c(1, n, 1, y_phys, n, 1, n(1), &
                                              y_spec, [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        
        ! We also need the eigenvalues of the matrix (FFT of the row)
        ! Re-using the forward plan logic for the row
        call fftw_execute_dft_r2c(forward_plan, row_phys, a_spec)
        
        ! Transform physical product 'y' to Fourier space
        call fftw_execute_dft_r2c(forward_plan, y_phys, y_spec)

        ! 5. Solve in Spectral Space: x_hat = y_hat / a_hat
        y_spec = y_spec / a_spec

        ! 6. Transform back to Physical Space
        backward_plan = fftw_plan_many_dft_c2r(1, n, 1, y_spec, [n(1)/2+1], 1, n(1)/2+1, &
                                               x_final, n, 1, n(1), FFTW_ESTIMATE)
        call fftw_execute_dft_c2r(backward_plan, y_spec, x_final)

        ! 7. Normalize (FFTW does not scale inverse) and Check Error
        x_final = x_final / real(N_size, DP)
        
        err = maxval(abs(x_orig - x_final))
        print *, "Test 3: 1D Periodic FE Solver (Matrix Inverse)"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        if (err < 1e-12) then
            print *, "  Result: SUCCESS"
        else
            print *, "  Result: FAILURE"
        end if

        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        deallocate(x_orig, y_phys, x_final, row_phys, y_spec, a_spec)
    end subroutine test_1d_periodic_fe_solver

subroutine test_1d_differentiation()
        integer, parameter :: DP = C_DOUBLE
        integer, parameter :: DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), rank, howmany
        real(DP), allocatable :: phys_in(:), phys_out(:), phys_ref(:)
        complex(DCP), allocatable :: spec_data(:)
        type(C_PTR) :: plan_fwd, plan_bwd
        integer :: i, N_size
        real(DP) :: L, dx, k_val, err

        real(DP), parameter :: pi = 4.0_DP * atan(1.0_DP)

        N_size = 16
        n(1) = N_size
        howmany = 1
        L = 2.0_DP * pi
        dx = L / N_size

        allocate(phys_in(N_size), phys_out(N_size), phys_ref(N_size))
        allocate(spec_data(N_size/2 + 1))

        ! 1. Initialize Signal: f(x) = sin(x), f'(x) = cos(x)
        do i = 1, N_size
            phys_in(i)  = sin((i-1) * dx)
            phys_ref(i) = cos((i-1) * dx)
        end do

        ! 2. Create Plans
        plan_fwd = fftw_plan_many_dft_r2c(1, n, 1, phys_in, n, 1, n(1), &
                                          spec_data, [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        plan_bwd = fftw_plan_many_dft_c2r(1, n, 1, spec_data, [n(1)/2+1], 1, n(1)/2+1, &
                                          phys_out, n, 1, n(1), FFTW_ESTIMATE)

        ! 3. Forward Transform
        call fftw_execute_dft_r2c(plan_fwd, phys_in, spec_data)

        ! 4. Spectral Differentiation: Multiply by i*k
        ! k = 2*pi*m / L, where m is the mode index
        do i = 1, N_size/2 + 1
            k_val = (2.0_DP * pi / L) * real(i-1, DP)
            spec_data(i) = spec_data(i) * (0.0_DP, 1.0_DP) * k_val
        end do

        ! Handle the Nyquist frequency if N is even (it must be purely real)
        spec_data(N_size/2 + 1) = cmplx(real(spec_data(N_size/2 + 1)), 0.0_DP, kind=DCP)

        ! 5. Backward Transform
        call fftw_execute_dft_c2r(plan_bwd, spec_data, phys_out)

        ! 6. Normalize and Compare
        phys_out = phys_out / real(N_size, DP)
        err = maxval(abs(phys_out - phys_ref))

        print *, "Test 1: 1D R2C Differentiation (sin(x) -> cos(x))"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        if (err < 1e-12) then
            print *, "  Result: SUCCESS"
        else
            print *, "  Result: FAILURE"
        end if

        call fftw_destroy_plan(plan_fwd)
        call fftw_destroy_plan(plan_bwd)
        deallocate(phys_in, phys_out, phys_ref, spec_data)
    end subroutine test_1d_differentiation

    ! ==========================================================================
    ! TEST 4: Cyclic Shift Property
    ! Verifies: f(x - a) <-> exp(-i * k * a) * F(k)
    ! ==========================================================================
    subroutine test_1d_cyclic_shift()
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), N_size
        real(DP), allocatable :: x_orig(:), x_shifted(:), x_recovered(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: p_fwd, p_bwd
        integer :: i, shift_bins
        real(DP) :: phase, arg, err

        N_size = 32
        n(1) = N_size
        shift_bins = 5 ! Shift the signal by 5 grid points

        allocate(x_orig(N_size), x_shifted(N_size), x_recovered(N_size))
        allocate(spec(N_size/2 + 1))

        ! Gaussian pulse
        do i = 1, N_size
            x_orig(i) = exp(-((real(i-N_size/2, DP)/4.0_DP)**2))
            x_shifted(i) = exp(-((real(mod(i-1-shift_bins+N_size, N_size)-N_size/2+1, DP)/4.0_DP)**2))
        end do

        p_fwd = fftw_plan_many_dft_r2c(1, n, 1, x_orig, n, 1, n(1), spec, [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        p_bwd = fftw_plan_many_dft_c2r(1, n, 1, spec, [n(1)/2+1], 1, n(1)/2+1, x_recovered, n, 1, n(1), FFTW_ESTIMATE)

        call fftw_execute_dft_r2c(p_fwd, x_orig, spec)

        ! Apply shift in Fourier Space: multiply by exp(-i * 2*pi * k * shift / N)
        do i = 1, N_size/2 + 1
            arg = -2.0_DP * 3.141592653589793_DP * real(i-1, DP) * real(shift_bins, DP) / real(N_size, DP)
            spec(i) = spec(i) * cmplx(cos(arg), sin(arg), kind=DCP)
        end do

        call fftw_execute_dft_c2r(p_bwd, spec, x_recovered)
        x_recovered = x_recovered / real(N_size, DP)

        err = maxval(abs(x_recovered - x_shifted))
        print *, "Test 4: 1D Cyclic Shift (Phase shift in Fourier Space)"
        print '(A, E12.4)', "  Max Absolute Error: ", err
        if (err < 1e-12) then
            print *, "  Result: SUCCESS"
        else
            print *, "  Result: FAILURE"
        end if

        call fftw_destroy_plan(p_fwd)
        call fftw_destroy_plan(p_bwd)
        deallocate(x_orig, x_shifted, x_recovered, spec)
    end subroutine test_1d_cyclic_shift


    ! ==========================================================================
    ! TEST 5: Parseval's Theorem
    ! Verifies: Sum(|x_n|^2) = (1/N) * Sum(|X_k|^2)
    ! ==========================================================================
    subroutine test_1d_parseval_energy()
        integer, parameter :: DP = C_DOUBLE, DCP = C_DOUBLE_COMPLEX
        integer(C_INT) :: n(1), N_size
        real(DP), allocatable :: x(:)
        complex(DCP), allocatable :: spec(:)
        type(C_PTR) :: plan
        real(DP) :: energy_phys, energy_spec, err
        integer :: i

        N_size = 64
        n(1) = N_size
        allocate(x(N_size), spec(N_size/2 + 1))
        call random_number(x)

        plan = fftw_plan_many_dft_r2c(1, n, 1, x, n, 1, n(1), spec, [n(1)/2+1], 1, n(1)/2+1, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan, x, spec)

        energy_phys = sum(x**2)

        ! In R2C, modes 1 and N/2+1 (if N is even) are unique.
        ! Intermediate modes represent two points in the full complex spectrum.
        energy_spec = abs(spec(1))**2  ! DC component
        do i = 2, N_size/2
            energy_spec = energy_spec + 2.0_DP * abs(spec(i))**2
        end do
        energy_spec = energy_spec + abs(spec(N_size/2 + 1))**2 ! Nyquist
        energy_spec = energy_spec / real(N_size, DP)

        err = abs(energy_phys - energy_spec)
        print *, "Test 5: Parseval's Theorem (Energy Conservation)"
        print '(A, E12.4)', "  Energy Difference: ", err
        if (err < 1e-12) then
            print *, "  Result: SUCCESS"
        else
            print *, "  Result: FAILURE"
        end if

        call fftw_destroy_plan(plan)
        deallocate(x, spec)
    end subroutine test_1d_parseval_energy

end program test_fftw_interface
