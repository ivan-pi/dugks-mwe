program test_periodic_solver
    use periodic_solver_mod
    implicit none

    call test_trivial
    call test_harmonic
    call test_stencil_integrity
    call test_batched_stencil_integrity

contains

    subroutine test_trivial()
    integer, parameter :: nx = 8, ny = 8
    type(PeriodicMassSolver) :: solver
    real(8) :: f(nx, ny), u(nx, ny)
    integer :: i, j
    logical :: success

    ! 1. Initialize the solver
    ! This sets up the FFTW plans and pre-calculates the stencil eigenvalues
    call solver%init(nx, ny)

    ! 2. Set RHS to a vector of ones (Constant field)
    f = 1.0d0

    ! 3. Solve M u = f
    call solver%solve(u, f)

    ! 4. Verify Results
    print *, "Testing Constant Field (All Ones)..."
    success = .true.
    check: do j = 1, ny
        do i = 1, nx
            if (abs(u(i,j) - 1.0d0) > 1.0d-12) then
                print *, "Error at ", i, j, " Value: ", u(i,j)
                success = .false.
                exit check
            end if
        end do
    end do check

    if (success) then
        print *, "Verification Successful! u is 1.0 everywhere."
    else
        print *, "Verification Failed."
    end if

    ! 5. Try a "Delta" function (Single Pulse)
    f = 0.0d0
    f(nx/2, ny/2) = 1.0d0
    call solver%solve(u, f)

    print *, ""
    print *, "Testing Delta Function (Pulse at Center)..."
    print *, "Solution at center (should be > 1 due to mass matrix inversion):"
    print '(F10.6)', u(nx/2, ny/2)
    print *, "Solution at immediate neighbor (should be negative/oscillatory):"
    print '(F10.6)', u(nx/2 + 1, ny/2)

    ! Clean up
    call solver%destroy()
    end subroutine

    subroutine test_harmonic()
    use periodic_solver_mod
    implicit none

    integer, parameter :: nx = 32, ny = 32
    type(PeriodicMassSolver) :: solver
    real(8) :: f(nx, ny), u(nx, ny), u_theory(nx, ny)
    real(8) :: x, y, pi, arg_x, arg_y, lambda_theory
    integer :: i, j, kx, ky

    pi = acos(-1.0d0)
    call solver%init(nx, ny)

    ! 1. Choose a test frequency (Mode 1 in X, Mode 2 in Y)
    kx = 1
    ky = 2

    ! 2. Calculate the theoretical eigenvalue for this 9-point stencil
    ! For a stencil with weights [w_center, w_adj, w_diag]:
    ! Lambda = w_c + 2*w_adj*(cos(tx) + cos(ty)) + 4*w_diag*cos(tx)*cos(ty)
    arg_x = 2.0d0 * pi * kx / nx
    arg_y = 2.0d0 * pi * ky / ny

    lambda_theory = (4.0d0/9.0d0) + &
                    (2.0d0/9.0d0) * (cos(arg_x) + cos(arg_y)) + &
                    (4.0d0/36.0d0) * (cos(arg_x) * cos(arg_y))

    ! 3. Set RHS as the harmonic function
    do j = 1, ny
        do i = 1, nx
            x = real(i-1, 8)
            y = real(j-1, 8)
            f(i, j) = cos(2.0d0 * pi * kx * x / nx) * cos(2.0d0 * pi * ky * y / ny)
        end do
    end do

    ! 4. Solve M u = f
    call solver%solve(u, f)

    ! 5. Compare against Theory: u_theory = f / lambda_theory
    u_theory = f / lambda_theory

    print *, "Testing Harmonic Mode (kx=1, ky=2)..."
    print *, "Theoretical Eigenvalue:", lambda_theory
    print *, "Max Absolute Error:    ", maxval(abs(u - u_theory))

    if (maxval(abs(u - u_theory)) < 1.0d-12) then
        print *, "SUCCESS: Solver perfectly matches spectral theory."
    else
        print *, "FAILURE: Discrepancy detected."
    end if

    call solver%destroy()

    end subroutine test_harmonic

    subroutine test_stencil_integrity()
    use periodic_solver_mod
    implicit none

    integer, parameter :: nx = 16, ny = 16
    type(PeriodicMassSolver) :: solver
    real(8) :: u_orig(nx, ny), f_manual(nx, ny), u_rec(nx, ny)
    real(8) :: pi, x, y, err
    integer :: i, j, ip, im, jp, jm
    real(8), parameter :: w_c = 4.0d0/9.0d0
    real(8), parameter :: w_a = 1.0d0/9.0d0
    real(8), parameter :: w_d = 1.0d0/36.0d0

    pi = acos(-1.0d0)
    call solver%init(nx, ny)

    ! 1. Create a non-symmetric test function u
    do j = 1, ny
        do i = 1, nx
            x = real(i-1, 8)/nx
            y = real(j-1, 8)/ny
            u_orig(i, j) = sin(2*pi*x) + 0.5*cos(4*pi*y) + 0.3*sin(2*pi*(x+y))
        end do
    end do

    ! 2. Manually apply the stencil (The "Definition" of M*u)
    ! This tests every direction: Up, Down, Left, Right, and all 4 Diagonals
    do j = 1, ny
        do i = 1, nx
            ! Periodic indices
            ip = mod(i, nx) + 1
            im = mod(i + nx - 2, nx) + 1
            jp = mod(j, ny) + 1
            jm = mod(j + ny - 2, ny) + 1

            f_manual(i, j) = w_c * u_orig(i, j) + &                          ! Center
                             w_a * (u_orig(ip, j) + u_orig(im, j) + &         ! Right, Left
                                    u_orig(i, jp) + u_orig(i, jm)) + &        ! Up, Down
                             w_d * (u_orig(ip, jp) + u_orig(im, jm) + &       ! UR, DL
                                    u_orig(im, jp) + u_orig(ip, jm))          ! UL, DR
        end do
    end do

    ! 3. Solve using the FFT solver: u_rec = M^-1 * f_manual
    call solver%solve(u_rec, f_manual)

    ! 4. Error Analysis
    err = maxval(abs(u_rec - u_orig))
    print *, "Rigorous Stencil Integrity Test"
    print *, "Max Difference (L-infinity):", err

    if (err < 1.0d-13) then
        print *, "PASSED: FFT solver is exactly consistent with the 9-point stencil."
    else
        print *, "FAILED: Potential indexing mismatch in stencil vs FFT kernel."
    end if

    call solver%destroy()
    end subroutine test_stencil_integrity


    subroutine test_batched_stencil_integrity()
    use batch_periodic_solver_mod
    implicit none

    integer, parameter :: nx = 16, ny = 16
    type(BatchMassSolver) :: solver
    real(8), dimension(0:nx+1,0:ny+1) :: u_orig, f_manual, u_rec
    real(8) :: pi, x, y, err
    integer :: i, j, ip, im, jp, jm
    real(8), parameter :: w_c = 4.0d0/9.0d0
    real(8), parameter :: w_a = 1.0d0/9.0d0
    real(8), parameter :: w_d = 1.0d0/36.0d0

    pi = acos(-1.0d0)
    call solver%init(nx, ny, 1, nhalo=1)

    ! 1. Create a non-symmetric test function u
    do j = 1, ny
        do i = 1, nx
            x = real(i-1, 8)/nx
            y = real(j-1, 8)/ny
            u_orig(i, j) = sin(2*pi*x) + 0.5*cos(4*pi*y) + 0.3*sin(2*pi*(x+y))
        end do
    end do

    ! 2. Manually apply the stencil (The "Definition" of M*u)
    ! This tests every direction: Up, Down, Left, Right, and all 4 Diagonals
    do j = 1, ny
        do i = 1, nx
            ! Periodic indices
            ip = mod(i, nx) + 1
            im = mod(i + nx - 2, nx) + 1
            jp = mod(j, ny) + 1
            jm = mod(j + ny - 2, ny) + 1

            f_manual(i, j) = w_c * u_orig(i, j) + &                          ! Center
                             w_a * (u_orig(ip, j) + u_orig(im, j) + &         ! Right, Left
                                    u_orig(i, jp) + u_orig(i, jm)) + &        ! Up, Down
                             w_d * (u_orig(ip, jp) + u_orig(im, jm) + &       ! UR, DL
                                    u_orig(im, jp) + u_orig(ip, jm))          ! UL, DR
        end do
    end do

    ! 3. Solve using the FFT solver: u_rec = M^-1 * f_manual
    call solver%solve(u_rec, f_manual)

    ! 4. Error Analysis
    err = maxval(abs(u_rec(1:nx,1:ny) - u_orig(1:nx,1:ny)))
    print *, "Rigorous Stencil Integrity Test"
    print *, "Max Difference (L-infinity):", err

    if (err < 1.0d-13) then
        print *, "PASSED: FFT solver is exactly consistent with the 9-point stencil."
    else
        print *, "FAILED: Potential indexing mismatch in stencil vs FFT kernel."
    end if

    call solver%destroy()
    end subroutine test_batched_stencil_integrity

end program test_periodic_solver
