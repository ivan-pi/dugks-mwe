module periodic_solver_mod
    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: PeriodicMassSolver

    type :: PeriodicMassSolver
        integer :: nx, ny
        real(8) :: scale
        type(C_PTR) :: plan_f = C_NULL_PTR
        type(C_PTR) :: plan_b = C_NULL_PTR
        complex(8), allocatable :: lambdas(:,:)
    contains
        procedure :: init => init_solver
        procedure :: solve
        procedure :: destroy
    end type

    include 'fftw3.f03'

contains

    subroutine init_solver(this, nx, ny)
        class(PeriodicMassSolver), intent(out) :: this
        integer, intent(in) :: nx, ny

        ! Buffers must match the R2C/C2R requirements exactly
        real(8), allocatable    :: work_real(:,:)
        complex(8), allocatable :: work_comp(:,:)

        this%nx = nx
        this%ny = ny
        this%scale = 1.0d0 / (real(nx, 8) * real(ny, 8))

        allocate(this%lambdas(nx/2 + 1, ny))
        allocate(work_real(nx, ny))
        allocate(work_comp(nx, ny))

        ! 1. Create Plans
        ! Use the distinct REAL and COMPLEX buffers
        this%plan_f = fftw_plan_dft_r2c_2d(nx, ny, work_real, work_comp, FFTW_MEASURE)
        this%plan_b = fftw_plan_dft_c2r_2d(nx, ny, work_comp, work_real, FFTW_MEASURE)

        ! 2. Pre-calculate Eigenvalues
        work_real = 0.0d0
        work_real(1, 1) = 4.0d0/9.0d0
        work_real(2, 1) = 1.0d0/9.0d0;  work_real(nx, 1) = 1.0d0/9.0d0
        work_real(1, 2) = 1.0d0/9.0d0;  work_real(1, ny) = 1.0d0/9.0d0
        work_real(2, 2) = 1.0d0/36.0d0; work_real(nx, ny) = 1.0d0/36.0d0
        work_real(nx, 2) = 1.0d0/36.0d0; work_real(2, ny) = 1.0d0/36.0d0

        ! Execute transform using the correct buffers
        call fftw_execute_dft_r2c(this%plan_f, work_real, this%lambdas)

    end subroutine

    subroutine solve(this, u, f)
        class(PeriodicMassSolver), intent(inout) :: this
        real(8), intent(in)  :: f(this%nx, this%ny)
        real(8), intent(out) :: u(this%nx, this%ny)

        complex(8) :: f_hat(this%nx/2 + 1, this%ny)
        real(8)    :: work_in(this%nx, this%ny)

        ! 1. Copy input to local work buffer (to preserve f)
        work_in = f

        ! 2. Forward Transform
        call fftw_execute_dft_r2c(this%plan_f, work_in, f_hat)

        ! 3. Point-wise division
        f_hat = f_hat / this%lambdas

        ! 4. Backward Transform (Inverse)
        ! Note: FFTW c2r transform overwrites the input buffer (f_hat)
        call fftw_execute_dft_c2r(this%plan_b, f_hat, work_in)

        ! 5. Scale and Output
        u = work_in * this%scale
    end subroutine

    subroutine destroy(this)
        class(PeriodicMassSolver), intent(inout) :: this
        if (c_associated(this%plan_f)) call fftw_destroy_plan(this%plan_f)
        if (c_associated(this%plan_b)) call fftw_destroy_plan(this%plan_b)
        if (allocated(this%lambdas)) deallocate(this%lambdas)
    end subroutine

end module


module batch_periodic_solver_mod
    use, intrinsic :: iso_c_binding
    implicit none
    private

    public :: fftw_init_threads
    public :: fftw_plan_with_nthreads
    public :: fftw_cleanup_threads

    public :: BatchMassSolver

    type :: BatchMassSolver
        integer :: nx, ny, k_batch, nhalo
        real(8) :: scale
        type(C_PTR) :: plan_f = C_NULL_PTR
        type(C_PTR) :: plan_b = C_NULL_PTR
        complex(8), allocatable :: lambdas(:,:)
    contains
        procedure :: init => init_batch_solver
        procedure :: solve => solve_batch
        procedure :: update_halos
        procedure :: destroy
    end type

    include 'fftw3.f03'

contains

    subroutine init_batch_solver(this, nx, ny, k_batch, nhalo)
        class(BatchMassSolver), intent(out) :: this
        integer, intent(in) :: nx, ny, k_batch, nhalo

        ! We swap the order for FFTW: [slowest, fastest]
        integer :: n(2), inembed(2), onembed(2)
        integer :: istride, idist, ostride, odist

        ! Temporary pointers for planning
        real(8), target :: dummy_r(1-nhalo:nx+nhalo, 1-nhalo:ny+nhalo, k_batch)
        complex(8), target :: dummy_c(nx/2 + 1, ny, k_batch)

        this%nx = nx
        this%ny = ny
        this%k_batch = k_batch
        this%nhalo = nhalo
        this%scale = 1.0d0 / (real(nx, 8) * real(ny, 8))

        ! --- CRITICAL CHANGE: SWAP DIMENSIONS FOR FORTRAN ---
        ! n(1) is slowest varying (ny), n(2) is fastest varying (nx)
        n = [ny, nx]

        ! inembed(1) is ny dimension (slowest), inembed(2) is nx dimension (fastest)
        ! Leading dimension (fastest) is nx + 2*nhalo
        inembed = [ny + 2*nhalo, nx + 2*nhalo]
        istride = 1
        idist   = (nx + 2*nhalo) * (ny + 2*nhalo)

        ! onembed(1) is ny dimension, onembed(2) is nx/2 + 1 dimension
        onembed = [ny, nx/2 + 1]
        ostride = 1
        odist   = (nx/2 + 1) * ny
        ! ----------------------------------------------------

        allocate(this%lambdas(nx/2 + 1, ny))

        ! 1. Forward Plan
        this%plan_f = fftw_plan_many_dft_r2c(2, n, k_batch, &
            dummy_r(1,1,1), inembed, istride, idist, &
            dummy_c(1,1,1), onembed, ostride, odist, &
            FFTW_MEASURE)

        ! 2. Backward Plan
        this%plan_b = fftw_plan_many_dft_c2r(2, n, k_batch, &
            dummy_c(1,1,1), onembed, ostride, odist, &
            dummy_r(1,1,1), inembed, istride, idist, &
            FFTW_MEASURE)

        call calculate_lambdas_r2c(this)

    contains

        subroutine calculate_lambdas_r2c(this)
            class(BatchMassSolver), intent(inout) :: this
            real(8) :: kernel(this%nx, this%ny)
            complex(8) :: out(this%nx/2 + 1, this%ny)
            type(C_PTR) :: p

            kernel = 0.0d0
            kernel(1, 1) = 4.0d0/9.0d0
            kernel(2, 1) = 1.0d0/9.0d0;  kernel(nx, 1) = 1.0d0/9.0d0
            kernel(1, 2) = 1.0d0/9.0d0;  kernel(1, ny) = 1.0d0/9.0d0
            kernel(2, 2) = 1.0d0/36.0d0; kernel(nx, ny) = 1.0d0/36.0d0
            kernel(nx, 2) = 1.0d0/36.0d0; kernel(2, ny) = 1.0d0/36.0d0

            p = fftw_plan_dft_r2c_2d(this%nx, this%ny, kernel, out, FFTW_ESTIMATE)
            call fftw_execute_dft_r2c(p, kernel, out)
            this%lambdas = out
            call fftw_destroy_plan(p)
        end subroutine

    end subroutine

    subroutine solve_batch(this, u, f)
        class(BatchMassSolver), intent(inout) :: this
        real(8), intent(inout) :: f(1-this%nhalo:this%nx+this%nhalo, 1-this%nhalo:this%ny+this%nhalo, this%k_batch)
        real(8), intent(out) :: u(1-this%nhalo:this%nx+this%nhalo, 1-this%nhalo:this%ny+this%nhalo, this%k_batch)

        ! Workspace for the spectral domain
        ! IMPORTANT: c2r is destructive. We need a fresh copy of f_hat
        ! for each transform if we were doing them sequentially,
        ! but here plan_many handles the whole block.
        complex(8), allocatable :: f_hat(:,:,:)
        allocate(f_hat(this%nx/2 + 1, this%ny, this%k_batch))

        ! 1. Forward Transform (Real -> Complex)
        call fftw_execute_dft_r2c(this%plan_f, f(1,1,1), f_hat(1,1,1))

        ! 2. Point-wise division
        ! (Vectorized across spatial dims, loop across batch)
        block
            integer :: k
            do k = 1, this%k_batch
                f_hat(:,:,k) = f_hat(:,:,k) / this%lambdas
            end do
        end block

        ! 3. Backward Transform (Complex -> Real)
        ! This will write directly into the interior of u
        call fftw_execute_dft_c2r(this%plan_b, f_hat(1,1,1), u(1,1,1))

        ! 4. Normalization
        u(1:this%nx, 1:this%ny, :) = u(1:this%nx, 1:this%ny, :) * this%scale

        deallocate(f_hat)
    end subroutine

    subroutine destroy(this)
        class(BatchMassSolver), intent(inout) :: this
        if (c_associated(this%plan_f)) call fftw_destroy_plan(this%plan_f)
        if (c_associated(this%plan_b)) call fftw_destroy_plan(this%plan_b)
        if (allocated(this%lambdas)) deallocate(this%lambdas)
    end subroutine


    subroutine update_halos(this, u)
        class(BatchMassSolver), intent(in) :: this
        real(8), intent(inout) :: u(1-this%nhalo:this%nx+this%nhalo, 1-this%nhalo:this%ny+this%nhalo,this%k_batch)
        call doit(u,this%nx,this%ny,this%nhalo,this%k_batch)
    contains
        subroutine doit(u, nx, ny, nhalo, k_batch)
            integer, intent(in) :: nx, ny, nhalo, k_batch
            real(8), intent(inout) :: u(1-nhalo:nx+nhalo, 1-nhalo:ny+nhalo, k_batch)

            ! 1. Wrap X-direction (Side strips only)
            ! Copy interior left to exterior right
            u(nx+1:nx+nhalo, 1:ny, :) = u(1:nhalo, 1:ny, :)
            ! Copy interior right to exterior left
            u(1-nhalo:0, 1:ny, :)     = u(nx-nhalo+1:nx, 1:ny, :)

            ! 2. Wrap Y-direction (Full width, including the X-halos just updated)
            ! This step fills the 4 corner blocks automatically
            u(1-nhalo:nx+nhalo, ny+1:ny+nhalo, :) = u(1-nhalo:nx+nhalo, 1:nhalo, :)
            u(1-nhalo:nx+nhalo, 1-nhalo:0, :)     = u(1-nhalo:nx+nhalo, ny-nhalo+1:ny, :)

        end subroutine
    end subroutine update_halos

end module
