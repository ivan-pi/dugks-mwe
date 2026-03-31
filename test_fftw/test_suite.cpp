#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

// Helper to compare values
bool is_close(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

/**
 * TEST 1: 2D Batch with Custom Padding (nembed)
 * Use std::complex<double> which is binary compatible with fftw_complex
 */
void test_2d_batch_embedded() {
    std::cout << "Running Test 1: 2D Batch Embedded (nembed)..." << std::endl;

    const int rank = 2;
    int n[] = {4, 4};    // Logical size
    int howmany = 2;     // Testing 2 slices

    // Physical allocation including 1 ghost cell on all sides (6x6)
    int inembed[] = {6, 6};
    int onembed[] = {6, 6/2 + 1}; // Rows=6, Complex Columns=4

    int idist = inembed[0] * inembed[1]; // 36 doubles
    int odist = onembed[0] * onembed[1]; // 24 complex numbers

    std::vector<double> in(idist * howmany, 0.0);
    std::vector<std::complex<double>> out(odist * howmany);

    // --- Slice 0 ---
    // Place 10.0 in the first logical element (row 1, col 1 of the 6x6)
    in[0 * idist + 1 * inembed[1] + 1] = 10.0;

    // --- Slice 1 ---
    // Place 20.0 in the first logical element of the second batch
    in[1 * idist + 1 * inembed[1] + 1] = 20.0;

    // Pointer to the first logical element of Slice 0 (Batch 0)
    // Offset only by the padding of the first slice
    double* in_ptr = &in[1 * inembed[1] + 1];

    // Similarly for the output, point to the first logical complex number
    // Offset by 1 row (onembed[1] wide)
    std::complex<double>* out_ptr = &out[1 * onembed[1] + 0];

    // We are passing a pointer to the first logical element (row 1, col 1).
    // The library is expected to use 'inembed' to determine the stride
    // between rows and 'idist' to jump between batch slices.

    fftw_plan plan = fftw_plan_many_dft_r2c(rank, n, howmany,
                                            in_ptr, inembed, 1, idist,
                                            reinterpret_cast<fftw_complex*>(out_ptr), onembed, 1, odist,
                                            FFTW_ESTIMATE);

    fftw_execute(plan);

    // Verify both batches
    // Note: Because we passed out_ptr as the start of the logical data,
    // the DC for slice 0 is at out_ptr[0] and slice 1 is at out_ptr[odist]
    bool slice0_ok = is_close(out_ptr[0].real(), 10.0);
    bool slice1_ok = is_close(out_ptr[odist].real(), 20.0);

    if (slice0_ok && slice1_ok) {
        std::cout << "  Result: SUCCESS" << std::endl;
        std::cout << "    Slice 0 DC: " << out_ptr[0].real() << std::endl;
        std::cout << "    Slice 1 DC: " << out_ptr[odist].real() << std::endl;
    } else {
        std::cout << "  Result: FAILURE" << std::endl;
        if (!slice0_ok) std::cout << "    Slice 0 failed! Got: " << out_ptr[0].real() << std::endl;
        if (!slice1_ok) std::cout << "    Slice 1 failed! Got: " << out_ptr[odist].real() << std::endl;
    }

    fftw_destroy_plan(plan);
}

/**
 * TEST 2: Interleaved Strides (istride)
 */
void test_stride_interleaved() {
    std::cout << "Running Test 2: Interleaved Strides (istride)..." << std::endl;

    int n[] = {8};
    int howmany = 1;
    int istride = 3;

    std::vector<double> in(n[0] * istride, 0.0);
    std::vector<std::complex<double>> out(n[0] / 2 + 1);

    // Fill "Channel 2" (index 1, 4, 7...)
    for (int i = 0; i < n[0]; ++i) in[i * istride + 1] = 1.0;

    fftw_plan plan = fftw_plan_many_dft_r2c(1, n, howmany,
                                            &in[1], NULL, istride, n[0],
                                            reinterpret_cast<fftw_complex*>(out.data()), NULL, 1, n[0]/2+1,
                                            FFTW_ESTIMATE);

    fftw_execute(plan);

    if (is_close(out[0].real(), 8.0)) {
        std::cout << "  Result: SUCCESS (DC component: " << out[0].real() << ")" << std::endl;
    } else {
        std::cout << "  Result: FAILURE" << std::endl;
    }

    fftw_destroy_plan(plan);
}

/**
 * TEST 3: In-Place with Distances
 */
void test_inplace_advanced() {
    std::cout << "Running Test 3: In-Place Advanced..." << std::endl;

    int n[] = {10};
    int howmany = 1;

    // In-place requires padding. N real entries need space for N/2+1 complex entries.
    // That is (10/2 + 1) * 2 = 12 doubles.
    int real_padded_size = 2 * (n[0] / 2 + 1);
    std::vector<double> data(real_padded_size, 0.0);

    for (int i = 0; i < n[0]; ++i) data[i] = 1.0;

    fftw_plan plan = fftw_plan_many_dft_r2c(1, n, howmany,
                                            data.data(), NULL, 1, real_padded_size,
                                            reinterpret_cast<fftw_complex*>(data.data()), NULL, 1, real_padded_size/2,
                                            FFTW_ESTIMATE);

    fftw_execute(plan);

    if (is_close(data[0], 10.0)) {
        std::cout << "  Result: SUCCESS (DC component: " << data[0] << ")" << std::endl;
    } else {
        std::cout << "  Result: FAILURE" << std::endl;
    }

    fftw_destroy_plan(plan);
}

/**
 * TEST 4: 2D Tight Batch
 */
void test_2d_tight_batch() {
    std::cout << "Running Test 4: 2D Tight Batch (odist validation)..." << std::endl;

    const int rank = 2;
    int n[] = {4, 4};    // Logical 4x4
    int howmany = 3;     // 3 separate slices

    // Logical output is 4 x (4/2 + 1) = 4 x 3
    int idist = n[0] * n[1];             // 16 doubles
    int odist = n[0] * (n[1] / 2 + 1);    // 12 complex numbers

    std::vector<double> in(idist * howmany, 1.0); // All 1s
    std::vector<std::complex<double>> out(odist * howmany, 0.0);

    fftw_plan plan = fftw_plan_many_dft_r2c(rank, n, howmany,
                                            in.data(), NULL, 1, idist,
                                            reinterpret_cast<fftw_complex*>(out.data()), NULL, 1, odist,
                                            FFTW_ESTIMATE);
    fftw_execute(plan);

    // If odist is correct, each batch's DC component (index 0, 12, 24) should be 16.0
    bool success = true;
    for (int i = 0; i < howmany; ++i) {
        if (!is_close(out[i * odist].real(), 16.0)) {
            std::cout << "  Batch " << i << " failed! Expected 16, got " << out[i * odist].real() << std::endl;
            success = false;
        }
    }

    if (success) std::cout << "  Result: SUCCESS (Batch distances verified)" << std::endl;
    fftw_destroy_plan(plan);
}

/**
 * TEST 5: Spectral differentation in 2-d
 */
void test_2d_advanced_differentiation() {
    std::cout << "Running Test 5: 2D Advanced many_dft (Spectral Deriv)..." << std::endl;

    const int nx = 8, ny = 8;
    int n[] = {nx, ny};
    int rank = 2;
    int howmany = 1;

    // Add padding (ghost cells) to the physical layout
    // Logical: 8x8, Physical: 10x10
    int pad = 1;
    int inembed[] = {nx + 2 * pad, ny + 2 * pad};
    int onembed[] = {nx + 2 * pad, ny / 2 + 1}; // Complex padding on rows

    size_t in_total = inembed[0] * inembed[1];
    size_t out_total = onembed[0] * onembed[1];

    std::vector<double> phys_in(in_total, -1.0); // Garbage in ghosts
    std::vector<double> phys_out(in_total, 0.0);
    std::vector<std::complex<double>> spec(out_total);

    double dx = (2.0 * M_PI) / nx;
    double dy = (2.0 * M_PI) / ny;

    // 1. Initialize: f(x,y) = sin(x)cos(y)
    // Offset pointers by 'pad' to start at logical index (0,0)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            phys_in[(i + pad) * inembed[1] + (j + pad)] = std::sin(i * dx) * std::cos(j * dy);
        }
    }

    // 2. Setup Plans using the advanced pointers
    // in: pointer to phys_in[pad][pad]
    double* in_ptr = &phys_in[pad * inembed[1] + pad];
    double* out_ptr = &phys_out[pad * inembed[1] + pad];
    fftw_complex* spec_ptr = reinterpret_cast<fftw_complex*>(&spec[pad * onembed[1]]);

    fftw_plan p_fwd = fftw_plan_many_dft_r2c(rank, n, howmany,
                                            in_ptr, inembed, 1, in_total,
                                            spec_ptr, onembed, 1, out_total,
                                            FFTW_ESTIMATE);

    fftw_plan p_bwd = fftw_plan_many_dft_c2r(rank, n, howmany,
                                            spec_ptr, onembed, 1, out_total,
                                            out_ptr, inembed, 1, in_total,
                                            FFTW_ESTIMATE);

    // 3. Execute Forward
    fftw_execute(p_fwd);

    // 4. Differentiate w.r.t X: Multiply by i*kx
    for (int i = 0; i < nx; ++i) {
        int kx = (i <= nx / 2) ? i : i - nx;
        for (int j = 0; j < (ny / 2 + 1); ++j) {
            // Apply multiplier to the spectral data inside the "onembed" bounds
            spec[(i + pad) * onembed[1] + j] *= std::complex<double>(0, kx);
        }
    }

    // 5. Execute Backward
    fftw_execute(p_bwd);

    // 6. Verify inner 8x8 result (f'(x) = cos(x)cos(y))
    double max_err = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double expected = std::cos(i * dx) * std::cos(j * dy);
            double actual = phys_out[(i + pad) * inembed[1] + (j + pad)] / (nx * ny);
            max_err = std::max(max_err, std::abs(actual - expected));
        }
    }

    std::cout << "  Max Error (2D many_dft d/dx): " << max_err << std::endl;
    if (max_err < 1e-12) std::cout << "  Result: SUCCESS" << std::endl;
    else std::cout << "  Result: FAILURE" << std::endl;

    fftw_destroy_plan(p_fwd);
    fftw_destroy_plan(p_bwd);
}

/**
 * TEST 6: Spectral laplacian operator in 2-d
 *
 * Inspired by the discussion at:
 *      https://fortran-lang.discourse.group/t/spectral-derivative-in-2-d/6978
 */
void test_2d_simple_laplacian() {
    std::cout << "Running Test 6 (Simplified): 2D Laplacian of sin(2pi*x)..." << std::endl;

    const int nx = 4, ny = 4;
    int n[] = {nx, ny};
    int pad = 1;
    int inembed[] = {nx + 2 * pad, ny + 2 * pad};
    int onembed[] = {nx + 2 * pad, ny / 2 + 1};

    std::vector<double> phys_in(inembed[0] * inembed[1], 0.0);
    std::vector<double> phys_out(inembed[0] * inembed[1], 0.0);
    std::vector<std::complex<double>> spec(onembed[0] * onembed[1], 0.0);

    double dx = 1.0 / nx;
    double TWO_PI = 2.0 * M_PI;

    // 1. Initialize f(x,y) = sin(2*pi*x)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            phys_in[(i + pad) * inembed[1] + (j + pad)] = std::sin(TWO_PI * i * dx);
        }
    }

    double* in_ptr = &phys_in[pad * inembed[1] + pad];
    double* out_ptr = &phys_out[pad * inembed[1] + pad];
    fftw_complex* spec_ptr = reinterpret_cast<fftw_complex*>(&spec[pad * onembed[1]]);

    fftw_plan p_fwd = fftw_plan_many_dft_r2c(2, n, 1, in_ptr, inembed, 1, 0, spec_ptr, onembed, 1, 0, FFTW_ESTIMATE);
    fftw_plan p_bwd = fftw_plan_many_dft_c2r(2, n, 1, spec_ptr, onembed, 1, 0, out_ptr, inembed, 1, 0, FFTW_ESTIMATE);

    fftw_execute(p_fwd);

    // 2. Apply Laplacian Filter
    for (int i = 0; i < nx; ++i) {
        // Correct kx handling for nx=4: indices 0, 1, 2, 3 -> kx = 0, 1, -2, -1
        int kx = i;
        if (kx > nx / 2) kx -= nx;

        for (int j = 0; j < (ny / 2 + 1); ++j) {
            int ky = j;
            double k2 = static_cast<double>(kx * kx + ky * ky);

            // CRITICAL: FFTW Nyquist Frequency (nx/2)
            // For a Laplacian (even derivative), kx = nx/2 is fine.
            // For odd derivatives, the Nyquist imaginary part must be zeroed.
            spec[(i + pad) * onembed[1] + j] *= -k2 * TWO_PI * TWO_PI;
        }
    }

    fftw_execute(p_bwd);

    double max_err = 0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double expected = -TWO_PI * TWO_PI * std::sin(TWO_PI * i * dx);
            double actual = phys_out[(i + pad) * inembed[1] + (j + pad)] / (nx * ny);
            max_err = std::max(max_err, std::abs(actual - expected));
        }
    }

    std::cout << "  Max Error (Simple 2D Laplacian): " << max_err << std::endl;
    if (max_err < 1e-12) std::cout << "  Result: SUCCESS" << std::endl;
    else std::cout << "  Result: FAILURE" << std::endl;

    fftw_destroy_plan(p_fwd); fftw_destroy_plan(p_bwd);
}


void test_2d_advanced_laplacian() {
    std::cout << "Running Test 7: 2D Advanced many_dft (sin*cos Laplacian)..." << std::endl;

    const int nx = 32, ny = 32;
    const int pad = 2;
    int n[] = {nx, ny};
    int inembed[] = {nx + 2 * pad, ny + 2 * pad};
    int onembed[] = {nx + 2 * pad, ny / 2 + 1};

    std::vector<double> phys_in(inembed[0] * inembed[1], 0.0);
    std::vector<double> phys_out(phys_in.size(), 0.0);
    std::vector<std::complex<double>> spec(onembed[0] * onembed[1], 0.0);

    const double TWO_PI = 8.0 * std::atan(1.0);

    // 1. Initialize f(x,y) = sin(2pi*x)*cos(2pi*y)
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = static_cast<double>(i) / nx;
            double y = static_cast<double>(j) / ny;
            phys_in[(i + pad) * inembed[1] + (j + pad)] = std::sin(TWO_PI * x) * std::cos(TWO_PI * y);
        }
    }

    double* in_ptr = &phys_in[pad * inembed[1] + pad];
    double* out_ptr = &phys_out[pad * inembed[1] + pad];
    fftw_complex* spec_ptr = reinterpret_cast<fftw_complex*>(&spec[pad * onembed[1]]);

    fftw_plan p_fwd = fftw_plan_many_dft_r2c(2, n, 1, in_ptr, inembed, 1, 0, spec_ptr, onembed, 1, 0, FFTW_ESTIMATE);
    fftw_plan p_bwd = fftw_plan_many_dft_c2r(2, n, 1, spec_ptr, onembed, 1, 0, out_ptr, inembed, 1, 0, FFTW_ESTIMATE);

    fftw_execute(p_fwd);

    // 2. Apply Laplacian Filter: multiplier = -(kx^2 + ky^2) * (2*pi)^2
    for (int i = 0; i < nx; ++i) {
        int kx = (i <= nx / 2) ? i : i - nx;
        for (int j = 0; j < (ny / 2 + 1); ++j) {
            int ky = j;
            double k2 = static_cast<double>(kx * kx + ky * ky);
            double multiplier = -k2 * TWO_PI * TWO_PI;

            // Maintain Hermitian symmetry by zeroing Nyquist
            if (i == nx/2 || j == ny/2) multiplier = 0.0;

            spec[(i + pad) * onembed[1] + j] *= multiplier;
        }
    }

    fftw_execute(p_bwd);

    // 3. Verification against SymPy: -8 * pi^2 * sin(2*pi*x) * cos(2*pi*y)
    double max_err = 0;
    double sum_sq_err = 0;
    const double total_norm = 1.0 / (static_cast<double>(nx) * ny);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = static_cast<double>(i) / nx;
            double y = static_cast<double>(j) / ny;

            double expected = -8.0 * std::pow(M_PI, 2) * std::sin(TWO_PI * x) * std::cos(TWO_PI * y);
            double actual = phys_out[(i + pad) * inembed[1] + (j + pad)] * total_norm;

            double err = actual - expected;
            sum_sq_err += err * err;
            max_err = std::max(max_err, std::abs(err));
        }
    }

    double l2_err = std::sqrt(sum_sq_err / (nx * ny));
    std::cout << "  L2 Error: " << l2_err << std::endl;
    std::cout << "  Max Error: " << max_err << std::endl;
    // Tolerance 1e-11 is reasonable for a 2nd derivative operator
    if (max_err < 1e-11) std::cout << "  Result: SUCCESS" << std::endl;
    else std::cout << "  Result: FAILURE" << std::endl;

    fftw_destroy_plan(p_fwd); fftw_destroy_plan(p_bwd);
}


void test_2d_advanced_laplacian_complex() {
    std::cout << "Running Test 8: 2D Advanced many_dft (Complex Nested Laplacian)..." << std::endl;

    const int nx = 64; // Increased resolution to better capture nested harmonics
    const int ny = 64;
    const int pad = 2;
    int n[] = {nx, ny};
    int inembed[] = {nx + 2 * pad, ny + 2 * pad};
    int onembed[] = {nx + 2 * pad, ny / 2 + 1};

    std::vector<double> phys_in(inembed[0] * inembed[1], 0.0);
    std::vector<double> phys_out(phys_in.size(), 0.0);
    std::vector<std::complex<double>> spec(onembed[0] * onembed[1], 0.0);

    const double TWO_PI = 8.0 * std::atan(1.0);

    // 1. Initialize f = sin(2*cos(2pi*x)) * cos(3*sin(2pi*y))
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = static_cast<double>(i) / nx;
            double y = static_cast<double>(j) / ny;
            phys_in[(i + pad) * inembed[1] + (j + pad)] =
                std::sin(2.0 * std::cos(TWO_PI * x)) * std::cos(3.0 * std::sin(TWO_PI * y));
        }
    }

    double* in_ptr = &phys_in[pad * inembed[1] + pad];
    double* out_ptr = &phys_out[pad * inembed[1] + pad];
    fftw_complex* spec_ptr = reinterpret_cast<fftw_complex*>(&spec[pad * onembed[1]]);

    fftw_plan p_fwd = fftw_plan_many_dft_r2c(2, n, 1, in_ptr, inembed, 1, 0, spec_ptr, onembed, 1, 0, FFTW_ESTIMATE);
    fftw_plan p_bwd = fftw_plan_many_dft_c2r(2, n, 1, spec_ptr, onembed, 1, 0, out_ptr, inembed, 1, 0, FFTW_ESTIMATE);

    fftw_execute(p_fwd);

    // 2. Apply Laplacian Filter
    for (int i = 0; i < nx; ++i) {
        int kx = (i <= nx / 2) ? i : i - nx;
        for (int j = 0; j < (ny / 2 + 1); ++j) {
            int ky = j;
            double k2 = static_cast<double>(kx * kx + ky * ky);
            double multiplier = -k2 * TWO_PI * TWO_PI;

            // Zero out Nyquist to ensure Hermitian symmetry for the c2r transform
            if (i == nx/2 || j == ny/2) multiplier = 0.0;

            spec[(i + pad) * onembed[1] + j] *= multiplier;
        }
    }

    fftw_execute(p_bwd);

    // 3. Verification using the EXACT SymPy Expression
    double max_err = 0;
    double sum_sq_err = 0;
    const double total_norm = 1.0 / (static_cast<double>(nx) * ny);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = static_cast<double>(i) / nx;
            double y = static_cast<double>(j) / ny;

            // Direct implementation of SymPy's C++ output
            double term1 = -8.0 * (2.0 * std::pow(std::sin(TWO_PI * x), 2) * std::sin(2.0 * std::cos(TWO_PI * x))
                           + std::cos(TWO_PI * x) * std::cos(2.0 * std::cos(TWO_PI * x))) * std::cos(3.0 * std::sin(TWO_PI * y));

            double term2 = 12.0 * (std::sin(TWO_PI * y) * std::sin(3.0 * std::sin(TWO_PI * y))
                           - 3.0 * std::pow(std::cos(TWO_PI * y), 2) * std::cos(3.0 * std::sin(TWO_PI * y))) * std::sin(2.0 * std::cos(TWO_PI * x));

            double expected = std::pow(M_PI, 2) * (term1 + term2);
            double actual = phys_out[(i + pad) * inembed[1] + (j + pad)] * total_norm;

            double err = actual - expected;
            sum_sq_err += err * err;
            max_err = std::max(max_err, std::abs(err));
        }
    }

    double l2_err = std::sqrt(sum_sq_err / (nx * ny));
    std::cout << "  L2 Error: " << l2_err << std::endl;
    std::cout << "  Max Error: " << max_err << std::endl;
    // Set tolerance to 5e-11 to account for spectral truncation of the nested terms
    if (max_err < 5e-11) std::cout << "  Result: SUCCESS" << std::endl;
    else std::cout << "  Result: FAILURE" << std::endl;

    fftw_destroy_plan(p_fwd); fftw_destroy_plan(p_bwd);
}


void test_round_trip_precision() {
    std::cout << "Test 9: Round-Trip Precision (f -> FFT -> IFFT -> f)..." << std::endl;

    const int n_size = 32;
    std::vector<double> in(n_size * n_size);
    std::vector<std::complex<double>> spec(n_size * (n_size/2 + 1));
    std::vector<double> out(n_size * n_size);

    for (int i = 0; i < n_size * n_size; ++i) in[i] = 1.0 + static_cast<double>(i) / 1000.0;

    fftw_plan p_fwd = fftw_plan_dft_r2c_2d(n_size, n_size, in.data(),
                                           reinterpret_cast<fftw_complex*>(spec.data()), FFTW_ESTIMATE);
    fftw_plan p_bwd = fftw_plan_dft_c2r_2d(n_size, n_size, reinterpret_cast<fftw_complex*>(spec.data()),
                                           out.data(), FFTW_ESTIMATE);

    fftw_execute(p_fwd);
    fftw_execute(p_bwd);

    double max_err = 0;
    for (int i = 0; i < n_size * n_size; ++i) {
        double actual = out[i] / (n_size * n_size);
        max_err = std::max(max_err, std::abs(actual - in[i]));
    }

    std::cout << "  Round-trip Max Error: " << max_err << std::endl;
    if (max_err > 1e-13) {
        std::cout << "  CONCLUSION: Library lacks true double-precision accuracy." << std::endl;
    }

    fftw_destroy_plan(p_fwd); fftw_destroy_plan(p_bwd);
}

int main() {
    test_2d_batch_embedded();
    test_stride_interleaved();
    test_inplace_advanced();
    test_2d_tight_batch();
    test_2d_advanced_differentiation();
    test_2d_simple_laplacian();
    test_2d_advanced_laplacian();
    test_2d_advanced_laplacian_complex();
    test_round_trip_precision();
    return 0;
}
