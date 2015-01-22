#include "gtest/gtest.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <numeric>

#include <fftw3.h>

TEST(fftw, one_dim) {
	int n = 10;
	std::valarray<double> in(n);
	std::valarray<std::complex<double>> tmp(n / 2 + 1);
	std::valarray<double> out(n);

	for (int i = 0; i < n; ++i)
		in[i] = i;

	fftw_plan forward = fftw_plan_dft_r2c_1d(n,
			reinterpret_cast<double*>(&in[0]),
			reinterpret_cast<fftw_complex*>(&tmp[0]), 0);
	fftw_plan backward = fftw_plan_dft_c2r_1d(n,
			reinterpret_cast<fftw_complex*>(&tmp[0]),
			reinterpret_cast<double*>(&out[0]), 0);

	fftw_execute(forward);
	tmp /= sqrt(n);

	fftw_execute(backward);
	out /= sqrt(n);

	for (int i = 0; i < in.size(); ++i)
		ASSERT_NEAR(in[i], out[i], 1e-12);
}

TEST(fftw, two_dim) {
	int n_x = 10;
	int n_y = 10;
	std::valarray<double> in(n_x * n_y);
	std::valarray<std::complex<double>> tmp(n_x * (n_y / 2 + 1));
	std::valarray<double> out(n_x * n_y);

	for (int i = 0; i < n_x * n_y; ++i)
		in[i] = i;

	fftw_plan forward = fftw_plan_dft_r2c_2d(n_x, n_y,
			reinterpret_cast<double*>(&in[0]),
			reinterpret_cast<fftw_complex*>(&tmp[0]), FFTW_ESTIMATE);

	fftw_plan backward = fftw_plan_dft_c2r_2d(n_x, n_y,
			reinterpret_cast<fftw_complex*>(&tmp[0]),
			reinterpret_cast<double*>(&out[0]), FFTW_ESTIMATE);

	fftw_execute(forward);
	tmp /= sqrt(n_x * n_y);

	fftw_execute(backward);
	out /= sqrt(n_x * n_y);

	for (int i = 0; i < in.size(); ++i)
		ASSERT_NEAR(in[i], out[i], 1e-12);
}

TEST(fftw, three_dim) {
	int n_x = 10;
	int n_y = 10;
	int n_z = 10;
	std::valarray<double> in(n_x * n_y * n_z);
	std::valarray<std::complex<double>> tmp(n_x * n_y * (n_z / 2 + 1));
	std::valarray<double> out(n_x * n_y * n_z);

	for (int i = 0; i < n_x * n_y * n_z; ++i)
		in[i] = i;

	fftw_plan forward = fftw_plan_dft_r2c_3d(n_x, n_y, n_z,
			reinterpret_cast<double*>(&in[0]),
			reinterpret_cast<fftw_complex*>(&tmp[0]), FFTW_ESTIMATE);

	fftw_plan backward = fftw_plan_dft_c2r_3d(n_x, n_y, n_z,
			reinterpret_cast<fftw_complex*>(&tmp[0]),
			reinterpret_cast<double*>(&out[0]), FFTW_ESTIMATE);

	fftw_execute(forward);
	tmp /= sqrt(n_x * n_y * n_z);

	fftw_execute(backward);
	out /= sqrt(n_x * n_y * n_z);

	for (int i = 0; i < in.size(); ++i)
		ASSERT_NEAR(in[i], out[i], 1e-12);
}
