#include <fftw3.h>

#include "scft/rk2.hpp"

#include <valarray>
#include <complex>
#include <vector>

std::vector<std::valarray<double>> rk2(const std::valarray<double> &q0,
		const std::valarray<double> &wA, const std::valarray<double> &wB,
		double f, int n_t) {
	std::vector<std::valarray<double>> result;

	result.push_back(q0);
	for (int iter = 1; iter <= n_t; ++iter) {
		const auto last = result.back();

		size_t n = last.size();
		std::valarray<std::complex<double>> tmp(n);
		std::valarray<double> out(n);

		double t = iter / double(n_t);
		double delta_t = 1.0 / double(n_t);
		double tmp1 = (1 - t <= f) ? 1 : 0;
		double tmp2 = (!(1 - t <= f)) ? 1 : 0;
		std::valarray<double> w = wA * tmp1 + wB * tmp2;
		std::valarray<double> W = std::exp(w);

		auto next = last;
		next = W * next;

		fftw_plan forward = fftw_plan_dft_r2c_1d(n,
				reinterpret_cast<double*>(&next[0]),
				reinterpret_cast<fftw_complex*>(&tmp[0]), 0);

		fftw_execute(forward);
		tmp /= sqrt(n);

//		tmp = std::exp(...);

		fftw_plan backward = fftw_plan_dft_c2r_1d(n,
				reinterpret_cast<fftw_complex*>(&tmp[0]),
				reinterpret_cast<double*>(&next[0]), 0);

		fftw_execute(backward);
		next /= sqrt(n);

		next = W * next;

		result.push_back(next);
	}

	return result;
}
