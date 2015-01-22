#pragma once

#include <vector>
#include <valarray>

std::vector<std::valarray<double>> rk2(const std::valarray<double> &q0,
		const std::valarray<double> &wA, const std::valarray<double> &wB,
		double f, int n_t);


