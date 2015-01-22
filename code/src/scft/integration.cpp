#include "scft/integration.hpp"

#include <cassert>

/**
 * Numerische Quadratur mittels summierter Trapezformel.
 * @param f Funktionswerte
 * @param from Linkes Ende des Intervalls
 * @param to Rechtes Ende des Intervalls
 * @return Integral von f über [from, to]
 */
double integration::trapez(const std::valarray<double> &f, double from,
		double to) {
	double result = 0;
	unsigned int n = f.size();
	unsigned int m = n - 1;

	assert(m % 2 == 0);
	double h = (to - from) / m;

	result += h * f[std::slice(1, m - 1, 1)].sum();
	result += (h / 2) * (f[0] + f[m]);

	return result;
}

/**
 * Numerische Quadratur mittels summierter Simpsonformel.
 * @param f Funktionswerte
 * @param from Linkes Ende des Intervalls
 * @param to Rechtes Ende des Intervalls
 * @return Integral von f über [from, to]
 */
double integration::simpson(const std::valarray<double> &f, double from,
		double to) {
	double result = 0;
	unsigned int n = f.size();
	unsigned int m = n - 1;

	assert(m % 2 == 0);
	double h = (to - from) / m;

	result += h * 2 * f[std::slice(1, m - 1, 1)].sum();
	result += h * 2 * f[std::slice(1, (m - 1) / 2.0 + 0.5, 2)].sum();
	result += h * (f[0] + f[m]);
	result /= 3;

	return result;
}
