#include "gtest/gtest.h"

#include <numeric>

#include "scft/integration.hpp"

TEST(trapez, constant) {
	int size = 100001;
	std::valarray<double> f(size);
	double c = 2.5;
	f = c;

	ASSERT_NEAR(c, integration::trapez(f, 0, 1), 1e-10);
}

TEST(trapez, linear) {
	int size = 100001;
	std::valarray<double> linear(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	ASSERT_NEAR(0.5, integration::trapez(linear, 0, 1), 1e-10);
}

TEST(trapez, quadratic) {
	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);
	f = linear * linear;

	ASSERT_NEAR(1 / 3.0, integration::trapez(f, 0, 1), 1e-4);
}

TEST(trapez, sin) {
	const double pi = std::acos(-1);

	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::sin(linear * pi);
	ASSERT_NEAR(2, integration::trapez(f, 0, pi), 1e-3);

	f = std::sin(linear * 2 * pi);
	ASSERT_NEAR(0, integration::trapez(f, 0, 2 * pi), 1e-10);
}

TEST(trapez, cos) {
	const double pi = std::acos(-1);

	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::cos(linear * pi);
	ASSERT_NEAR(0, integration::trapez(f, 0, pi), 1e-10);

	f = std::cos(linear * 2 * pi);
	ASSERT_NEAR(0, integration::trapez(f, 0, 2 * pi), 1e-10);
}

TEST(trapez, exp) {
	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::exp(linear);
	ASSERT_NEAR(std::exp(1) - 1, integration::trapez(f, 0, 1), 1e-4);
}

TEST(simpson, constant) {
	int size = 100001;
	std::valarray<double> f(size);
	double c = 2.5;
	f = c;

	ASSERT_NEAR(c, integration::simpson(f, 0, 1), 1e-10);
}

TEST(simpson, linear) {
	int size = 100001;
	std::valarray<double> linear(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	ASSERT_NEAR(0.5, integration::simpson(linear, 0, 1), 1e-10);
}

TEST(simpson, quadratic) {
	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);
	f = linear * linear;

	ASSERT_NEAR(1 / 3.0, integration::simpson(f, 0, 1), 1e-4);
}

TEST(simpson, sin) {
	const double pi = std::acos(-1);

	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::sin(linear * pi);
	ASSERT_NEAR(2, integration::simpson(f, 0, pi), 1e-3);

	f = std::sin(linear * 2 * pi);
	ASSERT_NEAR(0, integration::simpson(f, 0, 2 * pi), 1e-10);
}

TEST(simpson, cos) {
	const double pi = std::acos(-1);

	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::cos(linear * pi);
	ASSERT_NEAR(0, integration::simpson(f, 0, pi), 1e-10);

	f = std::cos(linear * 2 * pi);
	ASSERT_NEAR(0, integration::simpson(f, 0, 2 * pi), 1e-10);
}

TEST(simpson, exp) {
	int size = 100001;
	std::valarray<double> linear(size);
	std::valarray<double> f(size);
	for (int i = 0; i < size; ++i)
		linear[i] = i / double(size - 1);

	f = std::exp(linear);
	ASSERT_NEAR(std::exp(1) - 1, integration::simpson(f, 0, 1), 1e-4);
}

