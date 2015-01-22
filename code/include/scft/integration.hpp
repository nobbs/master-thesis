#pragma once

#include <valarray>

namespace integration {

double trapez(const std::valarray<double> &f, double from, double to);
double simpson(const std::valarray<double> &f, double from, double to);

}
