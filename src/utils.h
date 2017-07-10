#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp ;
using namespace arma ;
using namespace std  ;

/*double gamma(double x); */

double lgamma(double x);

double digamma(double x);

double trigamma(double x);

vec gamma_vec(vec x);

vec lgamma_vec(vec x);

vec digamma_vec(vec x);

vec trigamma_vec(vec x);

#endif
