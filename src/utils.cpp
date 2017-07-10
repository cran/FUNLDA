#include "utils.h"

/*double gamma(double x)
{
  double y = Rf_gammafn(x);
  return y;
}*/

double lgamma(double x)
{
  double y = Rf_lgammafn(x);
  return y;
}

double digamma(double x)
{
  double y = Rf_digamma(x);
  return y;
}

double trigamma(double x)
{
  double y = Rf_trigamma(x);
  return y;
}

vec gamma_vec(vec x)
{
  int n = x.n_elem;
  vec y = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    y(i) = Rf_gammafn(x(i));
  }
  return y;
}

vec lgamma_vec(vec x)
{
  int n = x.n_elem;
  vec y = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    y(i) = Rf_lgammafn(x(i));
  }
  return(y);
}

vec digamma_vec(vec x)
{
  int n = x.n_elem;
  vec y = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    y(i) = Rf_digamma(x(i));
  }
  return y;
}

vec trigamma_vec(vec x)
{
  int n = x.n_elem;
  vec y = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    y(i) = Rf_trigamma(x(i));
  }
  return y;
}

