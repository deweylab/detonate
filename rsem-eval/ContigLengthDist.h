/*
 * ContigLengthDist.h
 */

#ifndef CONTIGLENGTHDIST_H_
#define CONTIGLENGTHDIST_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<algorithm>

#include "boost/math/special_functions/log1p.hpp"
#include "boost/math/special_functions/expm1.hpp"
#include "boost/math/distributions/negative_binomial.hpp"


class ContigLengthDist {
public:
	ContigLengthDist(double, double, int, int, int);
	~ContigLengthDist() {}
	
	double calcLogPrior(double, int);
	double calcLogNorm(double, int);
	double calcLogF(double, int);

private:
	static const double a0;
	static const double x1;
	static const double INF; // infinite

	int L, w;
	double mu;
	boost::math::negative_binomial nbdist;
	std::vector<double> nbpdf, nbpsum; // negative binomial pdf and partial sum

	double log1pexp(double x) {
		return (x >= x1 ? x : boost::math::log1p(exp(x)));
	}

	// a should be smaller than 0
	double log1mexp(double a) {
		return (a >= -a0 ? log(-boost::math::expm1(a)) : boost::math::log1p(-exp(a)));

	}

	double log_add(double a, double b) {
		if (a < b) { double tmp = a; a = b; b = tmp; }
		return (b > -INF ? a + log1pexp(b - a) : a);
	}

	// exp(a) - exp(b),  a must > b
	double log_minus(double a, double b) {
		assert(a > b);
		return (b > -INF ? a + log1mexp(b - a) : a);
	}

	double calcLogNum(double, int); // numerator
	double calcLogDenom(double, int); // denominator

	bool calcLogFSlow(double, int, double&);
	bool calcLogFFast(double, int, double&);
	bool calcLogFFast2(double, int, double&);
};

const double ContigLengthDist::a0 = log(2.0);
const double ContigLengthDist::x1 = 33.3;
const double ContigLengthDist::INF = 1e300;

ContigLengthDist::ContigLengthDist(double nb_r, double nb_p, int L, int w, int maxl) : L(L), w(w), nbdist(nb_r, nb_p) {
	int upper = std::max(L + L - w , maxl + 2 * (L - w) - 2);
	nbpdf.assign(upper + 1, 0.0);
	nbpsum.assign(upper + 1, 0.0);
	mu = boost::math::mean(nbdist);
	for (int t = 0; t <= upper; t++) {
		nbpdf[t] = boost::math::pdf(nbdist, t);
		nbpsum[t] = (t > 0 ? nbpsum[t - 1] : 0.0) + t * nbpdf[t];
	}
}

double ContigLengthDist::calcLogPrior(double lambda, int l) {
	if (l <= L) return 0.0;
	double p = exp(-lambda);
	return calcLogNum(p, l) - calcLogDenom(p, l);
}

double ContigLengthDist::calcLogNorm(double lambda, int l) {
	return log1mexp(-lambda);
}

bool ContigLengthDist::calcLogFSlow(double lambda, int l, double& result) {
	assert(l >= L);
	double log_p = -lambda, log_q = log1mexp(-lambda);
	std::vector<double> f(l + 5, -INF);

	f[L] = 0.0;
	for (int i = L + 1; i <= l; i++)
		for (int j = std::max(w, L + L - i); j < L; j++)
			f[i] = log_add(f[i], f[i - (L - j)] + (L - j - 1) * log_p + log_q);

	result = f[l];
	return true;
}

bool ContigLengthDist::calcLogFFast(double lambda, int l, double& result) {
	assert(l >= L);
	double log_p = -lambda, log_q = log1mexp(-lambda);
	std::vector<double> f(l + 5, -INF); // l + 5 -> a bit extra so that even l == L, f[L + 1] still exists
	int id;
	double subtrahend;

	f[L] = 0.0;
	f[L + 1] = log_q;
	for (int i = L + 2; i <= l; i++) {
		id = i - 1 - (L - w);
		if (id >= L) {
		  subtrahend = f[id] + (L - w) * log_p + log_q;

		  if (f[i - 1] <= subtrahend) { 
		    result = -INF;
		    return false;
		  }

		  f[i] = log_minus(f[i - 1], subtrahend);
		}
		else f[i] = f[i - 1];

	}

	result = f[l];
	return true;
}

bool ContigLengthDist::calcLogFFast2(double lambda, int l, double& result) {
  std::vector<double> c;
  double p = exp(-lambda), q = 1.0 - p;
  double pprod, pconst; // pprod, partial product; pr, probability const in the formula
  double logf = 0.0;

  assert(l >= L);

  if (l == L) logf = 0.0;
  else if (l == L + 1) logf = log(q);
  else {
    c.assign(l + 1, 0.0);
    c[L] = 1.0; c[L + 1] = q;

    pprod = q; 
    pconst = pow(p, L - w) * q;
    
    logf = log(q);
    for (int i = L + 2; i <= l; i++) {
      if (i <= L + (L - w)) c[i] = 1.0;
      else {
	pprod /= c[i - 1 - (L - w)];

	if (pprod <= pconst) {
	  result = -INF;
	  return false;
	}

	c[i] = 1.0 - pconst / pprod;
	pprod *= c[i];
      }
      logf += log(c[i]);
    }
  }

  return logf;
}

double ContigLengthDist::calcLogF(double lambda, int l) {
  double result;

  //  if (!calcLogFFast(lambda, l, result)) assert(calcLogFSlow(lambda, l, result));
  assert(calcLogFSlow(lambda, l, result));

  return result;
}

double ContigLengthDist::calcLogNum(double p, int l) {
	double num = 0.0, value;

	for (int t = l; t <= l + L - w; t++) { 
	  value = nbpdf[t] * (t - l + 1) * pow(p, t - l);
	  num += (value < 0.0 ? 0.0 : value);
	}

	for (int t = l + L - w + 1; t <= l + 2 * (L - w) - 2; t++) {
	  value = nbpdf[t] * ((2 * (L - w) + 1 - (t - l)) * pow(p, t - l) + 2 * pow(p, L - w) * (1.0 - pow(p, (t - l) - (L - w))) / (1.0 - p));
	  num += (value < 0.0 ? 0.0 : value);
	}

	value = boost::math::cdf(complement(nbdist, l + 2 * (L - w) - 2)) * (2 * pow(p, L - w) * (1.0 - pow(p, L - w)) / (1.0 - p) - (l + 2 * (L - w) - 1) * pow(p, 2 * (L - w)));
	value += pow(p, 2 * (L - w)) * (mu - nbpsum[l + 2 * (L - w) - 2]);
	num += (value < 0.0 ? 0.0 : value);

	assert(num > 0.0);
	return log(num);
}

double ContigLengthDist::calcLogDenom(double p, int l) {
	double denom = 0.0, value;

	for (int t = L; t <= L + L - w; t++) {
	  value = nbpdf[t] * (1.0 - pow(p, t - L + 1)) / (1.0 - p);
	  denom += (value < 0.0 ? 0.0 : value);
	}

        value = boost::math::cdf(complement(nbdist, L + L - w)) * ((1.0 - pow(p, L - w + 1)) / (1.0 - p) - (L + L - w) * pow(p, L - w));
	value += (mu - nbpsum[L + L - w]) * pow(p, L - w);
	denom += (value < 0.0 ? 0.0 : value);

	assert(denom > 0.0);
	return log(denom);
}

#endif /* CONTIGLENGTHDIST_H_ */
