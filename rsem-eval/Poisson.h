#ifndef POISSON_H_
#define POISSON_H_

#include<cmath>
#include<cstdio>
#include<cassert>
#include<vector>

class Poisson {
 public:
  Poisson(double lambda) {
    assert(lambda > 0);
    this->lambda = lambda;
    log_factorial.clear();
    log_factorial.push_back(0); // log 0!
  }

  ~Poisson() {}

  void setLambda(double lambda) {
    assert(lambda > 0);
    this->lambda = lambda;
  }

  double getLogProb(int value) {
    return value * log(lambda) - lambda - getLogFactorial(value); 
  }

 private:
  double lambda;
  std::vector<double> log_factorial;

  double getLogFactorial(int k) {
    if ((int)log_factorial.size() <= k) {
      for (int i = log_factorial.size(); i <= k; i++) {
	log_factorial.push_back(log_factorial[i - 1] + log(1.0 * i));
      }
    }
    return log_factorial[k];
  }

};

#endif
