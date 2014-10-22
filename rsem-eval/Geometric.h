#ifndef GEOMETRIC_H_
#define GEOMETRIC_H_

#include<cmath>
#include<cstdio>
#include<cassert>

class Geometric {
 public:
  Geometric(double mean) {
    assert(mean > 0);
    p = 1.0 / mean;
  }

  ~Geometric() {}

  double getLogProb(int value) {
    return log(p) + (value - 1) * log(1.0 - p);
  }

 private:
  double p;
};

#endif
