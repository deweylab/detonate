#ifndef NOISEPROFILE_H_
#define NOISEPROFILE_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<string>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "simul.h"

class NoiseProfile {
public:
	NoiseProfile() {
		logp = 0.0;
		memset(c, 0, sizeof(c));
		memset(p, 0, sizeof(p));
		memset(log_probs, 0, sizeof(log_probs));
	}

	void init();
	void updateC(const std::string&);
	void update(const std::string&, double frac);
	void finish();
	void calcInitParams();

	double getLogProb(const std::string&);

	double getLogP() { return logp; }

	void collect(const NoiseProfile&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int);
	void finishSimulation();

private:
	static const int NCODES = 5;
	static const double prior[NCODES];

	double logp;
	double c[NCODES]; // counts in N0;
	double p[NCODES];
	
	double log_probs[NCODES]; // log_probs[i] = log(p[i])

	double *pc; // for simulation
};

void NoiseProfile::init() {
	memset(p, 0, sizeof(p));
}

void NoiseProfile::updateC(const std::string& readseq) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		++c[get_base_id(readseq[i])];
	}
}

void NoiseProfile::update(const std::string& readseq, double frac) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		p[get_base_id(readseq[i])] += frac;
	}
}

void NoiseProfile::finish() {
	double sum;

	logp = 0.0;
	sum = 0.0;
	for (int i = 0; i < NCODES; i++) {
	  p[i] += c[i];
	  sum += p[i];
	}
	if (sum > EPSILON) {
	  for (int i = 0; i < NCODES; i++) {
	    p[i] /= sum;
	    log_probs[i] = Log(p[i]);
	    if (c[i] > EPSILON) { logp += c[i] * log_probs[i]; }
	  }
	}
	else {
	  for (int i = 0; i < NCODES; i++) {
	    p[i] = 0.0;
	    log_probs[i] = LOGZERO; 
	  }
	}
}

void NoiseProfile::calcInitParams() {
	double sum;

	logp = 0.0;
	sum = 0.0;
	for (int i = 0; i < NCODES; i++) {
	  p[i] = c[i] + 1.0; // 1.0 pseudo count
	  sum += p[i];
	}
	for (int i = 0; i < NCODES; i++) {
		p[i] /= sum;
		log_probs[i] = Log(p[i]);
		if (c[i] > EPSILON) { logp += c[i] * log_probs[i]; }
	}
}

double NoiseProfile::getLogProb(const std::string& readseq) {
	double log_prob = 0.0;
	int len = readseq.size();

	for (int i = 0; i < len; i++) {
		log_prob += log_probs[get_base_id(readseq[i])];
	}

	return log_prob;
}

void NoiseProfile::collect(const NoiseProfile& o) {
	for (int i = 0; i < NCODES; i++)
		p[i] += o.p[i];
}

void NoiseProfile::read(FILE *fi) {
	int tmp_ncodes;

	memset(c, 0, sizeof(c));
	assert(fscanf(fi, "%d", &tmp_ncodes) == 1);
	assert(tmp_ncodes == NCODES);
	for (int i = 0; i < NCODES; i++)
	  assert(fscanf(fi, "%lf", &p[i]) == 1);
}

void NoiseProfile::write(FILE *fo) {
	fprintf(fo, "%d\n", NCODES);
	for (int i = 0; i < NCODES - 1; i++) {
		fprintf(fo, "%.10g ", p[i]);
	}
	fprintf(fo, "%.10g\n", p[NCODES - 1]);
}

void NoiseProfile::startSimulation() {
	pc = new double[NCODES];

	for (int i = 0; i < NCODES; i++) {
		pc[i] = p[i];
		if (i > 0) pc[i] += pc[i - 1];
	}
}

std::string NoiseProfile::simulate(simul* sampler, int len) {
	std::string readseq = "";

	for (int i = 0; i < len; i++) {
		readseq.push_back(getCharacter(sampler->sample(pc, NCODES)));
	}
	return readseq;
}

void NoiseProfile::finishSimulation() {
	delete[] pc;
}

#endif /* NOISEPROFILE_H_ */
