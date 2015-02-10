#ifndef NOISEQPROFILE_H_
#define NOISEQPROFILE_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<string>
#include<cassert>

#include "utils.h"
#include "RefSeq.h"
#include "simul.h"

class NoiseQProfile {
public:
	NoiseQProfile() {
		logp = 0.0;
		memset(c, 0, sizeof(c));
		memset(p, 0, sizeof(p));
		memset(log_probs, 0, sizeof(log_probs));
	}

	void init();
	void updateC(const std::string&, const std::string&);
	void update(const std::string&, const std::string&, double frac);
	void finish();
	void calcInitParams();

	double getLogProb(const std::string&, const std::string&);
	double getLogP() { return logp; }

	void collect(const NoiseQProfile&);

	void read(FILE*);
	void write(FILE*);

	void startSimulation();
	std::string simulate(simul*, int, const std::string&);
	void finishSimulation();

private:
	static const int NCODES = 5; // number of possible codes
	static const int SIZE = 100;
	static const double prior[NCODES]; // prior for each base

	double logp; //log prob;
	double c[SIZE][NCODES]; //counts in N0;
	double p[SIZE][NCODES]; //p[q][c] = p(c|q)
	double log_probs[SIZE][NCODES]; // log_probs[q]c[] = log(p(c|q))

	int c2q(char c) { assert(c >= 33 && c <= 126); return c - 33; }

	double (*pc)[NCODES]; // for simulation
};

void NoiseQProfile::init() {
	memset(p, 0, sizeof(p));
}

void NoiseQProfile::updateC(const std::string& readseq, const std::string& qual) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		++c[c2q(qual[i])][get_base_id(readseq[i])];
	}
}

void NoiseQProfile::update(const std::string& readseq, const std::string& qual, double frac) {
	int len = readseq.size();
	for (int i = 0; i < len; i++) {
		p[c2q(qual[i])][get_base_id(readseq[i])] += frac;
	}
}

void NoiseQProfile::finish() {
	double sum;

	logp = 0.0;
	for (int i = 0; i < SIZE; i++) {
		sum = 0.0;
		for (int j = 0; j < NCODES; j++) {
		  p[i][j] += c[i][j];
		  sum += p[i][j];
		}
		if (sum > EPSILON) {
		  for (int j = 0; j < NCODES; j++) {
		    p[i][j] /= sum;
		    log_probs[i][j] = Log(p[i][j]);
		    if (c[i][j] > EPSILON) { logp += c[i][j] * log_probs[i][j]; }
		  }
		}
		else {
		  for (int j = 0; j < NCODES; j++) {
		    p[i][j] = 0.0;
		    log_probs[i][j] = LOGZERO; 
		  }
		}
	}
}

void NoiseQProfile::calcInitParams() {
	double sum;

	logp = 0.0;
	for (int i = 0; i < SIZE; i++) {
		sum = 0.0;
		for (int j = 0; j < NCODES; j++) {
		  p[i][j] = c[i][j] + 1.0; // 1.0 pseudo-count
		  sum += p[i][j];
		}
		for (int j = 0; j < NCODES; j++) {
		  p[i][j] /= sum;
		  log_probs[i][j] = Log(p[i][j]);
		  if (c[i][j] > EPSILON) { logp += c[i][j] * log_probs[i][j]; }
		}
	}
}

double NoiseQProfile::getLogProb(const std::string& readseq, const std::string& qual) {
	double log_prob = 0.0;
	int len = readseq.size();

	for (int i = 0; i < len; i++) {
		log_prob += log_probs[c2q(qual[i])][get_base_id(readseq[i])];
	}

	return log_prob;
}

void NoiseQProfile::collect(const NoiseQProfile& o) {
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < NCODES; j++)
			p[i][j] += o.p[i][j];
	}
}

//If read from file, assume do not need to estimate from data
void NoiseQProfile::read(FILE *fi) {
	int tmp_size, tmp_ncodes;

	memset(c, 0, sizeof(c));

	assert(fscanf(fi, "%d %d", &tmp_size, &tmp_ncodes) == 2);
	assert(tmp_size == SIZE && tmp_ncodes == NCODES);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < NCODES; j++) 
		  assert(fscanf(fi, "%lf", &p[i][j]) == 1);
	}
}

void NoiseQProfile::write(FILE *fo) {
	fprintf(fo, "%d %d\n", SIZE, NCODES);
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < NCODES - 1; j++) { fprintf(fo, "%.10g ", p[i][j]); }
		fprintf(fo, "%.10g\n", p[i][NCODES - 1]);
	}
}

void NoiseQProfile::startSimulation() {
	pc = new double[SIZE][NCODES];

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < NCODES; j++) {
			pc[i][j] = p[i][j];
			if (j > 0) pc[i][j] += pc[i][j - 1];
		}
}

std::string NoiseQProfile::simulate(simul* sampler, int len, const std::string& qual) {
	std::string readseq = "";

	for (int i = 0; i < len; i++) {
		readseq.push_back(getCharacter(sampler->sample(pc[c2q(qual[i])], NCODES)));
	}
	return readseq;
}

void NoiseQProfile::finishSimulation() {
	delete[] pc;
}

#endif /* NOISEQPROFILE_H_ */
