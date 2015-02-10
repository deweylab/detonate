#ifndef SINGLEQMODEL_H_
#define SINGLEQMODEL_H_

#include<cmath>
#include<cstdio>
#include<cassert>
#include<cstring>
#include<string>
#include<algorithm>
#include<sstream>
#include<iostream>
#include<vector>

#include "utils.h"
#include "my_assert.h"
#include "Orientation.h"
#include "LenDist.h"
#include "RSPD.h"
#include "QualDist.h"
#include "QProfile.h"
#include "NoiseQProfile.h"

#include "ModelParams.h"
#include "RefSeq.h"
#include "Refs.h"
#include "SingleReadQ.h"
#include "SingleHit.h"
#include "ReadReader.h"

#include "simul.h"

class SingleQModel {
public:
	SingleQModel(Refs* refs = NULL) {
		this->refs = refs;
		M = (refs != NULL ? refs->getM() : 0);
		memset(N, 0, sizeof(N));
		estRSPD = false;
		needCalcConPrb = true;
		
		ori = new Orientation();
		gld = new LenDist();
		mld = NULL;
		rspd = new RSPD(estRSPD);
		qd = new QualDist();
		qpro = new QProfile();
		nqpro = new NoiseQProfile();

		mean = -1.0; sd = 0.0;

		seedLen = 0;

		loglik_mld_noise = 0.0;
	}

	//If it is not a master node, only init & update can be used!
	SingleQModel(ModelParams& params, bool isMaster = true) {
		M = params.M;
		memcpy(N, params.N, sizeof(params.N));
		refs = params.refs;
		estRSPD = params.estRSPD;
		mean = params.mean; sd = params.sd;
		seedLen = params.seedLen;
		needCalcConPrb = true;

		ori = NULL; gld = NULL; mld = NULL; rspd = NULL; qd = NULL; qpro = NULL; nqpro = NULL;

		loglik_mld_noise = 0.0;

		if (isMaster) {
			gld = new LenDist(params.minL, params.maxL);			
			if (mean >= EPSILON) {
				mld = new LenDist(params.mate_minL, params.mate_maxL);
			}
			if (!estRSPD) { rspd = new RSPD(estRSPD); }
			qd = new QualDist();
		}

		ori = new Orientation(params.probF);
		if (estRSPD) { rspd = new RSPD(estRSPD, params.B); }
		qpro = new QProfile();
		nqpro = new NoiseQProfile();
	}

	~SingleQModel() {
		refs = NULL;
		if (ori != NULL) delete ori;
		if (gld != NULL) delete gld;
		if (mld != NULL) delete mld;
		if (rspd != NULL) delete rspd;
		if (qd != NULL) delete qd;
		if (qpro != NULL) delete qpro;
		if (nqpro != NULL) delete nqpro;
	}

	void estimateFromReads(const char*);

	//if prob is too small, just make it 0
	double getLogConPrb(const SingleReadQ& read, const SingleHit& hit) const {
		double log_prob;
		int sid = hit.getSid();
		RefSeq &ref = refs->getRef(sid);
		int fullLen = ref.getFullLen();
		int totLen = ref.getTotLen();
		int dir = hit.getDir();
		int pos = hit.getPos();
		int readLen = read.getReadLength();
		int fpos = (dir == 0 ? pos : totLen - pos - readLen); // the aligned position reported in SAM file, should be a coordinate in forward strand

		general_assert(fpos >= 0, "The alignment of read " + read.getName() + " to transcript " + itos(sid) + " starts at " + itos(fpos) + \
				" from the forward direction, which should be a non-negative number! " + \
				"It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(fpos + readLen <= totLen,"Read " + read.getName() + " is hung over the end of transcript " + itos(sid) + "! " \
				+ "It is possible that the aligner you use gave different read lengths for a same read in SAM file.");
		general_assert(readLen <= totLen, "Read " + read.getName() + " has length " + itos(readLen) + ", but it is aligned to transcript " \
				+ itos(sid) + ", whose length (" + itos(totLen) + ") is shorter than the read's length!");

		int effL;
		double value;

		if (mld != NULL) {
			int minL = std::max(readLen, gld->getMinL());
			int maxL = std::min(totLen - pos, gld->getMaxL());
			int pfpos; // possible fpos for fragment
			value = 0.0;
			for (int fragLen = minL; fragLen <= maxL; fragLen++) {
				pfpos = (dir == 0 ? pos : totLen - pos - fragLen);
				effL = std::min(fullLen, totLen - fragLen + 1);
				value += gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * mld->getAdjustedProb(readLen, fragLen);
			}
		}
		else {
			effL = std::min(fullLen, totLen - readLen + 1);
			value = gld->getAdjustedProb(readLen, totLen) * rspd->getAdjustedProb(fpos, effL, fullLen);
		}

		log_prob = Log(ori->getProb(dir) * value) + qpro->getLogProb(read.getReadSeq(), read.getQScore(), ref, pos, dir);

		return log_prob;
	}

	double getNoiseLogConPrb(const SingleReadQ& read) {
		double log_prob = Log(mld != NULL ? mld->getProb(read.getReadLength()) : gld->getProb(read.getReadLength()));
		log_prob += nqpro->getLogProb(read.getReadSeq(), read.getQScore());

		return log_prob;
	}

	double getLogP() { return loglik_mld_noise + nqpro->getLogP(); }

	void init();

	void update(const SingleReadQ& read, const SingleHit& hit, double frac) {
		if (frac <= EPSILON) return;

		const RefSeq& ref = refs->getRef(hit.getSid());

		int dir = hit.getDir();
		int pos = hit.getPos();

		if (estRSPD) {
			int fullLen = ref.getFullLen();

			// Only use one strand to estimate RSPD
			if (ori->getProb(0) >= ORIVALVE && dir == 0) {
				rspd->update(pos, fullLen, frac);
			}

			if (ori->getProb(0) < ORIVALVE && dir == 1) {
				int totLen = ref.getTotLen();			  
				int readLen = read.getReadLength();
				
				int pfpos, effL; 

				if (mld != NULL) {
					int minL = std::max(readLen, gld->getMinL());
					int maxL = std::min(totLen - pos, gld->getMaxL());
					double sum = 0.0;
					assert(maxL >= minL);
					std::vector<double> frag_vec(maxL - minL + 1, 0.0);

					for (int fragLen = minL; fragLen <= maxL; fragLen++) {
						pfpos = totLen - pos - fragLen;
						effL = std::min(fullLen, totLen - fragLen + 1);
						frag_vec[fragLen - minL] = gld->getAdjustedProb(fragLen, totLen) * rspd->getAdjustedProb(pfpos, effL, fullLen) * mld->getAdjustedProb(readLen, fragLen);
						sum += frag_vec[fragLen - minL];
					}
					assert(sum >= EPSILON);
					for (int fragLen = minL; fragLen <= maxL; fragLen++) {
						pfpos = totLen - pos - fragLen;
						rspd->update(pfpos, fullLen, frac * (frag_vec[fragLen - minL] / sum));
					}
				}
				else {
					rspd->update(totLen - pos - readLen, fullLen, frac);
				}
			}
		}
		qpro->update(read.getReadSeq(), read.getQScore(), ref, pos, dir, frac);
	}

	void updateNoise(const SingleReadQ& read, double frac) {
		if (frac <= EPSILON) return;

		nqpro->update(read.getReadSeq(), read.getQScore(), frac);
	}

	void finish();

	void collect(const SingleQModel&);

	bool getNeedCalcConPrb() { return needCalcConPrb; }
	void setNeedCalcConPrb(bool value) { needCalcConPrb = value; }

	void read(const char*);
	void write(const char*);

	const LenDist& getGLD() { return *gld; }

	void startSimulation(simul*, const std::vector<double>&);
	bool simulate(READ_INT_TYPE, SingleReadQ&, int&);
	void finishSimulation();

	int getModelType() const { return model_type; }

private:
	static const int model_type = 1;
	static const int read_type = 1;

	int M;
        READ_INT_TYPE N[3];
	Refs *refs;
	double mean, sd;
	int seedLen;
	//double *p1, *p2; P_i' & P_i'';

	bool estRSPD; // true if estimate RSPD
	bool needCalcConPrb; //true need, false does not need

	Orientation *ori;
	LenDist *gld, *mld;
	RSPD *rspd;
	QualDist *qd;
	QProfile *qpro;
	NoiseQProfile *nqpro;

	simul *sampler; // for simulation
	double *theta_cdf; // for simulation

	double *mw; // for masking

	double loglik_mld_noise;
};

void SingleQModel::estimateFromReads(const char* readFN) {
	int s;
	char readFs[2][STRLEN];
	SingleReadQ read;

	assert(mld == NULL);
	gld->init();

	std::vector<READ_INT_TYPE> noise_len_stat(gld->getMaxL() + 1, 0);
	
	for (int i = 0; i < 3; i++)
		if (N[i] > 0) {
			genReadFileNames(readFN, i, read_type, s, readFs);
			ReadReader<SingleReadQ> reader(s, readFs, refs->hasPolyA(), seedLen); // allow calculation of calc_lq() function

			READ_INT_TYPE cnt = 0;
			while (reader.next(read)) {
				gld->update(read.getReadLength(), 1.0);
				qd->update(read.getQScore());
				if (i == 0) { 
				  nqpro->updateC(read.getReadSeq(), read.getQScore()); 
				  ++noise_len_stat[read.getReadLength()];
				}
				
				++cnt;
				if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " READS PROCESSED"<< std::endl; }
			}

			if (verbose) { std::cout<< "estimateFromReads, N"<< i<< " finished."<< std::endl; }
		}
	
        // Calculate mate length related log-likelihood                                                                                                                                                   
        gld->finish();
	loglik_mld_noise = 0.0;
        for (int len = gld->getMinL(); len <= gld->getMaxL(); ++len)
          if (noise_len_stat[len] > 0) loglik_mld_noise += noise_len_stat[len] * log(gld->getProb(len));

	qd->finish();
	nqpro->calcInitParams();
}

void SingleQModel::init() {
	if (estRSPD) rspd->init();
	qpro->init();
	nqpro->init();
}

void SingleQModel::finish() {
	if (estRSPD) rspd->finish();
	qpro->finish();
	nqpro->finish();
	needCalcConPrb = true;
}

void SingleQModel::collect(const SingleQModel& o) {
	if (estRSPD) rspd->collect(*(o.rspd));
	qpro->collect(*(o.qpro));
	nqpro->collect(*(o.nqpro));
}

//Only master node can call
void SingleQModel::read(const char* inpF) {
	int val;
	FILE *fi = fopen(inpF, "r");
	if (fi == NULL) { fprintf(stderr, "Cannot open %s! It may not exist.\n", inpF); exit(-1); }

	assert(fscanf(fi, "%d", &val) == 1);
	assert(val == model_type);

	ori->read(fi);
	gld->read(fi);
	assert(fscanf(fi, "%d", &val) == 1);
	if (val > 0) {
		if (mld == NULL) mld = new LenDist();
		mld->read(fi);
	}
	rspd->read(fi);
	qd->read(fi);
	qpro->read(fi);
	nqpro->read(fi);

	fclose(fi);
}

//Only master node can call. Only be called at EM.cpp
void SingleQModel::write(const char* outF) {
	FILE *fo = fopen(outF, "w");

	fprintf(fo, "%d\n", model_type);
	fprintf(fo, "\n");

	ori->write(fo);  fprintf(fo, "\n");
	gld->write(fo);  fprintf(fo, "\n");
	if (mld != NULL) {
		fprintf(fo, "1\n");
		mld->write(fo);
	}
	else { fprintf(fo, "0\n"); }
	fprintf(fo, "\n");
	rspd->write(fo); fprintf(fo, "\n");
	qd->write(fo);   fprintf(fo, "\n");
	qpro->write(fo); fprintf(fo, "\n");
	nqpro->write(fo);

	fclose(fo);
}

void SingleQModel::startSimulation(simul* sampler, const std::vector<double>& theta) {
	this->sampler = sampler;

	theta_cdf = new double[M + 1];
	for (int i = 0; i <= M; i++) {
		theta_cdf[i] = theta[i];
		if (i > 0) theta_cdf[i] += theta_cdf[i - 1];
	}

	rspd->startSimulation(M, refs);
	qd->startSimulation();
	qpro->startSimulation();
	nqpro->startSimulation();
}

bool SingleQModel::simulate(READ_INT_TYPE rid, SingleReadQ& read, int& sid) {
	int dir, pos, readLen, fragLen;
	std::string name;
	std::string qual, readseq;
	std::ostringstream strout;

	sid = sampler->sample(theta_cdf, M + 1);

	if (sid == 0) {
		dir = pos = 0;
		readLen = (mld != NULL ? mld->simulate(sampler, -1) : gld->simulate(sampler, -1));
		qual = qd->simulate(sampler, readLen);
		readseq = nqpro->simulate(sampler, readLen, qual);
	}
	else {
		RefSeq &ref = refs->getRef(sid);
		dir = ori->simulate(sampler);
		fragLen = gld->simulate(sampler, ref.getTotLen());
		if (fragLen < 0) return false;

		int effL = std::min(ref.getFullLen(), ref.getTotLen() - fragLen + 1);
		pos = rspd->simulate(sampler, sid, effL);
		if (pos < 0) return false;
		if (dir > 0) pos = ref.getTotLen() - pos - fragLen;

		if (mld != NULL) {
			readLen = mld->simulate(sampler, fragLen);
			if (readLen < 0) return false;
			qual = qd->simulate(sampler, readLen);
			readseq = qpro->simulate(sampler, readLen, pos, dir, qual, ref);
		}
		else {
			qual = qd->simulate(sampler, fragLen);
			readseq = qpro->simulate(sampler, fragLen, pos, dir, qual, ref);
		}
	}

	strout<<rid<<"_"<<dir<<"_"<<sid<<"_"<<pos;
	name = strout.str();

	read = SingleReadQ(name, readseq, qual);

	return true;
}

void SingleQModel::finishSimulation() {
	delete[] theta_cdf;

	rspd->finishSimulation();
	qd->finishSimulation();
	qpro->finishSimulation();
	nqpro->finishSimulation();
}

#endif /* SINGLEQMODEL_H_ */
