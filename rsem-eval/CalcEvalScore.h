/*
 * CalcEvalScore.h
 */

#ifndef CALCEVALSCORE_H_
#define CALCEVALSCORE_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<cassert>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<algorithm>

#include "utils.h"
#include "my_assert.h"

#include "Refs.h"
#include "GroupInfo.h"
#include "Transcript.h"
#include "Transcripts.h"
#include "ContigLengthDist.h"

class CalcEvalScore {
public:
	CalcEvalScore(Refs&, double, double, int, int, char*, char*, char* = NULL);
	~CalcEvalScore() { assert(cld != NULL); delete cld; }

	void writeScoresTo(char*);
	// must be called after writeScoresTo
	void generateExpressionFiles(GroupInfo&, Transcripts&, char*);

private:
	static const double logp_per_base;

	struct Item {
		int sid;
		double log_conprb;

		Item() : sid(0), log_conprb(0.0) {}
	};

	int M, L;
	READ_INT_TYPE N, N0, N1, N2;
	HIT_INT_TYPE nHits;
	double noiseLogP;

	Refs& refs;
	ContigLengthDist *cld;

	std::vector<int> lens;

	std::vector<HIT_INT_TYPE> s;
	std::vector<Item> hits;

	std::vector<double> theta, counts, logBayesFactors;

	void calcAssemblyPrior(std::vector<double>&, double&, double&, double&);
	void calcMaximumDataLikelihood(double&);
	void calcDataLikelihood(double&, double&, int&, int&);
	void calcCorrectionScore(std::vector<double>&, double&);

	std::string scoreF;
	void loadDataLikelihood(double&, double&, int&, int&, double&);
};

const double CalcEvalScore::logp_per_base = log(0.25);

CalcEvalScore::CalcEvalScore(Refs& refs, double nb_r, double nb_p, int L, int w, char *statName, char *imdName, char *score_file) : L(L), refs(refs)  {
	std::ifstream fin;
	std::string line;
	int tmpVal;

	if (score_file != NULL) scoreF = score_file;

	M = refs.getM();
	int maxl = 0;
	lens.assign(M + 1, 0);
	for (int i = 1; i <= M; i++) { 
		lens[i] = refs.getRef(i).getFullLen();
		maxl = std::max(maxl, lens[i]);
	}

	cld = new ContigLengthDist(nb_r, nb_p, L, w, maxl);

	// loading thetaF
	char thetaF[STRLEN];
	sprintf(thetaF, "%s.theta", statName);
	fin.open(thetaF);
	fin>> tmpVal;
	general_assert(tmpVal == M + 1, "Number of contigs does not match! You might provide a wrong " + cstrtos(thetaF) + "!");
	theta.assign(M + 1, 0.0);
	for (int i = 0; i <= M; i++) {
		fin>>theta[i];
		//if (i > 0 && lens[i] < L) assert(theta[i] == 0);
	}
	fin.close();


	// loading cntF
	char cntF[STRLEN];
	sprintf(cntF, "%s.cnt", statName);
	fin.open(cntF);
	fin>> N0>> N1>> N2>> N;
	READ_INT_TYPE tmp;
	fin>> tmp>> tmp>> tmp;
	fin>> nHits;
	fin.close();
	assert(N2 == 0 && N0 + N1 == N);

	if (scoreF == "") {
	  // loading conprbF
		char conprbF[STRLEN];
		sprintf(conprbF, "%s.conprb", imdName);
		fin.open(conprbF);
		general_assert(fin.is_open(), "Cannot open " + cstrtos(conprbF) + "!");

		READ_INT_TYPE tmpn0, tmpn1;
		HIT_INT_TYPE tmpnhits;
		fin>> tmpVal>> tmpn0>> tmpn1>> tmpnhits>> noiseLogP;
		assert(N0 == tmpn0 && N1 == tmpn1 && nHits == tmpnhits);

		general_assert(tmpVal == M, "M in " + cstrtos(conprbF) + " is not consistent with the reference!");
		getline(fin, line);

		READ_INT_TYPE i = 0; // point to current location of s
		HIT_INT_TYPE j = 0; // point to current location of hits

		s.assign(N1 + 1, 0);
		hits.assign(nHits, Item());

		while (getline(fin, line)) {
			std::istringstream strin(line);
			int sid;
			double log_conprb;

			while (strin>> sid>> log_conprb) {
				hits[j].sid = sid;
				hits[j++].log_conprb = log_conprb;
			}
			s[++i] = j;
		}
		assert(i == N1 && j == nHits);
		fin.close();
	}

	logBayesFactors.assign(M + 1, 0.0);
}

void CalcEvalScore::writeScoresTo(char *outF) {
	std::vector<double> lambda;
	double bic_term, prior_contig_lengths, prior_sequence_bases;
	double data_ll, correction_score;
	double total_score;
	double true_data_ll;

	// Sufficient statistics
	double numAlignedReads;
	int numEmptyContigs, numInvalidContigs;
	
	lambda.assign(M + 1, 0.0);
	lambda[0] = N * theta[0];
	for (int i = 1; i <= M; i++)
	  lambda[i] = std::max(1e-8, N * theta[i] / (lens[i] + L + 1)); // no matter if lens[i] >= L
	calcAssemblyPrior(lambda, bic_term, prior_contig_lengths, prior_sequence_bases);

	// calculate true maximum likelihood score
	if (scoreF == "") calcMaximumDataLikelihood(true_data_ll);
	
	// renomalize lambda to get adjusted theta vector
	double sum = 0.0;
	theta[0] = lambda[0];
	sum += theta[0];
	for (int i = 1; i <= M; i++)
	  if (theta[i] > EPSILON) {
	    theta[i] = lambda[i] * std::max(lens[i] - L + 1, 1);
	    sum += theta[i];
	  }
	  else theta[i] = 0.0;
	for (int i = 0; i <= M; i++) theta[i] /= sum;
	
	if (scoreF == "") calcDataLikelihood(data_ll, numAlignedReads, numEmptyContigs, numInvalidContigs);
	else loadDataLikelihood(data_ll, numAlignedReads, numEmptyContigs, numInvalidContigs, true_data_ll);

	calcCorrectionScore(lambda, correction_score);

	total_score = bic_term + prior_contig_lengths + prior_sequence_bases + data_ll - correction_score;

	FILE *fo = fopen(outF, "w");

	fprintf(fo, "Score\t%.2f\n", total_score);
	fprintf(fo, "BIC_penalty\t%.2f\n", bic_term);
	fprintf(fo, "Prior_score_on_contig_lengths_(f_function_canceled)\t%.2f\n", prior_contig_lengths);
	fprintf(fo, "Prior_score_on_contig_sequences\t%.2f\n", prior_sequence_bases);
	fprintf(fo, "Data_likelihood_in_log_space_without_correction\t%.2f\n", data_ll);
	fprintf(fo, "Correction_term_(f_function_canceled)\t%.2f\n", correction_score);
	fprintf(fo, "Number_of_contigs\t%d\n", M);
	fprintf(fo, "Expected_number_of_aligned_reads_given_the_data\t%.2f\n", numAlignedReads);
	fprintf(fo, "Number_of_contigs_smaller_than_expected_read/fragment_length\t%d\n", numInvalidContigs);
	fprintf(fo, "Number_of_contigs_with_no_read_aligned_to\t%d\n", numEmptyContigs);
	fprintf(fo, "Maximum_data_likelihood_in_log_space\t%.2f\n", true_data_ll);
	fprintf(fo, "Number_of_alignable_reads\t%llu\n", (unsigned long long)N1);
	fprintf(fo, "Number_of_alignments_in_total\t%llu\n", (unsigned long long)nHits); 
	fclose(fo);
}

void CalcEvalScore::calcAssemblyPrior(std::vector<double>& lambda, double& bic_term, double& prior_contig_lengths, double& prior_sequence_bases) {
	if (verbose) { printf("Calculating assembly priors.\n"); }

	double log_prior;

	bic_term = - 0.5 * (M + 1) * log(1.0 * N); // even if some contigs are less than L, still penalize them

	for (int i = 1; i <= M; i++) logBayesFactors[i] -= 0.5 * log(1.0 * N);

	prior_contig_lengths = 0.0;
	for (int i = 1; i <= M; i++) {
	  log_prior = cld->calcLogPrior(lambda[i], lens[i]);
	  prior_contig_lengths += log_prior;
	  logBayesFactors[i] += log_prior;
	}

	prior_sequence_bases = 0.0; // even if some contigs are less than L, still penalize each base of the contig
	for (int i = 1; i <= M; i++) {
		std::string seq = refs.getRef(i).getSeq();
		double psb = 0.0;
		for (int j = 0; j < lens[i]; j++) psb += (seq[j] != 'N');
		psb *= logp_per_base;
		
		prior_sequence_bases += psb;
		logBayesFactors[i] += psb;
	}

	if (verbose) { printf("Assembly priors are calculated!\n"); }
}

void CalcEvalScore::calcMaximumDataLikelihood(double& true_data_ll) {
	HIT_INT_TYPE fr, to;
	double sum, max_value;

	if (verbose) { printf("Calculating maximum data likelihood in log space.\n"); }

	true_data_ll = noiseLogP + (N0 > 0 ? N0 * log(theta[0]) : 0.0);
	for (READ_INT_TYPE i = 0; i < N1; i++) {
		fr = s[i]; to = s[i + 1];
		max_value = hits[fr].log_conprb;
		for (HIT_INT_TYPE j = fr + 1; j < to; j++)
		  if (max_value < hits[j].log_conprb) max_value = hits[j].log_conprb;
		sum = 0.0;
		for (HIT_INT_TYPE j = fr; j < to; j++) 
		  sum += theta[hits[j].sid] * exp(hits[j].log_conprb - max_value);
		assert(sum > EPSILON);
		true_data_ll += log(sum) + max_value;
	}

	if (verbose) { printf("Maximum data likelihood is calculated!\n"); }
}

void CalcEvalScore::calcDataLikelihood(double& data_ll, double& numAlignedReads, int& numEmptyContigs, int& numInvalidContigs) {
	HIT_INT_TYPE fr, to;
	double sum, max_value, frac;
	std::vector<double> fracs;

	if (verbose) { printf("Calculating data likelihood in log space.\n"); }

	counts.assign(M + 1, 0.0);
	counts[0] = N0;

	data_ll = noiseLogP + (N0 > 0 ? N0 * log(theta[0]) : 0.0);
	for (READ_INT_TYPE i = 0; i < N1; i++) {
		fr = s[i]; to = s[i + 1];
		max_value = hits[fr].log_conprb;
		for (HIT_INT_TYPE j = fr + 1; j < to; j++)
		  if (max_value < hits[j].log_conprb) max_value = hits[j].log_conprb;
		sum = 0.0;
		fracs.resize(to - fr);
		for (HIT_INT_TYPE j = fr; j < to; j++) {
		  fracs[j - fr] = theta[hits[j].sid] * exp(hits[j].log_conprb - max_value);
		  sum += fracs[j - fr];
		}
		assert(sum > EPSILON);
		data_ll += log(sum) + max_value;

		assert(hits[fr].sid == 0);
		frac = fracs[0] / sum;
		if (frac > EPSILON) counts[0] += frac;
		for (HIT_INT_TYPE j = fr + 1; j < to; j++) {
			frac = fracs[j - fr] / sum;
			if (frac > EPSILON) {
			  counts[hits[j].sid] += frac;
			  assert(theta[hits[j].sid] > EPSILON);
			  //slightly different from the supplementary text
			  logBayesFactors[hits[j].sid] += frac * (log(theta[hits[j].sid]) + hits[j].log_conprb - log(theta[hits[j].sid] + theta[0]) - hits[fr].log_conprb);
			}
		}
	}

	numAlignedReads = 0.0;
	numEmptyContigs = numInvalidContigs = 0;
	for (int i = 1; i <= M; i++) {
		numAlignedReads += counts[i];
		if (lens[i] < L) ++numInvalidContigs;
		if (counts[i] < 0.005) ++numEmptyContigs; // since expected counts are rouned to 2 digits after the decimal point
	}

	if (verbose) { printf("Data likelihood is calculated!\n"); }
}

void CalcEvalScore::calcCorrectionScore(std::vector<double>& lambda, double& correction_score) {
	double log_norm;

	if (verbose) { printf("Calculating correction score.\n"); }
	
	correction_score = 0.0;
	for (int i = 1; i <= M; i++) if (theta[i] > EPSILON) {
	    log_norm = cld->calcLogNorm(lambda[i], lens[i]);
	    correction_score += log_norm;
	    logBayesFactors[i] -= log_norm;
	}

	if (verbose) { printf("Correction score is calculated!\n"); }
}

void CalcEvalScore::generateExpressionFiles(GroupInfo& gi, Transcripts& transcripts, char* outName) {
	int m = gi.getm();
	double denom;
	std::vector<double> frac, fpkm, tpm, isopct;
	std::vector<double> glens, gene_eels, gene_counts, gene_tpm, gene_fpkm;
	char outF[STRLEN];
	FILE *fo;

	denom = 0.0;
	frac.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++) {
		frac[i] = theta[i];
		denom += frac[i];
	}
	general_assert(denom > EPSILON, "No alignable reads?!");
	for (int i = 1; i <= M; i++) frac[i] /= denom;

	// calculate FPKM
	fpkm.assign(M + 1, 0.0);
	for (int i = 1; i <= M; i++)
	  fpkm[i] = frac[i] * 1e9 / std::max(lens[i] - L + 1, 1);

	// calculate TPM
	tpm.assign(M + 1, 0.0);
	denom = 0.0;
	for (int i = 1; i <= M; i++) denom += fpkm[i];
	for (int i = 1; i <= M; i++) tpm[i] = fpkm[i] / denom * 1e6;

	//calculate IsoPct, etc.
	isopct.assign(M + 1, 0.0);

	glens.assign(m, 0.0); gene_eels.assign(m, 0.0);
	gene_counts.assign(m, 0.0); gene_tpm.assign(m, 0.0); gene_fpkm.assign(m, 0.0);

	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		for (int j = b; j < e; j++) {
			gene_counts[i] += counts[j];
			gene_tpm[i] += tpm[j];
			gene_fpkm[i] += fpkm[j];
		}

		if (gene_tpm[i] <= EPSILON) {
			double frac = 1.0 / (e - b);
			for (int j = b; j < e; j++) {
				glens[i] += lens[j] * frac;
				gene_eels[i] += std::max(0, (lens[j] - L + 1)) * frac;
			}
		}
		else {
			for (int j = b; j < e; j++) {
				isopct[j] = tpm[j] / gene_tpm[i];
				glens[i] += lens[j] * isopct[j];
				gene_eels[i] += std::max(0, (lens[j] - L + 1)) * isopct[j];
			}
		}
	}

	//isoform level results
	sprintf(outF, "%s.isoforms.results", outName);
	fo = fopen(outF, "w");
	fprintf(fo, "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tCPM\tFPKM\tIsoPct\tcontig_impact_score\n");
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		fprintf(fo, "%s\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str(),
			lens[i], std::max(0, lens[i] - L + 1), counts[i], tpm[i], fpkm[i], isopct[i] * 1e2, logBayesFactors[i]);
	}
	fclose(fo);

	//gene level results
	sprintf(outF, "%s.genes.results", outName);
	fo = fopen(outF, "w");
	fprintf(fo, "gene_id\ttranscript_id(s)\tlength\teffective_length\texpected_count\tCPM\tFPKM\n");
	for (int i = 0; i < m; i++) {
		int b = gi.spAt(i), e = gi.spAt(i + 1);
		fprintf(fo, "%s\t", transcripts.getTranscriptAt(b).getGeneID().c_str());
		for (int j = b; j < e; j++) {
			fprintf(fo, "%s%c", transcripts.getTranscriptAt(j).getTranscriptID().c_str(), (j < e - 1 ? ',' : '\t'));
		}
		fprintf(fo, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", glens[i], gene_eels[i], gene_counts[i], gene_tpm[i], gene_fpkm[i]);
	}
	fclose(fo);
}

void CalcEvalScore::loadDataLikelihood(double& data_ll, double& numAlignedReads, int& numEmptyContigs, int& numInvalidContigs, double& true_data_ll) {
  std::ifstream fin(scoreF.c_str());
  std::string line, name;

  for (int i = 0; i < 4; i++) getline(fin, line);
  fin>> name>> data_ll;
  for (int i = 0; i < 3; i++) getline(fin, line);
  fin>> name>> numAlignedReads;
  fin>> name>> numInvalidContigs;
  fin>> name>> numEmptyContigs;
  fin>> name>> true_data_ll;
  fin.close();
}

#endif /* CALCEVALSCORE_H_ */
