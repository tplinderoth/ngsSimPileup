// ngsSimPile.h

#ifndef NGSSIMPILE_H_
#define NGSSIMPILE_H_

#include <vector>
#include "generalUtils.h"
#include "SFS.h"

const char* version = "0.2.0"; // version 27 February 2018

// FUNCTION PROTOTYPES

int parseInputs (int argc, char** argv, std::vector<double>* altfreq, std::vector<double>* admix, double& minfreq,
		double& maxfreq,std::vector< Array<double> >* fitness, double* inbreed, double& minF, double& maxF, double& coverage, double& covsd,
		int& nind,double& avgqual, double& maxqual, double& betab, int& qcode, std::string& seqname, double& theta, unsigned int& nsites, std::string& outfile, bool& fold);
int setFitness (std::vector<const char*>* wchar, std::vector< Array<double> >* wvals);
int setInbreeding(const char* Fchar, double* Fval);
int setAltFreq (std::vector<const char*>* fchar, std::vector<double>* fvals);
int setAdmix (std::vector<const char*>* achar, std::vector<double>* avals);
int setDefaultVals (std::vector<double>* admix, std::vector<double>* altfreq, double* inbreed, std::vector< Array<double> >* fitness);
void simNormal (unsigned int nsites, double altfreq, double lbound, double ubound, std::vector<double>* mvec, const double inbreed,
		std::vector< Array<double> >* fitness, SiteData* dat, std::fstream& seqfile, std::fstream& parfile, bool foldsfs);
void simParalog (unsigned int nsites, std::vector<double>* fvec, double lbound, double ubound, std::vector<double>* mvec, const double inbreed,
		std::vector< Array<double> >* fitness, SiteData* dat, std::fstream& seqfile, std::fstream& parfile, bool foldsfs);
double alleleFreq (double freq, double lobound, double upbound, SFS* sfs); // draw allele frequency
void info (const double& covmean, const double& covsd, const double& minfreq, const double& maxfreq, const double& avgqual, const double& maxqual, const double& beta,
		const double& encode, const std::string& name, const double* F, const bool& foldsfs, int help = 0, const char* v = version);

#endif /* NGSSIMPILE_H_ */
