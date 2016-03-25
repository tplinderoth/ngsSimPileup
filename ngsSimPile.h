// ngsSimPile.h

#ifndef NGSSIMPILE_H_
#define NGSSIMPILE_H_

#include <vector>

const char* version = "0.0.6"; // version 9 Feb 2015

// FUNCTION PROTOTYPES

int parseInputs (int argc, char** argv, std::vector<double>& altfreq, std::vector<double>& admix, double& minfreq, double& maxfreq,
double& coverage, int& nind, double& avgqual, double& maxqual, double& betab, int& qcode,
std::string& seqname, int& paralog_sites, int& normal_sites, std::string& outfile);
void simNormal (int nsites, double lbound, double ubound, std::vector<double>* mvec, SiteData* dat, std::fstream& seqfile, std::fstream& parfile);
void simParalog (std::vector<double>* fvec, std::vector<double>* mvec, int nsites, double lbound, double ubound, SiteData* dat,
std::fstream& seqfile, std::fstream& parfile);
void info (const double& cov, const double& minfreq, const double&maxfreq, const double& avgqual, const double& maxqual,
		const double& beta, const double& encode, const std::string& name, int help = 0, const char* v = version);

#endif /* NGSSIMPILE_H_ */
