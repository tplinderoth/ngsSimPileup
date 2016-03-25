// simPileup.h

#ifndef SIMPILEUP_H_
#define SIMPILEUP_H_

#include <vector>
#include <map>
#include <string>
#include <fstream>

// TYPE DEFINITIONS

typedef std::multimap<int, int> rdat;

// CLASS DEFINITIONS

class SiteData
{
public:

	// PUBLIC FUNCTIONS

	SiteData (double depth, int n, double avgq, double maxq, double beta_b, int code); // constructor
	double getCoverage (double); // assigns coverage
	double getInds (int); // assigns number of individuals
	double getMaxQ (double); // assigns maximum quality score
	void getBetaPar (double, double, double []); // assigns beta distribution shape parameters alpha and beta
	int getCode (int); // assigns quality score encoding; offset 33 or 64
	double getDupProb (std::vector<double>* admixpro); // assigns the probability that a read comes from any duplicate copy of a locus
	void assignSeqData (const double); // assigns sequence data to member object seqdat
	std::vector<rdat> genSeqData(const double afreq, const double pdup = 1.0); // generates sequencing data for a site
	void getGeno (int [], const int, const double); // calculates genotype numbers based on HWE
	void doParalog (std::vector<double>* altf, std::vector<double>* admixpro); // generates reads and quality scores for a site comprised of paralogs
	void mergeLoci (std::vector<rdat>* locus2); // combines data from another site with seqdat member
	void printSite (std::fstream& out, int code); // dumps the pileup format line to output
	void printParam (std::fstream& parfile, std::vector<double>* f, std::vector<double>* m); // prints simulation parameters
	int coverage (std::vector<rdat>*); // returns the site coverage
//	void alleleSwap (); // simulates effects of multiplexed barcode swapping

	// PUBLIC OBJECTS

	std::string name; // chromosome name
	int pos; // site position
	int code; // quality score encoding
	double dupReadProb; // probability that a read comes from any duplicate copy of a locus
	double swapProb; // probability that allele is derived from another individual
	std::vector<rdat> seqdat; // reads and quality scores for site
	std::vector<int> realcov; // actual sequencing coverage for each individual

private:

	// PRIVATE FUNCTIONS

	void doError ( std::vector<rdat>* data ); // creates sequencing errors based on quality scores
	char randAllele (); // randomly picks a reference allele
	char pileRead (int obsread, char ref, char alt); // converts numeric read code to pileup symbol
	char pickStrand (int, char* = NULL); // decides what strand a read comes from

	// PRIVATE OBJECTS

	double cov; // average individual coverage
	int nind; // number of individuals
	double qmax; // max quality score
	double betap[2]; // beta shape params; betap[0] = alpha, betap[1] = beta
};

// FUNCTION PROTOTYPES

double unif ();
int poisson (double);
int maxIndex (const double arr [], const int size);

#endif /* SIMPILEUP */
