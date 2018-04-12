// simPileup.h

#ifndef SIMPILEUP_H_
#define SIMPILEUP_H_

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include "generalUtils.h"
#include "RandomGenerator.h"

// TYPE DEFINITIONS

typedef std::multimap<int, int> rdat;

// CLASS DEFINITIONS

class SiteData
{
public:

	// PUBLIC FUNCTIONS

	SiteData (double depth, double depthsd, int n, double avgq, double maxq, double beta_b, int code, unsigned int* randseed = NULL); // constructor
	void setCoverage (double depth, double depth_sddev); // assigns coverage
	double getInds (int); // assigns number of individuals
	double getMaxQ (double); // assigns maximum quality score
	void getBetaPar (double, double, double []); // assigns beta distribution shape parameters alpha and beta
	int getCode (int); // assigns quality score encoding; offset 33 or 64
	double getDupProb (std::vector<double>* admixpro); // assigns the probability that a read comes from any duplicate copy of a locus
	void assignSeqData (const double altfreq, const double inbreedcoef, const Array<double>* fitness, double cnvf=0); // assigns sequence data to member object seqdat
	std::vector<rdat> genSeqData (const double afreq, const double inbreedcoef, const Array<double>* fitness,
			const unsigned int copyn = 1, unsigned int copy = 0, std::vector<double>* admix = NULL, double cnvf=0); // generates sequencing data for a site
	void getGeno (int genoNum [], double n, const double p, const double F, const Array<double>* w, int ncnv [], double cnvf=0); // generate genotypes based on HWE, inbreeding, and fitness
	void doParalog (std::vector<double>* altf, const double inbreedcoef, std::vector< Array<double> >* fitness, std::vector<double>* admixpro); // generates reads and quality scores for a site comprised of paralogs
	void mergeLoci (std::vector<rdat>* locus2); // combines data from another site with seqdat member
	void printSite (std::fstream& out, int code); // dumps the pileup format line to output
	void printParam (std::fstream& parfile, std::vector<double>* f, std::vector<double>* m, const double F, std::vector< Array<double> >* w); // prints simulations parameters
	int coverage (std::vector<rdat>*); // returns the site coverage
	int getIndN() const; // return number of individuals
//	void alleleSwap (); // simulates effects of multiplexed barcode swapping

	// PUBLIC OBJECTS

	std::string name; // chromosome name
	int pos; // site position
	int code; // quality score encoding
	double swapProb; // probability that allele is derived from another individual
	std::vector<rdat> seqdat; // reads and quality scores for site
	std::mt19937::result_type seed; // random number generator seed
	std::mt19937 randomgen; // random number generator

private:

	// PRIVATE FUNCTIONS

	void doError ( std::vector<rdat>* data ); // creates sequencing errors based on quality scores
	char randAllele (); // randomly picks a reference allele
	char pileRead (int obsread, char ref, char alt); // converts numeric read code to pileup symbol
	char pickStrand (int, char* = NULL); // decides what strand a read comes from

	// PRIVATE OBJECTS

	double cov; // average individual coverage - treated as a normally distributed random variable
	double covsd; // standard deviation for average individual coverage
	std::vector<int> _covvec; // sequencing coverage for each individual
	int nind; // number of individuals
	double qmax; // max quality score
	double betap[2]; // beta shape params; betap[0] = alpha, betap[1] = beta
	Random rand_0_1; // random generator for values in range (0,1)
};

// FUNCTION PROTOTYPES

double unif ();
int poisson (double);
int maxIndex (const double arr [], const int size);

#endif /* SIMPILEUP */
