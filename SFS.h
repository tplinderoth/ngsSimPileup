// SFS.h

#include <vector>

class SFS {
public:
	SFS(unsigned int nind, bool fold = false);
	double rfreq(); // draw allele frequency from SFS distribution
	double getProb(unsigned int category); // get probability of specific SFS category
private:
	void setProbs (); // set probability of each SFS category
	unsigned int _nind; // number of individuals
	unsigned int _ncat; // number of SFS categories
	std::vector<double> _probs; // probabilities of each SFS category
	std::vector<double> _scaledprobs; // probabilities scaled for easy sampling using uniform random variates
};
