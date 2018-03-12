// simPiluep.cpp

#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include "simPileup.h"
#include "rbeta.h"

// FUNCTION DEFINITIONS

// SiteData constructor
SiteData::SiteData (double depth, double depthsd, int n, double avgq, double maxq, double beta_b, int qoffset, unsigned int* randseed)
	: name ("chr"),
	  pos (0),
	  swapProb(0.0),
	  seed(std::chrono::system_clock::now().time_since_epoch().count())

{
	if (randseed) seed = *randseed;
	randomgen.seed(seed);

	setCoverage(depth, depthsd);
	nind = getInds(n);
	qmax = getMaxQ(maxq);
	getBetaPar(avgq, beta_b, betap);
	code = getCode(qoffset);
	_covvec.resize(nind, 0);
}

// SiteData::setCoverage assigns coverage
void SiteData::setCoverage (double depth, double depth_sddev)
{
	if (depth >= 0)
		cov = depth;
	else {
		fprintf(stderr, "Attempt to set negative coverage in call to SiteData::setCoverage\n");
		exit(EXIT_FAILURE);
	}

	if (depth_sddev >= 0)
		covsd = depth_sddev;
	else {
		fprintf(stderr, "Attempt to set negative coverage standard deviation in call to SiteData::setCoverage\n");
		exit(EXIT_FAILURE);
	}
}

// SiteData::getInds assigns number of individuals
double SiteData::getInds (int n)
{
	if (n > 0)
		return n;
	else
	{
		fprintf(stderr, "Number of Individuals for SiteData object < 1\n");
		exit(EXIT_FAILURE);
	}
}

// SiteData::getMaxQ assigns maximum quality score
double SiteData::getMaxQ (double q)
{
	if (q > 0.0 && q <= 93.0) //can't be zero because of beta
		return q;
	else
	{
		fprintf(stderr, "Maximum quality for SiteData object is outside (0.0, 93.0]\n");
		exit(EXIT_FAILURE);
	}
}

// SiteData::GetBetaPar assigns beta distribution shape parameters alpha and beta
void SiteData::getBetaPar (double qbar, double b, double para [])
{
	if (qbar <= 0.0 || qbar > 93.0)
	{
		fprintf(stderr, "Average quality for SiteData object is outside (0.0, 93.0]\n");
		exit(EXIT_FAILURE);
	}

	if (b <= 0 )
	{
		fprintf(stderr, "Beta shape parameter for SiteData object <= 0\n");
		exit(EXIT_FAILURE);
	}

	double scaled_mu = qbar / qmax; // scale beta mean on [0,1]
	para[0] = (scaled_mu * b) / (1 - scaled_mu); // alpha shape parameter
	para[1] = b; // beta shape parameter
}

// SiteData::getCode assigns quality score encoding
int SiteData::getCode (int qualcode)
{
	if (qualcode == 33 || qualcode == 64)
		return qualcode;
	else
	{
		fprintf(stderr, "Quality score encoding can only have offset 33 or 64 -> exiting\n");
		exit(EXIT_FAILURE);
	}
}

// SiteData::getDupProb assigns the probability that a read comes from any duplicate copy of a locus
double SiteData::getDupProb (std::vector<double>* admixpro)
{
	double prob = 0;
	if (!admixpro->empty())
	{
		for (std::vector<double>::iterator i = admixpro->begin() + 1; i != admixpro->end(); ++i)
			prob += *i;
	}
	else
	{
		fprintf(stderr, "No read admixture proportions supplied to SiteData::getDupProb -> exiting\n");
		exit(EXIT_FAILURE);
	}
	return prob;
}

/*
// SiteData::genSeqData generates reads and quality scores for a site
std::vector<rdat> SiteData::genSeqData (const double afreq, const double inbreedcoef, const Array<double>* fitness, const unsigned int copyn, unsigned int copy, std::vector<double>* admix)
{
	// this function version draws a new poisson parameter value for each individual at a site

	int ind_geno [3] = {}; // # ref homo, # het, # alt homo
	getGeno(ind_geno, nind, afreq, inbreedcoef, fitness);
	int readn = 0; // number of reads for an individual
	int ind = 0;
	int q; // quality score
	double allcov = 0; // number of all reads from duplicates before mapping bias
	static double mappedcov = 0; // number of reads from duplicates that actually map
	double maxp = 1.0/copyn; // proportion of total reads coming from a single locus if contribution is equal
	double pdup = admix ? (*admix)[copy] : 1.0; // admix proportion for current locus
	std::vector<rdat> data (nind); // read data
	static std::vector<double>::iterator iter;
	static int prevpos = -1;

	double covmin = 0.0;
	double covmax = 1.0/0.0;
	static NormalGenerator normgen(cov, covsd, covmin, covmax);
	double poismean = 0.0;

	if (pdup > 1.0 || pdup < 0.0)
	{
		fprintf(stderr, "Invalid admixture proportion in SiteData::genSeqData -> exiting\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 3; i++) // go over all 3 genotype configurations
	{
		for (int j = 0; j < ind_geno[i]; j++) // go over number of individuals in current genotype
		{
			if (prevpos != pos) {
				poismean = 0.0;

				for (unsigned int k=0; k < copyn; ++k)
					poismean += normgen();

				std::poisson_distribution<int> poissdist(poismean);
				mappedcov = 0.0;
				allcov = poissdist(randomgen);
				if (copyn > 1) {
						for(iter = admix->begin(); iter != admix->end(); ++iter)
							mappedcov += (*iter < maxp) ? *iter * allcov : maxp * allcov;
					mappedcov = ceil(mappedcov);
				} else
					mappedcov = allcov;
			}
			if (pdup < 1.0)
			{
				readn = 0;
				for (int k = 0; k < mappedcov; ++k)
				{
					if (rand_0_1() < pdup) ++readn;
				}
			}
			else
				readn = mappedcov;
			if (readn == 0) // missing data for individual
				data[ind].insert( rdat::value_type ( '*', '*') );
			else
			{
				for (int k = 0; k < readn; k++) // assign reads and quality scores
				{
					q = floor(rbeta(betap[0], betap[1]) * qmax + 0.5);

					if (i == 0) // ref homo
						data[ind].insert( rdat::value_type ( q, 0 ) );
					else if (i == 2) // ref alt
						data[ind].insert( rdat::value_type ( q, 1 ) );
					else // het
					{
						if (rand_0_1() < 0.5)
							data[ind].insert( rdat::value_type ( q, 0 ) );
						else
							data[ind].insert( rdat::value_type ( q, 1 ) );
					}
				}
			}

			ind++;
		}
	}

	doError(&data); // introduce sequencing errors
	prevpos = pos;

	return data;
}
*/


// SiteData::genSeqData generates reads and quality scores for a site
std::vector<rdat> SiteData::genSeqData (const double afreq, const double inbreedcoef, const Array<double>* fitness, const unsigned int copyn, unsigned int copy, std::vector<double>* admix)
{
	//this function version use the same poisson parameter value to simulate the coverage for each individual at a site

	int ind_geno [3] = {}; // # ref homo, # het, # alt homo
	getGeno(ind_geno, nind, afreq, inbreedcoef, fitness);
	int readn = 0; // number of reads for an individual
	int ind = 0;
	int q; // quality score
	double allcov = 0; // number of all reads from duplicates before mapping bias
	static double mappedcov = 0; // number of reads from duplicates that actually map
	double maxp = 1.0/copyn; // proportion of total reads coming from a single locus if contribution is equal
	double pdup = admix ? (*admix)[copy] : 1.0; // admix proportion for current locus
	std::vector<rdat> data (nind); // read data
	static std::vector<double>::iterator iter;
	static int prevpos = -1;

	double covmin = 0.0;
	double covmax = 1.0/0.0;
	static NormalGenerator normgen(cov, covsd, covmin, covmax);
	double poismean = copyn*normgen();
	std::poisson_distribution<int> poissdist(poismean);

	if (pdup > 1.0 || pdup < 0.0)
	{
		fprintf(stderr, "Invalid admixture proportion in SiteData::genSeqData -> exiting\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 3; i++) // go over all 3 genotype configurations
	{
		for (int j = 0; j < ind_geno[i]; j++) // go over number of individuals in current genotype
		{
			if (prevpos != pos) {
				mappedcov = 0.0;
				allcov = poissdist(randomgen);
				if (copyn > 1) {
						for(iter = admix->begin(); iter != admix->end(); ++iter)
							mappedcov += (*iter < maxp) ? *iter * allcov : maxp * allcov;
					mappedcov = ceil(mappedcov);
				} else
					mappedcov = allcov;
			}
			if (pdup < 1.0)
			{
				readn = 0;
				for (int k = 0; k < mappedcov; ++k)
				{
					if (rand_0_1() < pdup) ++readn;
				}
			}
			else
				readn = mappedcov;
			if (readn == 0) // missing data for individual
				data[ind].insert( rdat::value_type ( '*', '*') );
			else
			{
				for (int k = 0; k < readn; k++) // assign reads and quality scores
				{
					q = floor(rbeta(betap[0], betap[1]) * qmax + 0.5);

					if (i == 0) // ref homo
						data[ind].insert( rdat::value_type ( q, 0 ) );
					else if (i == 2) // ref alt
						data[ind].insert( rdat::value_type ( q, 1 ) );
					else // het
					{
						if (rand_0_1() < 0.5)
							data[ind].insert( rdat::value_type ( q, 0 ) );
						else
							data[ind].insert( rdat::value_type ( q, 1 ) );
					}
				}
			}

			ind++;
		}
	}

	doError(&data); // introduce sequencing errors
	prevpos = pos;

	return data;
}


// SiteData::getGeno calculates genotype numbers based on HWE
void SiteData::getGeno (int genoNum [], double n, const double p, const double F, const Array<double>* w)
{
	// ADD EXCEPTION HANDLING FOR CHECKING ALLELE FREQ, INBREEDING COEF, FITNESS

	if (p < 0.0 || p > 1.0)
	{
		fprintf(stderr, "Alternate allele frequency out of [0,1] range in call to SiteData::genGeno -> exiting\n");
		exit(EXIT_FAILURE);
	}

	if (F < 0.0 || F > 1.0)
	{
		fprintf(stderr, "Inbreeding coefficient out of [0,1] range in call to SiteData::getGeno -> exiting\n");
		exit(EXIT_FAILURE);
	}

	if (w->size() == 3)
	{
		for(unsigned int i=0; i<3; ++i)
		{
			if ((*w)[i] < 0.0 || (*w)[i] > 1.0)
			{
				fprintf(stderr, "Fitness values out of [0,1] range in call to SiteData::getGeno -> exiting\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	else
	{
		fprintf(stderr, "Missing fitness values in call to SiteData::getGeno -> exiting\n");
		exit(EXIT_FAILURE);
	}


	double q = 1.0 - p;

	// genotype frequencies accounting for inbreeding and fitness
	double paa = (*w)[2] * ( q * q * (1.0 - F) + q * F );
	double pAa = (*w)[1] * ( 2.0 * p * q * (1.0 - F) );
	double pAA = (*w)[0] * ( p * p * (1.0 - F) + p * F );

	// calculate average fitness and scale the genotype frequencies
	double wbar = paa + pAa + pAA;
	paa /= wbar;
	pAa /= wbar;
	pAA /= wbar;

	// expected number of each genotype
	double rhomo = n * paa;
	double het = n * pAa;
	double ahomo = n * pAA;

	genoNum[0] = static_cast<int>( floor(rhomo + 0.5) );
	genoNum[1] = static_cast<int>( floor(het + 0.5) );
	genoNum[2] = static_cast<int>( floor(ahomo + 0.5) );

	// ensure that the sum of genotype numbers equals the sample size
	int genoSum = genoNum[0] + genoNum[1] + genoNum[2];

	if (genoSum == n)
		return;
	else
	{
		double diff [3] = {};
		diff[0] = abs(rhomo - (double)genoNum[0]);
		diff[1] = abs(het - (double)genoNum[1]);
		diff[2] = abs(ahomo - (double)genoNum[2]);
		int index = 0;
		int nchanges = 0;

		while (genoSum != n)
		{
			index = maxIndex(diff, 3);

			if (genoSum < n)
			{
				if (index == 0)
				{
					genoNum[0]++;
					diff[0] = 0;
				}
				else if (index == 1)
				{
					genoNum[1]++;
					diff[1] = 0;
				}
				else
				{
					genoNum[2]++;
					diff[2] = 0;
				}
			}
			else // genoSum > n
			{
				if (index == 0)
				{
					genoNum[0]--;
					diff[0] = 0;
				}
				else if (index == 1)
				{
					genoNum[1]--;
					diff[1] = 0;
				}
				else
				{
					genoNum[2]--;
					diff[2] = 0;
				}
			}

			nchanges++;
			if (nchanges == 3) // reset the real differences
			{
				diff[0] = abs(rhomo - (double)genoNum[0]);
				diff[1] = abs(het - (double)genoNum[1]);
				diff[2] = abs(ahomo - (double)genoNum[2]);
				nchanges = 0;
			}
			genoSum = genoNum[0] + genoNum[1] + genoNum[2];
		}

	}
}

// SiteData::assignSeqData assigns sequence data to member object seqdat
void SiteData::assignSeqData (const double altfreq, const double inbreedcoef, const Array<double>* fitness)
{
	seqdat.resize(nind);
	seqdat = genSeqData(altfreq, inbreedcoef, fitness);
}

// SiteData::doParalog generates reads and quality scores for a site comprised of paralogs - see simPileup.h
void SiteData::doParalog (std::vector<double>* altf, const double inbreedcoef, std::vector< Array<double> >* fitness, std::vector<double>* admixpro)
{
	if (altf->size() < 2)
	{
		fprintf(stderr, "At least 2 alternate allele frequencies required for SiteData::doParalog -> exiting\n");
		exit(EXIT_FAILURE);
	}
	if (admixpro->size() < 2)
	{
		fprintf(stderr, "Too few admixture proportions supplied to SiteData::doParalog -> exiting\n");
		exit(EXIT_FAILURE);
	}

	std::vector<rdat> dup;

	//generate data for duplicate loci
	for (unsigned int i = 0; i < altf->size(); i++)
	{
		if (i == 0)
			seqdat = genSeqData((*altf)[i], inbreedcoef, &(*fitness)[i], altf->size(), i, admixpro); // get seq data for locus 1
		else
		{
			dup = genSeqData((*altf)[i], inbreedcoef, &(*fitness)[i], altf->size(), i, admixpro); // get seq data for duplicates
			mergeLoci(&dup); //combine data with seqdat member
		}
	}
}

// SiteData::mergeLoci combines data from another site with seqdat member
void SiteData::mergeLoci (std::vector<rdat>* locus2)
{
	std::vector<rdat>::iterator ind2;
	rdat::const_iterator dat1;
	rdat::const_iterator dat2;
	int i = 0;

	if (!seqdat.empty())
	{
		for (ind2 = locus2->begin(); ind2 != locus2->end(); ++ind2)
		{
			dat1 = seqdat[i].begin();
			dat2 = ind2->begin();

			if (dat1->second == '*' || dat2->second == '*') // missing data
			{
				if (dat2->second != '*') // locus 1 only missing data
					seqdat[i].swap(*ind2);
				else // locus1 and locus 2 missing data OR locus 2 only missing data
				{
					i++;
					continue;
				}
			}
			else
				seqdat[i].insert(ind2->begin(), ind2->end());

			i++;
		}
	}
	else
	{
		fprintf(stderr, "Pileup member seqdat is empty in call to SiteData::mergeLoci -> exiting\n");
		exit(EXIT_FAILURE);
	}
}

/*
// SiteData::alleleSwap simulates the effects of multiplexed barcode swapping
void SiteData::alleleSwap ( std::vector<rdat>* data, double pswap, double ubswap, double lbswap)
{

	std::vector<rad>::iterator;
	rdat::const_iterator dat;

	if (static_cast<int> (pswap) == 82 || static_cast<int> (pswap) == 114)
	{
		swapProb = unif();
	}
	else if (pswap >= 0.0 && pswap <= 1.0)
	{
		swapProb = pswap;
	}
	else
	{
		fprintf(stderr, "Invalid allele swap probability in SiteData::alleleSwap -> exiting\n");
		exit(EXIT_FAILURE);
	}

	if ( !data.empty() )
	{

	}


}
*/

// SiteData::doError creates sequencing errors based on quality scores
void SiteData::doError ( std::vector<rdat>* data )
{
    if (data->empty())
    {
            fprintf(stderr, "warning: attempted SiteData::doError on empty container\n");
            return;
    }

	double epsilon;
	double x, err_readp;

	for ( std::vector<rdat>::iterator ind = data->begin(); ind != data->end(); ++ind)
	{
		for (rdat::iterator dat = ind->begin(); dat != ind->end(); ++dat)
		{
			if (dat->second == '*') // missing data for individual
				continue;

			x = rand_0_1(); // draw random number
			epsilon = pow(10.0, -static_cast<double>(dat->first)/10.0);

			if (x <= epsilon)
			{
				err_readp = rand_0_1(); // 1/3 probability of picking either of the wrong reads

				if (dat->second == 0) // ref change
				{
					if (err_readp <= static_cast<double>(1)/3)
						dat->second = 2; // ref error read1
					else if (err_readp > static_cast<double>(2)/3)
						dat->second = 3; // ref error read2
					else
						dat->second = 4; // ref error read3
				}
				else // alt change
				{
					if (err_readp <= static_cast<double>(1)/3)
						dat->second = 5; // alt error read1
					else if (err_readp > static_cast<double>(2)/3)
						dat->second = 6; // alt error read2
					else
						dat->second = 7; // alt error read3
				}

			}

		}
	}
}

// SiteData::coverage returns the site coverage
int SiteData::coverage (std::vector<rdat>* site)
{
	int depth = 0;

	for (std::vector<rdat>::const_iterator ind = site->begin(); ind != site->end(); ++ind)
		depth += ind->size();

	return depth;
}

// SiteData::printSite dumps the pileup format line to output
void SiteData::printSite (std::fstream& out, int code)
{

	std::string reads; // reads for an individual
	std::string scores; // quality scores for an individual
	int i = 0;
	char ref = randAllele(); // get reference allele
	char alt; // alt allele
	unsigned int depth; // coverage for site

	do
		alt = randAllele();
	while (alt == ref);

	if (!seqdat.empty())
	{
		out << name << "\t" << pos << "\t" << ref;

		for (std::vector<rdat>::const_iterator ind = seqdat.begin(); ind != seqdat.end(); ++ind)
		{
			reads.resize(ind->size());
			scores.resize(ind->size());
			depth = ind->size();

			for (rdat::const_iterator readdat = ind->begin(); readdat != ind->end(); ++ readdat)
			{

				if (readdat->second == 0) //ref
				{
					reads[i] = pickStrand(0);
					scores[i] = static_cast<char> (readdat->first + code);
				}
				else if (readdat->second == 1) //alt
				{
					reads[i] = pickStrand(1, &alt); // reverse strand
					scores[i] = static_cast<char> (readdat->first + code);
				}
				else if (readdat->second > 1 && readdat->second < 8) //error
				{
					reads[i] = pileRead(readdat->second, ref, alt);
					scores[i] = static_cast<char> (readdat->first + code);
				}
				else // missing data for individual
				{
					reads[i] = '*';
					scores[i] = '*';
					depth = 0;
				}

				i++;
			}

			out << "\t" << depth << "\t" << reads << "\t" << scores;

			reads.clear();
			scores.clear();
			i = 0;
		}
	}

	out << "\n";
}

// SiteData::printParam prints simulation parameters for site
void SiteData::printParam (std::fstream& parfile, std::vector<double>* f, std::vector<double>* m, const double F, std::vector< Array<double> >* w)
{
	// prints (1) seqid, (2) position, (3) admixture proportions, (4) allele frequency, (5) fitness values, (6) inbreeding coefficient

	parfile.width(12);
	parfile << std::left << name << "\t";
	parfile.width(12);
	parfile << std::left << pos << "\t";
	unsigned int i = 0;

	// print admixture proportion
	if (f->size() == 1)
	{
		parfile.width(12);
		parfile << std::left << 1.0 << "\t";
	}
	else if (f->size() > 1)
	{
		for (i = 0; i < m->size(); i++)
		{
			if (i == m->size() - 1)
			{
				parfile.width(12);
				parfile << std::left << (*m)[i] << "\t";
			}
			else
				parfile << (*m)[i] << ";";
		}
	}
	else
	{
		fprintf(stderr, "No allele frequency provided to SiteData::printParam -> exiting\n");
		exit(EXIT_FAILURE);
	}

	// print allele frequency
	for (i = 0; i < f->size(); i++)
	{
		if (i == f->size() - 1)
		{
			parfile.width(12);
			parfile << (*f)[i] << "\t";
		}
		else
			parfile << (*f)[i] << ";";
	}

	// print fitness
	for (i=0; i < w->size(); ++i)
	{
		for (int j=0; j<3; ++j)
		{
			if (i == w->size()-1 && j == 2)
			{
				parfile.width(12);
				parfile << (*w)[i][j] << "\t";
			}
			else
				parfile << (*w)[i][j] << ";";
		}
	}

	// print inbreeding coefficient
	parfile << F << "\n";

}

// SiteData::randAllele randomly picks a reference allele
char SiteData::randAllele ()
{
	double rval = rand_0_1();
	if (rval <= 0.25)
		return('A');
	else if (rval <= 0.5)
		return('C');
	else if (rval <= 0.75)
		return('G');
	else
		return('T');
}

// SiteData::pickStrand decides what strand a read comes from
char SiteData::pickStrand (int readcode, char* allele)
{
	char read;
	double strand = rand_0_1();
	if (readcode == 0)
		if (strand <= 0.5)
			read = '.';
		else
			read = ',';
	else
		if (strand <= 0.5)
			read = *allele;
		else
			read = tolower(*allele);

	return read;
}


// SiteData::pileRead converts numeric read code to pileup symbol
char SiteData::pileRead (int obs, char ref, char alt)
{
	static const char a [] = "CGT";
	static const char c [] = "AGT";
	static const char g [] = "ACT";
	static const char t [] = "ACG";
	char err;

	if (obs > 1 && obs < 5) // ref change
	{
		if (ref == 'A')
			err = a[obs - 2];
		else if (ref == 'C')
			err = c[obs - 2];
		else if (ref == 'G')
			err = g[obs - 2];
		else
			err = t[obs - 2];

		if (rand_0_1() > 0.5) // decide strand
			err = tolower(err);

	}
	else if (obs > 4 && obs < 8) // alt change
	{
		if (alt == 'A')
			err = a[obs - 5];
		else if (alt == 'C')
			err = c[obs - 5];
		else if (alt == 'G')
			err = g[obs - 5];
		else
			err = t[obs - 5];

		if (err == ref)
			err = pickStrand(0);
		else
		{
			if (rand_0_1() > 0.5)
				err = tolower(err);
		}
	}
	else
	{
		fprintf(stderr, "Invalid error type in SiteData::pileRead -> exiting\n");
		exit(EXIT_FAILURE);
	}

	return err;
}

int SiteData::getIndN() const {return nind;}

// unif draws uniform random number from [0,1]
double unif ()
{
	return (static_cast <double> (rand()) / (RAND_MAX) );
}


// poisson draws random number from Poisson(lambda)
/* I don't like the way this one is working - switching to std::poisson_distribution with better random number generation
int poisson (double lamda)
{
	double u = unif();
	double p = exp(-lamda);
	int i = 0;
	int x = i;

	while (u > p)
	{
		p += (lamda * p) / (i + 1);
		i++;
		x = i;
	}

	return x;
}
*/

// maxIndex returns the array index corresponding to the first occurrence of the largest element in the array
int maxIndex (const double arr [], const int size)
{
	double maxval = arr[0];
	int i = 0;

	for (i = 1; i < size; i++)
	{
		if (maxval < arr[i])
			maxval = arr[i];
	}

	for (i = 0; i < size; i++)
	{
		if (arr[i] == maxval)
		break;
	}

	return i;
}
