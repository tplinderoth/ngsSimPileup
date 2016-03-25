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
SiteData::SiteData (double depth, int n, double avgq, double maxq, double beta_b, int qoffset)
	: name ("chr"),
	  pos (0),
	  dupReadProb(0.0),
	  swapProb(0.0)

{
	cov = getCoverage(depth);
	nind = getInds(n);
	qmax = getMaxQ(maxq);
	getBetaPar(avgq, beta_b, betap);
	code = getCode(qoffset);
	realcov.resize(nind);
}

// SiteData::getCoverage assigns coverage
double SiteData::getCoverage (double depth)
{
	if (depth >= 0)
		return depth;
	else
	{
		fprintf(stderr, "Coverage for SiteData object <= 0\n");
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

// SiteData::genSeqData generates reads and quality scores for a site
std::vector<rdat> SiteData::genSeqData (const double afreq, const double pdup)
{
	if (afreq < 0.0 || afreq > 1.0)
	{
		fprintf(stderr, "Alternate allele frequency out of bounds in call to SiteData::genSeqData -> exiting\n");
		exit(EXIT_FAILURE);
	}

	int ind_geno [3] = {}; // # ref homo, # het, # alt homo
	getGeno(ind_geno, nind, afreq);
	int readn = 0; // number of reads for an individual
	int ind = 0;
	int q; // quality score
	std::vector<rdat> data (nind); // read data

	for (int i = 0; i < 3; i++) // go over all 3 genotype configurations
		for (int j = 0; j < ind_geno[i]; j++) // go over number of individuals in current genotype
		{
			if (pdup == 1.0)
			{
				readn = static_cast<int>(poisson(cov));
				realcov[ind] = readn;
			}
			else if (pdup < 1.0 && pdup >= 0)
				readn = floor(pdup * (static_cast<double>(realcov[ind]) + (dupReadProb * static_cast<double>(realcov[ind])) / (1.0 - dupReadProb)) + 0.5);
			else
			{
				fprintf(stderr, "Invalid admixture proportion in SiteData::genSeqData -> exiting\n");
				exit(EXIT_FAILURE);
			}
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
						if (unif() <= 0.5)
							data[ind].insert( rdat::value_type ( q, 0 ) );
						else
							data[ind].insert( rdat::value_type ( q, 1 ) );
					}
				}
			}

			ind++;
		}

	doError(&data); // introduce sequencing errors

	return data;
}

// SiteData::getGeno calculates genotype numbers based on HWE
void SiteData::getGeno (int genoNum [], const int n, const double f)
{
	double rhomo = (double)n * (1 - f) * (1 - f);
	double het = (double)n * 2 * f * (1 - f);
	double ahomo = (double)n * f * f;

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
void SiteData::assignSeqData (const double altfreq)
{
	seqdat.resize(nind);
	seqdat = genSeqData(altfreq);
}

// SiteData::doParalog generates reads and quality scores for a site comprised of paralogs - see simPileup.h
void SiteData::doParalog (std::vector<double>* altf, std::vector<double>* admixpro)
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

	seqdat = genSeqData((*altf)[0]); // get seq data for locus 1

	std::vector<rdat> dup;
	dupReadProb = getDupProb(admixpro);

	//generate data for duplicate loci
	for (unsigned int i = 1; i < altf->size(); i++)
	{
		dup = genSeqData((*altf)[i], (*admixpro)[i]); // get seq data for duplicate
		mergeLoci(&dup); //combine data with seqdat member
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
			{
				for (dat2 = ind2->begin(); dat2 != ind2->end(); ++dat2)
					seqdat[i].insert(*dat2);
			}

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

			x = unif(); // draw random number
			epsilon = pow(10.0, -1 * static_cast<double>(dat->first) / 10.0);

			if (x <= epsilon)
			{
				err_readp = unif(); // 1/3 probability of picking either of the wrong reads

				if (dat->second == 0) // ref change
				{
					if (err_readp <= (double)1/3)
						dat->second = 2; // ref error read1
					else if (err_readp > (double)2/3)
						dat->second = 3; // ref error read2
					else
						dat->second = 4; // ref error read3
				}
				else // alt change
				{
					if (err_readp <= (double)1/3)
						dat->second = 5; // alt error read1
					else if (err_readp > (double)2/3)
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
void SiteData::printParam (std::fstream& parfile, std::vector<double>* f, std::vector<double>* m)
{
	parfile.width(12);
	parfile << std::left << name << "\t";
	parfile.width(12);
	parfile << std::left << pos << "\t";
	unsigned int i = 0;
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
	for (i = 0; i < f->size(); i++)
	{
		if (i == f->size() - 1)
			parfile << (*f)[i] << "\n";
		else
			parfile << (*f)[i] << ";";
	}
}

// SiteData::randAllele randomly picks a reference allele
char SiteData::randAllele ()
{
	double rval = unif();
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
	double strand = unif();
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
	static char a [] = "CGT";
	static char c [] = "AGT";
	static char g [] = "ACT";
	static char t [] = "ACG";
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

		if (unif() > 0.5) // decide strand
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
			if (unif() > 0.5)
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

// unif draws uniform random number from [0,1]
double unif ()
{
	return (static_cast <double> (rand()) / (RAND_MAX) );
}


// poisson draws random number from Poisson(lambda)
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
