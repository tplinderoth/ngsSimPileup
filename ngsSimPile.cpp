// ngsSimPileup.cpp

/*
 * TO DO:
 * (1) add exception handling
 * (2) add theta parameter to number SNPs (maybe...?)
 * (3) Random fitness and inbreeding values - and bounds for random variables
 *
 * SFS sampling method is untested
 */

#include <cstdlib>
#include <cstring>
#include <time.h>
#include <iostream>
#include "simPileup.h"
#include "ngsSimPile.h"

int main (int argc, char** argv)
{

	// declare and initialize variables
	std::vector<double> altfreq; // alternate allele frequencies for loci
	std::vector<double> admix; // probability that read comes from locus
	double deletionf = 0.0; // CNV delection frequency
	double inbreed = 0.0; // inbreeding coefficient for locus
	std::vector< Array<double> > fitness; // relative fitness values for [AA, Aa, aa] genotypes
	double minfreq = 0.0; // minimum random alternate allele frequency
	double maxfreq = 1.0; // maximum random alternate allele frequency
	double minF = 0.0; // minimum inbreeding coefficient
	double maxF = 1.0; // maximum inbreeding coefficient
	double covmean = 10.0; // mean for normally distributed average individual coverage
	double covstddev = 1.0; // standard deviation for normally distributed average individual coverage
	int nind = 0; // number diploid individuals
	double avgqual = 38.0; // average Phred quality score
	double maxqual = 41.0; // maximum Phred quality score
	double betab = 0.3; // beta shape parameter for beta distribution on quality scores
	int qcode = 33; // quality score offset (added to quality score in pileup format)
	std::string seqname("chr1"); // chromosome name
	double theta = 0.01; // population-size scaled mutation rate
	unsigned int nsites = 0; // number of sites to simulate
	std::string outfile; // file to dump results to
	bool foldsfs = false; // use the folded SFS to draw allele frequencies from

	if (argc==1)
	{
		info (covmean, covstddev, minfreq, maxfreq, avgqual, maxqual, betab, qcode, seqname, &inbreed, foldsfs, deletionf);
		return 0;
	}

	// read-in, check, and assign inputs
	int validin = parseInputs(argc, argv, &altfreq, &admix, minfreq, maxfreq, &fitness, &inbreed, minF, maxF, covmean, covstddev, nind,
			avgqual, maxqual, betab, qcode, seqname, theta, nsites, outfile, foldsfs, deletionf);
	if (validin > 1)
		return 0;
	else if (validin)
	{
		fprintf(stderr, "--> exiting\n");
		return 1; // terminate
	}

	// open output files
	std::fstream pileout; // create fstream object for pileup results
	getFILE(pileout, outfile.c_str(), "out"); // open pileup file
	std::string paramsf = outfile;
	paramsf.append(".param");
	std::fstream parout; // create fstream object for parameters
	getFILE(parout, paramsf.c_str(), "out"); // open parameter file
	std::cerr << "Dumping results to:\n" << outfile << "\n" << paramsf << "\n";

	// initialize pileup object
	SiteData sdat(covmean, covstddev, nind, avgqual, maxqual, betab, qcode);
	sdat.name = seqname;
	sdat.pos = 1;

	// generate sequencing data

	if (altfreq.size() > 1) // simulate paralogous sites
	{
		fprintf(stderr, "simulating %u paralogous sites\n", nsites);
		simParalog(nsites, &altfreq, minfreq, maxfreq, &admix, inbreed, &fitness, &sdat, pileout, parout, foldsfs);
	}
	else // simulate normal sites
	{
		fprintf(stderr, "simulating %u non-paralogous sites\n", nsites);
		simNormal(nsites, altfreq[0], minfreq, maxfreq, &admix, inbreed, &fitness, &sdat, pileout, parout, foldsfs, deletionf);
	}

	// prepare to exit program

	pileout.close();
	parout.close();

	std::cerr << "Finished!\n";

	return 0;
}

int parseInputs (int argc, char** argv, std::vector<double>* altfreq, std::vector<double>* admix, double& minfreq, double& maxfreq,
		std::vector< Array<double> >* fitness, double* inbreed, double& minF, double& maxF, double& coverage, double& covsd,
		int& nind,double& avgqual, double& maxqual, double& betab, int& qcode, std::string& seqname, double& theta, unsigned int& nsites, std::string& outfile, bool& fold, double& deletionf)
{
	// read and assign input

	int argpos = 1;

	// set up container to temporary hold parameter values
	std::vector<char const*> tmpcontainer;
	tmpcontainer.reserve(12);

		while (argpos < argc)
		{
			if ( strcmp(argv[argpos], "-help") == 0 )
			{
				info(coverage, covsd, minfreq, maxfreq, avgqual, maxqual, betab, qcode, seqname, inbreed, fold, 1);
				return 2;
			}
			else if ( strcmp(argv[argpos], "-altfreq") == 0 ) // parse alternate allele frequencies
			{
				++argpos;
				while (argpos < argc)
				{
					if (argv[argpos][0] == '-' && static_cast<int>(argv[argpos][0]) + static_cast<int>(argv[argpos][1]) > 102) break;
					tmpcontainer.push_back(argv[argpos]);
					++argpos;
				}
				if(setAltFreq(&tmpcontainer, altfreq))
					return 1;
				tmpcontainer.clear();
				continue;
			}
			else if ( strcmp(argv[argpos], "-admix") == 0) // parse read admixture proportions
			{
				++argpos;
				while (argpos < argc)
				{
					if (argv[argpos][0] == '-' && static_cast<int>(argv[argpos][0]) + static_cast<int>(argv[argpos][1]) > 102) break;
					tmpcontainer.push_back(argv[argpos]);
					++argpos;
				}
				if (setAdmix(&tmpcontainer, admix))
					return 1;
				tmpcontainer.clear();
				continue;
			}
			else if (strcmp(argv[argpos], "-fitness") == 0) // parse relative fitnesses
			{
				++argpos;
				while (argpos < argc)
				{
						if (argv[argpos][0] == '-' && static_cast<int>(argv[argpos][0]) + static_cast<int>(argv[argpos][1]) > 102) break;
						tmpcontainer.push_back(argv[argpos]);
						++argpos;
				}
				if(setFitness(&tmpcontainer, fitness))
					return 1;
				tmpcontainer.clear();
				continue;
			}
			else if (strcmp(argv[argpos], "-inbreed") == 0) // parse inbreeding coefficient
			{
				if(setInbreeding(argv[argpos+1], inbreed))
					return 1;
			}
			else if ( strcmp(argv[argpos], "-minfreq") == 0 )
				minfreq = atof(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-maxfreq") == 0 )
				maxfreq = atof(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-minF") == 0 )
				minF = atof(argv[argpos+1]);
			else if ( strcmp(argv[argpos], "-maxF") == 0)
				maxF = atof(argv[argpos+1]);
			else if ( strcmp(argv[argpos], "-coverage") == 0 )
				coverage = atof(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-covstddev") == 0 )
				covsd = atof(argv[argpos+1]);
			else if ( strcmp(argv[argpos], "-nind") == 0 )
				nind = atoi(argv[argpos + 1] );
			else if ( strcmp( argv[argpos], "-avgqual") == 0 )
				avgqual = atof(argv[argpos + 1]);
			else if ( strcmp( argv[argpos], "-maxqual") == 0 )
				maxqual = atof(argv[argpos + 1]);
			else if ( strcmp( argv[argpos], "-betab") == 0 )
				betab = atof(argv[argpos + 1]);
			else if ( strcmp( argv[argpos], "-qcode") == 0 )
				qcode = atoi(argv[argpos + 1]);
			else if ( strcmp( argv[argpos], "-nsites") == 0)
				nsites = atoi(argv[argpos + 1]);
			else if ( strcmp ( argv[argpos], "-seqname") == 0 )
				seqname = argv[argpos + 1];
			else if ( strcmp(argv[argpos], "-theta") == 0)
				theta = atof(argv[argpos+1]);
			else if ( strcmp( argv[argpos], "-outfile") == 0 )
				outfile = argv[argpos + 1];
			else if ( strcmp(argv[argpos], "-foldsfs") == 0)
				fold = atoi(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-cnvdel") == 0)
				deletionf = atof(argv[argpos+1]);
			else
			{
				fprintf(stderr, "\nUnknown argument: %s\n", argv[argpos]);
				return 1;
			}
			argpos += 2;
		}

		// set default parameter values
		if (setDefaultVals(admix, altfreq, inbreed, fitness))
			return 1;


		// check user inputs
		if (minfreq < 0.0 || minfreq > 1.0)
		{
			fprintf(stderr, "\n-minfreq out of boundsn");
			return 1;
		}
		if (maxfreq < 0.0 || maxfreq > 1.0)
		{
			fprintf(stderr, "\n-maxfreq out of bounds\n");
			return 1;
		}
		if (minF < 0.0 || minF > 1.0)
		{
			fprintf(stderr, "\n-minF out of bounds\n");
			return 1;
		}
		if (maxF < 0.0 || maxF > 1.0)
		{
			fprintf(stderr, "\n-maxF out of bounds\n");
			return 1;
		}
		if (minfreq > maxfreq)
		{
			fprintf(stderr, "\n-minfreq greater than -maxfreq\n");
			return 1;
		}
		if (coverage < 0)
		{
			fprintf(stderr, "\n-coverage cannot be negative\n");
			return 1;
		}
		if (covsd < 0)
		{
			fprintf(stderr, "\n-covstddev cannot be negative\n");
			return 1;
		}
		if (nind <= 0 )
		{
			fprintf(stderr, "\n-nind must be a positive integer\n");
			return 1;
		}
		if (avgqual <= 0)
		{
			fprintf(stderr, "\n-avgqual must be a positive integer\n");
			return 1;
		}
		if ( qcode != 33)
			if ( qcode != 64)
			{
				fprintf(stderr, "\n-qcode must be 33 or 64\n");
				return 1;
			}
		if (maxqual <= 0 || maxqual < avgqual)
		{
			fprintf(stderr, "\n-maxqual must be a positive integer at least as great as avgqual\n");
			return 1;
		}
		if ( qcode + maxqual > 126 )
		{
			fprintf(stderr, "\noffset maximum quality score exceeds 126\n");
			return 1;
		}
		if ( betab <= 0 )
		{
			fprintf(stderr, "\n-betab must be greater than zero\n");
			return 1;
		}
		if ( nsites < 1)
		{
			fprintf(stderr, "\n-nsites should be > 0\n");
			return 1;
		}
		if (theta <= 0)
		{
			fprintf(stderr, "\n-theta must be > 0\n");
			return 1;
		}
		if ( outfile.empty())
		{
			fprintf(stderr, "\nMust supply -outfile\n");
			return 1;
		}
		if (fold != 0)
		{
			if (fold != 1)
			{
				fprintf(stderr, "\n-fold must be 0 or 1\n");
				return 1;
			}
		}
		if (deletionf > 1.0 || deletionf < 0.0) {
			fprintf(stderr,"\n-cnvdel out of range [0,1]\n");
			return 1;
		}
		if (deletionf + minfreq > 1) {
			fprintf(stderr, "Sum of allele frequencies > 1: lower -cnvdel or -minfreq\n");
			return 1;
		}
		if (deletionf + maxfreq > 1) {
			fprintf(stderr, "Sum of deletion and maximum SNP frequency > 1: consider lowering -cnvdel or -maxfreq\n");
		}

		return 0;
}


int setAltFreq (std::vector<const char*>* fchar, std::vector<double>* fvals)
{
	if (fchar->size() < 1)
	{
		fprintf(stderr, "No allele frequencies specified in call to setAltFreq\n");
		return 1;
	}

	double f = 0.0;
	fvals->resize(fchar->size());
	for (unsigned int i = 0; i < fvals->size(); ++i)
	{
		if ( isdigit((*fchar)[i][0]) )
		{
			f = atof((*fchar)[i]);
			if (f >= 0.0 && f <= 1.0)
				(*fvals)[i] = f;
			else
			{
				fprintf(stderr, "\nNumeric argument to -altfreq not in [0,1]\n");
				return 1;
			}
		}
		else if ( strcmp((*fchar)[i], "R") == 0 || strcmp((*fchar)[i], "r") == 0)
			(*fvals)[i] = 82.0; // ascii R = 82
		else if ( strcmp((*fchar)[i], "S") == 0 || strcmp((*fchar)[i], "s") == 0)
			(*fvals)[i] = 83.0; //ascii S = 83
		else
		{
			fprintf(stderr, "\nInvalid argument to -altfreq\n");
			return 1;
		}
	}

	return 0;
}

int setAdmix (std::vector<const char*>* achar, std::vector<double>* avals)
{
	if (achar->size() < 1)
	{
		fprintf(stderr, "No admixture proportions specified in call to setAdmix\n");
		return 1;
	}

	double m = 0.0;
	avals->resize(achar->size());
	for (unsigned int i = 0; i < avals->size(); ++i)
	{
		if ( isdigit((*achar)[i][0]) )
		{
			m = atof((*achar)[i]);
			if (m >= 0.0 && m <= 1.0)
				(*avals)[i] = m;
			else
			{
				fprintf(stderr, "Values passed to -admix must be in [0,1]\n");
				return 1;
			}
		}
	}

	return 0;
}

int setFitness (std::vector<const char*>* wchar, std::vector< Array<double> >* wvals)
{
	if (wchar->size() % 3)
	{
		fprintf(stderr, "Number of supplied fitness values not a multiple of 3 in call to setFitness\n");
		return 1;
	}

	wvals->resize(wchar->size()/3);
	double w = 0.0;
	unsigned int k = 0;
	for (unsigned int i = 0; i < wvals->size(); ++i)
	{
		(*wvals)[i].setSize(3);
		for (unsigned int j = 0; j < 3; ++j)
		{
			if ( isdigit((*wchar)[k][0]) )
			{
				w = atof((*wchar)[k]);
				if (w >= 0.0 && w <= 1.0)
					(*wvals)[i][j] = w;
				else
				{
					fprintf(stderr, "\nfitness values must be in range [0,1]\n");
					return 1;
				}
			}
			else
			{
				fprintf(stderr, "\nfitness values must be numeric\n");
				return 1;
			}
			++k;
		}
	}

	return 0;
}

int setInbreeding(const char* Fchar, double* Fval)
{
	if (Fchar)
	{
		if ( isdigit(Fchar[0]))
		{
			*Fval = atof(Fchar);
			if (*Fval > 1.0 || *Fval < 0.0)
			{
				fprintf(stderr, "\ninbreeding coefficient outside of range [0,1]\n");
				return 1;
			}
		}
		else if (strcmp(Fchar, "R") == 0 || strcmp(Fchar, "r") == 0)
		{
			*Fval = 82.0; // ascii R = 82
		}
		else
		{
			fprintf(stderr, "\nInvalid argument to -inbreed\n");
			return 1;
		}
	}
	else
	{
		fprintf(stderr, "No inbreeding value provided in call to setInbreeding\n");
		return 1;
	}

	return 0;
}

int setDefaultVals (std::vector<double>* admix, std::vector<double>* altfreq, double* inbreed, std::vector< Array<double> >* fitness)
{
	unsigned int i = 0;

	// default alternative frequency = random
	if (altfreq->size() < 1)
	{
		altfreq->push_back(82.0);
		fprintf(stderr, "Assuming random allele frequencies for single copy sites\n");
	}

	// default admixture proportion = 1/(number loci)
	if (admix->size() < 1)
	{
		admix->resize(altfreq->size());
		if (altfreq->size() > 1)
			fprintf(stderr, "Assuming equal contribution of reads from each copy\n");

		for (i=0; i < admix->size(); ++i)
			(*admix)[i] = 1.0/static_cast<double>(admix->size());
	}

	if (admix->size() != altfreq->size())
	{
		fprintf(stderr, "Number of admixture proportions and loci differ\n");
		return 1;
	}

	if (vecsum(*admix) != 1.0)
	{
		fprintf(stderr, "Admixture proportions must sum to 1\n");
		return 1;
	}

	// default F = 0
	if (!inbreed)
		*inbreed = 0.0;

	// default fitness is uniform (i.e. no differential fitness)
	if (fitness->size() < 1)
		fitness->resize(altfreq->size());

	if (fitness->size() != altfreq->size())
	{
		fprintf(stderr, "(Number of fitness values)/3 does not equal the number of loci\n");
		return 1;
	}

	for (i=0; i < altfreq->size(); ++i)
	{
		if ((*fitness)[i].size() < 3)
		{
			(*fitness)[i].setSize(3);
			for (unsigned int j=0; j < 3; ++j)
				(*fitness)[i][j] = 1.0;
		}
	}

	return 0;
}

void simNormal (unsigned int nsites, double altfreq, double lbound, double ubound, std::vector<double>* mvec, const double inbreed,
		std::vector< Array<double> >* fitness, SiteData* dat, std::fstream& seqfile, std::fstream& parfile, bool foldsfs, double delcnvf)
{
	// ADD EXCEPTION HANDLING
	std::vector<double> f (1);
	SFS sfs(static_cast<unsigned int>(dat->getIndN()), foldsfs);

	for (unsigned int i = 0; i < nsites; ++i)
	{
		f[0] = alleleFreq(altfreq, lbound, ubound, &sfs, delcnvf);
		dat->assignSeqData(f[0], inbreed, &(*fitness)[0], delcnvf);
		dat->printSite(seqfile, dat->code);
		dat->printParam(parfile, &f, mvec, inbreed, fitness); // print to parameter file
		dat->pos++;
	}
}

void simParalog (unsigned int nsites, std::vector<double>* fvec, double lbound, double ubound, std::vector<double>* mvec, const double inbreed,
		std::vector< Array<double> >* fitness, SiteData* dat, std::fstream& seqfile, std::fstream& parfile, bool foldsfs)
{
	std::vector<double> infreq (fvec->size());
	unsigned int locus = 0;
	SFS sfs(static_cast<unsigned int>(dat->getIndN()), foldsfs);

	for (unsigned int i = 0; i < nsites; ++i)
	{
		for (locus = 0; locus < fvec->size(); ++locus)
			infreq[locus] = alleleFreq((*fvec)[locus], lbound, ubound, &sfs);

		dat->doParalog(&infreq, inbreed, fitness, mvec);
		dat->printSite(seqfile, dat->code); // print pileup file
		dat->printParam(parfile, &infreq, mvec, inbreed, fitness); // print to parameter file
		dat->pos++;
	}
}

double alleleFreq (double freq, double lobound, double upbound, SFS* sfs, double cnvf)
{
	// ADD EXCEPTION HANDLING
	double f;
	int sample_thresh = 1000;
	int n = 0;

	if (freq == 82.0) {
		do {
			f = decimalUnifBound(lobound, upbound);
		} while (f+cnvf > 1);
	}
	else if (freq == 83.0)
	{
		if (upbound == lobound)
		{
			fprintf(stderr, "Allele frequency range too narrow for SFS sampling method\n");
			exit(EXIT_FAILURE);
		}
		n = 0;
		do
		{
			f = sfs->rfreq();
			++n;
			if (n > sample_thresh)
				fprintf(stderr, "SFS sampling attempts exceed %i -> consider widening allele frequency bounds\n", sample_thresh);
		} while (f > upbound || f < lobound);
	}
	else
	{
		if (freq >= 0.0 && freq <= 1.0)
			f = freq;
		else
		{
			fprintf(stderr, "Allele frequency %f out of [0,1] range\n", freq);
			exit(EXIT_FAILURE);
		}
	}

	return f;
}

void info (const double& covmean, const double& covsd, const double& minfreq, const double& maxfreq, const double& avgqual, const double& maxqual, const double& beta,
		const double& encode, const std::string& name, const double* F, const bool& foldsfs, double cnvdel, int help, const char* v)
{
	std::cerr << "\nngsSimPileup version " << v << "\n\nInput:\n"
	<< "\n-nind INT number of diploid individuals"
	<< "\n-coverage DOUBLE mean average per-site sequencing depth [" << covmean << "]"
	<< "\n-covstddev DOUBLE standard deviation of the average per-site sequencing depth [" << covsd << "]"
	<< "\n-nsites INT number of sites to simulate"
	<< "\n-minfreq DOUBLE minimum alternate allele frequency [" << minfreq << "]"
	<< "\n-maxfreq DOUBLE maximum alternate allele frequency [" << maxfreq << "]"
	<< "\n-altfreq DOUBLE|R|S vector of alternate allele frequencies for paralogous sites ('R': random frequency, F: draw frequency from SFS)"
	<< "\n-cnvdel DOUBLE frequency of deletion CNV [" << cnvdel << "]"
	<< "\n-foldsfs 0|1 draw allele frequencies from unfolded (0) or folded (1) SFS [" << foldsfs << "]"
	<< "\n-admix DOUBLE vector of probabilities that a read comes from the locus with given alternate allele frequency [1/number_copies]"
	<< "\n-inbreed DOUBLE inbreeding coefficient (F) [" << *F << "]"
	<< "\n-fitness DOUBLE vector of fitness values ranging [0,1] for reference homozygotes, heterozygotes, and alternate homozygotes for each locus [equal fitness]"
	<< "\n-avgqual DOUBLE average Phred read quality [" << avgqual << "]"
	<< "\n-maxqual DOUBLE maximum Phred read quality [" << maxqual << "]"
	<< "\n-betab DOUBLE shape parameter > 0 for beta-distributed quality scores [" << beta << "]"
	<< "\n-qcode 33|64 quality score encoding; 64 for illumina v1.3+ or 33 for illumina v1.8+ [" << encode << "]"
	<< "\n-seqname STRING name of chromosome/contig [" << name << "]"
	<< "\n-outfile FILE output file"
	<< "\n\noutput:\n"
	<< "\n.param file fields:"
	<< "\n(1) sequence name (2) position (3) probability read comes from locus (4) alternate allele frequency (5) genotype fitnesses (6) inbreeding coefficient\n"
	<< "\n";

	if (help)
	{
	std::cerr << "Extended help:\n\n"
	<< "For -altfreq an 'R' or 'r' argument will generate a random allele frequency\n"
	<< "ex: -altfreq 0 1 R 0.35 <= paralogous site comprised of loci having alternative allele frequencies 0.0, 1.0, random, and 0.35\n\n"
	<<	"For -admix the order of read probabilities corresponds to the loci in the -altfreq vector\n"
	<< "-altfreq 0 1 R 0.35 -admix 0.5 0.2 0.1 0.2 <= paralogous site comprised of:\n"
	<< "locus1: alt allele freq of 0.0 and P(read comes from locus 1) = 0.5\n"
	<< "locus2: alt allele freq of 1.0 and P(read comes from locus 2) = 0.2\n"
	<< "locus3: random alt allele freq and P(read comes from locus 3) = 0.1\n"
	<< "locus4: alt allele freq of 0.35 and P(read comes from locus 4) = 0.2\n\n"
	<< "-betab provides the beta shape parameter for a beta distribution which models the quality score distribution:\n"
	<< "a smaller value skews the distribution towards the maximum quality score\n"
	<< "a larger value makes the distribution more 'normal' around the average quality score\n\n"
	<< "-fitness 1.0 1.0 1.0 0.4 1.0 0.6 would indicate that all genotypes for locus 1 have equal fitness and heterozygote advantage at locus 2\n";
	}
}
