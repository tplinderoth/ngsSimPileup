// ngsSimPileup.cpp

#include <cstdlib>
#include <cstring>
#include <time.h>
#include <iostream>
#include "simPileup.h"
#include "generalUtils.h"
#include "ngsSimPile.h"

int main (int argc, char** argv)
{
	srand ((unsigned)time(NULL));

	// declare and initialize variables
	std::vector<double> altfreq; // alternate allele frequencies for loci
	std::vector<double> admix; // probability that read comes from locus
	double minfreq = 0.0; // minimum random alternate allele frequency
	double maxfreq = 1.0; // maximum random alternate allele frequency
	double coverage = 10.0; // average individual sequencing depth
	int nind = 0; // number diploid individuals
	double avgqual = 38.0; // average Phred quality score
	double maxqual = 41.0; // maximum Phred quality score
	double betab = 0.3; // beta shape parameter for beta distribution on quality scores
	int qcode = 33; // quality score offset (added to quality score in pileup format)
	std::string seqname("chr1"); // chromosome name
	int paralog_sites = 0; // number of paralogous sites to simulate
	int normal_sites = 0; // number of non-paralogous sites to simulate
	std::string outfile; // file to dump results to

	if (argc==1)
	{
		info (coverage, minfreq, maxfreq, avgqual, maxqual, betab, qcode, seqname);
		return 0;
	}

	// read-in, check, and assign inputs
	int validin = parseInputs (argc, argv, altfreq, admix, minfreq, maxfreq, coverage, nind, avgqual, maxqual, betab, qcode,
			seqname, paralog_sites, normal_sites, outfile);
	if (validin)
		return 0; // terminate

	// open output files
	std::fstream pileout; // create fstream object for pileup results
	getFILE(pileout, outfile.c_str(), "out"); // open pileup file
	std::string paramsf = outfile;
	paramsf.append(".param");
	std::fstream parout; // create fstream object for parameters
	getFILE(parout, paramsf.c_str(), "out"); // open parameter file
	std::cerr << "Dumping results to:\n" << outfile << "\n" << paramsf << "\n";

	// initialize pileup object
	SiteData sdat(coverage, nind, avgqual, maxqual, betab, qcode);
	sdat.name = seqname;
	sdat.pos = 1;

	// generate sequencing data

	if (normal_sites > 0) // simulate normal sites
		simNormal (normal_sites, minfreq, maxfreq, &admix, &sdat, pileout, parout);

	if (paralog_sites > 0) // simulate paralogous sites
		simParalog (&altfreq, &admix, paralog_sites, minfreq, maxfreq, &sdat, pileout, parout);

	// prepare to exit program

	pileout.close();
	parout.close();

	std::cerr << "Finished!\n";

	return 0;
}

int parseInputs (int argc, char** argv, std::vector<double>& altfreq, std::vector<double>& admix, double& minfreq, double& maxfreq,
	double& coverage, int& nind, double& avgqual, double& maxqual, double& betab, int& qcode,
	std::string& seqname, int& paralog_sites, int& normal_sites, std::string& outfile)
{
	// read and assign input

	int argpos = 1;
	double f = -1.0;
	double m = -1.0;

		while (argpos < argc)
		{
			if ( strcmp(argv[argpos], "-help") == 0 )
			{
				info (coverage, minfreq, maxfreq, avgqual, maxqual, betab, qcode, seqname, 1);
				return 1;
			}
			else if ( strcmp(argv[argpos], "-altfreq") == 0 ) // parse alternate allele frequencies
			{
				argpos++;
				while ( *argv[argpos] != '-')
				{
					if ( isdigit (argv[argpos][0]) ) //*argv[argpos]
					{
						f = atof(argv[argpos]);
						if (f >= 0.0 && f <= 1.0)
							altfreq.push_back(f);
						else
						{
							fprintf(stderr, "\nNumeric argument to -altfreq not in [0,1] -> exiting\n");
							return 1; // terminate
						}
					}
					else if ( strcmp(argv[argpos], "R") == 0 || strcmp(argv[argpos], "r") == 0)
						altfreq.push_back( 82.0 ); // ascii R = 82
					else
					{
						fprintf(stderr, "\nInvalid argument to -altfreq -> exiting\n");
						return 1; // terminate
					}
					argpos++;
					if (argpos >= argc)
						break;
				}
				altfreq.resize(altfreq.size());
				continue;
			}
			if ( strcmp(argv[argpos], "-admix") == 0) // parse read admixture proportions
			{
				argpos++;
				while ( *argv[argpos] != '-')
				{
					if ( isdigit (*argv[argpos]) )
					{
						m = atof(argv[argpos]);
						if (m >= 0.0 && m <= 1.0)
							admix.push_back(m);
						else
						{
							fprintf(stderr, "\nArgument to -admix not in [0,1] -> exiting\n");
							return 1; // terminate
						}
					}
					else
					{
						fprintf(stderr, "\nInvalid argument to -admix -> exiting\n");
						return 1; // terminate
					}
				argpos++;
				if (argpos >= argc)
					break;
				}
				admix.resize(admix.size());
				continue;
			}
			else if ( strcmp(argv[argpos], "-minfreq") == 0 )
				minfreq = atof(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-maxfreq") == 0 )
				maxfreq = atof(argv[argpos + 1]);
			else if ( strcmp(argv[argpos], "-coverage") == 0 )
				coverage = atof( argv[argpos + 1]);
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
			else if ( strcmp( argv[argpos], "-paralog_sites") == 0)
				paralog_sites = atoi(argv[argpos + 1]);
			else if ( strcmp( argv[argpos], "-normal_sites") == 0 )
				normal_sites = atoi(argv[argpos + 1]);
			else if ( strcmp ( argv[argpos], "-seqname") == 0 )
				seqname = argv[argpos + 1];
			else if ( strcmp( argv[argpos], "-outfile") == 0 )
				outfile = argv[argpos + 1];
			else
			{
				fprintf(stderr, "\nUnknown argument: %s\n", argv[argpos]);
				return 1;
			}
			argpos += 2;
		}

		// check user inputs

		if (minfreq < 0.0 || minfreq > 1.0)
		{
			fprintf(stderr, "\n-minfreq out of bounds -> exiting\n");
			return 1;
		}
		if (maxfreq < 0.0 || maxfreq > 1.0)
		{
			fprintf(stderr, "\n-maxfreq out of bounds -> exiting\n");
			return 1;
		}
		if (minfreq > maxfreq)
		{
			fprintf(stderr, "\n-minfreq greater than -maxfreq -> exiting\n");
			return 1;
		}
		if (coverage <= 0)
		{
			fprintf(stderr, "\n-coverage must greater than zero -> exiting\n");
			return 1;
		}
		if (nind <= 0 )
		{
			fprintf(stderr, "\n-nind must be a positive integer -> exiting\n");
			return 1;
		}
		if (avgqual <= 0)
		{
			fprintf(stderr, "\n-avgqual must be a positive integer -> exiting\n");
			return 1;
		}
		if ( qcode != 33)
			if ( qcode != 64)
			{
				fprintf(stderr, "\n-qcode must be 33 or 64 -> exiting\n");
				return 1;
			}
		if (maxqual <= 0 || maxqual < avgqual)
		{
			fprintf(stderr, "\n-maxqual must be a positive integer at least as great as avgqual -> exiting\n");
			return 1;
		}
		if ( qcode + maxqual > 126 )
		{
			fprintf(stderr, "\noffset maximum quality score exceeds 126 -> exiting\n");
			return 1;
		}
		if ( betab <= 0 )
		{
			fprintf(stderr, "\n-betab must be greater than zero -> exiting\n");
			return 1;
		}
		if ( paralog_sites < 0)
		{
			fprintf(stderr, "\n-paralog_sites less than zero -> exiting\n");
			return 1;
		}
		if ( normal_sites < 0)
		{
			fprintf(stderr, "\n-normal_sites less than zero -> exiting\n");
			return 1;
		}
		if ( paralog_sites + normal_sites <= 0)
		{
			fprintf(stderr, "\nNo sites to simulate -> exiting\n");
			return 1;
		}
		if ( paralog_sites > 0 && altfreq.empty())
		{
			fprintf(stderr, "\nMust supply -altfreq for paralogous sites -> exiting\n");
			return 1;
		}
		if (!admix.empty())
		{
			if (!altfreq.empty() && admix.size() != altfreq.size())
			{
				fprintf(stderr, "\nNumber of arguments for -altfreq and -admix differ -> exiting\n");
				return 1;
			}
			if (altfreq.size() == 1 || paralog_sites == 0)
			{
				if (admix[0] != 1.0)
				{
					admix[0] = 1.0;
					fprintf(stderr, "\nNo duplicates, so admixture proportion set to 1.0\n");
				}
			}
			double admixsum = 0;
			for (std::vector<double>::const_iterator i = admix.begin(); i != admix.end(); ++i)
				admixsum += *i;
			if (admixsum != 1.0)
			{
				fprintf(stderr, "-admix arguments must sum to 1 -> exiting\n");
				return 1;
			}
		}
		else
		{
			fprintf(stderr, "Must supply admixture proportions with -admix -> exiting\n");
			return 1;
		}
		if ( outfile.empty())
		{
			fprintf(stderr, "\nMust supply -outfile -> exiting\n");
			return 1;
		}

		return 0;
}

void simNormal (int nsites, double lbound, double ubound, std::vector<double>* mvec, SiteData* dat, std::fstream& seqfile, std::fstream& parfile)
{
	std::vector<double> f (1);

	for (int i = 0; i < nsites; i++)
	{
		f[0] = decimalUnifBound(lbound, ubound);
		dat->assignSeqData(f[0]);
		dat->printSite(seqfile, dat->code);
		dat->printParam(parfile, &f, mvec); // print to parameter file
		dat->pos++;
	}
}

void simParalog (std::vector<double>* fvec, std::vector<double>* mvec, int nsites, double lbound, double ubound, SiteData* dat, std::fstream& seqfile, std::fstream& parfile)
{
	std::vector<double> infreq (fvec->size());
	unsigned int locus = 0;

	for (int i = 0; i < nsites; i++)
	{
		for (locus = 0; locus < fvec->size(); locus++)
		{
			if ((*fvec)[locus] == 82.0) // random
				infreq[locus] = decimalUnifBound(lbound, ubound);
			else
				infreq[locus] = (*fvec)[locus];
		}

		dat->doParalog(&infreq, mvec);
		dat->printSite(seqfile, dat->code); // print pileup file
		dat->printParam(parfile, &infreq, mvec); // print to parameter file
		dat->pos++;
	}
}

void info (const double& cov, const double& minfreq, const double&maxfreq, const double& avgqual, const double& maxqual,
		const double& beta, const double& encode, const std::string& name, int help, const char* v)
{
	std::cerr << "\nngsSimPileup version " << v << "\n\nInput:\n"
	<< "\n-nind INT number of diploid individuals"
	<< "\n-coverage DOUBLE average per-site sequencing depth [10.0]"
	<< "\n-normal_sites INT number of normal sites to simulate"
	<< "\n-paralog_sites INT number of paralogous sites to simulate"
	<< "\n-minfreq DOUBLE minimum alternate allele frequency [0]"
	<< "\n-maxfreq DOUBLE maximum alternate allele frequency [1.0]"
	<< "\n-altfreq DOUBLE|R vector of alternate allele frequencies for paralogous sites (use 'R' for random frequency)"
	<< "\n-admix DOUBLE vector of probabilities that a read comes from the locus with given alternate allele frequency"
	<< "\n-avgqual DOUBLE average Phred read quality [38.0]"
	<< "\n-maxqual DOUBLE maximum Phred read quality [41.0]"
	<< "\n-betab DOUBLE shape parameter > 0 for quality score distribution (beta dist) [0.3]"
	<< "\n-qcode 33|64 quality score encoding; 64 for illumina v1.3+ or 33 for illumina v1.8+ [33]"
	<< "\n-seqname STRING name of chromosome/contig [chr1]"
	<< "\n-outfile FILE output file"
	<< "\n\noutput:\n"
	<< "\n.param file fields:"
	<< "\n(1) sequence name (2) position (3) probability read comes from locus (4) alternate allele frequency for locus\n"
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
	<< "a larger value makes the distribution more 'normal' around the average quality score\n\n";
	}
}
