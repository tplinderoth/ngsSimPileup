// SFS.cpp

#include "SFS.h"
#include "generalUtils.h"

SFS::SFS (unsigned int n, bool fold)
	: _nind(n),
	  _ncat(fold ? _nind : 2*_nind-1)
{
	setProbs();
}

void SFS::setProbs ()
{
	unsigned int i=0;
	double denom = 0.0;
	_probs.resize(_ncat);

	for (double i = 1.0; i <= _ncat; ++i)
		denom += 1.0/i;

	for (i=0; i < _probs.size(); ++i)
		_probs[i] = 1.0/(i*denom);

	_scaledprobs.resize(_ncat);
	for (i=0; i < _scaledprobs.size(); ++i)
		_scaledprobs[i] = _probs[i]/_probs[0];
}

double SFS::rfreq ()
{
	double univar = decimalUnifBound(0.0, 1.0);
	unsigned int sfscat = _scaledprobs.size()-1;

	while (univar > _scaledprobs[sfscat])
		--sfscat;

	return((sfscat+1)/static_cast<double>((_ncat+1)));
}

double SFS::getProb (unsigned int category)
{
	if (category == 0 || category == 2*_nind+1)
	{
		fprintf(stderr, "SFS does not include fixed categories\n");
		return 0.0;
	}
	return _probs[category];
}
