// generalUtils.cpp
// version 0.0.5

#include <cstdio>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <sstream>
#include "generalUtils.h"

// getCString converts string to C-Style pointer string
char * getCString (std::string s)
{
 	int slen = s.length();
	char * cstring = new char[slen + 1];
	s.copy(cstring, slen, 0);
	cstring[ slen ] = '\0';

	return cstring;
}

// getFILE is a wrapper for getting files
void getFILE(std::fstream &fp, const char* fname, const char* mode)
{
	int writeFile = 0;
	if (strcmp(mode, "out") == 0)
	{
		writeFile = 1;

		if(writeFile && fexists(fname))
		{
			fprintf(stderr,"File already exists: %s -> exiting...\n",fname);
			exit(0);
		}

		fp.open(fname, std::ios::out);
	}
	else if (strcmp(mode, "app") == 0)
		fp.open(fname, std::ios::app);
	else if (strcmp(mode, "in") == 0)
		fp.open(fname, std::ios::in);

	if( !fp )
	{
		fprintf(stderr,"Error opening FILE handle for file: %s -> exiting...\n",fname);
		fp.close();
		exit(0);
	}
}

// fexists finds out if a file exists
int fexists(const char* str)
{
	struct stat buffer;
 	return (stat(str, &buffer )==0 );
}

// readFileChunk chunks in file
void readFileChunk (std::fstream &infile, std::vector<std::string>& datavec, unsigned int* chunk, int* end)
{

	unsigned int i = 0;
	std::string(line);
	int cap = 0;

	if (!datavec.empty())
	{
		cap = datavec.capacity();
		datavec.clear();
		datavec.resize(cap);
	}

	while (i < *chunk)
	{
		if (getline(infile, line))
		{
			if (i < datavec.capacity())
				datavec[i] = line;
			else
				datavec.push_back(line);
			i++;
		}
		else
		{
			datavec.resize(i); // remove empty elements
			*end = 1;
			*chunk = i; // number of lines processed upon returning
			return;
		}
	}
}

// split splits a string based on a delimiter
std::vector<std::string> split (const std::string& s, char delim)
{
        std::vector<std::string> elems;
        std::stringstream ss(s);
        std::string sholder;
        while (std::getline(ss, sholder, delim))
        {
                if (!sholder.empty())
                        elems.push_back(sholder);
        }
        return elems;
}

// draws random, bounded, decimal number
double decimalUnifBound (double min, double max )
{
        return rand() / (static_cast<double>(RAND_MAX) + 1) * (max - min) + min;
}

AssertStyleException::AssertStyleException(const char* err)
        : _error(err)
{}

const char* AssertStyleException::what() const throw()
{
        std::stringstream message;
        message << "Assert type exception occurred:\n";
        if (_error) message << _error;
        return message.str().c_str();
}

PreConditionException::PreConditionException(const char* err)
	: AssertStyleException(err)
{}

const char* PreConditionException::what() const throw()
{
        std::stringstream message;
        message << "Precondition exception occurred:\n";
        if (_error) message << _error;
        return message.str().c_str();
}
