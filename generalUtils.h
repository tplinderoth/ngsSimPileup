// generalUtils.h
// version 0.0.5

#ifndef GENERALUTILS_H_
#define GENERALUTILS_H_

#include <string>
#include <fstream>
#include <vector>

// TEMPLATES

template <class T>
class matrix
{
public:
	matrix (int rown, int coln); // constructor
	~matrix(); // destructor
	T max (); // returns maximum value in matrix
	T min (); // returns minimum value in matrix
	int ncol;
	int nrow;
	T** data;
private:
	int setCol(int); // sets number columns
	int setRow(int); // sets number rows
	void allocate (int row, int col); // reserves memory for matrix
};

// MATRIX CLASS TEMPLATE FUNCTIONS

template<class T> matrix<T>::matrix (int rown, int coln)
{
	int sx = setCol(coln);
	int sy = setRow(rown);
	if (sx != 0 || sy != 0)
	{
		fprintf(stderr, "Could not allocate memory for matrix object -> exiting\n");
		exit(EXIT_FAILURE);
	}
	else
		allocate(rown, coln);
}

template<class T> matrix<T>::~matrix ()
{
	for (int i = 0; i < nrow; i++)
		delete [] data[i];
	delete [] data;
}

template<class T> int matrix<T>::setCol(int c)
{
	if (c <= 0)
	{
		fprintf(stderr, "Number of rows not a positive integer in matrix object\n");
		return (1);
	}
	else
		ncol = c;

	return 0;
}

template<class T> int matrix<T>::setRow(int r)
{
	if (r <= 0)
	{
		fprintf(stderr, "Number of columns not a positive integer in matrix object\n");
		return (1);
	}
	else
		nrow = r;

	return 0;
}

template<class T> void matrix<T>::allocate (int row, int col)
{
	data = new T* [row];
	int i = 0;
	int j = 0;

	for (i = 0; i < row; i++)
		data[i] = new T [col];

	// initialize values
	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
			data[i][j] = 0;
}

template<class T> T matrix<T>::max ()
{
	T max = data[0][0];

	for (int i = 0; i < nrow; i++)
		for (int j = 0; j < ncol; j++)
		{
			if (data[i][j] > max)
				max = data[i][j];
		}

	return max;
}

template<class T> T matrix<T>::min ()
{
	T min = data[0][0];

	for (int i = 0; i < nrow; i++)
		for (int j = 0; j < ncol; j++)
		{
			if (data[i][j] < min)
				min = data[i][j];
		}

	return min;
}

// END MATRIX CLASS TEMPLATE FUNCTIONS

// FUNCTION PROTOTYPES

char * getCString (std::string);
void getFILE (std::fstream &, const char*, const char*);
int fexists (const char*);
void readFileChunk (std::fstream &infile, std::vector<std::string>& datavec, unsigned int* chunk, int* end);
std::vector<std::string> split (const std::string&, char);
double decimalUnifBound (double min, double max );

#endif /* GENERALUTILS_H_ */
