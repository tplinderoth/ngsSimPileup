// generalUtils.h
// version 0.0.5

#ifndef GENERALUTILS_H_
#define GENERALUTILS_H_

#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>

// TEMPLATES

// ARRAY TEMPLATE
template <class T>
class Array
{
public:
		T& operator[] (size_t i)
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range\n", i);
					exit(1);
				}
                return data[i];
        }

		T operator[] (size_t i) const
        {
				if (i >= size())
				{
					fprintf(stderr, "ERROR: subscript %lu out of range\n", i);
					exit(1);
				}
                return data[i];
        }

        void setSize(size_t size)
        {
				sz = size;
                data = new T[size];
                for(unsigned int long i = 0; i < size; ++i)
                	data[i] = 0;
        }

         size_t size() const
        {
        	 return sz;
        }

         Array ()
			 : data(0),
			   sz(0)
         { }

         Array ( const Array& oldarr)
			 : sz( oldarr.sz )
         {
        	 data = new T[sz];
        	 for( size_t i = 0; i < sz; ++i)
        		 data[i]= oldarr.data[i];
         }

        ~Array ()
        {
        	delete [] data;
        	data = 0;
        }


private:
        T* data;
        size_t sz;
};

// MATRIX TEMPLATE

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

// EXCEPTION CLASSES

class AssertStyleException : public std::exception
{
public:
        AssertStyleException (const char* err = NULL);
        virtual const char* what () const throw();
protected:
        const char* _error;
};

class PreConditionException : public AssertStyleException
{
public:
	PreConditionException (const char* err = NULL);
	virtual const char* what () const throw();
};

// EXCEPTION MESSAGE FORMATTER CLASS

class ExceptionFormatter
{
public:
	ExceptionFormatter() {}
	~ExceptionFormatter() {}

	template <typename T> ExceptionFormatter& operator<< (const T& value)
	{
		stream_ << value;
		return *this;
	}

	// for converting to std::string

	std::string str() const {return stream_.str();}
	operator std::string () const {return stream_.str();} // cast operator so that anything requiring std::string can take type ExceptionFormatter

	enum ConvertToString {to_str}; // ConvertToString is a variable within the ExceptionFormatter class that can only take 'to_str' as its value
	std::string operator>> (ConvertToString) {return stream_.str();} // overloaded >> so that ExceptionFormater::to_str will return stream_ as a string

	// for converting to C style strings
    const char* cstr() const {return stream_.str().c_str();}
    operator const char* () const {return stream_.str().c_str();}

    enum ConvertToCString {to_cstr};
    const char* operator>> (ConvertToCString) {return stream_.str().c_str();}

private:
	std::stringstream stream_;
	//ExceptionFormatter(const ExceptionFormatter&);
	//ExceptionFormatter& operator= (ExceptionFormatter&);
};

// TEMPLATE FUNCTIONS

template<class T> T vecsum (const std::vector<T>& v)
{
	// sum of all values in a vector
	T sum = 0.0;
	for (unsigned int i = 0; i < v.size(); ++i)
		sum += v[i];

	return sum;
}

// FUNCTION PROTOTYPES

char * getCString (std::string);
void getFILE (std::fstream &, const char*, const char*);
int fexists (const char*);
void readFileChunk (std::fstream &infile, std::vector<std::string>& datavec, unsigned int* chunk, int* end);
std::vector<std::string> split (const std::string&, char);
double decimalUnifBound (double min, double max);


#endif /* GENERALUTILS_H_ */
