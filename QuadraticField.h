//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>

/** 
    class Rational
    This class is for representing fractions that involve precision arithmetic. That is, not floating point
    arithmetic. We leave every operation as a fraction.
*/

class Rational {
public:
	Rational(){num = 1; den = 1;};
	Rational(int p, int q);

	int GetP() const {return num;};
	int GetQ() const {return den;};
	std::string print();
	const Rational operator*(const Rational &other) const;
	const Rational operator+(const Rational &other) const;
	const Rational operator-(const Rational &other) const;

	Rational Inverse();

	bool IsZero();
	bool FindSqrt();
	int GetSqrtP(){return sqrtp;}
	int GetSqrtQ(){return sqrtq;}
	Rational GetSqrt(){return Rational(sqrtp, sqrtq);}
	double ToFloat();

private:

	unsigned int gcd(int p, int q) const;
	unsigned int abs(int p) const;

	int num, den;
	int sqrtp, sqrtq;
};

/** 
    struct Coordinates
    
    Simple wrapper for the type of coordinates, which is std::vector<Rational>.
*/

struct Coordinates {

    Coordinates(const std::vector<Rational> &values_)
        : values(values_) {}
    Coordinates(int size)
        : values(size) {}
    
    std::string print();

    std::vector<Rational> values;

};

/**
    struct CoordinateRange

    Holds iterators for a range of coordinates. Makes it easier to pass around subsets of coordinates
    for operations inside class FieldTower without using copying.

*/

struct CoordinateRange {

    CoordinateRange(std::vector<Rational>::iterator begin_, std::vector<Rational>::iterator end_, int size_)
        : begin(begin_), end(end_), size(size_) {}

    CoordinateRange(Coordinates &x)
        : begin(x.values.begin()), end(x.values.end()), size(x.values.size()) {}

    CoordinateRange firstHalf();
    CoordinateRange secondHalf();

    std::vector<Rational>::iterator begin, end;
    int size;

};

/**
    class QuadraticFieldTower

*/

class QuadraticFieldTower {

public:
    QuadraticFieldTower()
        : topCoordLength(1) {}
        
    void addIfNoSqrRoot(Coordinates x); 
    bool hasSqrRoot(Coordinates x);
    Coordinates multiply(Coordinates &x, Coordinates &y);
    Coordinates add(Coordinates &x, Coordinates &y);
    int getNLevels() {return squaresOfRoots.size();}
    std::string print();

private:

    void add(CoordinateRange x, CoordinateRange y, CoordinateRange result);
    void multiply(int level, CoordinateRange x, CoordinateRange y, CoordinateRange result);

    std::vector<Coordinates> squaresOfRoots;
    int topCoordLength;

};
