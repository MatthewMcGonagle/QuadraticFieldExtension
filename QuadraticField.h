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

    /**
        member function Rational
        Defaults to giving representation of the number 1.
    */
	Rational(){num = 1; den = 1;};

    /**
        member function Ratoinal
        Initializes the number as the fraction p/q.
        @param p Numerator of number.
        @param q Denominator of number.
    */
	Rational(int p, int q);

    /**
        member function GetP
        Gets numerator.
        @return The numerator of the fraction.
    */
	int GetP() const {return num;};

    /**
        member function GetQ
        Return denominator.
        @return The denominator of the fraction.
    */
	int GetQ() const {return den;};

    /**
        member function print
        Returns a string representation of the fraction. When denominator is 1, this is simply the 
            represetation of the numerator. So, the representation is as simple as possble.
        @return String representation of the fraction.
    */
	std::string print();

    /**
        member function operator*
        Allows multiplication of fractions. Note, this will try to remove common factors before
        multiplying inorder to reduce the chance of integer overflow.
        @param other A reference to the other fraction in the multiplication.
        @return The result of the multiplication.
    */
	const Rational operator*(const Rational &other) const;

    /**
        member function operator+
        Allows addition of fractions. Note, this will try to remove common factors before
        performing the addition. This will reduce the chance of integer overflow error.
        @param other A reference to the other fraction in the addition.
        @return The result of the addition.
    */
	const Rational operator+(const Rational &other) const;

    /** 
        member function opertor-
        Allows the subtraction of fractions. Note, this will try to remove common factors before
        subtracting. This will reduce the chance of integer overflow error.
        @param other A reference to the fraction being subtracted.
        @return The result of the subtraction.
    */
	const Rational operator-(const Rational &other) const;

    /**
        member function Inverse
        Gets the inverse of the fraction.
        @return Returns the fraction representing q/p when this fraction is p/q.
    */
	Rational Inverse();

    /**
        member function IsZero
        Tests if fraction is 0.
        @return True if the fraction is 0.
    */
	bool IsZero();
    
    /**
        member function FindSqrt
        Finds the square root of the fraction IF the square root is also a fraction. Use this with the 
        member function GetSqrt.  
        @return True if the square root of the fraction is also a fraction.
    */
	bool FindSqrt();

    /**
        member function GetSqrtP
        Gets the numerator of the square root when it exists.
        @return The numerator of the square root.
    */
	int GetSqrtP(){return sqrtp;}

    /**
        member function GetSqrtQ
        Gets the denominator of the square root when it exists.
        @return The denominator of the square root.
    */
	int GetSqrtQ(){return sqrtq;}

    /**
        member function GetSqrt
        After testing if the square root exists as a rational number, this returns the square root
        if it is also a rational number. Should be used with the member function FindSqrt.
        @return The fraction that is the square root of this number (if it exists).
    */
	Rational GetSqrt(){return Rational(sqrtp, sqrtq);}

    /**
        member function toFloat
        Gets the floating point approximation of the rational number. This is mostly for testing purposes.
        @return The floating point number representing (approximating) the fraction.
    */
	double toFloat();

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
    std::complex<float> convertToComplex(Coordinates &x);

private:

    void add(CoordinateRange x, CoordinateRange y, CoordinateRange result);
    void multiply(int level, CoordinateRange x, CoordinateRange y, CoordinateRange result);

    std::vector<Coordinates> squaresOfRoots;
    std::vector<std::complex<float> > complexRoots;
    int topCoordLength;

};
