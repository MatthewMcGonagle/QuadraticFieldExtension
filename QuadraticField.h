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

    /** 
        member function gcd
        Gets the greatest common divisor of two integers. Useful for reducing by common factors to keep
        integers from overflowing.
        @param p First integer.
        @param q Second integer.
        @return The greatest common divisor.
    */
	unsigned int gcd(int p, int q) const;

    /**
        member function abs
        Returns the absolute values of an integer.
        @param p The integer.
        @return The absolute value of p.
    */ 
	unsigned int abs(int p) const;

	int num, den;
	int sqrtp, sqrtq;
};

/** 
    struct Coordinates
    
    Simple wrapper for the type of coordinates, which is std::vector<Rational>.
*/

struct Coordinates {

    /**
        initializer Coordinates
        Initializes the values based on a pre-assigned vector<Rational> giving the coordinate values.   
    */
    Coordinates(const std::vector<Rational> &values_)
        : values(values_) {}

    /** 
        initializer Coordinates
        Initializes the values based on their total length. Uses the defualt values for the Rational class.
    */
    Coordinates(int size)
        : values(size) {}

    /**
        member function print
        Gives the string reperentation of the coordinates. This is meant to be human readable.
        @return The string representation. 
    */
    std::string print();

    std::vector<Rational> values;

};

/**
    struct CoordinateRange

    Holds iterators for a range of coordinates. Makes it easier to pass around subsets of coordinates
    for operations inside class FieldTower without using copying.

*/

struct CoordinateRange {

    /** 
        initializer CoordinateRange
        Constructs the coordinates range using iterators for the beginning of the range and the end of the range.
        We also include the size so we don't need to waste time calculating it from the iterators.
        @param begin_ An iterator pointing to the beginning of this coordinate range.
        @param end_ An iterator pointing to the end of this coordinate range.
        @param size The size of the coordinate range. It should match the number of elements from begin_ to end_.
            The caller is responsible for making sure this is true.
    */
    CoordinateRange(std::vector<Rational>::iterator begin_, std::vector<Rational>::iterator end_, int size_)
        : begin(begin_), end(end_), size(size_) {}

    /**
        initializer CoordinateRange
        Constructs the coordinates range from a reference to a complete set of coordinates. So the coordinate
        range points to the entirety of the coordinates.
        @param x Reference to coordinates range to point to.
    */
    CoordinateRange(Coordinates &x)
        : begin(x.values.begin()), end(x.values.end()), size(x.values.size()) {}

    /**
        member function firstHalf
        Constructs a new coordinate range pointing to the first half of this coordinate range.
        @return The coordinate range for the first half of this coordinate range.
    */
    CoordinateRange firstHalf();

    /**
        member function secondHalf
        Constructs a new coordinate range pointing to the second half of this coordinate range.
        @return The coordinate range for the second half of this coorindate range.
    */
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
