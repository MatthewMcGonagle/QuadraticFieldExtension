//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>
#include <exception>

/** 
    class Rational
    This class is for representing fractions that involve precision arithmetic. That is, not floating point
    arithmetic. We leave every operation as a fraction.
*/

class Rational {
public:

    /**
        member function Rational
        Defaults to giving representation of the number 0.
    */
	Rational(){num = 0; den = 1;};

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
	Rational inverse();

    /**
        member function IsZero
        Tests if fraction is 0.
        @return True if the fraction is 0.
    */
	bool isZero();
    
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
    function printCoords
    
    Prints coordinates in nice human readable format.
    @param x Reference to coordinate vector.
    @return A string holding the human readable format.
**/

std::string printCoords(std::vector<Rational> &x);

/** 
    function minDegree
    
    Gets the minimum level of Quadratic Field Tower necessary to give meaning to the coordinates.
    @param x Reference to coordinate vector.
    @return The minimum level.
**/

    unsigned int minDegree(std::vector<Rational> &x);

/**
    function padCoordsToLevel
    
    Pads a coordinate vector with zeroes to make it match a given Quadratic Field 
    Tower level n (i.e. length = 2**n).

    @param x Reference to coordinate vector.
    @param level The level to pad the vector to.
**/

    void padCoordsToLevel(std::vector<Rational> &x);

/**
    function padCoordsToSize

    Pads a coordinate vector with zeroes to make it a certain size.
    
    @param x Reference to coordinate vector.
    @param size The size to pad the vector to.
**/

    void padCoordsToSize(std::vector<Rational> &x, int size);
/**
    struct CoordinateRange

    Holds iterators for a range of coordinates. Makes it easier to pass around subsets of coordinates
    for operations inside class FieldTower without using copying.

*/

struct CoordinateRange {

    CoordinateRange() : size(0) {}
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
        Constructs the coordinates range from a reference to an entire std::vector<Rational>.
        @param x Reference to the coordinates vector.
    **/
    CoordinateRange(std::vector<Rational> &x)
        : begin(x.begin()), end(x.end()), size(x.size()) {}

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

    /**
        member function isZero
        Tests if all Rational numbers in coordinate range are zero.
        @return True if all zero
    */
    bool isZero();

    /**
        member function scale
        Scales each number by scalar. Does so in place.
        @param scalar The value to scale by.
    */
    void scale(Rational scalar);

    /**
        member function set
        Sets each number to prescribed Rational number. Does so in place.
        @param x The Rational to set each number to.
    */
    void set(Rational x);

    /**
        member function copyVals
        Copy the values from the coordinates referenced by another CoordinateRange instance. If lengths
        don't match, then it only copies for smaller of two sizes.
        @param x The CoordinateRange to copy from.
    */
    void copyVals(CoordinateRange x);

    std::vector<Rational>::iterator begin, end;
    int size;

};

/**
    class QuadraticFieldTower
    This class is responsible for holding the chain of squares of roots that define the entire field. It
    also handles the operations on coordinates of elements in the field, such as multiplication and finding
    square roots.
*/

class QuadraticFieldTower {

public:
    /**
        initializer QuadraticFieldTower
        No roots have been added yet. So it starts as only working with rational numbers, that is coordinates
        of length 1.
    */
    QuadraticFieldTower()
        : topCoordLength(1) {}
        
    /**
        member function addIfNoSqrRoot
        Tests if number has square root in current state of field. If not, then this number is added
        to the list of square roots, and we expand the field to handle coordinates that are now twice as
        large as before.
        @param x The coordiantes of the number to check.
    */
    void addIfNoSqrRoot(std::vector<Rational> x); 

    /**
        member function multiply
        Multiplies two numbers together. Uses their coordinates.
        @param x Coordinates of first number.
        @param y Coordinates of second number.
        @return The coordinates of the result of the multiplication.
    */
    std::vector<Rational> multiply(std::vector<Rational> &x, std::vector<Rational> &y);

    /**
        member function add
        Adds two numbers together using their coordinate representation.
        @param x Coordinates of the first number.
        @param y Coordinates of the second number.
        @return The coordinates of the result of the addition.
    */
    std::vector<Rational> add(std::vector<Rational> &x, std::vector<Rational> &y);

    /**
        member function getNLevels
        The returns the number of square roots that are used in the current coordinate system. Note,
        the size of the coordinate system should be 2**getNLevels.
        @return The number of square roots.
    */ 
    int getNLevels() {return squaresOfRoots.size();}

    /**
        member function print
        Prints out the square roots that are used in the current coordinate system. This completely 
        describes the Quadratic Field.
        @return String containing the printing.
    */
    std::string print();

    /** 
        member function convertToComplex
        Finds the complex approximation of a number based on the coordinates representing it. This is
        based on the complex approximation of all of the square roots.
        @return The complex approximation. 
    */
    std::complex<float> convertToComplex(std::vector<Rational> &x);

    /**
        memberfunction hasSqrt
        Finds if the number represented by the coordinates x has a square root in the current number
        system that can be handled by the instance of QuadraticFieldTower. If the square root exists,
        then it is stored in solution.

        @param x Reference to the coordinates representing the number to put under the square root. 
        @param sol Reference to coordinates to hold the square root if it exists.
        @return True if the square root exists in the current number system used by this
            instance of QuadraticFieldTower. 
    **/
    bool hasSqrt(std::vector<Rational> &x, std::vector<Rational> &sol);

private:

    /**
        memberfunction add

        Adds values in two coordinate ranges x and y, and stores result in values in third coordinate range
        result. It is fine to do something like add(x, y, x). Note, guaranteed to work properly only when
        x, y, and result have the same size (as ranges).

        @param x First coordinate range to add.
        @param y Second coordinate range to add.
        @param result Where to store result of adding ranges.
    **/
    void add(CoordinateRange x, CoordinateRange y, CoordinateRange result);
    void multiply(int level, CoordinateRange x, CoordinateRange y, CoordinateRange result);
    void inverse(int level, CoordinateRange x, CoordinateRange sol);
    bool hasSqrt(Rational &x, Rational &sol);
    bool hasSqrt(CoordinateRange x, int level, CoordinateRange sol);

    std::vector<std::vector<Rational> > squaresOfRoots;
    std::vector<std::complex<float> > complexRoots;
    int topCoordLength;

};
