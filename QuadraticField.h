//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>

class Rational {
public:
	Rational(){num = 1; den = 1; }; 
	Rational(int p, int q);

	int GetP() const {return num;};
	int GetQ() const {return den;};
	std::string print();
    Rational& operator *= (const Rational& rhs);
	Rational operator*(const Rational &other) const;
    Rational& operator+=(const Rational &rhs);
	Rational operator+(const Rational &other) const;
    Rational& operator-=(const Rational &rhs);
	Rational operator-(const Rational &other) const;

	Rational Inverse();

	bool IsZero() const;
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

/////////////////////////////////

class Coordinates;

class CoordinatesRange {
    public:
    CoordinatesRange(){length = -1;};
    CoordinatesRange( std::vector<Rational>::iterator beginIt_, int length_);
    CoordinatesRange( std::vector<Rational> coord);
    CoordinatesRange( Coordinates & coord);

    CoordinatesRange& operator += (const CoordinatesRange& rhs);
    CoordinatesRange& operator -= (const CoordinatesRange& rhs);
    CoordinatesRange& operator *= (const Rational& rhs);
    std::vector<Rational>::iterator begin() {return beginIt;};
    std::vector<Rational>::iterator middle() {return beginIt + length / 2;};
    std::pair<CoordinatesRange, CoordinatesRange> splitToPair();
    int size() { return length;};

    std::vector<Rational>::iterator beginIt;
    int length;
};

////////////////////////////////

class Coordinates {

    public:
        Coordinates() {coordinates = std::vector<Rational>(0); level = -1;};
        Coordinates(std::vector<Rational> coord, int level_){ coordinates = coord; level = level_;};

        Coordinates& operator+=(const Coordinates& rhs);
        Coordinates operator+(const Coordinates& rhs) const;
        Coordinates& operator-=(const Coordinates& rhs);
        Coordinates operator-(const Coordinates& rhs) const;
        Coordinates& operator*=(const Rational& scaling);
        Coordinates operator*(const Rational& scaling) const;

        std::vector<Rational>::iterator begin() {return coordinates.begin(); };
        std::vector<Rational>::const_iterator begin() const {return coordinates.begin();};
        std::vector<Rational>::iterator end() {return coordinates.end(); }
        std::vector<Rational>::const_iterator end() const {return coordinates.end();};
        std::vector<Rational>::iterator middle() {return coordinates.begin() + size() / 2;};
        std::vector<Rational>::const_iterator middle() const {coordinates.begin() + size() / 2;};

        int size() const {return coordinates.size();};
        bool isZero() const;
        std::string print() ;

        std::vector<Rational> coordinates;        
        int level;

}; 

//////////////////////////////////

// For quick computation, should only overload *= and += ?

class QuadraticFieldTower;

class FieldElement {

public:
	FieldElement();
	FieldElement(QuadraticFieldTower* field, Coordinates &coordinates_); 
	
	std::string Print();
	FieldElement& operator+=(const FieldElement& rhs);
    FieldElement& operator-=(const FieldElement& rhs);
	// FieldElement& operator*=(const FieldElement& rhs);

    FieldElement operator+(const FieldElement& rhs);
    FieldElement operator-(const FieldElement& rhs);
    // FieldElement operator*(const FieldElement& rhs);

    Coordinates coordinates;
    QuadraticFieldTower* field;

private:

};

//////////////////////////////////
class QuadraticFieldTower {
	public:
		QuadraticFieldTower(Rational Square_);
		std::string Print();
		std::string PrintRootList();
		std::complex<double> toComplex(Coordinates& x);
        std::complex<double> toComplex(FieldElement &x);

        bool sqrtExists(std::vector<Rational> &coords);
        // void multiplyInPlace( const std::vector<Rational>& lhs, const std::vector<Rational>& rhs
        //                     , int lhsLevel, int rhsLevel, std::vector<Rational>& solution);
        void multiplyInPlace( const Coordinates& lhs, const Coordinates& rhs, Coordinates& solution); 
        // std::vector<Rational> multiply( const std::vector<Rational>& lhs, const std::vector<Rational>& rhs
        //                               , int lhsLevel, int rhsLevel);
        Coordinates multiply( const Coordinates& lhs, const Coordinates& rhs);
        Coordinates inverse(Coordinates& x, int level);
        void AddSquare(Coordinates &square);

        int getNumSquares() {return numsquares;}
        int getDegree() {return degree;}
        int getLevel() {return numsquares - 1;}

	private:
        bool sqrtExistsResult, inverseExists;
		int degree, numsquares;
        std::vector<Rational> sqrtResult;
		std::vector< std::vector<Rational> > squares;
        std::vector<Coordinates> rootsSquared;
		std::vector< std::complex<double> > complexroots;

        // void multiplyLhsLargest( std::vector<Rational>::const_iterator lhs
        //                        , std::vector<Rational>::const_iterator rhs
        //                        , std::vector<Rational>::iterator solutionIt
        //                        , const int lhsLength, const int rhsLength, const int lhsLevel) const;
        void multiplyLhsLargest( CoordinatesRange lhs, CoordinatesRange rhs, CoordinatesRange solution, int level);
        Coordinates getSqrt(Coordinates& x, int level);
};

