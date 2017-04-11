//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>

class Rational {
public:
	Rational(){num = 1; den = 1; name = std::string();};
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
	std::string name;
};

/////////////////////////////////

class CoordinatesRange {
    public:
    CoordinatesRange(){length = 0;};
    CoordinatesRange( std::vector<Rational>::iterator beginIt_
                    , std::vector<Rational>::iterator endIt_
                    , int length_);

    CoordinatesRange& operator += (const CoordinatesRange& rhs);
    std::vector<Rational> operator + (const CoordinatesRange& rhs) const;
    CoordinatesRange& operator -= (const CoordinatesRange& rhs);
    std::vector<Rational> operator - (const CoordinatesRange& rhs) const;


    std::vector<Rational>::iterator beginIt, endIt;
    int length;
};

////////////////////////////////

class Coordinates {

    public:
        Coordinates() {coordinates = std::vector<Rational>(0);};
        Coordinates(std::vector<Rational> coord){ coordinates = coord;};

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

        std::vector<Rational> coordinates;        

}; 

//////////////////////////////////

// For quick computation, should only overload *= and += ?

class QuadraticFieldTower;

class FieldElement {

public:
	FieldElement();
	FieldElement(QuadraticFieldTower* field, int level, std::vector<Rational> coords_);
	
	std::string Print();
	FieldElement& operator+=(const FieldElement& rhs);
    FieldElement& operator-=(const FieldElement& rhs);
	FieldElement& operator*=(const FieldElement& rhs);

    FieldElement operator+(const FieldElement& rhs);
    FieldElement operator-(const FieldElement& rhs);
    FieldElement operator*(const FieldElement& rhs);

	std::vector<Rational> coords;
    QuadraticFieldTower* field;
    int level;

private:

};

//////////////////////////////////
class QuadraticFieldTower {
	public:
		QuadraticFieldTower(Rational Square_);
		void AddSquare(std::vector<Rational> coords);
		std::string Print();
		std::string PrintRootList();
		std::complex<double> toComplex(std::vector<Rational> & coords);
        std::complex<double> toComplex(FieldElement &x);

        bool sqrtExists(std::vector<Rational> &coords);
        void multiplyInPlace( const std::vector<Rational>& lhs, const std::vector<Rational>& rhs
                            , int lhsLevel, int rhsLevel, std::vector<Rational>& solution);
        std::vector<Rational> multiply( const std::vector<Rational>& lhs, const std::vector<Rational>& rhs
                                      , int lhsLevel, int rhsLevel);

	private:
        bool sqrtExistsResult;
		int degree, numsquares;
        std::vector<Rational> sqrtResult;
		std::vector< std::vector<Rational> > squares;
		std::vector< std::complex<double> > complexroots;

        void multiplyLhsLargest(std::vector<Rational>::const_iterator lhs, std::vector<Rational>::const_iterator rhs, std::vector<Rational>::iterator solutionIt, int lhsLength, int rhsLength, int lhsLevel);
        std::vector<Rational> getSqrt(std::vector<Rational> coords, int level);

};

