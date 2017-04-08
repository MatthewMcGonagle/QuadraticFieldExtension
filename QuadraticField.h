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
	std::string name;
};


//////////////////////////////////

// class QuadraticField;

// For quick computation, should only overload *= and += ?

class QuadraticFieldTower;

class FieldElement {
public:
	FieldElement();
	FieldElement(QuadraticFieldTower* field, int level, std::vector<Rational> coords_);
	
	std::string Print();
	FieldElement& operator+=(const FieldElement& rhs);
	FieldElement& operator*=(const FieldElement& rhs);

    FieldElement operator+(const FieldElement& rhs);
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
		std::complex<double> CoordsToComplex(std::vector<Rational> & coords);

        FieldElement multiply(const std::vector<Rational>& lhs, const std::vector<Rational>& rhs);
        void multiply(std::vector<Rational>::const_iterator lhs, std::vector<Rational>::const_iterator rhs, std::vector<Rational>::iterator solutionIt, int length, int level);

	private:
		int degree, numsquares;
		std::vector< std::vector<Rational> > squares;
		std::vector< std::complex<double> > complexroots;
};

