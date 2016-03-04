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
	float ToFloat();

private:

	unsigned int gcd(int p, int q) const;
	unsigned int abs(int p) const;

	int num, den;
	int sqrtp, sqrtq;
	std::string name;
};


//////////////////////////////////

class QuadraticField;

// For quick computation, should only overload *= and += ?
class FieldElement {
public:
	FieldElement();
	FieldElement(QuadraticField* field_, std::vector<Rational> coords_);
	
	std::string Print();
	FieldElement& operator+=(const FieldElement& rhs);
	FieldElement& operator*=(const FieldElement& rhs);

	std::vector<Rational> coords;
	QuadraticField* field;

private:

	FieldElement subfieldelem(int n);
	void multiply(std::vector<Rational>::iterator lhs, std::vector<Rational>::const_iterator rhs);
};

//////////////////////////////////
class QuadraticFieldTower {
	public:
		QuadraticFieldTower(Rational Square_);
		void AddSquare(std::vector<Rational> coords);
		void Add(std::vector<Rational> & lhs, std::vector<Rational> & rhs);
		void Product(std::vector<Rational> & lhs, std::vector<Rational> & rhs);
		std::string Print();
		std::string PrintRootList();
		std::string PrintCoords(std::vector<Rational> & coords);
		std::complex<float> CoordsToComplex(std::vector<Rational> & coords);

	private:
		int degree, numsquares;
		std::vector< std::vector<Rational> > squares;
		std::vector< std::complex<float> > complexroots;
};
//////////////////////////////////
class QuadraticField {

public:
	QuadraticField();
	QuadraticField(Rational root);
	QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_);
	~QuadraticField();

	int GetDegree(){return degree;}
	int GetExtensionIndex(){return extensioni;}
	QuadraticField* GetBaseField() { return basefield;}
	FieldElement GetRoot();
	std::string Print();

	static const int numscratch = 3;	
	std::vector<Rational> Scratch[numscratch];
	std::vector<Rational> Result;

private:
	int degree, extensioni;
	std::vector<Rational> Root;
	QuadraticField* basefield;
	std::string name;

	std::string PrintRootList();
};
