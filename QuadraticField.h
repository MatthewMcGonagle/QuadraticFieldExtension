//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <vector>

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

private:

	unsigned int gcd(int p, int q) const;
	unsigned int abs(int p) const;

	int num, den;
	int sqrtp, sqrtq;
	std::string name;
};


//////////////////////////////////

class QuadraticField;

class FieldElement {
public:
	FieldElement();
	FieldElement(QuadraticField* field_, std::vector<Rational>::iterator coords_);
	
	std::string Print();

	std::vector<Rational>::iterator coords;
	QuadraticField* field;

};

//////////////////////////////////

class QuadraticField {

public:
	QuadraticField();
	QuadraticField(Rational root);
	QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_);
	~QuadraticField();

	int GetDegree(){return degree;}

	static const int numscratch = 3;	
	std::vector<Rational> Scratch[numscratch];
	std::vector<Rational> Result;

private:
	int degree;
	std::vector<Rational> Root;
	QuadraticField* basefield;
	std::string name;
};
