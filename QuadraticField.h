//Author: Matthew McGonagle
#include <string>
#include <sstream>

class Rational {
public:
	Rational(){num = 1; den = 1; name = std::string();};
	Rational(int p, int q);

	int GetP() {return num;};
	int GetQ() {return den;};
	std::string print();
	Rational Product(Rational r);
	Rational Sum(Rational r);
	bool FindSqrt();
	int GetSqrtP(){return sqrtp;}
	int GetSqrtQ(){return sqrtq;}

private:

	unsigned int gcd(int p, int q);
	unsigned int abs(int p);

	int num, den;
	int sqrtp, sqrtq;
	std::string name;
};

class QuadraticField {

public:
	QuadraticField();
	QuadraticField(Rational root);
	QuadraticField(QuadraticField* basefield, int degbase, Rational* root);
	~QuadraticField();
	bool FindSqrt(Rational* r, int n);
	void Product(Rational* a, Rational *b);
	std::string Print(Rational* element);
private:
	int degree;
	Rational *root, *current;
	Rational SqrtResult;
	QuadraticField* basefield;
	std::string name;
};