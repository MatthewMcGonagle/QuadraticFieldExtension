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

private:

	unsigned int gcd(int p, int q);
	unsigned int abs(int p);
	int num, den;
	std::string name;
};

class QuadraticField {

public:
	QuadraticField();
	QuadraticField(QuadraticField* basefield, int degbase, Rational* root);
	~QuadraticField();
	bool FindSqrt(Rational r);
	bool FindSqrt(Rational* r, int n);
	Rational GetSqrtResult(){return SqrtResult;}
private:
	int degree;
	Rational* root;
	Rational SqrtResult;
	QuadraticField* basefield;
};