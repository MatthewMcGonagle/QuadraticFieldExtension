//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <vector>

class Rational {
public:
	Rational(){num = 1; den = 1; name = std::string();};
	Rational(int p, int q);

	int GetP() {return num;};
	int GetQ() {return den;};
	std::string print();
	Rational Product(Rational r);
	Rational Sum(Rational r);
	Rational Minus(Rational r);
	Rational Inverse();

	bool IsZero();
	bool FindSqrt();
	int GetSqrtP(){return sqrtp;}
	int GetSqrtQ(){return sqrtq;}
	Rational GetSqrt(){return Rational(sqrtp, sqrtq);}

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
	QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_);
	~QuadraticField();
	bool FindSqrt(std::vector<Rational> root_);
	void Product(std::vector<Rational> a, std::vector<Rational> b);
	std::string Print(Rational* element);
	std::string Print(std::vector<Rational> coordinates);

	int GetDegree() {return degree;}
	std::vector<Rational> GetSqrt() {return Result;}
private:
	int degree;
	std::vector<Rational> Result, Root;
	QuadraticField* basefield;
	std::string name;
};