//Author: Matthew McGonagle

class Rational {
public:
	Rational(){num = 1; den = 1;};
	Rational(int p, int q);
	int GetP() {return num;};
	int GetQ() {return den;};
private:

	unsigned int gcd(int p, int q);
	unsigned int abs(int p);
	int num, den;
};

class QuadraticField {

public:
	QuadraticField();
	bool FindSqrt(Rational r);
	Rational GetSqrtResult(){return SqrtResult;}
private:
	int degree;
	Rational SqrtResult;
};