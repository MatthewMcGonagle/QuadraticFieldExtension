// Author: Matthew McGonagle
#include "QuadraticField.h"

Rational::Rational(int p, int q) {
	if (q == 0)
		q = 1;

	if( q <0) {
		p = -p;
		q = -q;
	}

	unsigned g = gcd(p, q);
	if (g>1) {
		p /= g;
		q /=g;
	}
	num = p;
	den = q;

	name = std::string();
}

unsigned int Rational::gcd(int p, int q) {
	unsigned int a = abs(p), b = abs(q), r;
	if(b == 0)
		return 1;
	do {
		r = a % b;
		a = b;
		b = r;
	} while( r >0);
	return a;
}

unsigned int Rational::abs(int p) {
	if (p < 0)
		return (unsigned)(-p);
	else
		return (unsigned)p;
}

std::string Rational::print() {

	std::ostringstream convert;

	convert << num;
	name = convert.str();
	if(den == 1)
		return name;
	name += "/";
	convert.str(std::string());
	convert << den;
	name += convert.str();

	return name;
}

QuadraticField::QuadraticField() {
	degree = 1;
	SqrtResult = Rational(1,1);
	root = NULL;
	basefield= NULL;
}

QuadraticField::QuadraticField(QuadraticField* basefield_, int degbase, Rational* root_){
	degree = 2*degbase;
	root = new Rational[degree];
	for (int i=0; i<degbase; i++)
		root[i] = root_[i];
	basefield = basefield_;
}

QuadraticField::~QuadraticField() {
	delete root;

}

bool QuadraticField::FindSqrt(Rational r) {
	int p = r.GetP(), q = r.GetQ();
	int sqrtp, sqrtq;
	if (p>0 && q>0){
		sqrtp = sqrtq = 0;
		do
			sqrtp++;
		while(sqrtp*sqrtp < p);
		do
			sqrtq++;
		while(sqrtq*sqrtq < q);
		if(sqrtp*sqrtp == p && sqrtq*sqrtq ==q)
		{
			SqrtResult = Rational(sqrtp, sqrtq);
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

bool QuadraticField::FindSqrt(Rational* r, int n) {
	if(n != degree)
		return false;

}

