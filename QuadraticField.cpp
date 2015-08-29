// Author: Matthew McGonagle
#include "QuadraticField.h"

Rational::Rational(int p, int q) {
	if (q == 0)
		q = 1;

	unsigned g = gcd(p, q);
	if (g>1) {
		p /= g;
		q /=g;
	}
	num = p;
	den = q;
}

unsigned int Rational::gcd(int p, int q) {
	unsigned int a = p, b = q, r;
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

QuadraticField::QuadraticField() {
	degree = 1;
	SqrtResult = Rational(1,1);
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

