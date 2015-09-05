// Author: Matthew McGonagle
#include "QuadraticField.h"

Rational::Rational(int p, int q) {
	if (q == 0)
		q = 1;

	if( q <0) {
		p = -p;
		q = -q;
	}

	unsigned int g = gcd(p, q);
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

Rational Rational::Product(Rational r) {
	return Rational(num*r.GetP(), den*r.GetQ());
}

Rational Rational::Sum(Rational r) {
	int lcm = den*r.GetQ()/gcd(den, r.GetQ());
	int pnew = num*lcm/den + r.GetP()*lcm/r.GetQ();
	return Rational(pnew, lcm);
}

bool Rational::FindSqrt() {

	if (num>0 && den>0){
		sqrtp = sqrtq = 0;
		do
			sqrtp++;
		while(sqrtp*sqrtp < num);

		do
			sqrtq++;
		while(sqrtq*sqrtq < den);

		if(sqrtp*sqrtp == num && sqrtq*sqrtq ==den)
			return true;
		else
			return false;
	}
	else
		return false;
}

QuadraticField::QuadraticField() {
	degree = 1;
	SqrtResult = Rational(1,1);
	root = NULL;
	current = NULL;
	basefield= NULL;

	name = std::string();
}

QuadraticField::QuadraticField(Rational root_) {
	degree = 2;
	root  = new Rational[1];
	current = new Rational[2];
	basefield = NULL;

	root[0] = root_;
	name = std::string();
}

QuadraticField::QuadraticField(QuadraticField* basefield_, int degbase, Rational* root_){
	degree = 2*degbase;
	root = new Rational[degree];
	current = new Rational[degree];
	for (int i=0; i<degbase; i++)
		root[i] = root_[i];
	basefield = basefield_;
}


QuadraticField::~QuadraticField() {
	if(root != nullptr)
		delete[] root;
	if(current != nullptr)
		delete[] current;

}

bool QuadraticField::FindSqrt(Rational* r, int n) {
	if(n > degree)
		return false;
	

}

void QuadraticField::Product(Rational *a, Rational *b) {
	if(degree == 2) {
		current[0] = a[1].Product(b[1]);
		current[0] = root[0].Product(current[0]);
		current[0] = current[0].Sum(a[0].Product(b[0]));

		current[1] = a[0].Product(b[1]);
		current[1] = current[1].Sum(a[1].Product(b[0]));
	}
		
}

std::string QuadraticField::Print(Rational *element){
	if (degree == 2) {
		name = element[0].print();
		name += std::string(" + ");
		name += element[1].print();
		name += std::string(" (");
		name += root[0].print();
		name += std::string(")^1/2");
		return name;
	}
	else
		return std::string("Missed it?");
}
