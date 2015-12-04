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

unsigned int Rational::gcd(int p, int q) const {
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

unsigned int Rational::abs(int p) const {
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

const Rational Rational::operator*(const Rational &r) const {
	// First remove common divisors of opposing numerators and denominators! Reduce the chance of overflow error.
	// Multiplication is (a/b) * (p/q)
	
	int a=this->num, b=this->den, p=r.GetP(), q=r.GetQ(), g;
	g = gcd(a,q);
	a/=g;
	q/=g;

	g = gcd(b, p);
	b/=g;
	p/=g;
	return Rational(a*p, b*q);
}

const Rational Rational::operator+(const Rational &r) const {
	// Do division first to reduce chance of integer overflow errors.
	int lcm, pnew;
	lcm = r.GetQ()/gcd(this->den, r.GetQ());
	lcm *= this->den;
	pnew = this->num*(lcm/this->den);
	pnew += r.GetP()*(lcm/r.GetQ());
	return Rational(pnew, lcm);
}

const Rational Rational::operator-(const Rational &r) const {
	// Do division first to reduce chance of integer overflow errors.
	int lcm, pnew;
	lcm = r.GetQ()/gcd(this->den, r.GetQ());
	lcm *= this->den;
	pnew = this->num*(lcm/this->den);
	pnew -= r.GetP()*(lcm/r.GetQ());
	return Rational(pnew, lcm);
}

Rational Rational::Inverse() {
	if(num == 0)
		return Rational(num, den);
	else
		return Rational(den, num);
}

bool Rational::IsZero() {
	if(num== 0)
		return true;
	else
		return false;
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


//////////////////////////////////////////////
FieldElement::FieldElement() {

}

FieldElement::FieldElement(QuadraticField* field_, std::vector<Rational>::iterator coords_) {
	field = field_;
	coords = coords_;
}

std::string FieldElement::Print() {
	std::string name = std::string();
	for(int i=0; i < field->GetDegree(); i++) {
		name += coords[i].print();
		name += " ";
	}
	return name;
}

//////////////////////////////////////////////

QuadraticField::QuadraticField() {
	degree = 1;
	Result = std::vector<Rational>(1);
	basefield= NULL;

	name = std::string();
	Result = std::vector<Rational>(1);

	for(int i=0; i<numscratch; i++) 
		Scratch[i] = std::vector<Rational>(1);	
}

QuadraticField::QuadraticField(Rational root_) {
	degree = 2;
	Root = std::vector<Rational>(1);
	basefield = NULL;

	Root[0] = root_;
	name = std::string();

	Result = std::vector<Rational>(degree);
	for(int i=0; i<numscratch; i++)
		Scratch[i] = std::vector<Rational>(degree);
}

QuadraticField::QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_) {
	basefield = basefield_;
	degree = 2*basefield->GetDegree();
	if(basefield->GetDegree() == root_.size())
		Root = root_;
	else
		Root = std::vector<Rational>(basefield->GetDegree());

	Result = std::vector<Rational>(degree);
	for(int i=0; i<numscratch; i++)
		Scratch[i] = std::vector<Rational>(degree);
}


QuadraticField::~QuadraticField() {

}

