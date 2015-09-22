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

Rational Rational::Product(Rational r) {
	// First remove common divisors of opposing numerators and denominators! Reduce the chance of overflow error.
	// Multiplication is (a/b) * (p/q)
	int a=num, b=den, p=r.GetP(), q=r.GetQ(), g;
	g = gcd(a,q);
	a/=g;
	q/=g;

	g = gcd(b, p);
	b/=g;
	p/=g;
	return Rational(a*p, b*q);
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

Rational Rational::Sum(Rational r) {
	int lcm = den*r.GetQ()/gcd(den, r.GetQ());
	int pnew = num*(lcm/den) + r.GetP()*(lcm/r.GetQ());
	return Rational(pnew, lcm);
}

Rational Rational::Minus(Rational r) {
	int lcm = den*r.GetQ()/gcd(den, r.GetQ());
	int pnew = num*(lcm/den) - r.GetP()*(lcm/r.GetQ());
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

///////////////////////////////////////////////

QuadraticField::QuadraticField() {
	degree = 1;
	Result = std::vector<Rational>(1);
	basefield= NULL;

	name = std::string();
	Result = std::vector<Rational>(1);
}

QuadraticField::QuadraticField(Rational root_) {
	degree = 2;
	Root = std::vector<Rational>(1);
	basefield = NULL;

	Root[0] = root_;
	name = std::string();

	Result = std::vector<Rational>(degree);
}

QuadraticField::QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_) {
	degree = 2*basefield_->GetDegree();
	if(basefield->GetDegree() == root_.size())
		Root = root_;
	else
		Root = std::vector<Rational>(basefield_->GetDegree());
}


QuadraticField::~QuadraticField() {

}

bool QuadraticField::FindSqrt(std::vector<Rational> root_) {
	if(root_.size() != degree)
		return false;
	else if(degree==2) {
		// We look for a square root of the form a+b*R^(1/2) where a,b rational.
		// For the notes, we are finding the square root of X+Y*R^(1/2)
		// so that X = root_[0] and Y = root_[1].
		Rational newroot(1,1), newrootB(1,1); // Extra variables that will be necessary later
		if(root_[1].IsZero()) {
			// The case of finding square root of X
			// Solve a^2 + R b^2 = X and 2ab = 0

			// The case that a=0
			// Look at b = (X/R)^(1/2)
			newroot = root_[0]*Root[0].Inverse();
			if(newroot.FindSqrt()) { 
				Result[0] = Rational(0,1);
				Result[1] = newroot.GetSqrt();
				return true;
			}
			
			// The case that b=0
			// Look at a = X^(1/2)
			newroot = root_[0];
			if(newroot.FindSqrt()) {
				Result[0] = newroot.GetSqrt();
				Result[1] = Rational(0,1);
				return true;
			}

			// Neither Case Works, so there is no square root
			return false;
		}
		else { 
			// Looking at the hardest case of Y non-zero
			// Solve a^2 + R b^2 = X and 2ab = Y
			// Now know that a and b both are non-zero,
			// So it is safe to substitute b = Y/2a into the 
			// first equation.
			// So we must solve a^2 + (R/4) Y^2/a^2 = X
			// equivalent to a^4 - X a^2 + (R/4) Y^2 = 0.

			// Need to find a in Q satisfying quartic polynomial.
			// Quadratic Formula tells us a^2 is in Q iff
			// (X^2 - R Y^2)^(1/2) is in Q. First check this.

			newroot = root_[0]*root_[0] - Root[0]*root_[1]*root_[1];
			if(newroot.FindSqrt()) {

				// Know that a^2 is in the Q. Need to check square root
				// of either positive or negative root to see if any a
				// is in the field

				newroot = newroot.GetSqrt();

				// Check the plus root
				newrootB = (root_[0]+newroot)*Rational(1,2);
				if(newrootB.FindSqrt()) {
					Result[0] = newrootB.GetSqrt();
					Result[1] = root_[1]*Rational(1,2)*Result[0].Inverse();
					return true;
				}

				//Check the minus root
				newrootB = (root_[0]-newroot)*Rational(1,2);
				if(newrootB.FindSqrt()) {
					Result[0] = newrootB.GetSqrt();
					Result[1] = root_[1]*Rational(1,2)*Result[0].Inverse();
					return true;
				}

				// Neither plus/minus root has a root in Q. So there is no rational root
				return false;

			}
			else
				return false;
		}
	}
	else {
		return false;
	}

}


void QuadraticField::Product(std::vector<Rational> a, std::vector<Rational> b) {
	if(a.size() < degree || b.size() < degree)
		return;
	else {

	}
}

std::string QuadraticField::Print(Rational *element){
	if (degree == 2) {
		name = element[0].print();
		name += std::string(" + ");
		name += element[1].print();
		name += std::string(" (");
		name += Root[0].print();
		name += std::string(")^1/2");
		return name;
	}
	else
		return std::string("Missed it?");
}

std::string QuadraticField::Print(std::vector<Rational> element){
	if (degree == 2) {
		name = element[0].print();
		name += std::string(" + ");
		name += element[1].print();
		name += std::string(" (");
		name += Root[0].print();
		name += std::string(")^1/2");
		return name;
	}
	else
		return std::string("Missed it?");
}
