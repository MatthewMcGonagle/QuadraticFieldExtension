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

Rational Rational::Minus(Rational r) {
	int lcm = den*r.GetQ()/gcd(den, r.GetQ());
	int pnew = num*lcm/den - r.GetP()*lcm/r.GetQ();
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
		// We look for a square root of the form a+R^(1/2) b where a,b rational
		Rational newroot(1,1), newrootB(1,1); // Extra variables that will be necessary later
		if(root_[0].IsZero()) {
			newroot = root_[0].Product(Root[0].Inverse());
			//Check is a==0 works.
			if(newroot.FindSqrt()) { 
				Result[0] = Rational(0,1);
				Result[1] = newroot.GetSqrt();
				return true;
			}
			
			// Check b==0 works
			newroot = root_[1];
			if(newroot.FindSqrt()) {
				Result[0] = newroot.GetSqrt();
				Result[1] = Rational(0,1);
				return true;
			}

			// Neither Case Works, so there is no square root
			return false;
		}
		else { 
			// We are the in the non-zero case and is the hardest case
			newroot = root_[0].Product(root_[0]).Minus( Root[0].Product(root_[1].Product(root_[1])));
			if(newroot.FindSqrt()) {
				newroot = newroot.GetSqrt();

				// Check the plus root
				newrootB = root_[0].Sum(newroot);
				newrootB = newrootB.Product(Root[0].Product(Rational(2,1)).Inverse());
				if(newrootB.FindSqrt()) {
					Result[0] = root_[1].Product(Rational(2,1).Product(newrootB.GetSqrt()).Inverse());
					Result[1] = newrootB.GetSqrt();
					return true;
				}

				//Check the minus root
				newrootB = root_[0].Minus(newroot);
				newrootB = newrootB.Product(Root[0].Product(Rational(2,1)).Inverse());
				if(newrootB.FindSqrt()) {
					Result[0] = root_[1].Product(Rational(2,1).Product(newrootB.GetSqrt()).Inverse());
					Result[1] = newrootB.GetSqrt();
					return true;
				}

				// Neither plus/minus root has a root. So there is no rational root
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
