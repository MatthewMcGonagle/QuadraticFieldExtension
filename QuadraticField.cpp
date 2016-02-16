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
	int whichroots, currentroot;
	std::ostringstream ss;
	std::string name = std::string();
	for(int i=0; i < field->GetDegree(); i++) {
		if(i > 0)
			name += " + ";
		name += coords[i].print();
		whichroots = i;	
		currentroot = 1;
		if(i > 0)
			name += " ";
		while(whichroots > 0) {
			if(whichroots % 2 == 1) { 
				ss.str(std::string());	
				ss << currentroot;
				name += "r" + ss.str();
			}
			whichroots /= 2;
			currentroot++;
		}
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
	extensioni = 1;
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
	extensioni = 1+basefield->GetExtensionIndex();
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

FieldElement QuadraticField::GetRoot() {
	return FieldElement(this, Root.begin()); 
}

std::string QuadraticField::Print() {
	std::ostringstream oss;
	std::string name = std::string();

	name += "Q( ";
	int idegree = degree, rootnum = 1;
	while (idegree > 1) {	
		name += "r";
		oss.str(std::string());
		oss << rootnum;
		name += oss.str();
		if (idegree > 2)
			name += ", ";
		idegree /= 2;
		rootnum++;
	}
	name += ") where\n";
	name += PrintRootList();	
	return name;	
}

std::string QuadraticField::PrintRootList() {
	std::ostringstream oss;
	std::string name = std::string();
	if (degree > 2) {
		name += basefield->PrintRootList();
		FieldElement myroot = FieldElement(basefield, Root.begin());
		oss.str(std::string());
		oss << extensioni;
		name += "r";
		name += oss.str();
		name += " = Sqrt of ";
		name += myroot.Print();
		name += "\n";
	}
	else {
		name += "r1 = Sqrt of ";
		name += Root[0].print();
		name += "\n";	
	}

	return name;
}
