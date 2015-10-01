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

QuadraticField::QuadraticField() {
	degree = 1;
	Result = std::vector<Rational>(1);
	basefield= NULL;

	name = std::string();
	Result = Scratch = std::vector<Rational>(1);
}

QuadraticField::QuadraticField(Rational root_) {
	degree = 2;
	Root = std::vector<Rational>(1);
	basefield = NULL;

	Root[0] = root_;
	name = std::string();

	Result = Scratch = std::vector<Rational>(degree);
}

QuadraticField::QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_) {
	basefield = basefield_;
	degree = 2*basefield->GetDegree();
	if(basefield->GetDegree() == root_.size())
		Root = root_;
	else
		Root = std::vector<Rational>(basefield->GetDegree());

	Result = Scratch = std::vector<Rational>(degree);
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

std::string QuadraticField::Print(std::vector<Rational> element){
	if (degree == 2) {
		name = element[0].print();
		name += std::string(" + ");
		name += element[1].print();
		name += std::string(" (");
		name += Root[0].print();
		name += std::string(")^0.5");
		return name;
	}
	else {
		CoordinateChunk x = CoordinateChunk(&element, 0, degree);
		return Print(x);

	}
}

std::string QuadraticField::PrintR(std::vector<Rational> coords) {
	std::string name=std::string();
	std::ostringstream convert;

	for(int i=0; i<degree; i++) {
		if(i != 0)
			name += std::string(" + ");
		name += coords[i].print();
		if(i != 0)
			name += std::string(" ");
		for(int j=1, k=1; j<degree; j*=2, k++) {
			if((i/j)%2==1) {
				convert << k;
				name += std::string("R") += std::string(convert.str());
				convert.str(std::string());
			}
		}
	}
	return name;
}

std::string QuadraticField::Print(CoordinateChunk x){
	if (degree == 2) {
		name = x.Get(0).print();
		name += std::string(" + ");
		name += x.Get(1).print();
		name += std::string(" (");
		name += Root[0].print();
		name += std::string(")^0.5");
		return name;
	}
	else {
		name = basefield->Print(x);
		name += std::string(" + (");
		name += basefield->Print(CoordinateChunk(x.coords, degree/2, degree/2));
		name += std::string(") (");
		name += basefield->Print(Root);
		name += std::string(")^0.5");
		return name;

	}
}

std::string QuadraticField::PrintName() {
	std::string name = std::string("Q(");
	name += PrintRootList();
	name += std::string(")");
	return name;
}

std::string QuadraticField::PrintNameR() {
	std::string name = std::string("Q(");
	name += PrintRootListR();
	name += std::string(")");
	return name;

}

std::string QuadraticField::PrintRootList() {
	std::string name = std::string("");
	if(degree==2) {
		name += std::string("(");
		name += Root[0].print();
		name += std::string(")^0.5");
	}
	else {
		name += std::string("(");
		name += basefield->Print(Root);
		name += std::string(")^0.5, ");
		name += basefield->PrintRootList();
	}
	return name;
}

std::string QuadraticField::PrintRootListR() {
	std::string name = std::string();
	std::ostringstream convert;
	int rootnum=0, power=1;

	if(degree==2) {
		name += std::string("R1 = Sqrt(");
		name += Root[0].print();
		name += std::string(")");
	}
	else {
		//compute the logarithm of degree
		do {
			rootnum++;
			power*=2;
		}while(power < degree);
		convert << rootnum;
		name += std::string("R") += convert.str() += std::string(" = Sqrt(");
		name += basefield->PrintR(Root);
		name += std::string("), ");
		name += basefield->PrintRootListR();
	}
	return name;
}

void QuadraticField::ZeroResult() {
	for(int i =0; i<degree; i++)
		Result[i] = Rational(0,1);
}

void QuadraticField::SetResult(std::vector<Rational>& Result_) {
	for(int i=0; i<degree; i++)
		Result[i] = Result_[i];
}

QuadraticField QuadraticField::operator+ (std::vector<Rational> & other) {
	QuadraticField dummyfield(this, Root);
	for(int i=0; i<degree; i++)
		dummyfield.Result[i] = Result[i]+other[i];
	return dummyfield;
}

void QuadraticField::Sum(std::vector<Rational> &a, std::vector<Rational> &b){
	for(int i=0; i<degree; i++)
		Result[i] = a[i] + b[i];
}

void QuadraticField::Product(std::vector<Rational> &a, std::vector<Rational> &b) {
	if(degree == 2) {
		Result[0] = a[0]*b[0] + Root[0]*a[1]*b[1];
		Result[1] = a[0]*b[1] + a[1]*b[0];
	} else {
		// Make CoordinateChunk 's and then call Product(CoordinateChunk, CoordinateChunk, CoordinateChunk)
		CoordinateChunk x,y,z;
		x = CoordinateChunk(&a, 0, degree);
		y = CoordinateChunk(&b, 0, degree);
		z = CoordinateChunk(&Result, 0, degree);
		Product(x,y,z);
	}
}

void QuadraticField::Product(CoordinateChunk a, CoordinateChunk b, CoordinateChunk Result_) {
	//All CoordinateChunk.size should match degree
	if(degree==2) {
		Result_.Get(0) = a.Get(0)*b.Get(0) + Root[0]*a.Get(1)*b.Get(1);
		Result_.Get(1) = a.Get(1)*b.Get(0) + a.Get(0)*b.Get(1);
	} else {

		// Containers for recursion
		CoordinateChunk x,y,z;

		// First multiply the like terms and the root

		x = CoordinateChunk(a.coords, 0, degree/2);
		y = CoordinateChunk(b.coords, 0, degree/2);
		z = CoordinateChunk(Result_.coords, 0, degree/2);
		basefield->Product(x,y,z);

		x = CoordinateChunk(a.coords, degree/2, degree/2);
		y = CoordinateChunk(b.coords, degree/2, degree/2);
		z = CoordinateChunk(Result_.coords, degree/2, degree/2);
		basefield->Product(x,y,z);

		x = z;
		y = CoordinateChunk(&Root, 0, degree/2);
		z = CoordinateChunk(&Scratch, 0, degree/2);
		basefield->Product(x,y,z);

		for(int i=0; i<degree/2; i++)
			Result_.Get(i) = Result_.Get(i) + Scratch[i];

		// Now the cross terms
		x = CoordinateChunk(a.coords, 0, degree/2);
		y = CoordinateChunk(b.coords, degree/2, degree/2);
		z = CoordinateChunk(Result_.coords, degree/2, degree/2);
		basefield->Product(x,y,z);

		x = CoordinateChunk(a.coords, degree/2, degree/2);
		y = CoordinateChunk(b.coords, 0, degree/2);
		z = CoordinateChunk(&Scratch, 0, degree/2);
		basefield->Product(x,y,z);

		for(int i=0; i<degree/2; i++)
			Result_.Get(degree/2+i) = Result_.Get(degree/2+i) + Scratch[i];
	}
}

////////////////////////////////////////

CoordinateChunk::CoordinateChunk() {
	size = 0;
	offset = 0;
	coords = nullptr;
}

CoordinateChunk::CoordinateChunk(std::vector<Rational>* coords_, int offset_, int size_) {
	size_ = size;
	offset = offset_;
	coords = coords_;
}