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

double Rational::ToFloat() {
	double result = num;
        result /= den;
	return result;	
}

//////////////////////////////////////////////
/// Functions for struct CoordinateRange
//////////////////////////////////////////////

CoordinateRange CoordinateRange::firstHalf() {

    std::vector<Rational>::iterator newEnd;
    int newSize;

    newSize = size / 2;
    newEnd = begin + newSize;

    return CoordinateRange(begin, newEnd, newSize);
}

CoordinateRange CoordinateRange::secondHalf() {

    std::vector<Rational>::iterator newBegin;
    int newSize;

    newSize = size / 2;
    newBegin = begin + newSize;
    
    return CoordinateRange(newBegin, end, newSize);
}

//////////////////////////////////////////////
/// Functions for class QuadraticFieldTower
//////////////////////////////////////////////

void QuadraticFieldTower::addIfNoSqrRoot(Coordinates x) {

    if (x.values.size() != topCoordLength)
        return;

    if (hasSqrRoot(x)) {
        return;
    }
   
    squaresOfRoots.push_back(x); 
    topCoordLength *= 2;
}

bool QuadraticFieldTower::hasSqrRoot(Coordinates x) {

    return false;
}

void QuadraticFieldTower::add(CoordinateRange x, CoordinateRange y, CoordinateRange result) {

    std::vector<Rational>::iterator xIt, yIt, rIt;
    for(xIt = x.begin, yIt = y.begin, rIt = result.begin;
        xIt != x.end && yIt != y.end && rIt != result.end;
        xIt++, yIt++, rIt++) {
        *rIt = *xIt + *yIt;
    } 
}

void QuadraticFieldTower::multiply(int level, CoordinateRange x, CoordinateRange y, CoordinateRange result) {

    std::vector<Rational> scratch(result.size / 2);
    CoordinateRange x1 = x.firstHalf(),
                    x2 = x.secondHalf(),
                    y1 = y.firstHalf(),
                    y2 = y.secondHalf(),
                    r1 = result.firstHalf(),
                    r2 = result.secondHalf(),
                    square = CoordinateRange(squaresOfRoots[level].values.begin(),
                                             squaresOfRoots[level].values.end(),
                                             squaresOfRoots[level].values.size()),
                    sc(scratch.begin(), scratch.end(), scratch.size());
    int newLevel = level - 1;

    // Multiply the non-cross terms. Need to multiply in the square of the root for this level.

    multiply(newLevel, x2, y2, r1);
    multiply(newLevel, r1, square, r2);
    multiply(newLevel, x1, y1, r1);
    add(r1, r2, r1);

    // Now multiply cross terms. Need to use scratch since result of non-cross multiply is stored in
    // r1.
    
    multiply(newLevel, x1, y2, r2);
    multiply(newLevel, x2, y1, sc);
    add(r2, sc, r2);
      
}
