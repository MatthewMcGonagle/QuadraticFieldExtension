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

	std::ostringstream s;
    std::string name = std::string();

	s << num;
	if(den == 1)
		return s.str();
	s << "/";
	s << den;

	return s.str();
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

double Rational::toFloat() {
	double result = num;
        result /= den;
	return result;	
}

//////////////////////////////////////////////
/// Functions for Coordinates
//////////////////////////////////////////////

std::string Coordinates::print() {

	std::ostringstream convert;
    std::vector<std::string> names;
    std::vector<std::string>::iterator nameIt;
    std::string finalName, rootName;
    std::vector<Rational>::iterator coordIt;
    int level = -1, rootIncr;
    unsigned int logcompute = values.size();
   
    // First check trivial cases.

    if(values.size() == 0)
        return std::string();

    if(values.size() == 1)
        return values.begin() -> print();
 
    // Compute the log base 2 of size, rounded up to largest integer.
    while (logcompute > 0) {
        level++;
        logcompute = logcompute >> 1;  
    }

    // Get the print out of the Rational values for each coordinate.

    names.reserve(values.size());
    for(coordIt = values.begin(); coordIt != values.end(); coordIt++)
        names.push_back(coordIt -> print());

    // Put in Root Info. Outer iteration is on root, and then inner loop puts it into
    // correct coordinate position.

    rootIncr = 2;
    for(int rootNum = 0; rootNum < level; rootNum++, rootIncr*=2) {

        // Clear the output string stream.

        convert.clear();
        convert.str("");
        convert << "r" << rootNum;
        rootName = convert.str();

        nameIt = names.begin() + rootIncr - 1;
        for(int j = rootIncr - 1; j < names.size(); nameIt += rootIncr, j += rootIncr) 
            nameIt -> append(rootName);
    }
    // Concatenate all of the names for each coordinate into one single name.

    for(nameIt = names.begin(); nameIt != names.end(); nameIt++) {
        if (nameIt != names.begin())
            finalName.append(" + ");
        finalName.append(*nameIt);
    }

    return finalName;
    
    return std::string("Coordinates::Print Needs Implementation");

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

    std::complex<float> newConversion;

    if (x.values.size() != topCoordLength)
        return;

    if (hasSqrRoot(x)) {
        return;
    }
  
    newConversion = convertToComplex(x); 
    newConversion = std::sqrt(newConversion);
    complexRoots.push_back(newConversion);
    squaresOfRoots.push_back(x); 
    topCoordLength *= 2;
}

bool QuadraticFieldTower::hasSqrRoot(Coordinates x) {

    return false;
}

Coordinates QuadraticFieldTower::add(Coordinates &x, Coordinates &y) {

    Coordinates result(topCoordLength);
    CoordinateRange xRange(x), yRange(y), rRange(result);

    if(xRange.size != topCoordLength || yRange.size != topCoordLength)
        return Coordinates(std::vector<Rational>());

    add(xRange, yRange, rRange);
    
    return result;
}

Coordinates QuadraticFieldTower::multiply(Coordinates &x, Coordinates &y) {

    Coordinates result(topCoordLength);
    CoordinateRange xRange(x), yRange(y), rRange(result);

    if (x.values.size() != topCoordLength || y.values.size() != topCoordLength)
        return Coordinates(std::vector<Rational>());

    multiply(squaresOfRoots.size() - 1, xRange, yRange, rRange);

    return result;

}

void QuadraticFieldTower::add(CoordinateRange x, CoordinateRange y, CoordinateRange result) {

    std::vector<Rational>::iterator xIt, yIt, rIt;
    for(xIt = x.begin, yIt = y.begin, rIt = result.begin;
        xIt != x.end && yIt != y.end && rIt != result.end;
        xIt++, yIt++, rIt++) {
        *rIt = *xIt + *yIt;
    } 
}

std::string QuadraticFieldTower::print() {

    std::ostringstream s;

    for(int i = 0; i < squaresOfRoots.size(); i++) {
        s << "r" << i << " = "
          << squaresOfRoots[i].print()
          << std::endl;
    } 

    return s.str();
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

    if (level == -1) {
        *result.begin = *x.begin * (*y.begin); 
        return;
    }
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

std::complex<float> QuadraticFieldTower::convertToComplex(Coordinates &x) {

    std::complex<float> half1, half2;
    std::vector<std::complex<float> > complexCoord;
    std::complex<float> result(0.0, 0.0);

    // If the coordinates are too large, then we can't compute the conversion.

    if (x.values.size() > topCoordLength) 
        return std::complex<float>(0.0, 0.0);

    // If x is just rational, then we can just directly compute the floating point
    // conversion.

    if(x.values.size() == 1) 
        return std::complex<float>(x.values[0].toFloat(), 0.0);
   

    // For more coordinates, we figure out the complex values of each coordinate and then add together.

    // First get complex representation of just Rational part of each coordinate.

    for (int i = 0; i < x.values.size(); i++) 
        complexCoord.push_back(std::complex<float>(x.values[i].toFloat(), 0.0));


    // Now loop over each root square and multiply it into the appropriate coordinate.        

    for(int rootN = 0, coordInc = 2; rootN < complexRoots.size(); rootN++, coordInc *= 2) {

        for(int i = coordInc - 1; i < complexCoord.size(); i += coordInc) 
            complexCoord[i] *= complexRoots[rootN]; 
    }  

    // Now add all of the coordinate conversions together.

    for(int i  = 0; i < complexCoord.size(); i++)
        result += complexCoord[i];

    return result;
}


