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

Rational Rational::inverse() {
	if(num == 0)
		return Rational(num, den);
	else
		return Rational(den, num);
}

bool Rational::isZero() {
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

std::string printCoords(std::vector<Rational> &x) {

 	std::ostringstream convert;
    std::vector<std::string> names;
    std::vector<std::string>::iterator nameIt;
    std::string finalName, rootName;
    std::vector<Rational>::iterator coordIt;
    int level = -1, rootIncr;
    unsigned int logcompute = x.size();
   
    // First check trivial cases.

    if(x.size() == 0)
        return std::string();

    if(x.size() == 1)
        return x.begin() -> print();
 
    // Compute the log base 2 of size, rounded up to largest integer.
    while (logcompute > 0) {
        level++;
        logcompute = logcompute >> 1;  
    }

    // Get the print out of the Rational values for each coordinate.

    names.reserve(x.size());
    for(coordIt = x.begin(); coordIt != x.end(); coordIt++)
        names.push_back(coordIt -> print());

    // Put in Root Info. Outer iteration is on root, and then inner loop puts it into
    // correct coordinate position.

    rootIncr = 1;
    for(int rootNum = 0; rootNum < level; rootNum++, rootIncr*=2) {

        // Clear the output string stream.

        convert.clear();
        convert.str("");
        convert << "r" << rootNum;
        rootName = convert.str();

        nameIt = names.begin() + rootIncr;
        for(int j = rootIncr; j < names.size(); nameIt += rootIncr, j += rootIncr) 
            for(int k = 0; k < rootIncr; k++, nameIt++, j++)
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

unsigned int minDegree(std::vector<Rational> &x) {

    unsigned int max = 1, level = 0; 

    if (x.size() == 0)
        return 1; 

    while(max < x.size()) {
        max <<= 1;
        level++;
    }

    return level;
}

void padCoordsToLevel(std::vector<Rational> &x, unsigned int level) {

    unsigned int minLevel = minDegree(x);
    int missing;

    if (minLevel > level)
        return;
   
    // Remember that 1 << minLevel == 2**minLevel, a power of 2.
    missing = (1 << minLevel) - x.size();

    if(missing < 1)
        return;

    x.insert(x.end(), missing, Rational(0,1));
}

void padCoordsToSize(std::vector<Rational> &x, int size) {

    int missingSize = size - x.size();

    if(missingSize > 0) 
        x.insert(x.end(), missingSize, Rational(0,1));
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

bool CoordinateRange::isZero() {

    std::vector<Rational>::iterator numIt;
    for(numIt = begin; numIt != end; numIt++) 
        if (! numIt -> isZero())
            return false;

    return true;
}

void CoordinateRange::scale(Rational scalar) {


    std::vector<Rational>::iterator numIt;
    for(numIt = begin; numIt != end; numIt++) 
        *numIt = (*numIt) * scalar;

}

void CoordinateRange::set(Rational x) {


    std::vector<Rational>::iterator numIt;
    for(numIt = begin; numIt != end; numIt++)
        *numIt = x;
}

void CoordinateRange::copyVals(CoordinateRange x) {

    std::vector<Rational>::iterator myIt, xIt;
    for (myIt = begin, xIt = x.begin; myIt != end, xIt != x.end; myIt++, xIt++)
        *myIt = *xIt; 
}

//////////////////////////////////////////////
/// Functions for class QuadraticFieldTower
//////////////////////////////////////////////

void QuadraticFieldTower::addIfNoSqrRoot(std::vector<Rational> x) {

    std::complex<float> newConversion;
    unsigned int minLevel = minDegree(x);
    std::vector<Rational> sqrt(topCoordLength);

    // First handle size comparisons.

    if (x.size() < topCoordLength)
        padCoordsToSize(x, topCoordLength);
    else if (x.size() > topCoordLength)
        return;

    // Now check if there is a square root.

    if (hasSqrt(x, sqrt)) {
        return;
    }

    // x is the right size and has no square root. Therefore, we may add it.  

    newConversion = convertToComplex(x); 
    newConversion = std::sqrt(newConversion);
    complexRoots.push_back(newConversion);
    squaresOfRoots.push_back(x); 
    topCoordLength *= 2;
}

bool QuadraticFieldTower::hasSqrRoot(std::vector<Rational> x) {

    return false;
}

std::vector<Rational> QuadraticFieldTower::add(std::vector<Rational> &x, std::vector<Rational> &y) {

    std::vector<Rational> result;
    std::vector<Rational>::iterator addBegin, addEnd,
                                    resultIt, addIt;
    
    if (x.size() > y.size()) {
        result = x;
        addBegin = y.begin();
        addEnd = y.end();
    }
    else {
        result = y;
        addBegin = x.begin();
        addEnd = x.end();
    }

    for ( resultIt = result.begin(), addIt = addBegin; 
          resultIt != result.end() && addIt != addEnd;
          resultIt++, addIt++ ) {

        *resultIt = *resultIt + *addIt;
    } 
   
    return result;
}

std::vector<Rational> QuadraticFieldTower::multiply(std::vector<Rational> &x, std::vector<Rational> &y) {

    std::vector<Rational> result(topCoordLength);
    CoordinateRange xRange(x), yRange(y), rRange(result);

    // Check lengths of coordinates to multiply.

    if (x.size() > topCoordLength || y.size() > topCoordLength)
        throw std::invalid_argument("Coordinate size is too large.");

    if (x.size() < topCoordLength) {
        padCoordsToSize(x, topCoordLength);
        xRange = CoordinateRange(x);
    }

    if (y.size() < topCoordLength) {
        padCoordsToSize(y, topCoordLength);
    }

    // Now multiply. 

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
          << printCoords(squaresOfRoots[i])
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
                    square = CoordinateRange(squaresOfRoots[level].begin(),
                                             squaresOfRoots[level].end(),
                                             squaresOfRoots[level].size()),
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

std::complex<float> QuadraticFieldTower::convertToComplex(std::vector<Rational> &x) {

    std::complex<float> half1, half2;
    std::vector<std::complex<float> > complexCoord;
    std::complex<float> result(0.0, 0.0);

    // If the coordinates are too large, then we can't compute the conversion.

    if (x.size() > topCoordLength) 
        return std::complex<float>(0.0, 0.0);

    // If x is just rational, then we can just directly compute the floating point
    // conversion.

    if(x.size() == 1) 
        return std::complex<float>(x[0].toFloat(), 0.0);
   

    // For more coordinates, we figure out the complex values of each coordinate and then add together.

    // First get complex representation of just Rational part of each coordinate.

    for (int i = 0; i < x.size(); i++) 
        complexCoord.push_back(std::complex<float>(x[i].toFloat(), 0.0));


    // Now loop over each root square and multiply it into the appropriate coordinate.        

    for(int rootN = 0, coordInc = 1; rootN < complexRoots.size(); rootN++, coordInc *= 2) {

        for(int i = coordInc; i < complexCoord.size(); i += 2 * coordInc) 
            for(int k = i; k < i + coordInc && k < complexCoord.size(); k++)
                complexCoord[k] *= complexRoots[rootN]; 
    }  

    // Now add all of the coordinate conversions together.

    for(int i  = 0; i < complexCoord.size(); i++)
        result += complexCoord[i];

    return result;
}

bool QuadraticFieldTower::hasSqrt(std::vector<Rational> &x, std::vector<Rational> &sol) {

    CoordinateRange xRange(x), solRange(sol);

    // First deal with sizes of x and solution.

    if (x.size() > topCoordLength)
        return false;
    else if (x.size() < topCoordLength) {
        padCoordsToSize(x, topCoordLength);
        xRange = CoordinateRange(x);
    }

    if (sol.size() > topCoordLength)
        solRange = CoordinateRange(sol.begin(), sol.begin() + topCoordLength, topCoordLength);

    else if(sol.size() < topCoordLength) {
        padCoordsToSize(sol, topCoordLength);
        solRange = CoordinateRange(sol); 
    }

    // Now use internal function to check for square root using the coordinate ranges.

    return hasSqrt(xRange, squaresOfRoots.size() - 1, solRange);
}

void QuadraticFieldTower::inverse(int level, CoordinateRange x, CoordinateRange sol) {

    std::vector<Rational> scratch1(x.size / 2), scratch2(x.size / 2);
    CoordinateRange scratchR1(scratch1), scratchR2(scratch2),
    // Form of x = a + r b.
                    a = x.firstHalf(), b = x.secondHalf(), 
                    solR1 = sol.firstHalf(), solR2 = sol.secondHalf(),
                    rSquare = CoordinateRange(squaresOfRoots[level]);

    if(level < 0) {
        *(sol.begin) = x.begin -> inverse(); 
        return;
    }

    // Store inverse of a**2 - r**2 * b**2 into scratch. Note, that at the end,
    // second half of scratch should be all 0.

    multiply(level - 1, rSquare, b, scratchR1); 
    multiply(level - 1, scratchR1, b, scratchR2);
    scratchR2.scale(Rational(-1, 1));
    multiply(level - 1, a, a, scratchR1);
    add(scratchR1, scratchR2, scratchR2);
    inverse(level - 1, scratchR2, scratchR1);
   
    // Now multiply scratch against a - r * b 
    multiply(level - 1, scratchR1, a, solR1); 
    multiply(level - 1, scratchR1, b, solR2); 
    solR2.scale(Rational(-1,1));
}

bool QuadraticFieldTower::hasSqrt(Rational &x, Rational &sol) {

    int num = x.GetP(), den = x.GetQ(), sqrtp = 0, sqrtq = 0;

	if (num > 0 && den > 0){

		do
			sqrtp++;

		while ( sqrtp*sqrtp < num );

		do
			sqrtq++;

		while ( sqrtq*sqrtq < den );

		if( sqrtp*sqrtp == num && sqrtq*sqrtq ==den ) {
            sol = Rational(sqrtp, sqrtq);
			return true;
        }
		else
			return false;
	}
	else
		return false;
}

bool QuadraticFieldTower::hasSqrt(CoordinateRange x, int level, CoordinateRange sol) {

    // Use x = a + r * b
    std::vector<Rational> scratch1(x.size/2), scratch2(x.size/2), scratch3(x.size / 2);
    CoordinateRange a = x.firstHalf(), b = x.secondHalf(),
                    scratch1R(scratch1), scratch2R(scratch2), scratch3R(scratch3),
                    solR1 = sol.firstHalf(), solR2 = sol.secondHalf(),
                    r; 
    bool result;
  
    if (level < 0) {

       result = hasSqrt(*x.begin, *sol.begin); 
       return result;

    } 

    r = CoordinateRange(squaresOfRoots[level]);

    if (b.isZero()) {
        result = hasSqrt(a, level - 1, sol.firstHalf());
        return result;
    }     

    else {

        // Need to test sqrt of a**2 - r**2 * b**2

        multiply(level - 1, b, b, scratch1R); 
        multiply(level - 1, scratch1R, r, scratch2R);
        scratch2R.scale(Rational(-1, 1));
        multiply(level - 1, a, a, scratch1R);
        add(scratch1R, scratch2R, scratch3R);
        
        // Now test square root.

        result = hasSqrt(scratch3R, level - 1, scratch1R);

        if (!result)
            return result;

        // Now test c = sqrt[(a + sqrt(a**2 - r**2 * b**2)) / 2 ].

        add(scratch1R, a, scratch2R);
        scratch2R.scale(Rational(1,2)); 
        result = hasSqrt(scratch2R, level - 1, scratch3R);
    
        // When square root found, get other half of square root and store in solution.
        if (result) {

            // Set c is scratch3R.
            solR1.copyVals(scratch3R);

            // Second half root, d, is b / (2c). 
            scratch3R.scale(Rational(2,1));
            inverse(level - 1, scratch3R, scratch2R);
            multiply(level - 1, b, scratch2R, solR2);
             
            return result;
        }

        // Now test c = sqrt[(a - sqrt(a**2 - r**2 * b**2)) / 2 ].

        scratch1R.scale(Rational(-1, 1));
        add(scratch1R, a, scratch2R);
        scratch2R.scale(Rational(1,2));
        result = hasSqrt(scratch2R, level - 1, scratch3R);

        // When square root found, get the other half of the square root and store in solution.
        if (result) {

            // Set c is scratch3R.
            solR1.copyVals(scratch3R);

            // Second half root, d, is b / (2c). 
            scratch3R.scale(Rational(2,1));
            inverse(level - 1, scratch3R, scratch2R);
            multiply(level - 1, b, scratch2R, solR2);

            return result;
        }
        
    }

    // Didn't return when checking square roots in last steps, so it must not exist.
    return false;
} 

