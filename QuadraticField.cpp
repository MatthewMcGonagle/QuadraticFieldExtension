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

Rational& Rational::operator *= (const Rational& rhs) {
	// First remove common divisors of opposing numerators and denominators! Reduce the chance of overflow error.
	// Multiplication is (a/b) * (p/q)
	
	int a=this->num, b=this->den, p=rhs.GetP(), q=rhs.GetQ(), g;
	g = gcd(a,q);
	a/=g;
	q/=g;

	g = gcd(b, p);
	b/=g;
	p/=g;

    num = a * p;
    den = b * q;

    return *this;
}

Rational Rational::operator*(const Rational &rhs) const {
	// First remove common divisors of opposing numerators and denominators! Reduce the chance of overflow error.
	// Multiplication is (a/b) * (p/q)
	
// 	int a=this->num, b=this->den, p=r.GetP(), q=r.GetQ(), g;
// 	g = gcd(a,q);
// 	a/=g;
// 	q/=g;
// 
// 	g = gcd(b, p);
// 	b/=g;
// 	p/=g;
// 	return Rational(a*p, b*q);
    Rational result = *this;
    return result *= rhs; 
}

Rational& Rational::operator += (const Rational &rhs) {
    // Do division first to reduce chance of integer overflow errors.
	int lcm, pnew;
	lcm = rhs.GetQ()/gcd(this->den, rhs.GetQ());
	lcm *= this->den;
	pnew = this->num*(lcm/this->den);
	pnew += rhs.GetP()*(lcm/rhs.GetQ());

    num = pnew;
    den = lcm;
    return *this;
}

Rational Rational::operator+(const Rational &rhs) const {
    Rational result = *this;
    return result += rhs;
}

Rational& Rational::operator-=(const Rational &rhs) {
	// Do division first to reduce chance of integer overflow errors.
	int lcm, pnew;
	lcm = rhs.GetQ()/gcd(this->den, rhs.GetQ());
	lcm *= this->den;
	pnew = this->num*(lcm/this->den);
	pnew -= rhs.GetP()*(lcm/rhs.GetQ());

    num = pnew;
    den = lcm;

	return *this; 

}

Rational Rational::operator-(const Rational &rhs) const {
    Rational result = *this;
    return result -= rhs;

}

Rational Rational::Inverse() {
	if(num == 0)
		return Rational(num, den);
	else
		return Rational(den, num);
}

bool Rational::IsZero() const {
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

CoordinatesRange::CoordinatesRange( std::vector<Rational>::iterator beginIt_, int length_ ) {
    beginIt = beginIt_;
    length = length_;
}

CoordinatesRange::CoordinatesRange( std::vector<Rational> coord) {
    beginIt = coord.begin();
    length = coord.size();
}

CoordinatesRange& CoordinatesRange::operator += ( const CoordinatesRange& rhs) {
    std::vector<Rational>::iterator iIt = beginIt, jIt = rhs.beginIt;
    int numAddable;
    if (length < rhs.length)
        numAddable = length;
    else
        numAddable = rhs.length;

    for ( int i = 0; i < numAddable; i++, iIt++, jIt++)  
        *iIt += *jIt; 

    return *this;

}

CoordinatesRange& CoordinatesRange::operator -= ( const CoordinatesRange& rhs) {
    std::vector<Rational>::iterator iIt = beginIt, jIt = rhs.beginIt;
    int numAddable;
    
    if (length < rhs.length)
        numAddable = length;
    else
        numAddable = rhs.length;

    for ( int i = 0; i < numAddable; i++, iIt++, jIt++)  
        *iIt -= *jIt; 

    return *this;

}

CoordinatesRange& CoordinatesRange::operator *= (const Rational& rhs) {
    std::vector<Rational>::iterator iIt = beginIt;

    for (int i = 0; i < length; i++, iIt++)
        *iIt *= rhs;

    return *this;
}

std::pair<CoordinatesRange, CoordinatesRange> CoordinatesRange::splitToPair() {
    CoordinatesRange first(begin(), length / 2),
                     second(middle(), length / 2);
    
    return std::pair<CoordinatesRange, CoordinatesRange>(first, second);
}

//////////////////////////////////////////////

Coordinates& Coordinates::operator += (const Coordinates& rhs) {
    std::vector<Rational>::iterator iIt;
    std::vector<Rational>::const_iterator jIt;

    for( iIt = coordinates.begin(), jIt = rhs.coordinates.begin()
       ; iIt < coordinates.end() && jIt < rhs.coordinates.end()
       ; iIt++, jIt++)
        
        *iIt = *iIt + *jIt; 

    return *this;
         
}

Coordinates Coordinates::operator + (const Coordinates& rhs) const {
    Coordinates result = *this;
    result += rhs;
    return result;
} 

Coordinates& Coordinates::operator -= (const Coordinates& rhs) {
    std::vector<Rational>::iterator iIt;
    std::vector<Rational>::const_iterator jIt;

    for( iIt = coordinates.begin(), jIt = rhs.coordinates.begin()
       ; iIt < coordinates.end() && jIt < rhs.coordinates.end()
       ; iIt++, jIt++)
        
        *iIt = *iIt - *jIt; 

    return *this;

} 

Coordinates Coordinates::operator - (const Coordinates& rhs) const {
    Coordinates result = *this;
    result -= rhs;
    return result;
}

Coordinates& Coordinates::operator *= (const Rational& scaling) {
    for ( std::vector<Rational>::iterator iIt = coordinates.begin()
        ; iIt != coordinates.end()
        ; iIt++)
        *iIt *= scaling; 

    return *this;
}

Coordinates Coordinates::operator * (const Rational& scaling) const {
   Coordinates result = *this;
   return result; 
}

bool Coordinates::isZero() const {
    bool stillAllZero = true;
    std::vector<Rational>::const_iterator coordIt, 
                                          end = coordinates.end();

    for (coordIt = coordinates.begin(); coordIt != end && stillAllZero; coordIt++) {
        if (!(*coordIt).IsZero())
            stillAllZero = false;
    }

    return stillAllZero;

}
/////////////////////////////////////////////

FieldElement::FieldElement() {
    level = -1;
    field = NULL;
}

FieldElement::FieldElement(QuadraticFieldTower* field_, int level_, std::vector<Rational> coords_) {
	coords = coords_;
    level = level_; 
    field = field_;
}

std::string FieldElement::Print() {
	int whichroots, currentroot;
	std::ostringstream ss;
	std::string name = std::string();
	for(int i=0; i < coords.size(); i++) {
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


FieldElement& FieldElement :: operator += (const FieldElement& rhs) {

	int maxindex = coords.size(); 

	if(rhs.coords.size() < maxindex)
		maxindex = rhs.coords.size();

	for(int i = 0; i < maxindex; i++) 
		coords[i] = coords[i] + rhs.coords[i];	

	return *this;
}

FieldElement& FieldElement :: operator -= (const FieldElement& rhs) {

	int maxindex = coords.size(); 

	if(rhs.coords.size() < maxindex)
		maxindex = rhs.coords.size();

	for(int i = 0; i < maxindex; i++) 
		coords[i] = coords[i] - rhs.coords[i];	

	return *this;
}

FieldElement& FieldElement :: operator *= (const FieldElement& rhs) {
	coords = field -> multiply(coords, rhs.coords, level, rhs.level);
    return *this;
}

FieldElement FieldElement :: operator + (const FieldElement& rhs) {
    std::vector<Rational> result = coords;
    FieldElement resultField(field, level, result);

    return resultField += rhs; 
}

FieldElement FieldElement :: operator - (const FieldElement &rhs) {
    std::vector<Rational> result = coords;
    FieldElement resultField(field, level, result);

    return resultField -= rhs;
}

FieldElement FieldElement :: operator * (const FieldElement& rhs) {
    FieldElement solution = *this;
    return solution *= rhs;
}

/////////////////////////////////////////////

QuadraticFieldTower::QuadraticFieldTower(Rational square_) {
	std::complex<double> newcomplexroot = std::complex<double>(square_.ToFloat(), 0);
	newcomplexroot = std::sqrt(newcomplexroot);

	squares = std::vector< std::vector<Rational> >(1, std::vector<Rational>(1, square_) );
	complexroots = std::vector< std::complex<double> >(1, newcomplexroot); 
	degree = 2;
	numsquares = 1; 
}

std::string QuadraticFieldTower::Print() {
	std::ostringstream oss;
	std::string name = std::string();

	name += "Q( ";
	for(int i=0; i<numsquares; i++){
		name += "r";
		oss.str(std::string());
		oss << i+1;
		name += oss.str();
		if(i < numsquares-1)
			name += ", ";		
	}
	
	name += " ) where\n";
	name += PrintRootList();
	return name;	

}

std::string QuadraticFieldTower::PrintRootList() {
	std::ostringstream oss;
	std::string name = std::string();
	int numcoords = 1, whichroot, reducedindex;

	for(int i = 0; i < numsquares; i++, numcoords *= 2) {
		name += "r";
		oss.str(std::string());
		oss << i+1;
		name += oss.str();
		name += " = Sqrt( ";
		name += FieldElement(this, i, squares[i]).Print(); //PrintCoords(squares[i]);	
		name += " )\n";		
	}
	return name;
}

std::complex<double> QuadraticFieldTower::toComplex(std::vector<Rational> & coords) {
	std::complex<double> result(0.0,0.0), temp;
	int whichroot, mask;
	
	for(int i=0; i<coords.size(); i++) {
		temp = std::complex<double>(1.0, 0.0);
		whichroot = 0;
		mask = 1;
		while(mask < coords.size()) {
			if((i & mask) != 0) {
				temp *= complexroots[whichroot];
			}
			mask *= 2;
			whichroot++;
		}
		temp *= std::complex<double>( coords[i].ToFloat(), 0.0);
		result += temp;
	}
	return result;
}

std::complex<double> QuadraticFieldTower::toComplex(FieldElement &x) {
    return toComplex(x.coords);
}

bool QuadraticFieldTower::sqrtExists(std::vector<Rational> &coords) {
    if (coords.size() < 2) {
       sqrtExistsResult = coords[0].FindSqrt(); 
       sqrtResult[0] = coords[0].GetSqrt(); 
       return sqrtExistsResult;
    }

    return true;
}

void QuadraticFieldTower::multiplyInPlace(const std::vector<Rational>& lhs, const std::vector<Rational>& rhs, int lhsLevel, int rhsLevel, std::vector<Rational>& solution) {

    std::vector<Rational>::const_iterator lhsIt, rhsIt;
    int newLhsSize, newRhsSize, newLhsLevel;

    if (lhs.size() > rhs.size()) {
        lhsIt = lhs.begin();
        rhsIt = rhs.begin();
        newLhsLevel = lhsLevel;
        newLhsSize = lhs.size();
        newRhsSize = rhs.size();
    }
    else {
        lhsIt = rhs.begin();
        rhsIt = lhs.begin();
        newLhsLevel = rhsLevel;
        newLhsSize = rhs.size();
        newRhsSize = lhs.size();
    }

    solution = std::vector<Rational> (newLhsSize); 

    multiplyLhsLargest(lhsIt, rhsIt, solution.begin(), newLhsSize, newRhsSize, newLhsLevel);

}

std::vector<Rational> QuadraticFieldTower::multiply (const std::vector<Rational>& lhs, const std::vector<Rational>& rhs, int lhsLevel, int rhsLevel) {

    std::vector<Rational> solution;

    multiplyInPlace(lhs, rhs, lhsLevel, rhsLevel, solution); 
    return solution;
 
}

Coordinates QuadraticFieldTower::multiply (const Coordinates &lhs, const Coordinates& rhs, int lhsLevel, int rhsLevel) {
    std::vector<Rational> result;

    multiplyInPlace(lhs.coordinates, rhs.coordinates, lhsLevel, rhsLevel, result);
    return Coordinates(result, lhsLevel);
}

Coordinates QuadraticFieldTower::inverse(Coordinates& x, int level) {
    // MISSING non-zero check
    Coordinates a, b, scratch; // x = a + r * b. Think of a and b as existing in the sub-field.
    std::vector<Rational> result(x.size());
    int middleoffset, sublevel = level - 1;

    if(level < 1) {
       result[0] = x.coordinates[0].Inverse(); 
       return Coordinates(result, level);
    }
    
    middleoffset = x.size() / 2;
    a = Coordinates( std::vector<Rational> ( x.coordinates.begin()
                                           , x.coordinates.begin() + middleoffset
                                           ) 
                   , sublevel );
    b = Coordinates( std::vector<Rational> ( x.coordinates.begin() + middleoffset
                                           , x.coordinates.end()
                                           ) 
                   , sublevel );
    

    // a^2 - r^2 * b^2.
    scratch = Coordinates(squares[level - 1], sublevel);
    scratch = multiply(b, scratch, level - 1, level - 1);
    b *= Rational(-1, 1);
    scratch = multiply(b, scratch, level - 1, level - 1);
    scratch += multiply(a, a, level - 1, level - 1);
    scratch = inverse(scratch, level - 1);
    
    // Need to concatenate a and b. Then multiply by scratch to get the full inverse. 

    return Coordinates(std::vector<Rational>(0), -1);
}

void QuadraticFieldTower::multiplyLhsLargest (std::vector<Rational>::const_iterator lhsIt, std::vector<Rational>::const_iterator rhsIt, std::vector<Rational>::iterator solutionIt, const int lhsLength, const int rhsLength, const int lhsLevel) const {

    std::vector<Rational> scratch (lhsLength);
    int sublength = lhsLength / 2, sublevel = lhsLevel - 1;
    std::vector<Rational>::const_iterator lhsMiddleIt = lhsIt + lhsLength / 2, 
                                          rhsMiddleIt = rhsIt + rhsLength / 2; 
    std::vector<Rational>::iterator solMiddleIt = solutionIt + lhsLength / 2, 
                                    scratchMiddleIt = scratch.begin() + lhsLength / 2; 

    if (lhsLength < 2) {
        *solutionIt = *lhsIt * *rhsIt;
        return;
    }
   
    else if( lhsLength > rhsLength) {

        // First handle multiplication of cross terms. No need to use root here. 

        multiplyLhsLargest(lhsMiddleIt, rhsIt, scratchMiddleIt, sublength, rhsLength, sublevel);

        // Now handle the non-cross terms.

        multiplyLhsLargest(lhsIt, rhsIt, scratch.begin(), sublength, rhsLength, sublevel); 

        for( std::vector<Rational>::iterator iIt = scratch.begin(), solIt = solutionIt 
           ; iIt != scratch.end() 
           ; iIt++, solIt++
           ) 
           
           *solIt = *iIt;

    }

    else { 

        // First handle multiplication of cross terms. No need to use root here. 

        multiplyLhsLargest(lhsIt, rhsMiddleIt, scratch.begin(), sublength, sublength, sublevel);
        multiplyLhsLargest(lhsMiddleIt, rhsIt, scratchMiddleIt, sublength, sublength, sublevel);

        for( std::vector<Rational>::iterator iIt = scratch.begin(), jIt = scratchMiddleIt, solIt = solMiddleIt
           ; iIt < scratchMiddleIt 
           ; iIt++, jIt++, solIt++
           ) 

            *solIt = *iIt + *jIt;
        
        // Now do non-cross terms. Need to use the root here.

        multiplyLhsLargest(lhsMiddleIt, rhsMiddleIt, scratch.begin(), sublength, sublength, sublevel);
        multiplyLhsLargest(scratch.begin(), squares[lhsLevel-1].begin(), scratchMiddleIt, sublength, sublength, sublevel); 
        multiplyLhsLargest(lhsIt, rhsIt, scratch.begin(), sublength, sublength, sublevel);

        for( std::vector<Rational>::iterator iIt = scratch.begin(), jIt = scratchMiddleIt, solIt = solutionIt
           ; iIt < scratchMiddleIt
           ; iIt++, jIt++, solIt++) 
            
            *solIt = *iIt + *jIt;

    }
         
}

void QuadraticFieldTower::multiplyLhsLargest(CoordinatesRange lhs, CoordinatesRange rhs, CoordinatesRange solution, int lhsLevel) {
    std::vector<Rational> scratch1, scratch2, scratch;
    CoordinatesRange scratchRange;

    int sublength = lhs.size() / 2, sublevel = lhsLevel - 1;
    std::pair<CoordinatesRange, CoordinatesRange> lhsPair, rhsPair, solPair, scratchPair;

    if (lhsLevel < 1) {
        *solution.begin() = *lhs.begin() * *rhs.begin();
        return;
    }

    else if( lhs.size() > rhs.size() ) {
        lhsPair = lhs.splitToPair();
        solPair = solution.splitToPair();

        // First handle multiplication of cross terms. No need to use root here. 

        multiplyLhsLargest(lhsPair.second, rhs, solPair.second, sublevel);

        // Now handle the non-cross terms.

        multiplyLhsLargest(lhsPair.first, rhs, solPair.first, sublevel); 

    }

    else { 
        lhsPair = lhs.splitToPair();
        rhsPair = rhs.splitToPair();
        solPair = solution.splitToPair();
        scratch = std::vector<Rational>( sublength );
        scratchRange = CoordinatesRange(scratch);
 
        // First handle multiplication of cross terms. No need to use root here. 

        multiplyLhsLargest(lhsPair.first, rhsPair.second, solPair.second, sublevel);
        multiplyLhsLargest(lhsPair.second, rhsPair.first, scratchRange, sublevel);

        solPair.second += scratchRange; 

        // Now do non-cross terms. Need to use the root here.

        multiplyLhsLargest(lhsPair.second, rhsPair.second, scratchRange, sublevel);
        multiplyLhsLargest(scratchRange, CoordinatesRange(squares[lhsLevel-1]), solPair.first, sublevel); 
        multiplyLhsLargest(lhsPair.first, rhsPair.first, scratchRange, sublevel);
        
        solPair.first += scratchRange;

    } 
}

std::vector<Rational> QuadraticFieldTower::getSqrt(std::vector<Rational> coords, int level) {

    FieldElement a, b; // Think of coordinates as a + r * b where a and b are in the sub-field. 
    std::vector<Rational> subSqrt, subSqrt2;
    int middleoffset;

    if (level < 1) {
       sqrtExistsResult = coords[0].FindSqrt(); 
       if(sqrtExistsResult)
           return std::vector<Rational> (1, coords[0].GetSqrt()); 
       else
           return std::vector<Rational> (0);
    }

    middleoffset = coords.size() / 2;
    a = FieldElement( this
                    , level - 1
                    , std::vector<Rational>(coords.begin(), coords.begin() + middleoffset
                    ));    
    b = FieldElement( this
                    , level - 1
                    , std::vector<Rational>(coords.begin() + middleoffset, coords.begin()
                    ));

   subSqrt = getSqrt( (a * a - b * b).coords, level - 1); 
   if(!sqrtExistsResult)
        return std::vector<Rational>(0);

   subSqrt2 = getSqrt( (a + FieldElement(this, level - 1, subSqrt)).coords
                     , level - 1); 

   return subSqrt2;
}


Coordinates QuadraticFieldTower::getSqrt(Coordinates& x, int level) {
    Coordinates a, b; // Think of x = a + r * b where a and b are coordinates in the subfield.
    Coordinates subSqrt, subSqrt2, scratch; 
    std::vector <Rational> result;
    int middleoffset, sublevel;

    if (level < 1) {
        sqrtExistsResult = x.coordinates[0].FindSqrt(); 
           if(sqrtExistsResult) 
               result = std::vector<Rational> (1, x.coordinates[0].GetSqrt()); 
           else 
               result = std::vector<Rational> (0);
           return Coordinates(result, level);
    }

    middleoffset = x.coordinates.size() / 2;
    a = Coordinates( std::vector<Rational>( x.coordinates.begin()
                                          , x.coordinates.begin() + middleoffset
                                          ) 
                   , sublevel );
    b = Coordinates( std::vector<Rational>( x.coordinates.begin() + middleoffset
                                          , x.coordinates.begin() + middleoffset
                                          ) 
                   , sublevel );

    if(!b.isZero()) {

        subSqrt = multiply(a, a, sublevel, sublevel);
        subSqrt2 = multiply(b, b, sublevel, sublevel);
        subSqrt2 = multiply(Coordinates(squares[level], sublevel), subSqrt2, sublevel, sublevel);
        subSqrt -= subSqrt2;
        subSqrt = getSqrt(subSqrt, level - 1);
   
        if(! sqrtExistsResult) 
           return Coordinates(std::vector<Rational>(0), -1); 

        subSqrt2 = (a + subSqrt) * Rational(1,2);
        scratch = getSqrt(subSqrt2, sublevel);

        if(sqrtExistsResult && !scratch.isZero()) {

           subSqrt = b * Rational(1,2);
           subSqrt = multiply(subSqrt, inverse(scratch, sublevel), sublevel, sublevel);
           //return (scratch, subSqrt); 
        }
        
        subSqrt2 = (a - subSqrt) * Rational(1,2);
        scratch = getSqrt(subSqrt2, sublevel);

        if(sqrtExistsResult && !scratch.isZero()) {

            subSqrt = b * Rational(1, 2);
            subSqrt = multiply(subSqrt, inverse(scratch, sublevel), sublevel, sublevel); 
            // return (scratch, subSqrt);
        }

        sqrtExistsResult = false;
        return Coordinates(std::vector<Rational>(0), -1);

    }
    else { // now  in the case that b = 0
        subSqrt = getSqrt(a, sublevel);
        if(sqrtExistsResult) {
            // return (subSqrt, 0);
        }
    
        //scratch = a * inverse( Coordinates(squares[sublevel]), sublevel);
        subSqrt = getSqrt(scratch, sublevel);
        if(sqrtExistsResult) {
            // return (0, subSqrt);
        }
        
        sqrtExistsResult = false;
        return Coordinates(std::vector<Rational>(0), -1);    
    }

}
//////////////////////////////////////////////


