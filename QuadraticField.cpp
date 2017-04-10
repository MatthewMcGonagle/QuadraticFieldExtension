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

Rational Rational::operator+(const Rational &r) const {
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

Coordinates& Coordinates::operator += (const Coordinates& rhs) {
    std::vector<Rational>::iterator iIt;
    std::vector<Rational>::const_iterator jIt;

    for( iIt = coordinates.begin(), jIt = rhs.coordinates.begin()
       ; iIt < coordinates.end() && jIt < rhs.coordinates.end()
       ; iIt++, jIt++)
        
        *iIt = *iIt + *jIt; 

    return *this;
         
}

Coordinates Coordinates::operator + (const Coordinates& rhs) {
    Coordinates result(coordinates);
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

Coordinates Coordinates::operator - (const Coordinates& rhs) {
    Coordinates result(coordinates);
    result -= rhs;
    return result;
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

void QuadraticFieldTower::AddSquare(std::vector<Rational> coords) {
	std::complex<double> result(0.0, 0.0), temp(0.0, 0.0);
	int whichsquare, mask;

	if (coords.size() == degree) {
		squares.push_back(coords);
		for(int i=0; i<degree; i++) {
			temp = std::complex<double>(1.0,0.0);
			whichsquare = 0;
			mask = 1;
			while(mask < degree) {
				if( (i & mask) != 0) {
					temp *= complexroots[whichsquare];
				}	
				mask *= 2;
				whichsquare++;
			}
			temp *= std::complex<double>( coords[i].ToFloat(), 0.0);
			result += temp;	
		}
		result = std::sqrt(result);
		complexroots.push_back(result);
		degree *= 2;
		numsquares++;	
	}

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

std::vector<Rational> QuadraticFieldTower::multiply (const std::vector<Rational>& lhs, const std::vector<Rational>& rhs, int lhsLevel, int rhsLevel) {

    std::vector<Rational> solution;
    std::vector<Rational>::const_iterator lhsIt, rhsIt;
    int newLhsSize, newRhsSize, level;

    if (lhs.size() > rhs.size()) {
        lhsIt = lhs.begin();
        rhsIt = rhs.begin();
        level = lhsLevel;
        newLhsSize = lhs.size();
        newRhsSize = rhs.size();
    }
    else {
        lhsIt = rhs.begin();
        rhsIt = lhs.begin();
        level = rhsLevel;
        newLhsSize = rhs.size();
        newRhsSize = lhs.size();
    }

    solution = std::vector<Rational> (lhs.size()); 

    multiplyLhsLargest(lhs.begin(), rhs.begin(), solution.begin(), newLhsSize, newRhsSize, level);

    return solution;
 
}

void QuadraticFieldTower::multiplyLhsLargest (std::vector<Rational>::const_iterator lhsIt, std::vector<Rational>::const_iterator rhsIt, std::vector<Rational>::iterator solutionIt, int lhsLength, int rhsLength, int lhsLevel) {
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
   
    if( lhsLength > rhsLength) {

        // First handle multiplication of cross terms. No need to use root here. 

        multiplyLhsLargest(lhsMiddleIt, rhsIt, scratch.begin(), sublength, rhsLength, sublevel);

        // Now handle the non-cross terms.

        multiplyLhsLargest(lhsIt, rhsIt, scratchMiddleIt, sublength, rhsLength, sublevel); 

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
         
        // Now do non-cross terms. Need to use the root here.

        multiplyLhsLargest(lhsMiddleIt, rhsMiddleIt, scratch.begin(), sublength, sublength, sublevel);
        multiplyLhsLargest(scratch.begin(), squares[lhsLevel-1].begin(), scratchMiddleIt, sublength, sublength, sublevel); 
        multiplyLhsLargest(lhsIt, rhsIt, scratch.begin(), sublength, sublength, sublevel);

        for(std::vector<Rational>::iterator iIt = scratch.begin(), jIt = scratchMiddleIt, solIt = solutionIt;
            iIt < scratchMiddleIt;
            iIt++, jIt++, solIt++) {
            
            *solIt = *iIt + *jIt;
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
//////////////////////////////////////////////


