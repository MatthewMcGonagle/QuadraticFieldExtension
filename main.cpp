//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"
#include <math.h>


int main() {
	
	Rational q = Rational(1,2);
	QuadraticFieldTower fieldtower(q);
	std::vector<Rational> root(2, Rational(1,2)), coords(4), coords2(4);
	std::complex<float> check1, check2;

    // Add in two roots to build tower of fields.

	fieldtower.AddSquare(root);

	for(int i=0; i < 4; i++) {
		coords[i] = Rational(1,i+1);
		coords2[i] = Rational(pow(-1, i)*(i+1),i+2);
	}
	fieldtower.AddSquare(coords);

    // Check printing functions.
	
	std::cout << "Rational q = " << q.print() << std::endl;
	std::cout << "fieldtower = " << fieldtower.Print() << std::endl;

    // Check FieldElement.

    Rational myRationals[] = {Rational(1, 2), Rational(2, 3), Rational(3, 4), Rational(4, 5)},
             myRationals2[] = {Rational(1, 2), Rational(3, 4), Rational(5, 6), Rational(7, 8)};

    std::vector<Rational> coords3(myRationals, myRationals + sizeof(myRationals) / sizeof(Rational)),
                          coords4(myRationals2, myRationals2 + sizeof(myRationals2) / sizeof(Rational));
    FieldElement myElement(&fieldtower, 2, coords3),
                 myElement2(&fieldtower, 2, coords4),
                 mySum, 
                 myProduct; 

    std::cout << "myElement = " << myElement.Print() << std::endl
              << "myElement2 = " << myElement2.Print() << std::endl;

    mySum = myElement + myElement2;
    myProduct = myElement * myElement2;

    std::cout << "Their sum = " << mySum.Print() << std::endl
              << "Their product = " << myProduct.Print() << std::endl << std::endl;

    // Double check sums and products using conversion to complex.

    std::complex<double> myComplex = fieldtower.toComplex(myElement),
                         myComplex2 = fieldtower.toComplex(myElement2);

    std::cout << "Double check of sums and products using complex number conversion." << std::endl
              << "As a complex number, myElement = " << myComplex << std::endl
              << "As a complex number, myElement2 = " << myComplex2 << std::endl
              << "Sum of complex = " << myComplex + myComplex2 << std::endl
              << "Complex conversion of sum = " << fieldtower.toComplex(mySum) << std::endl
              << "Product of complex = " << myComplex * myComplex2 << std::endl
              << "Complex conversion of product = " << fieldtower.toComplex(myProduct) << std::endl;
                 
    return 0;

}
