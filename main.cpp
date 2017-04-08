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

    Rational myRationals[] = {Rational(0, 0), Rational(1, 1), Rational(1, 1), Rational(1, 1)},
             myRationals2[] = {Rational(0, 1), Rational(1, 1), Rational(0, 2), Rational(0, 1)};

    std::vector<Rational> coords3(myRationals, myRationals + sizeof(myRationals) / sizeof(Rational)),
                          coords4(myRationals2, myRationals2 + sizeof(myRationals2) / sizeof(Rational));
    FieldElement myElement(&fieldtower, 2, coords3),
                 myElement2(&fieldtower, 2, coords4);

    std::cout << "myElement = " << myElement.Print() << std::endl
              << "myElement2 = " << myElement2.Print() << std::endl
              << "Their sum = " << (myElement + myElement2).Print() << std::endl
              << "Their product = " << (myElement * myElement2).Print() << std::endl << std::endl;

    // // Check the adding function.

	// std::cout << "   ( " << fieldtower.PrintCoords(coords) << " )"  << std::endl
	// 	  << " + ( " << fieldtower.PrintCoords(coords2) << " )" << std::endl;

	// fieldtower.Add(coords, coords2);
	// std::cout << " = " << fieldtower.PrintCoords(coords) << std::endl << std::endl;

	// //coords.resize(2);
	// coords[0] = Rational(1,1);
	// coords[1] = Rational(0,1);
	// coords[2] = Rational(0,1);
	// coords[3] = Rational(1,1);	
	// std::cout << "   (" << fieldtower.PrintCoords(coords) << ")" << std::endl
	// 	  << " * (" << fieldtower.PrintCoords(coords) << ")" << std::endl;
	// check1 = fieldtower.CoordsToComplex(coords);
	// check1 *= check1;

	// fieldtower.Product(coords, coords);
	// std::cout << " =  " << fieldtower.PrintCoords(coords) << std::endl;
	// check2 = fieldtower.CoordsToComplex(coords);

	// std::cout << "Numerical double check" << std::endl
	// 	  << "check1 = " << check1 << std::endl
	// 	  << "check2 = " << check2 << std::endl;

	return 0;

}
