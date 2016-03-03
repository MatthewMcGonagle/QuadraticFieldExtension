//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"
#include <math.h>


int main() {
	
	Rational q = Rational(1,2);
	QuadraticFieldTower fieldtower(q);
	std::vector<Rational> root(2, Rational(1,2)), coords(4), coords2(4);
	
	fieldtower.AddSquare(root);

	for(int i=0; i < 4; i++) {
		coords[i] = Rational(1,i+1);
		coords2[i] = Rational(pow(-1, i)*(i+1),i+2);
	}
	fieldtower.AddSquare(coords);

	
	std::cout << "Rational q = " << q.print() << std::endl;
	std::cout << "fieldtower = " << fieldtower.Print() << std::endl;
	std::cout << "   ( " << fieldtower.PrintCoords(coords) << " )"  << std::endl
		  << " + ( " << fieldtower.PrintCoords(coords2) << " )" << std::endl;
	fieldtower.Add(coords, coords2);
	std::cout << " = " << fieldtower.PrintCoords(coords) << std::endl << std::endl;

	coords.resize(2);
	coords[0] = Rational(1,1);
	coords[1] = Rational(1,1);
	std::cout << "   (" << fieldtower.PrintCoords(coords) << ")" << std::endl
		  << " * (" << fieldtower.PrintCoords(coords) << ")" << std::endl;
	fieldtower.Product(coords, coords);
	std::cout << " =  " << fieldtower.PrintCoords(coords);

	return 0;

}
