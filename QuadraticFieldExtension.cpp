//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"


int main() {
	
	Rational q = Rational(1,2);
	QuadraticFieldTower fieldtower(q);
	std::vector<Rational> root(2, Rational(1,2)), coords(4);
	
	fieldtower.AddSquare(root);

	coords[0] = Rational(1,1);
	coords[1] = Rational(1,2);
	coords[2] = Rational(2,3);
	coords[3] = Rational(3,4);
	fieldtower.AddSquare(coords);

	std::cout << "Rational q = " << q.print() << std::endl;
	std::cout << "fieldtower = " << fieldtower.Print() << std::endl;
	return 0;

}
