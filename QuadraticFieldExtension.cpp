//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"


int main() {
	
	Rational q = Rational(1,2);
	QuadraticField F1=QuadraticField(q);
	std::vector<Rational> root(2), coords(4);
	root[0] = Rational(1,2);
	root[1] = Rational(1,2);
	FieldElement elem1 = FieldElement(&F1, root.begin());

	QuadraticField F2 = QuadraticField(&F1, root);
	coords[0] = Rational(1,1);
	coords[1] = Rational(1,2);
	coords[2] = Rational(2,3);
	coords[3] = Rational(3,4);
	FieldElement elem2 = FieldElement(&F2, coords.begin());

	QuadraticField F3 = QuadraticField(&F2, coords);

	std::cout << "Rational q = " << q.print() << std::endl;
	std::cout << "Field F1 = " << F1.Print() << std::endl;
	std::cout << "elem1 = " << elem1.Print() << std::endl;
	std::cout << "Field F2 = " << F2.Print() << std::endl;
	std::cout << "elem2 = " << elem2.Print() << std::endl;
	std::cout << "Field F3 = " << F3.Print() << std::endl;	
	return 0;

}
