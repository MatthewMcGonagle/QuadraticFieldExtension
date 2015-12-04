//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"


int main() {
	
	Rational q = Rational(1,2);
	QuadraticField F1=QuadraticField(q);
	std::vector<Rational> root(2);
	root[0] = Rational(1,2);
	root[1] = Rational(1,2);
	FieldElement elem1 = FieldElement(&F1, root.begin());

	QuadraticField F2 = QuadraticField(&F1, root);
	
	std::cout << "Rational q = " << q.print() << std::endl;
	std::cout << "Field F1 = Q( " << elem1.Print() << " )" << std::endl;
	std::cout << "Field F2 = " << std::endl;
	return 0;

}
