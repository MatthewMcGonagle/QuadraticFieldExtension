//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"
#include <math.h>


int main() {
	
    Rational testRational1(1,2), testRational2(2,3), testOp;
    Rational coordInit1[] = {Rational(2, 1), Rational(3, 1)},
             coordInit2[] = {Rational(1, 1), Rational(1, 1)},
             coordInit3[] = {Rational(1, 1), Rational(1, 1), Rational(1,1), Rational(1,1)};
    std::vector<Rational> coords1V(coordInit1, coordInit1 + sizeof(coordInit1) / sizeof(Rational)),
                          coords2V(coordInit2, coordInit2 + sizeof(coordInit2) / sizeof(Rational)),
                          coords3V(coordInit3, coordInit3 + sizeof(coordInit3) / sizeof(Rational)),
                          result_init(2);
    Coordinates rootToAdd(std::vector<Rational>(1, Rational(2,1))),
                rootToAdd2(std::vector<Rational>(2, Rational(1,1))),
                result(std::vector<Rational>(2)),
                coords1(coords1V),
                coords2(coords2V),
                coords3(coords3V);
    QuadraticFieldTower tower;

    std::cout << std::endl
              << "Test of Rational operations" << std::endl
              << testRational1.print() << " + " << testRational2.print() << " = ";
    testOp = testRational1 + testRational2;
    std::cout << testOp.print() << std::endl
              << testRational1.print() << " * " << testRational2.print() << " = ";
    testOp = testRational1 * testRational2;
    std::cout << testOp.print() << std::endl << std::endl;

    std::cout << "Test of Coordinates printing" << std::endl
              << "coords1 = " << coords1.print() << std::endl
              << "coords2 = " << coords2.print() << std::endl; 

    std::cout << std::endl 
              << "Test Field Tower" << std::endl
              << "Adding square root of " << rootToAdd.print() << std::endl; 
    tower.addIfNoSqrRoot(rootToAdd);
    std::cout << "Field Tower Squares of Roots are:" << std::endl
              << tower.print();
    std::cout << "Testing Addition" << std::endl
              << coords1.print() << "  +  " << coords2.print() << "  =  ";
    result = tower.add(coords1, coords2);
    std::cout << result.print() << std::endl;
    std::cout << "Testing Multiplication" << std::endl
              << coords1.print() << "  *  " << coords2.print() << "  =   "; 
    result = tower.multiply(coords1, coords2);
    std::cout << result.print() << std::endl;

    std::cout << "Now do a Numerical Test" << std::endl
              << "Complex Conversion of " << coords1.print() << " = "
              << tower.convertToComplex(coords1) << std::endl;
    std::cout << "Complex Conversion of " << coords2.print() << " = "
              << tower.convertToComplex(coords2) << std::endl;
    std::cout << "Complex Conversion of Product " << result.print()
              << " = " << tower.convertToComplex(result) << std::endl;
    std::cout << "Pure Complex Multiplication is = " 
              << tower.convertToComplex(coords1) * tower.convertToComplex(coords2) << std::endl;


    std::cout << std::endl
              << "Adding Square Root of " << rootToAdd2.print() << std::endl;
    tower.addIfNoSqrRoot(rootToAdd2);
    std::cout << "Field Tower Squares of Roots are:" << std::endl
              << tower.print();
    std::cout << "Testing Multiplication" << std::endl
              << "(" << coords3.print() << ") * (" << coords3.print() << ") = ";
    result = tower.multiply(coords3, coords3);
    std::cout << result.print() << std::endl; 
  
    std::cout << "Do Numerical Tests" << std::endl
              << "Complex Conversion of " << coords3.print() << " = "
              << tower.convertToComplex(coords3) << std::endl; 
    std::cout << "Complex Conversion of Product Result = "
              << tower.convertToComplex(result) << std::endl;
    std::cout << "Pure Complex Multiplication = " 
              << tower.convertToComplex(coords3) * tower.convertToComplex(coords3) << std::endl;

    std::cout << std::endl
              << "Finished" << std::endl;

	return 0;

}
