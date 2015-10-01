//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"

using namespace std;

void main() {

	QuadraticField Q;
	Rational r, sum, product;
	std::vector<Rational> elements(2,Rational(1,1));
	std::vector<Rational> elementsB(4); 
	elementsB[0] = Rational(1,2); 
	elementsB[1] = Rational(2,1); 
	elementsB[2] = Rational(3,1); 
	elementsB[3] = Rational(1,4);

	int p, q;
	cout << "Let us check the rationality of the Sqrt(p/q)" << endl
			<< "Please input integer p = " << endl;
	cin >> p;
	cout << "Please input integer q = " << endl;
	cin >> q;
	r = Rational(p,q);
	cout << endl;
	
	if(r.FindSqrt()) {
		Rational SqrtResult(r.GetSqrtP(), r.GetSqrtQ());

		cout << "Square root of "<< r.print()
			<< " is = " << SqrtResult.print() << endl;  

		sum = r+SqrtResult;
		cout << r.print() << " + " << SqrtResult.print() << " = " 
			<< sum.print() << endl;

		product = r*SqrtResult;
		cout << r.print() << " * " << SqrtResult.print() << " = "
			<< product.print() << endl;
	}
	else {

		cout << "There is no rational square root of " << r.print() << endl;
		QuadraticField QF(r);
		cout << "Constructed Quadratic Extension " << QF.PrintName() << endl;
		Rational coordinates[4] = {Rational(1,2), Rational(2,3), Rational(3,4), Rational(4,5)};
		cout << "Has elements such as " << QF.Print(elements) << endl << endl;

		elements[0] = Rational(1,1);
		elements[1] = Rational(1,1);
		QuadraticField QFB(&QF, elements);
		cout << "Constructed Quadratic Extension " << QFB.PrintName() << endl
				<< "Has elements such as " << QFB.Print(elementsB) << endl
				<< "This is messy. To clean up, we may use variables" << endl << endl
				<< "Can be rewritten as extension " << QFB.PrintNameR() << endl
				<< "Has elements such as " << QFB.PrintR(elementsB) << endl;
		QFB.Product(elementsB, elementsB);
		cout << "The square is " << QFB.PrintR(QFB.Result) << endl << endl;

		cout << "Now let us compute some square roots in " << QF.PrintNameR() << endl << endl;
		for(int i=0; i<100; i++){
			for(int j=0; j<100; j++) {
				elements[0] = Rational(i,1);
				elements[1] = Rational(j,1);
				if(QF.FindSqrt(elements)) {
					cout << "Sqrt(" << QF.PrintR(elements) << ") = ";
					elements = QF.GetSqrt();
					cout << QF.PrintR(elements) << endl;
					QF.Product(elements, elements);
				}
				
			}
		}
	}
	cout << "A pause. Enter anything.";
	cin >> p;
}