//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"

using namespace std;

void main() {

	QuadraticField Q;
	Rational r, sum, product;
	std::vector<Rational> elements(2,Rational(1,1));

	int p = 5, q=2;
	cout << "Please input integer p = " << endl;
	cin >> p;
	cout << "Please input integer q = " << endl;
	cin >> q;
	r = Rational(p,q);
	
	if(r.FindSqrt()) {
		Rational SqrtResult(r.GetSqrtP(), r.GetSqrtQ());

		cout << "Square root of "<< r.print()
			<< " is = " << SqrtResult.print() << endl;  

		sum = r.Sum(SqrtResult);
		cout << r.print() << " + " << SqrtResult.print() << " = " 
			<< sum.print() << endl;

		product = r.Product(SqrtResult);
		cout << r.print() << " * " << SqrtResult.print() << " = "
			<< product.print() << endl;
	}
	else {

		cout << "There is no rational square root of " << r.print() << endl;
		QuadraticField QF(r);// = QuadraticField(r);
		cout << "Constructed Quadratic Extension" << endl;
		Rational coordinates[2] = {Rational(1,2), Rational(2,3)};
		cout << "Has elements such as " << QF.Print(coordinates) << endl;

		for(int i=-10; i<10; i++){
			for(int j=-10; j<10; j++) {
				elements[0] = Rational(i,1);
				elements[1] = Rational(j,1);
				if(QF.FindSqrt(elements)) {
					cout << "Sqrt " << QF.Print(elements) << " = ";
					elements = QF.GetSqrt();
					cout << QF.Print(elements) << endl;
				}
				//else 
					//cout << QF.Print(elements) << " does not have a square root" << endl;
				
			}
		}
	}
	cout << "A pause. Enter anything.";
	cin >> p;
}