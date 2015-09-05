//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"

using namespace std;

void main() {

	QuadraticField Q;
	Rational r, sum, product;

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
	}
	cout << "A pause. Enter anything.";
	cin >> p;
}