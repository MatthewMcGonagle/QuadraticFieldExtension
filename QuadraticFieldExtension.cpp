//Author: Matthew McGonagle

#include <iostream>
#include "QuadraticField.h"

using namespace std;

void main() {

	QuadraticField Q;
	Rational r;

	int p = 5, q=2;
	cout << "Please input integer p = " << endl;
	cin >> p;
	cout << "Please input integer q = " << endl;
	cin >> q;
	r = Rational(p,q);
	
	if(Q.FindSqrt(r))
		cout << "Square root of "<< r.print()
			<< " is = " << Q.GetSqrtResult().print() << endl;  
	else
		cout << "There is no rational square root of " << r.print() << endl;
	cout << "A pause. Enter anything.";
	cin >> p;
}