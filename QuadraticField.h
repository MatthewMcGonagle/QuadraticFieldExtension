//Author: Matthew McGonagle
#include <string>
#include <sstream>
#include <vector>

class Rational {
public:
	Rational(){num = 1; den = 1; name = std::string();};
	Rational(int p, int q);

	int GetP() const {return num;};
	int GetQ() const {return den;};
	std::string print();
	const Rational operator*(const Rational &other) const;
	const Rational operator+(const Rational &other) const;
	const Rational operator-(const Rational &other) const;

	Rational Inverse();

	bool IsZero();
	bool FindSqrt();
	int GetSqrtP(){return sqrtp;}
	int GetSqrtQ(){return sqrtq;}
	Rational GetSqrt(){return Rational(sqrtp, sqrtq);}

private:

	unsigned int gcd(int p, int q) const;
	unsigned int abs(int p) const;

	int num, den;
	int sqrtp, sqrtq;
	std::string name;
};


//////////////////////////////////

class CoordinateChunk {
public:
	CoordinateChunk();
	CoordinateChunk(std::vector<Rational>* coords_, int offset_, int size_);
	Rational& Get(int i) {return (*coords)[offset+i];}

	std::vector<Rational>* coords;
	int offset, size;
};

////////////////////////////////////////


class QuadraticField {

public:
	QuadraticField();
	QuadraticField(Rational root);
	QuadraticField(QuadraticField* basefield_, std::vector<Rational> root_);
	~QuadraticField();
	bool FindSqrt(std::vector<Rational> root_);
	std::string Print(std::vector<Rational> coordinates);
	std::string Print(CoordinateChunk x);

	int GetDegree() {return degree;}
	std::vector<Rational> GetSqrt() {return Result;}

	void ZeroResult();
	void SetResult(std::vector<Rational>& Result_);
	QuadraticField operator+(std::vector<Rational> &other);
	void Sum(std::vector<Rational> &a, std::vector<Rational> &b);
	void Product(std::vector<Rational> &a, std::vector<Rational> &b);
	void Product(CoordinateChunk a, CoordinateChunk b, CoordinateChunk Result_);
	//QuadraticField operator*(Coordinates& other);

	std::vector<Rational> Result, Scratch;

private:
	int degree;
	std::vector<Rational> Root;
	QuadraticField* basefield;
	std::string name;
};