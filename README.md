# QuadraticFieldExtension

The purpose of this project is to create a C++ class library for dealing with field extensions of the rational numbers **Q** by the addition of square roots.

## Some Mathematical Explanation

Let us now take a moment to explain what is meant by field extensions of the rational numbers **Q**. First, the rational numbers **Q** consist of all fractions **p/q** where **p** and **q** are both integers. Furthermore, **q** is non-zero.

As a bit of notation let us use **Sqrt(x)** to be some choice of square root of the number **x**. Note that either choice differs by only a minus sign. Now, not every **Sqrt(p/q)** is a rational number. For example, **Sqrt(2)** is not a rational number.

We consider extending the rational numbers **Q** to include **Sqrt(2)** and all numbers that are generated from **Q** and **Sqrt(2)** by multiplication, division, addition, and subtraction. This is an object in algebra called a *Field* (not to be confused with a vector field), and we name it **Q(Sqrt(2))**. 

It is a mathematical fact that all of the numbers in **Q(Sqrt(2))** may be expressed as **f1 + f2 Sqrt(2)** where **f1, f2** are rational numbers. Here, we are using **f** to stand for "fraction". For example, this includes numbers such as **1 + 2 Sqrt(2)** and **1/2 + 2/3 Sqrt(2)**. An example of how division is included in **Q(Sqrt(2))** is given by **1/(1 + Sqrt(2)) = 1/(1 + Sqrt(2)) \* (1 - Sqrt(2))/(1 - Sqrt(2)) = -1 + Sqrt(2)**

Now, for ease of reading we make a notational change. Let **R1 = Sqrt(2)** denote our first square root. So now elements in **Q(R1) = Q(Sqrt(2))** are given by **f1 + f2 R1**.

Some square roots of elements in **Q(R1)** are suprisingly still in **Q(R1)**. For example. **Sqrt(3 + 2 R1) = 1 + R1**. However, similar to the rational case, not every **Sqrt(f1 + f2 R1)** is in the field **Q(R1)**. For example, **R2 = Sqrt( 1 + 2 R1)** is not in **Q(R1)**. 

So, we extend to a new field **Q(R1, R2)** generated by the rational numbers **Q**, **R1**, and **R2** under multiplication, division, addition, and subtraction. It is a mathematical fact that these numbers are all represented by **f1 + f2 R1 + f3 R2 + f4 R1R2**, where **f1, f2, f3, f4** are in **Q**. So, by extending using a square root, we have doubled the number of necessary fractions to represent our numbers.

Now, in this new field **Q(R1, R2)** we can again consider trying to find square roots of elements in **Q(R1,R2)** that are themselves also in **Q(R1, R2)**. For example, **Sqrt(14 + 10 R1 + 6 R2 + 4 R1R2) = 1 + 1 R1 + R2 + R1R2**. However, just like before, not every **Sqrt(f1 + f2 R1 + f3 R2 + f4 R1R2)** is in **Q(R1, R2)**. For example, **R3 = Sqrt(1 + R1 + R2 + R1R2)** is NOT in **Q(R1, R2)**.

We may then create another extension **Q(R1, R2, R3)** with elements represented using 8 fractions **f1, f2, f3, f4, f5, f6, f7, f8** by **f1 + f2 R1 + f3 R2 + f4 R1R2 + f5 R3 + f6 R1R3 + f7 R2R3 + f8 R1R2R3**. Furthermore, we can keep repeating this process for more square roots that aren't expressible in the current field extensions to create **Q(R1, R2, ..., RN)**. Elements in **Q(R1, R2, ..., RN)** are represented using **2<sup>N</sup>** fractions **(f1, f2, ...., f2<sup>N</sup>)**.

## Explanation of the Class Library

The purpose of the class library is to implement this mathematical process explicitly in C++. The first class **``Rational``** implements the rational numbers **Q** and comes with a method **``bool Rational::FindSqrt()``** to determine if a particular instance of **``Rational``** has a square root in **Q**.

The class **``CoordinateChunk``** represents segments of the **``Rational``** coordinates **(f1, f2, ..., f2<sup>N</sup>)**. The purpose of the **``CoordinateChunk``** class is to facilitate in calculations in field extensions by referring to parts of a **``std::vector<Rational>``** using references. This allows one to do these calculations without the need to copy values of **``std::vector<Rational>``**.

Instances of the class **``QuadraticField``** are extensions of other field by square roots. Instances of **``QuadraticField``** are used to build the tower of fields **Q(R1), Q(R1,R2), Q(R1, R2, R3), ...,** and **Q(R1, R2, ..., RN)**. The function **``bool QuadraticField::FindSqrt(std::vector<Rational> coordinates_element)``** computes if the element of field extension represented by the instance **``*this``** is in the field extension represented by **``*this``**. If not, then a new extension may be instancialized using the constructor **``QuadraticField::QuadraticField(QuadraticField* basefield, std::vector<Rational> coordinates_root)``**. 