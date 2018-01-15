---
layout: default
title: Quadratic Field Extension
---

The purpose of this project is a small library for doing infinite precision arithmetic for numbers made 
up of fractions and any number of square roots. 
By infinite precision, we mean that there is no rounding error. 
For those familiar enough with abstract algebra, they will recognize this as an implementation of 
quadratic field extensions of the fractions.

What do we mean by any number of square roots? First as a matter of notation, let <code> sqrt() </code> be the square root function. 
We mean numbers such as:

* `2/3`, 
* `1 + sqrt(2)`, 
* `1 + 2/3 * sqrt(2) + 7/5 * sqrt(1/2 + 4/5 * sqrt(2))`.

Let us next review how to use the classes in <code>QuadraticField.cpp</code> allow us to do 
arithmetic with such numbers.

## Adding a Square Root 

The operations are handled by the class `QuadraticFieldTower`.
First, let's create an instance using: 
```cpp
QuadraticFieldTower tower;
```

When first instantiated, an instance can only deal with fractions. 
We need to start by giving it a square root of a fraction to add to the list of numbers it can handle. 
The class will only add the square root if it can't be expressed in terms of fractions (or integers). 
For example, `sqrt(4) == 2`, so the class won't do anything if you try to add `sqrt(4)`. 
However, it is a mathematical fact that `sqrt(2)` is NOT a fraction, so we may add it.
So how would we add `sqrt(2)`?
 
First, we construct a vector with just one `Rational` element representing 2.
```cpp
std::vector<Rational> square1(1, Rational(2,1));
```

Then, we add it to `tower`.
```cpp
tower.addIfNoSqrRoot(square1);
```
`QuadraticFieldTower::addIfNoSqrRoot` will first test if the number has a square root existing in the current number system. If not, then it will extend the number system to handle the new square root.

Now, `tower` can handle numbers such as `1/2 + 2/3 * sqrt(2)`. Such numbers are represented using vectors of length 2. For example, `1/2 + 2/3 * sqrt(2)` is represented using a vector containing `{Rational(1,2), Rational(2,3)}`. For any such vector, we can see which number it represents using the `printCoords` function. Let's try this out.
```cpp
// Array holding coordinates of number.
Rational coordArray1[] = {Rational(1,2), Rational(2,3)}; 
// Put coordinates in a vector.
std::vector<Rational> coord1(coordArray1, coordArray1 + sizeof(coordArray1) / sizeof(Rational)); 

std::cout << printCoords(coord1);
```

We get
```
1/2 + 2/3r0
```
Here `r0` represents the first square root added to `tower`. To see a list of the square roots added to `tower`, we can just run 
```cpp
std::cout << tower.print();
```
We get
```
r0 = sqrt( 2 )
```

## Doing Some Operations

Now let's do some basic operations with our numbers. Let's add `1 + sqrt(2)` and `2 + 3 * sqrt(2)`. 
```cpp
// First set up the coordinates of the numbers 1 + sqrt(2), 2 + 3 * sqrt(2), and coordinates of result.
Rational coordArray1 = {Rational(1,1), Rational(1,1)},
         coordArray2 = {Rational(2,1), Rational(3,1)};
std::vector<Rational> coord1(coordArray1, coordArray1 + sizeof(coordArray1) / sizeof(Rational)),
                      coord2(coordArray2, coordArray2 + sizeof(coordArray2) / sizeof(Rational)),
                      result(2);

// Now use tower to add.
result = tower.add(coord1, coord2);
std::cout << printCoords(sum);
```

We get
```
3 + 4r0
```

Now, let's try multiplication.
```
result = tower.multiply(coord1, coord2);
std::cout << printCoords(result);
```
We get
```
8 + 5r0
```
which is correct.

Now, let's look at finding a square root.
``` cpp
Rational coordArray3[] = {Rational(22,1), Rational(12,1)};
std::vector<Rational> coord3(coordArray3, coordArray3 + sizeof(coordArray3) / sizeof(Rational)),
                      result(2);
bool sqrtExists;

std::cout << "coord3 = " << printCoords(coord3) << std::endl;
sqrtExists = tower.hasSqrt(coord3, result);
std::cout << "sqrtExists = " << sqrtExists << std::endl
          << "square root = " << printCoords(result) << std::endl;
```
We get
```
coord3 = 22 + 12r0
sqrtExists = 1
square root = 2 + 3r0 
```
Which is correct, because the square of `2 + 3 * sqrt(2)` is `22 + 12 * sqrt(2)`. Also, note that `bool` return value of `QuadraticFieldTower::hasSqrt` indicate whether the square root exists in the current number system.
It is not always true such a square root exists. For an example, see the next section.

## Adding More Square Roots

As we said before, not every number in our current system has a square root in the current system, i.e. of the form `fraction + fraction * sqrt(2)`. For example the square root of `1 + sqrt(2)` can not be found in the current number system. However, the number system can be extended again.
```cpp
Rational newRootArray[] = {Rational(1,1), Rational(2,1)};
std::vector<Rational> newRoot(newRootArray, newRootArray + sizeof(newRootArray) / sizeof(Rational)};

tower.addIfNoSqrRoot(newRoot);
cout << tower.print();
```
Now we have that the root list in `tower` is
```
r0 = sqrt( 2 )
r1 = sqrt( 1 + 1r0 )
```
Note that the function `QuadraticFieldTower::addIfNowSqrRoot` first tests to see if the square root exists in the current number system and only extends the number system if it does NOT exist. 

Now, numbers represented in our system require four fractions.
``` cpp
Rational coord4Array[] = {Rational(1,1), Rational(2,1), Rational(3,1), Rational(4,1)};
std::vector<Rational> coord4(coord4Array, coord4Array + sizeof(coord4Array) / sizeof(Rational));

std::cout << printCoords(coord4);
```
We get
```
1 + 2r0 + 3r1 + 4r0r1
```
We can similarly do operations for these new numbers.

Note that everytime we extend the number system by adding a square root, we double the number of coordinates necessary.
