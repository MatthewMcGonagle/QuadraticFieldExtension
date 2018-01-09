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

## How to Use Classes

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
 
