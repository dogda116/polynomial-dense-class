# polynomial-dense-class
C++ polynomial class. All coefficients are saved in vector. 

**Features:**
  * Template type of coefficients
  * Constructs new polynomial from vector, coefficient (zero-degree polynomial), pair of iterators
  * Comparison operators == and !=
  * Arithmetic operators +, -, * and respective +=, -=, *=
  * Operator [] returns coefficient of this degree or -1 if it is zero polynomial 
  * Operator () evaluates a polynomial at specific values. Implementation using Horner's method for more effectiveness
  * Operator << prints polynomial in output stream
  * begin() and end() to iterate through coefficients from lowest to highest degree
  * Operator & returns a composition f(g(x))
  * Operators / and % return quotient and remainder respectively
  * Operator , returns PolynomialGCD (Greatest common divisor)
 
