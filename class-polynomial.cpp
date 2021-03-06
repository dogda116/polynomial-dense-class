#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

template <typename T>
class Polynomial {
private:
    std::vector<T> coef;

    // delete leading non-significant zeros
    void Delete_Front_Zeros() {
        while (!coef.empty() && coef.back() == T())
            coef.pop_back();
    }
public:
    // initialize polynomial with vector of coefficients
    Polynomial(const std::vector<T>& v) : coef(v) {
        Delete_Front_Zeros();
    }

    // initialize constant polynomial
    Polynomial(T c = T()) {
        if (c != T())
            coef.push_back(c);
    }

    // initialize polynomial with the pair of iterators
    template <typename Iter>
    Polynomial(Iter first, Iter last) {
        while (first != last) {
            coef.push_back(*first++);
        }
        Delete_Front_Zeros();
    }
    
    // returns degree of polynomial or -1 if it is zero polynomial
    int Degree() const {
        for (int i = coef.size() - 1; i >= 0; --i) {
            if (coef[i] != T())
                return i;
        }
        return -1;
    }

    // returns coefficient before this degree
    T operator[] (size_t degree) const {
        if (coef.size() > degree)
            return coef[degree];
        else
            return T();  // returns default value if requested degree is higher than max degree
    }

    // returns true if two polynomials are equal, false otherwise
    bool operator == (const Polynomial<T>& other) const {
        int fd = Degree();
        int sd = other.Degree();
        if (fd != sd) {
            return false;
        } else {
            for (; fd >= 0; --fd) {
                if ((*this)[fd] != other[fd]) {
                    return false;
                }
            }
            return true;
        }
    }

    // returns true if two polynomials are not equal, false otherwise
    bool operator != (const Polynomial<T>& other) const {
        return !(*this == other);
    }

    // sum of two polynomials
    Polynomial<T> operator + (const Polynomial<T>& other) const {
        size_t pol_size = std::max(coef.size(), other.coef.size());
        std::vector<T> pol(pol_size, T());  // resizing with default values
        for (size_t i = 0; i != pol_size; ++i)
            pol[i] = (*this)[i] + other[i];
        return Polynomial<T> {pol};
    }

    // difference of two polynomials
    Polynomial<T> operator - (const Polynomial<T>& other) const {  
        size_t pol_size = std::max(coef.size(), other.coef.size());
        std::vector<T> pol(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            pol[i] = (*this)[i] - other[i];
        return Polynomial<T> {pol};
    }

    // product of two polynomials
    Polynomial<T> operator * (const Polynomial<T>& other) const {  
        size_t pol_size = coef.size() + other.coef.size();
        std::vector<T> pol(pol_size, T());
        for (size_t i = 0; i != coef.size(); ++i) {
            for (size_t j = 0; j != other.coef.size(); ++j)
                pol[i + j] += (*this)[i] * other[j];
        }
        return Polynomial<T> {pol};
    }

    // sum of two polynomials
    Polynomial<T>& operator += (const Polynomial<T>& other) {
        size_t pol_size = std::max(coef.size(), other.coef.size());
        coef.resize(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            coef[i] += other[i];
        Delete_Front_Zeros();
        return *this;
    }

    // difference of two polynomials
    Polynomial<T>& operator -= (const Polynomial<T>& other) {
        size_t pol_size = std::max(coef.size(), other.coef.size());
        coef.resize(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            coef[i] -= other[i];
        Delete_Front_Zeros();
        return *this;
    }

    // product of two polynomials
    Polynomial& operator *= (const Polynomial<T>& other) {
        Polynomial poly = *this * other;
        *this = poly;
        return *this;
    }

    // calculates f(value)
    T operator() (T value) const {  
        T ans = T();
        for (int i = coef.size() - 1; i >= 0; --i)
            ans = coef[i] + ans * value;  // Horner's method
        return ans;
    }

    typename std::vector<T>::const_iterator begin() const {
        return coef.begin();
    }

    typename std::vector<T>::const_iterator end() const {
        return coef.end();
    }

    // returns composition
    Polynomial<T> operator & (const Polynomial<T>& other) const {  
        Polynomial<T> composition;
        for (size_t i = 0; i != coef.size(); ++i) {
            if (coef[i] != T()) {
                Polynomial<T> tmp{T(coef[i])};
                for (size_t j = 1; j != i + 1; ++j)
                    tmp *= other;
                composition += tmp;
            }
        }
        return composition;
    }

    // divides one polynomial by another
    Polynomial<T> operator / (const Polynomial<T>& other) const {
        Polynomial<T> first = *this;
        std::vector<T> result(coef.size());
        int fd = first.Degree();
        while (fd >= other.Degree()) {
            T t = first[fd] / other[other.Degree()];
            std::vector<T> tmp;
            for (int i = 0; i != fd - other.Degree(); ++i)
                tmp.push_back(T());
            tmp.push_back(t);
            Polynomial<T> p {tmp};
            first -= other * p;
            result[fd - other.Degree()] += t;
            fd = first.Degree();
        }
        return Polynomial<T> {result};
    }

    // returns remainder
    Polynomial<T> operator % (const Polynomial<T>& other) const {  
        return *this - (*this / other) * other;
    }

    // returns PolynomialGCD (greatest common divisor)
    Polynomial<T> operator , (const Polynomial<T>& other) const {  
        Polynomial<T> first = *this, second = other;
        if (first.Degree() < second.Degree()) {
            std::swap(first, second);
        }
        while (second.Degree() > 0) {
            Polynomial<T> tmp = second;
            first = first % second;
            second = first;
            first = tmp;
        }
        if (second != T(0))
            return Polynomial(T(1));
        first = first / first[first.Degree()];
        return first;
    }
};

// overload "<<" operator to print polynomials as: std::cout << polynomial;
// example: x^3+2*x^2-x+3
template <typename T>
std::ostream& operator << (std::ostream& out, const Polynomial<T>& pol) {
    int deg = pol.Degree();
    if (deg == -1) {
        // zero polynomial
        out << '0';  
    } else {
        for (int i = deg; i != -1; --i) {
            if (pol[i] != T(0)) {
                // printing only degree with non-zero coefficient
                if (pol[i] != T(1) && pol[i] != T(-1)) {
                    if (i != deg && pol[i] > T(0)) {
                        out << "+";
                    }
                    out << pol[i];
                    if (i > 1)
                        out << "*x^" << i;
                    if (i == 1)
                        out << "*x";
                } else if (pol[i] == T(-1)) {
                    // coefficient -1 for non-zero degree prints as -
                    if (i > 1)
                        out << "-x^" << i;
                    if (i == 1)
                        out << "-x";
                    if (i == 0)
                        out << pol[i];
                } else if (pol[i] == T(1)) {
                    if (i != deg)
                        out << "+";
                    if (i > 1)
                        out << "x^" << i;
                    if (i == 1)
                        out << "x";
                    if (i == 0)
                        out << pol[i];
                }
            }
        }
    }
    return out;
}
