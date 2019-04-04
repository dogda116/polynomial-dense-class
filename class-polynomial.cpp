#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

template <typename T>
class Polynomial {
private:
    std::vector<T> coef;

    void Delete_Front_Zeros() {
        while (!coef.empty() && coef.back() == T())
            coef.pop_back();
    }
public:
    Polynomial(const std::vector<T>& v) : coef(v) {
        Delete_Front_Zeros();
    }

    Polynomial(T c = T()) {
        if (c != T())
            coef.push_back(c);
    }

    template <typename Iter>
    Polynomial(Iter first, Iter last) {
        while (first != last) {
            coef.push_back(*first++);
        }
        Delete_Front_Zeros();
    }

    int Degree() const {  // returns degree of polynomial or -1 if it is zero polynomial
        for (int i = coef.size() - 1; i >= 0; --i) {
            if (coef[i] != T())
                return i;
        }
        return -1;
    }

    T operator[] (size_t degree) const {  // returns coefficient before this degree
        if (coef.size() > degree)
            return coef[degree];
        else
            return T();  // returns default value if requested degree is higher than max degree
    }

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

    bool operator != (const Polynomial<T>& other) const {
        return !(*this == other);
    }

    Polynomial<T> operator + (const Polynomial<T>& other) const {  // sum of two polynomials
        size_t pol_size = std::max(coef.size(), other.coef.size());
        std::vector<T> pol(pol_size, T());  // resizing with default values
        for (size_t i = 0; i != pol_size; ++i)
            pol[i] = (*this)[i] + other[i];
        return Polynomial<T> {pol};
    }

    Polynomial<T> operator - (const Polynomial<T>& other) const {  // difference of two polynomials
        size_t pol_size = std::max(coef.size(), other.coef.size());
        std::vector<T> pol(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            pol[i] = (*this)[i] - other[i];
        return Polynomial<T> {pol};
    }

    Polynomial<T> operator * (const Polynomial<T>& other) const {  // product of two polynomials
        size_t pol_size = coef.size() + other.coef.size();
        std::vector<T> pol(pol_size, T());
        for (size_t i = 0; i != coef.size(); ++i) {
            for (size_t j = 0; j != other.coef.size(); ++j)
                pol[i + j] += (*this)[i] * other[j];
        }
        return Polynomial<T> {pol};
    }

    Polynomial<T>& operator += (const Polynomial<T>& other) {
        size_t pol_size = std::max(coef.size(), other.coef.size());
        coef.resize(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            coef[i] += other[i];
        Delete_Front_Zeros();
        return *this;
    }

    Polynomial<T>& operator -= (const Polynomial<T>& other) {
        size_t pol_size = std::max(coef.size(), other.coef.size());
        coef.resize(pol_size, T());
        for (size_t i = 0; i != pol_size; ++i)
            coef[i] -= other[i];
        Delete_Front_Zeros();
        return *this;
    }

    Polynomial& operator *= (const Polynomial<T>& other) {
        Polynomial poly = *this * other;
        *this = poly;
        return *this;
    }

    T operator() (T value) const {  // calculates f(value)
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

    Polynomial<T> operator & (const Polynomial<T>& other) const {  // returns composition
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

    Polynomial<T> operator % (const Polynomial<T>& other) const {  // returns remainder
        return *this - (*this / other) * other;
    }

    Polynomial<T> operator , (const Polynomial<T>& other) const {  // returns PolynomialGCD
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

template <typename T>
std::ostream& operator << (std::ostream& out, const Polynomial<T>& pol) {
    int deg = pol.Degree();
    if (deg == -1) {
        out << '0';  // zero polynomial
    } else {
        for (int i = deg; i != -1; --i) {
            if (pol[i] != T(0)) {  // printing only degree with non-zero coefficient
                if (pol[i] != T(1) && pol[i] != T(-1)) {
                    if (i != deg && pol[i] > T(0)) {
                        out << "+";
                    }
                    out << pol[i];
                    if (i > 1)
                        out << "*x^" << i;
                    if (i == 1)
                        out << "*x";
                } else if (pol[i] == T(-1)) {  // coefficient -1 for non-zero degree prints as -
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
