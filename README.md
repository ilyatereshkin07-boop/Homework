# Homework
#include <iostream>
#include <cmath>
#include <memory>
#include <stdexcept>

constexpr double PI = 3.14159265358979323846;
constexpr double EPS = 1e-9;


class Complex {
public:
    virtual ~Complex() = default;

   
    virtual double re() const = 0;   
    virtual double im() const = 0;  
    virtual double mod() const = 0;
    virtual double arg() const = 0;  


    static std::shared_ptr<Complex> make(double re, double im);

    
    virtual std::shared_ptr<Complex> clone() const = 0;

    
    virtual void print(std::ostream& os) const = 0;
};


class AlgebraicComplex : public Complex {
    double a_, b_;
public:
    AlgebraicComplex(double a = 0.0, double b = 0.0) : a_(a), b_(b) {}

    double re()  const override { return a_; }
    double im()  const override { return b_; }
    double mod() const override { return std::sqrt(a_ * a_ + b_ * b_); }
    double arg() const override { return std::atan2(b_, a_); }

    std::shared_ptr<Complex> clone() const override {
        return std::make_shared<AlgebraicComplex>(*this);
    }

    void print(std::ostream& os) const override {
        os << a_;
        if (b_ >= 0) os << " + " << b_ << "i";
        else         os << " - " << -b_ << "i";
    }
};


class TrigonometricComplex : public Complex {
    double r_, phi_;
public:
    TrigonometricComplex(double r = 0.0, double phi = 0.0) : r_(r), phi_(phi) {}

    double re()  const override { return r_ * std::cos(phi_); }
    double im()  const override { return r_ * std::sin(phi_); }
    double mod() const override { return r_; }
    double arg() const override { return phi_; }

    std::shared_ptr<Complex> clone() const override {
        return std::make_shared<TrigonometricComplex>(*this);
    }

    void print(std::ostream& os) const override {
        os << r_ << " * (cos(" << phi_ << ") + i*sin(" << phi_ << "))";
    }
};


std::shared_ptr<Complex> Complex::make(double re, double im) {
    return std::make_shared<AlgebraicComplex>(re, im);
}


std::shared_ptr<Complex> operator+(const Complex& x, const Complex& y) {
    return Complex::make(x.re() + y.re(), x.im() + y.im());
}
std::shared_ptr<Complex> operator-(const Complex& x, const Complex& y) {
    return Complex::make(x.re() - y.re(), x.im() - y.im());
}
std::shared_ptr<Complex> operator*(const Complex& x, const Complex& y) {
    return Complex::make(x.re() * y.re() - x.im() * y.im(),
        x.re() * y.im() + x.im() * y.re());
}
std::shared_ptr<Complex> operator/(const Complex& x, const Complex& y) {
    double denom = y.re() * y.re() + y.im() * y.im();
    if (denom < EPS) throw std::runtime_error("Division by zero");
    return Complex::make((x.re() * y.re() + x.im() * y.im()) / denom,
        (x.im() * y.re() - x.re() * y.im()) / denom);
}

bool operator==(const Complex& x, const Complex& y) {
    return std::abs(x.re() - y.re()) < EPS &&
        std::abs(x.im() - y.im()) < EPS;
}
bool operator!=(const Complex& x, const Complex& y) { return !(x == y); }

std::ostream& operator<<(std::ostream& os, const Complex& z) {
    z.print(os);
    return os;
}


std::shared_ptr<Complex> sqrt(const Complex& z) {
    double r = std::sqrt(z.mod());
    double phi = z.arg() / 2.0;
    return std::make_shared<TrigonometricComplex>(r, phi);
}


struct Roots {
    std::shared_ptr<Complex> x1;
    std::shared_ptr<Complex> x2;
};

Roots solveQuadratic(const Complex& a, const Complex& b, const Complex& c) {
    
    auto b2 = (*b.clone()) * (*b.clone());
    auto four_ac = Complex::make(4.0, 0.0) * a * c;
    auto D = *b2 - *four_ac;

    auto sqrtD = sqrt(*D);
    auto minusB = Complex::make(-b.re(), -b.im());
    auto twoA = Complex::make(2.0, 0.0) * a;

    Roots r;
    r.x1 = (*minusB + *sqrtD) / *twoA;
    r.x2 = (*minusB - *sqrtD) / *twoA;
    return r;
}


int main() {
    
    {
        AlgebraicComplex a(1, 0), b(-2, 0), c(5, 0);
        auto r = solveQuadratic(a, b, c);
        std::cout << "x^2 - 2x + 5 = 0\n";
        std::cout << "x1 = " << *r.x1 << "\n";
        std::cout << "x2 = " << *r.x2 << "\n\n";
    }

   
    {
        AlgebraicComplex a(1, 0), b(0, 2), c(-1, 0);  
        auto r = solveQuadratic(a, b, c);
        std::cout << "x^2 + 2i*x - 1 = 0\n";
        std::cout << "x1 = " << *r.x1 << "\n";
        std::cout << "x2 = " << *r.x2 << "\n\n";
    }

    
    {
        AlgebraicComplex z1(3, 4);              
        TrigonometricComplex z2(2, PI / 3);     
        std::cout << "z1 = " << z1 << "\n";
        std::cout << "z2 = " << z2 << "\n";
        std::cout << "z1 * z2 = " << *(z1 * z2) << "\n";
        std::cout << "|z1| = " << z1.mod() << "\n";
    }
    return 0;
}
