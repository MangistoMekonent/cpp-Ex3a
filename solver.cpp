//
// Created by mangisto on 24/04/2020.
//

using namespace std;
#include <iostream>
#include <complex>
#include "solver.hpp"
#include <cmath>
#include <bits/exception.h>
namespace solver {

    RealVariable operator==(RealVariable a, RealVariable b) {
        return a;
    };

    RealVariable operator+(const RealVariable a, const RealVariable b) {
        RealVariable z(b._re + a._re);
        return z;
    }

    RealVariable operator-(const RealVariable a, const RealVariable b) {
        RealVariable z(a._re - b._re);
        return z;
    };

    RealVariable operator*(const RealVariable a, const RealVariable b) {
        RealVariable z(a._re * b._re);
        return z;
    }

    const RealVariable operator/(const RealVariable a, const RealVariable b) {
        if (b._re != 0) {
            RealVariable z(a._re / b._re);
            return z;
        } else
            throw std::exception();

    }

    RealVariable operator+(int a, RealVariable b) {
        RealVariable z(b._re + a);
        return z;
    };

    RealVariable operator-(int a, RealVariable b) {
        RealVariable z(b._re - a);
        return z;
    }

    RealVariable operator*(const int a, RealVariable &b) {
        RealVariable z(b._re * a);

        return z;
    }

    RealVariable operator^(const RealVariable b,int a ) {
        RealVariable z(pow(b._re, a));

        return z;
    }

    RealVariable operator/(int a, const RealVariable b) {
        if (b._re != 0) {
            RealVariable z(a / b._re);
            return z;
        } else
            throw std::exception();
    }


    RealVariable operator+(double a, RealVariable b) {
        RealVariable z(b._re + a);
        return z;
    };

    RealVariable operator-(double a, RealVariable b) {
        RealVariable z(b._re - a);
        return z;
    }

    RealVariable operator*(double a, RealVariable b) {
        RealVariable z(b._re * a);

        return z;
    }

    RealVariable operator^(double a, const RealVariable b) {
        RealVariable z(pow(b._re, a));

        return z;
    }

    RealVariable operator/(double a, const RealVariable b) {
        if (b._re != 0) {
            RealVariable z(a / b._re);
            return z;
        } else
            throw std::exception();
    }

    ostream &operator<<(ostream &os, const RealVariable &c) {
        return (os << c._re);
    }

//*****************************

    ComplexVariable operator==(ComplexVariable a, ComplexVariable b) {
        return a;
    };


    const ComplexVariable operator+(const ComplexVariable a, const ComplexVariable b) {
        ComplexVariable z(b._re + a._re, b._im + a._im);
        return z;
    }

    const ComplexVariable operator-(const ComplexVariable a, const ComplexVariable b) {
       ComplexVariable z(a._re - b._re, a._im - b._im);
        return z;
    };

    const ComplexVariable operator*(const ComplexVariable a, const ComplexVariable b) {
        ComplexVariable z((a._re * b._re) - (a._im * b._im), (a._re * b._im) + (b._re * a._im));
        return z;
    }


    const ComplexVariable operator/(const ComplexVariable a, const ComplexVariable b) {
        if (b._re != 0) {
            ComplexVariable z((a._re * b._re + a._im * b._im) / (pow(b._re, 2) + pow(b._im, 2)),
                              (a._im * b._re - a._re * b._im) / (pow(b._re, 2) + pow(b._im, 2)));
            return z;
        } else
            throw std::exception();

    }

    ComplexVariable operator+(int a, ComplexVariable b) {
        ComplexVariable z(b._re + a, b._im);
        return z;
    };

    ComplexVariable operator-(int a, ComplexVariable b) {
        ComplexVariable z(b._re - a, b._im);
        return z;
    }

    ComplexVariable operator*(int a, ComplexVariable b) {
        ComplexVariable z(b._re * a, b._im * a);

        return z;
    }

    ComplexVariable operator^(ComplexVariable a, ComplexVariable b) {
        if (a._re == 0) {
            ComplexVariable z(1, 0);
            return z;
        } else if (a._re == 1)
            return b;
        else if (a._re == 2) {
            ComplexVariable z((b._re * b._re) - (b._im * b._im), (b._re * b._im) + (b._re * b._im));
            return z;
        } else
            throw std::exception();
    }

    ComplexVariable operator/(int a, ComplexVariable b) {
        if (b._re != 0) {
            ComplexVariable z((a * b._re + 0 * b._im) / (pow(b._re, 2) + pow(b._im, 2)),
                              (0 * b._re - a * b._im) / (pow(b._re, 2) + pow(b._im, 2)));
            return z;
        } else
            throw std::exception();
    }


    ComplexVariable operator+(double a, ComplexVariable b) {
        ComplexVariable z(b._re + a, b._im);
        return z;
    };

    ComplexVariable operator-(double a, ComplexVariable b) {
        ComplexVariable z(b._re - a, b._im);
        return z;
    }

    ComplexVariable operator*(double a, ComplexVariable b) {
        ComplexVariable z(b._re * a, b._im * a);

        return z;
    }

    ComplexVariable operator^(ComplexVariable b,double a ) {
        if (a == 0) {
            ComplexVariable z(1, 0);
            return z;
        } else if (a == 1)
            return b;
        else if (a == 2) {
            ComplexVariable z((b._re * b._re) - (b._im * b._im), (b._re * b._im) + (b._re * b._im));
            return z;
        } else
            throw std::exception();
    }

    ComplexVariable operator/(double a, ComplexVariable b) {
        if (b._re != 0) {
            ComplexVariable z((a * b._re + 0 * b._im) / (pow(b._re, 2) + pow(b._im, 2)),
                              (0 * b._re - a * b._im) / (pow(b._re, 2) + pow(b._im, 2)));
            return z;
        } else
            throw std::exception();
    }

    ostream &operator<<(ostream &os, const ComplexVariable &c) {
        return (os << c._re << '+' << c._im << 'i');
    }

    ostream &operator<<(ostream &os, const std::complex<double> &c) {
        return (os);
    }

    double solve(RealVariable a) {
        return a._re;
    }

    std::complex<double> solve(ComplexVariable y) {
        std::complex<double> x;
        return x;
    }

}

