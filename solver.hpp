//
// Created by mangisto on 24/04/2020.
//

#pragma once
using namespace std;
#include <ostream>
#include <complex>
namespace solver
{
    class RealVariable
    {
    public:
        double _re;

        RealVariable(double re =0):_re(re){}

        friend RealVariable operator == (RealVariable a , RealVariable b);
        friend  RealVariable operator + ( const RealVariable a , const RealVariable b );
        friend  RealVariable operator - ( const RealVariable a , const  RealVariable b );
        friend  RealVariable operator * ( const RealVariable a , const RealVariable b );
        friend const RealVariable operator ^ ( const RealVariable a , const RealVariable b );
        friend const RealVariable operator / ( const RealVariable a , const RealVariable b );
        friend RealVariable operator + (int a , solver::RealVariable b );
        friend RealVariable operator * ( const int a , solver::RealVariable& b );
        friend RealVariable operator - ( int a , RealVariable b );
        friend RealVariable operator^(const RealVariable b,int a );
        friend RealVariable operator / ( int a , const RealVariable b );

        friend RealVariable operator + (double a , RealVariable b );
        friend RealVariable operator * ( double a , RealVariable b );
        friend RealVariable operator - ( double a , RealVariable b );
        friend RealVariable operator ^ ( double a , const RealVariable b );
        friend RealVariable operator / ( double a , const RealVariable b );

        friend ostream& operator<< (ostream& os, const RealVariable& c);

    };

    class ComplexVariable
    {
    public:


        double _re;
        double _im;
        ComplexVariable() = default;
        ComplexVariable(double re, double im) : _re(re), _im(im) {};
        ComplexVariable(int re) : _re(re), _im(0) {};
        ComplexVariable(complex<double> y): _re(y.real()),_im(y.imag()){};
        friend ComplexVariable operator == (ComplexVariable a , ComplexVariable b);
        friend const ComplexVariable operator + ( const ComplexVariable a , const ComplexVariable b );
        friend const ComplexVariable operator - ( const ComplexVariable a , const  ComplexVariable b );
        friend const ComplexVariable operator * ( const ComplexVariable a ,const ComplexVariable b );
        friend const ComplexVariable operator / (  const ComplexVariable a , const ComplexVariable b );
        friend ComplexVariable operator + (int a , ComplexVariable b );
        friend ComplexVariable operator * ( int a , ComplexVariable b );
        friend ComplexVariable operator - ( int a , ComplexVariable b );
        friend ComplexVariable operator ^ (ComplexVariable a , ComplexVariable b );
        friend ComplexVariable operator / ( int a ,  ComplexVariable b );

        friend ComplexVariable operator + (double a , ComplexVariable b );
        friend ComplexVariable operator * ( double a , ComplexVariable b );
        friend ComplexVariable operator - ( double a , ComplexVariable b );
        friend ComplexVariable operator^(ComplexVariable b,double a );
        friend ComplexVariable operator / ( double a ,  ComplexVariable b );

        friend ostream& operator<< (ostream& os, const ComplexVariable& c);
        friend ostream& operator<< (ostream& os, const std::complex<double>& c);

    };

    double  solve( RealVariable a);

    std::complex<double>  solve(ComplexVariable y);

};

