//
// Created by mangisto on 27/04/2020.
//
/**
 * A demonstration test program for Gimatriya calculation.
 *
 * Use this as a basis to build a more complete Test.cpp file.
 *
 * Date: 2020-02
 */

#include "doctest.h"
#include <iostream>
#include <complex>
#include "solver.hpp"

using namespace std;
using solver::solve, solver::RealVariable, solver::ComplexVariable;

TEST_CASE("RealVariable Test case") {
    RealVariable x;
    ComplexVariable y;
//

            CHECK(solve(2*x-4 == 10)  == 7);
            CHECK(solve(x-4 == 10)  == 14);
            CHECK(solve(2*x == 10)  == 5);
            CHECK_THROWS(solve((x^2) == -16));
            CHECK_THROWS(solve((x^3) == 0));
            CHECK_THROWS(solve((x^-1) == 0));
            CHECK_THROWS(solve((x^4) == 0));
            CHECK_THROWS(solve((x^(0.5)) == 0));
            CHECK_THROWS(solve((x^-0.5) == 0));
            CHECK(solve((x^0) == 1));
            CHECK(solve((x^1) == 1));
            CHECK_UNARY(solve(6*(x^2)+11*x-35 == 0)==5/3||solve(6*(x^2)+11*x-35 == 0)==-3.5);
            CHECK_UNARY(solve(2*(x^2)-4*x-2 == 0)==1-sqrt(2)||solve(2*(x^2)-4*x-2 == 0)==1+sqrt(2));
            CHECK_UNARY(solve(-4*(x^2)-7*x==-12)==(-7-sqrt(241))/8||solve(-4*(x^2)-7*x == -12)==(-7+sqrt(241))/8);
            CHECK_UNARY(solve(20*(x^2)-15*x== 10)==(3+sqrt(41))/8||solve(20*(x^2)-15*x== 10)==(3-sqrt(41))/8);
            CHECK_UNARY(solve((x^2)-x==3)==(1+sqrt(13))/2||solve((x^2)-x==3)==(1-sqrt(13))/2);
            CHECK_UNARY(solve(5*(x^2)-5*x==9)==(5+sqrt(205))/10||solve(5*(x^2)-5*x==9)==(5-sqrt(205))/10);
            CHECK_UNARY(solve(3*(y^2)+4*y==-2)==-0.6666i||solve(3*(y^2)+4*y==-2)==+0.4714i);
            CHECK_UNARY(solve(-1*(x^2)+6*x==-18)==3+3*sqrt(3) || solve(-1*(x^2)+6*x==-18)==3-3*sqrt(3));
            CHECK_UNARY(solve(2*(x^2)==64)==4*sqrt(2) || solve(2*(x^2)==64)==-4*sqrt(2));
            CHECK_UNARY(solve(9*(y^2)==-49)==2.3333i || solve(9*(y^2)==-49)==-2.3333i);
            CHECK_UNARY(solve(-2*(y^2)==4)==1.4142i || solve(-2*(y^2)==4)==-1.4142i);
            CHECK_UNARY(solve(4*(x^2)==-81)==9/2 || solve(4*(x^2)==-81)==-9/2);
            CHECK_UNARY(solve(-1*(y^2)==9)==3i || solve(-1*(y^2)==9)==-3i);
            CHECK_UNARY(solve(3*(x^2)==36)==2*sqrt(3) || solve(3*(x^2)==36)==-2*sqrt(3));
            CHECK_UNARY(solve(6*(x^2)==-144)==2*sqrt(6) || solve(6*(x^2)==-144)==-2*sqrt(6));
            CHECK_UNARY(solve((x^2)-7*x==0)==7/2 || solve((x^2)-7*x==0)==0);
            CHECK_UNARY(solve(2*(x^2)+8*x==0)==-4 || solve(2*(x^2)+8*x==0)==0);
            CHECK_UNARY(solve(-1*(y^2)-9*y==0)==3i || solve(-1*(y^2)-9*y==0)==-3i);
            CHECK_UNARY(solve((x^2)+2*x==0)==-2 || solve((x^2)+2*x==0)==0);
            CHECK_UNARY(solve(-6*(x^2)-3*x==0)==-1/2 || solve(-6*(x^2)-3*x==0)==0);
            CHECK_UNARY(solve(-5*(x^2)+x==0)==sqrt(5)/5 || solve(-5*(x^2)+x==0)==-sqrt(5)/5);
            CHECK_UNARY(solve(-12*(x^2)+13*x==0)==13/12 || solve(-12*(x^2)+13*x==0)==0);
            CHECK_UNARY(solve(11*(x^2)-27*x==0)==27/11 || solve(11*(x^2)-27*x==0)==0);

            CHECK(solve(2*x-4 == 10)  == 7);
            CHECK(solve(x-4 == 10)  == 14);
            CHECK(solve(2*x == 10)  == 5);
            CHECK_THROWS(solve((x^2) == -16));
            CHECK_THROWS(solve((x^3) == 0));
            CHECK_THROWS(solve((x^-1) == 0));
            CHECK_THROWS(solve((x^4) == 0));
            CHECK_THROWS(solve((x^(0.5)) == 0));
            CHECK_THROWS(solve((x^-0.5) == 0));
            CHECK(solve((x^0) == 1));
            CHECK(solve((x^1) == 1));
            CHECK_UNARY(solve(6*(x^2)+11*x-35 == 0)==5/3||solve(6*(x^2)+11*x-35 == 0)==-3.5);
            CHECK_UNARY(solve(2*(x^2)-4*x-2 == 0)==1-sqrt(2)||solve(2*(x^2)-4*x-2 == 0)==1+sqrt(2));
            CHECK_UNARY(solve(-4*(x^2)-7*x==-12)==(-7-sqrt(241))/8||solve(-4*(x^2)-7*x == -12)==(-7+sqrt(241))/8);
            CHECK_UNARY(solve(20*(x^2)-15*x== 10)==(3+sqrt(41))/8||solve(20*(x^2)-15*x== 10)==(3-sqrt(41))/8);
            CHECK_UNARY(solve((x^2)-x==3)==(1+sqrt(13))/2||solve((x^2)-x==3)==(1-sqrt(13))/2);
            CHECK_UNARY(solve(5*(x^2)-5*x==9)==(5+sqrt(205))/10||solve(5*(x^2)-5*x==9)==(5-sqrt(205))/10);
            CHECK_UNARY(solve(3*(y^2)+4*y==-2)==-0.6666i||solve(3*(y^2)+4*y==-2)==+0.4714i);
            CHECK_UNARY(solve(-1*(x^2)+6*x==-18)==3+3*sqrt(3) || solve(-1*(x^2)+6*x==-18)==3-3*sqrt(3));
            CHECK_UNARY(solve(2*(x^2)==64)==4*sqrt(2) || solve(2*(x^2)==64)==-4*sqrt(2));
            CHECK_UNARY(solve(9*(y^2)==-49)==2.3333i || solve(9*(y^2)==-49)==-2.3333i);
            CHECK_UNARY(solve(-2*(y^2)==4)==1.4142i || solve(-2*(y^2)==4)==-1.4142i);
            CHECK_UNARY(solve(4*(x^2)==-81)==9/2 || solve(4*(x^2)==-81)==-9/2);
            CHECK_UNARY(solve(-1*(y^2)==9)==3i || solve(-1*(y^2)==9)==-3i);
            CHECK_UNARY(solve(3*(x^2)==36)==2*sqrt(3) || solve(3*(x^2)==36)==-2*sqrt(3));
            CHECK_UNARY(solve(6*(x^2)==-144)==2*sqrt(6) || solve(6*(x^2)==-144)==-2*sqrt(6));
            CHECK_UNARY(solve((x^2)-7*x==0)==7/2 || solve((x^2)-7*x==0)==0);
            CHECK_UNARY(solve(2*(x^2)+8*x==0)==-4 || solve(2*(x^2)+8*x==0)==0);
            CHECK_UNARY(solve(-1*(y^2)-9*y==0)==3i || solve(-1*(y^2)-9*y==0)==-3i);
            CHECK_UNARY(solve((x^2)+2*x==0)==-2 || solve((x^2)+2*x==0)==0);
            CHECK_UNARY(solve(-6*(x^2)-3*x==0)==-1/2 || solve(-6*(x^2)-3*x==0)==0);
            CHECK_UNARY(solve(-5*(x^2)+x==0)==sqrt(5)/5 || solve(-5*(x^2)+x==0)==-sqrt(5)/5);
            CHECK_UNARY(solve(-12*(x^2)+13*x==0)==13/12 || solve(-12*(x^2)+13*x==0)==0);
            CHECK_UNARY(solve(11*(x^2)-27*x==0)==27/11 || solve(11*(x^2)-27*x==0)==0);
            CHECK(solve(2*x-4 == 10)  == 7);
            CHECK(solve(x-4 == 10)  == 14);
            CHECK(solve(2*x == 10)  == 5);
            CHECK_THROWS(solve((x^2) == -16));
            CHECK_THROWS(solve((x^3) == 0));
            CHECK_THROWS(solve((x^-1) == 0));
            CHECK_THROWS(solve((x^4) == 0));
            CHECK_THROWS(solve((x^(0.5)) == 0));
            CHECK_THROWS(solve((x^-0.5) == 0));
            CHECK(solve((x^0) == 1));
            CHECK(solve((x^1) == 1));
            CHECK_UNARY(solve(6*(x^2)+11*x-35 == 0)==5/3||solve(6*(x^2)+11*x-35 == 0)==-3.5);
            CHECK_UNARY(solve(2*(x^2)-4*x-2 == 0)==1-sqrt(2)||solve(2*(x^2)-4*x-2 == 0)==1+sqrt(2));
            CHECK_UNARY(solve(-4*(x^2)-7*x==-12)==(-7-sqrt(241))/8||solve(-4*(x^2)-7*x == -12)==(-7+sqrt(241))/8);
            CHECK_UNARY(solve(20*(x^2)-15*x== 10)==(3+sqrt(41))/8||solve(20*(x^2)-15*x== 10)==(3-sqrt(41))/8);
            CHECK_UNARY(solve((x^2)-x==3)==(1+sqrt(13))/2||solve((x^2)-x==3)==(1-sqrt(13))/2);
            CHECK_UNARY(solve(5*(x^2)-5*x==9)==(5+sqrt(205))/10||solve(5*(x^2)-5*x==9)==(5-sqrt(205))/10);
            CHECK_UNARY(solve(3*(y^2)+4*y==-2)==-0.6666i||solve(3*(y^2)+4*y==-2)==+0.4714i);
            CHECK_UNARY(solve(-1*(x^2)+6*x==-18)==3+3*sqrt(3) || solve(-1*(x^2)+6*x==-18)==3-3*sqrt(3));
            CHECK_UNARY(solve(2*(x^2)==64)==4*sqrt(2) || solve(2*(x^2)==64)==-4*sqrt(2));
            CHECK_UNARY(solve(9*(y^2)==-49)==2.3333i || solve(9*(y^2)==-49)==-2.3333i);
            CHECK_UNARY(solve(-2*(y^2)==4)==1.4142i || solve(-2*(y^2)==4)==-1.4142i);
            CHECK_UNARY(solve(4*(x^2)==-81)==9/2 || solve(4*(x^2)==-81)==-9/2);
            CHECK_UNARY(solve(-1*(y^2)==9)==3i || solve(-1*(y^2)==9)==-3i);
            CHECK_UNARY(solve(3*(x^2)==36)==2*sqrt(3) || solve(3*(x^2)==36)==-2*sqrt(3));
            CHECK_UNARY(solve(6*(x^2)==-144)==2*sqrt(6) || solve(6*(x^2)==-144)==-2*sqrt(6));
            CHECK_UNARY(solve((x^2)-7*x==0)==7/2 || solve((x^2)-7*x==0)==0);
            CHECK_UNARY(solve(2*(x^2)+8*x==0)==-4 || solve(2*(x^2)+8*x==0)==0);
            CHECK_UNARY(solve(-1*(y^2)-9*y==0)==3i || solve(-1*(y^2)-9*y==0)==-3i);
            CHECK_UNARY(solve((x^2)+2*x==0)==-2 || solve((x^2)+2*x==0)==0);
            CHECK_UNARY(solve(-6*(x^2)-3*x==0)==-1/2 || solve(-6*(x^2)-3*x==0)==0);
            CHECK_UNARY(solve(-5*(x^2)+x==0)==sqrt(5)/5 || solve(-5*(x^2)+x==0)==-sqrt(5)/5);
            CHECK_UNARY(solve(-12*(x^2)+13*x==0)==13/12 || solve(-12*(x^2)+13*x==0)==0);
            CHECK_UNARY(solve(11*(x^2)-27*x==0)==27/11 || solve(11*(x^2)-27*x==0)==0);

}
TEST_CASE("ComplexVariable Test case ") {
    ComplexVariable a(3,4);
    ComplexVariable b(3,4);

            CHECK(a._re== b._re);
            CHECK(a._im== b._im);



            RealVariable x(1);

            CHECK((2*x+1)._re== 3);
}
//© 2020 GitHub,