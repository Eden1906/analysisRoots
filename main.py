import decimal
import math
from sympy.utilities.lambdify import lambdify
import sympy as sp
from mpmath import ln


#Eden levy 208376095 and Alex perednya 321223075

def float_range(a,b,section):
    while a<b:
        yield float(a)
        a += decimal.Decimal(section)


def bisctionMethod(func,a,b):
    eps = 0.0000000001
    fa = func(a)
    fb = func(b)
    counter = 0
    error = (-(ln((eps) / (b-a)) / ln(2)))
    if fa * fb > 0:
        return None
    while (b-a) > epsilon:
        c = (a+b) / 2
        fc = func(c)
        if fc == 0:
            print("Root:")
            return c
        if fa*fc > 0:
            a,fa = c,fc
        if fb*fc > 0 and counter<error:
            b,fb = c,fc
            counter = counter + 1
            print("itr: ",counter, "a: ",a,"b:",b,"f(a):",fa,"f(b):",fb,"f(c):",fc)
    if abs(f(c)) < 0.000001:
            print("Root:")
            return c
    print("Position c in the original function does not reset it, so there is no root!")


def newtonRaphson(f,a,b):
    eps = 0.00000001
    itr = 0
    xr = (a+b)/2
    next_xr = b
    error = (-(ln((eps) / (b-a)) / ln(2)))
    while abs(next_xr-xr)>eps and itr < error:
        if f(xr) == 0:
            print("cant divided by zero")
            return
        temp = next_xr
        next_xr = xr - (f(xr)/f_prime(xr))
        xr = temp
        itr = itr + 1
        print("itr:",itr,"xr:",xr,"f(x):",f(xr),"f'(x):",f_prime(xr))
    print("root:",xr)
    itr = 0

def secantMethod(f,a,b):
    eps = 0.0000000001
    itr = 0
    prev_xr = (a+b)/2
    xr = a + 0.1
    next_xr = xr
    error = (-(ln((eps) / (b-a)) / ln(2)))
    while abs(xr - prev_xr) > eps and itr < error:
        if f(xr) - f(prev_xr) == 0:
            print("cant divided by zero")
            return
        temp = next_xr
        next_xr = (prev_xr*f(xr)-xr*f(prev_xr)) / (f(xr)-f(prev_xr))
        xr = next_xr
        prev_xr = temp
        itr+=1
        print("itr:",itr,"xr:",next_xr)
    print("root:",next_xr)
    itr = 0



def program(f,a,b):
    print("choose the method\n for bisection press 1\n for newton-raphson press 2\n for secant press 3")
    flag = input()
    if flag == '1':
        iterationRange = 0.1 #sections
        iterationRange = float_range(a,b,iterationRange)
        if(f(0) == 0):
            print("0 is solution,therefor\nRoot:\n 0")
        for i in iterationRange:
            if f(a) * f(a+0.1) < 0:
                print(bisctionMethod(f,i,i+0.1))
            if f_prime(a)*f_prime(a+0.1)<0:
                print(bisctionMethod(f_prime,i,i+0.1))
            a = a + 0.1
    elif flag == '2':
        iterationRange1 = 0.1
        iterationRange1 = float_range(a,b,iterationRange1)
        for i in iterationRange1:
            if f(a) * f(a+0.1) < 0:
                newtonRaphson(f,i,i+0.1)
            a = a + 0.1
    elif flag == '3':
        iterationRange2 = 0.1
        iterationRange2 = float_range(a,b,iterationRange2)
        for i in iterationRange2:
            if f(a) * f(a+0.1) <0:
                secantMethod(f,i,i+0.1)
            a = a + 0.1
    else:
        program(f,a,b)


# ------------------------------------

x = sp.symbols('x')

#accuracy of result
epsilon = 0.0000000001

#Our function(f_prime is our derivative function)
f = x**4 + x**3 - 3*x**2
f_prime = f.diff(x)
f = lambdify(x,f)
f_prime = lambdify(x,f_prime)

#main method
program(f,-3,2)



