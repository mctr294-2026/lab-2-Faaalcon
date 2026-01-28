#include "roots.hpp"
#include <iostream>
#include <string>
#include <cmath>

double TOLERANCE = 1e-6;
int MAX_ITER = 1000000;

bool bisection(std::function<double(double)> f, double a, double b, double *root) {
    //bounds setup
    double fa = f(a);
    double fb = f(b);
    
//opposite signs check
if (fa * fb >= 0) {return false;}

for(int i = 0; i < MAX_ITER; i++){
    double mid = (a + b)/2;
    double fmid = f(mid);

 // Check for convergence
    if (std::abs(fmid) < TOLERANCE || std::abs((b - a) / 2.0) < TOLERANCE) {
        //found the root
        *root = mid;
        return true;
    } else {
        if(fa * fmid < 0){
        b = mid; fb= fmid;
        } else {
        a = mid; fa = fmid;}
    }
}
return false;
}

bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {
    double fa = f(a);
    double fb = f(b);
//opposite signs check
if (fa * fb >= 0) {return false;}

for(int i = 0; i < MAX_ITER; i++){
    double mid =  a - (fa*(b - a)) / (fb - fa);
    double fmid = f(mid);
    //if we found the root
    if(std::abs(fmid) < TOLERANCE ){*root = mid; return true;}
    //If f(a) and f(c) have opposite signs, set b = mid
    if (fa * fmid < 0) { b = mid; fb = fmid;
    //If f(b) and f(c) have opposite signs, set a = mid.
    } else {a = mid; fa = fmid;}    
}
return false;
}

bool newton_raphson(std::function<double(double)> f, std::function<double(double)>
                     g, double a, double b, double c, double *root) {
    double fa = f(a);
    double fb = f(b);
//opposite signs check
if (fa * fb >= 0) {return false;}

double x = c; //guess
for (int i = 0; i < MAX_ITER; ++i) { //run the loop like normal
//first check if we hit the right one
        if (x < a || x > b) {return false;}//is the guess to far out
        double fx = f(x);
        double gx = g(x); // Derivatives
        if (std::abs(fx) < TOLERANCE){*root = x; return true;}
        if (std::abs(gx) < 1e-12){return false;}
        //formula: x_new = x - f(x)/f'(x)
        double x_new = x - (fx / gx);
        // convergence? 
        if (std::abs(x_new - x) < TOLERANCE) {*root = x_new; return true;}
        //run it back if that wasn't it
        x = x_new;
}
return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    // Secant requires two points to start. 
    // We use the guess 'c' and a slightly perturbed point to approximate the tangent.
    double x0 = c;
    double x1 = c + 0.001; 
    
    // Ensure the second point is within bounds, if not, go the other way
    if (x1 > b) x1 = c - 0.001;

    for (int i = 0; i < MAX_ITER; ++i) {
        // Check bounds
        if (x1 < a || x1 > b) {
            return false;
        }

        double f0 = f(x0);
        double f1 = f(x1);

        // Check success
        if (std::abs(f1) < TOLERANCE) {
            *root = x1;
            return true;
        }

        // Avoid division by zero
        if (std::abs(f1 - f0) < 1e-12) {
            return false;
        }

        // Secant formula: x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        double x_new = x1 - f1 * (x1 - x0) / (f1 - f0);

        // Check step size convergence
        if (std::abs(x_new - x1) < TOLERANCE) {
            *root = x_new;
            return true;
        }

        // Update points
        x0 = x1;
        x1 = x_new;
    }

    return false;
}