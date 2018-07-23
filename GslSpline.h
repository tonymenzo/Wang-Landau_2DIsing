//  file: GslSpline.h
//
//  Header file for the GslSpline C++ class.
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02/10/09  original version
//
//  Notes:
//   * We encapsulate GSL spline functions in this class
//   * The GSL header file is included in GslSpline.cpp and
//      anywhere we want to create GslSpline objects.
//   * Based on the documentation for the GSL library under
//      "Interpolation" and on the "Using GSL Interpolation Functions"
//      780.20 handout.
//   * See the GSL documentation for error handling details.
//
//  To do:
//   * add other splining methods
//
//*****************************************************************

// The ifndef/define macro ensures that the header is only included once
#ifndef GSLSPLINE_H
#define GSLSPLINE_H

// include files
#include <string>     // C++ strings                                 
#include <gsl/gsl_spline.h>	  // header for gsl spline routines 

class Spline
{ 
  public:
    Spline (double* x_array, double* y_array, int num_pts, 
            std::string type);  // constructor
    ~Spline ();  // destructor

    // accessor functions
    double y (const double x);    // find y(x)
    double yp (const double x);   // find y'(x)
    double ypp (const double x);  // find y''(x)

  private:    
    std::string spline_type;       // type of spline 
    gsl_interp_accel *accel_ptr;   // accelerator pointer
    gsl_spline *spline_ptr;        // spline pointer

};

#endif
