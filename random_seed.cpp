//*********************************************************************
//
//  random_seed returns a starting seed for random number generators
//
//  Based on an contribution to the gsl-discuss@sources.redhat.com
//   discussion list by Robert Brown of Duke (rgb@phy.duke.edu)
//   on 01-Mar-2004.
//
//  Using /dev/urandom to avoid blocking, but the random numbers
//   may not be secure.
//
//  Revision history:
//   05-Mar-2004  original version based on Robert Brown's code snippet
//   16-Feb-2007  replaced Linux version with crude Windows version
//                  NOT A GOOD FIX!
//
//**********************************************************************   

// include files
#include <ctime>		// for the random number seed

// function prototypes
unsigned long int random_seed ();	// routine to generate a seed

//********************************************************************
unsigned long int
random_seed ()
{
  unsigned long int seed = time(0)*time(0);  // use the current time (crude!!)  
  return (seed);		// send back the seed!
}

