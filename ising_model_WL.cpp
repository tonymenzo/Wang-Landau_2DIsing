/***************************************************************************
//Programmer: Tony Menzo
//
//Summary: This program will be implementing Wang-Landau Sampling as an alternative to the Metropolis Algorithm in calculating thermodynamic variables  for the 2D Ising model 
//
//       Sources: (1)R.H. Landau and M.J. Paez,Computational Physics:  Problem Solving with Computers(Wiley-Interscience, 1997)
//       (2)Landau, D. P., Shan-Ho Tsai, and M. Exler. "A new approach to Monte Carlo simulations in statistical physics: Wang-Landau sampling." American Journal of Physics 72.10 (2004): 1294-1302.
//
//
//   See the readme.txt for description
// 
*****************************************************************************/

// include files
#include <iostream>		// cout and cin
#include <iomanip>		// manipulators like setprecision
#include <fstream>		// file input and output
#include <cmath>
#include <string>               // C++ strings                                 
#include <sstream>              // C++ stringstream class (could omit iostream)  
using namespace std;		// we need this when .h is omitted

#include <gsl/gsl_rng.h>	// GSL random number generators
#include <gsl/gsl_randist.h>	// GSL random distributions
#include "GslSpline.h"  // Header file for the GSL Spline class
#include <omp.h>     // This is what includes the OpenMP directives
#include <time.h>

// function prototypes
extern unsigned long int random_seed ();	// routine to generate a seed
int calculate_energy (int configuration[]);
double calculate_partition_Z(double g[], double energy[], int size,double E_high, double kT);

//Global constants
const int J_ising = 1;
const int dimension = 2;
const int L = 15;                 //number of spins in a given row
const int num_sites = L * L;  //total number of spins

//Inline functions
inline double sqr (double z) {return z*z;}  // inline function for z^2
/*******************************************************************/
int main()
{
   
  double f_min;          //minimum modification factor 
  double f_initial;      //initial mod factor
  double pow_decrease;   //power to raise f after each run must be [0.,1.]
  double flat_criterion; //how flat do you want the histogram
  int mc_steps;          //number of monte carlo simulations in energy space
  
  int E;                 //Initial Energy
  int E_new;             //New energy after random spin flip
  int E_initial;
  double prob;           //probability that a given energy will become new ene.
  int min_steps = 10000; //minimum number of mcs
  string type = "cubic"; //string for cubic spline
  
  double start,end;
  
  int config_WL[num_sites]; //Lattice of +/- spins

  
  
/*
  ofstream fit_Z ("g_E_fit_Z.dat");
  fit_Z << "# Z(T)         cubic spline fit "<<endl;
  ofstream timing_dat ("L_time.dat");
  timing_dat << "# L           time"<< endl;

  ofstream Z ("partition.dat");
  Z << "#  kT         Z(T)   "<<endl;
 */

  f_initial = 2.71;
  f_min = exp(pow(10,-8));
  pow_decrease = 0.5;
  flat_criterion = 0.05;
  

  //Calculating maximum energy  (1)
  int E_high;
  for (int i = 0; i < num_sites; i++)
  {
    config_WL[i] = 1;
    E_high = calculate_energy(config_WL);
  }
  
  E_high = abs(E_high);

  cout << "Max Energy: "<<E_high << endl;
  double g[2*E_high];  //array for energy density values
  int hist[2*E_high];  //array for histogram values
  

  //  Set up the GSL random number generators (rng's)
  gsl_rng *rng_ptr = gsl_rng_alloc (gsl_rng_taus);   // allocate an rng 
  gsl_rng_set (rng_ptr, random_seed());	             // seed the rng 


  // generate a random starting configuration (2)
  for (int i = 0; i < num_sites; i++)
  {
    double random = gsl_ran_flat (rng_ptr, 0., 1.);
    if (random > 0.5)
    {
      config_WL[i] = -1;  // spin down
    } 
    else
    {
      config_WL[i] = +1;  // spin up
    }
     E_initial = calculate_energy(config_WL);
     
  } 
  
  E = E_initial;
//initialize density of states g(E) for earch E(using log(g(E))!!! So g(E)=1...
//(3)
  for (int i=0; i<2*E_high;i++)
  { 
    g[i] = 0;
  }

//For a total count of the histogram values
  int histTotal[2*E_high];
  for (int i =0;i< 2*E_high;i++)
  { 
     histTotal[i] = 0;
  }

//Phew...let's do the Monte Carlo now... (4)

//Let's see how increasing L scales with time
start = omp_get_wtime();
   

while (f_initial> f_min){
  
//Only take ln of f and g!!
  double lnf= log(f_initial);
  int cont = 0;   //determines if we stop the do-while
  
//Start with flat histogram every loop
    for( int i = 0; i< 2*E_high; i++)
    {
       hist[i] = 0;
    }
     
//Monte Carlo
  do { 
        for (mc_steps = 0; mc_steps < min_steps; mc_steps++)
        {
           for (int i = 0; i < num_sites; i++)  // Entire loop is only one mcs  
           {
              // pick a random lattice site
              double random = gsl_ran_flat (rng_ptr, 0., 1.);
              int id = int(random * (num_sites));  // from 0 to num_sites

              // flip that spin (i.e., if +/- 1, change to -/+ 1)
              config_WL[id] *= -1; 

              E_new = calculate_energy( config_WL);  // new energy
        
           
              if (g[E_high - E_new] <= g[E_high-E])
              {
                 E=E_new;
              }
              else
              {
                prob = exp(g[E_high-E]-g[E_high-E_new]);
           
                //Accept if proposed move has lower g (prob >=1
                // or if random number [0:1] is less than prob
                double random1 = gsl_ran_flat(rng_ptr,0.,1.);
                if(prob >= random1) 
                {
                  E=E_new;
                }
                //Otherwise reject the change and revert back
                else    
                {
                   config_WL[id] *= -1;
                }
           

              }
           
              //Update g(i)
              g[E_high-E]+=lnf;

              //Update histograms 
              hist[E_high-E]++;
              histTotal[E_high-E]++;
         
           }
        }
      
      //lets find maximum and minimum values in histogram H-max and H-min
      int hmax = 0;
      for (int i =0; i< 2*E_high;i++)
      {
         if(hist[i] > hmax) 
         {
           hmax = hist[i];
         }
      }
       int hmin = 1000000000;
      for (int i =0; i< 2*E_high;i++)
      {
         if(hist[i] < hmin) 
         {
           hmin = hist[i];
         }
      }
      //quantifying flattness (average value/height)
      double flatness = ((double)hmax - (double)hmin)/((double)hmax 
                        +(double)hmin);
      //if the histogram is flat enough we can stop doing monte carlo loops
      if (flatness < flat_criterion)
      {
         cont = 1;
      }
      
    }while (cont);
    
   //Tell me what my mod factor is 
   cout << fixed << setprecision(9)<< "Modification factor f is currently " 
    << f_initial << endl;

   //descrease f and start the whole process over again
   f_initial=pow(f_initial,pow_decrease);
  
 }

   end = omp_get_wtime();     // get the end time

   cout << "Time to calculate was " << (end - start) <<" seconds."<< endl;


// Normalize the density of states, knowing that the lowest energy state is double degenerate lnCorrection = log( (exp(lngE[0])+exp(lngE[-1]))/4. ) 
//(5)
  
    double lgC;
    if (g[E_high]<g[0]){
        lgC = g[0] + log(1+ exp(g[E_high]-g[0])) - log(2.);
    }

    else{
        lgC = g[E_high] + log(1+ exp(g[0]-g[E_high])) - log(2.);
    }
    for (int i = 0;i<2*E_high;i++) 
    {
      g[i] -= lgC;
      
    }
    g[0]=0.0;
  

    ofstream out ("g_E.dat"); 
    out << "#   E              g(E)          Total Samplings Histogram     "
    <<   endl;

    //print out my values for g(E) and histTotal
    int n=0;
    for (int i = 0;i< 2*E_high;i++)
    {  
        if(g[i]>0.){
        out << fixed << setprecision(12)<<" " 
          << (double)(i-E_high) << "      " 
          << g[i]<< "  " <<"  "<<  histTotal[i]<<endl;
          n++;
        }
     }
    out.close();
    //I'm only interesting in non-zero energies, let's get them
    double reduced_g_array[n];
  
    double E_g_value = 0.0;
    int t=0;
    for (int i =0; i<2*E_high;i++)
    {
      if(g[i] > 0.)
      {
         E_g_value = g[i];
         t++;
         reduced_g_array[t-1] = E_g_value;
      } 
    
   
    }

    //need the values of the energies
    double counting_array[n];
    int index=0;
    for (int i =0; i<2*E_high;i++)
    {
      if(g[i]>0.0)
      {
         index++;
         counting_array[index-1] = (double)i-E_high;
      } 
      
    }


    ofstream fit ("g_E_fit.dat");
    fit << "# E          cubic spline fit "<<endl;
    //Lets fit the data and print to "partition.dat"
    //(6)
    Spline my_cubic_spline_g (counting_array, reduced_g_array, n, type);
    for(int i = counting_array[0]; i<=counting_array[n-1];i+=2)
    {
      
      // Evaluate the spline
      double yval = my_cubic_spline_g.y (double(i));   
      fit << " "<<i<< "      " << yval<<endl;
      
    }

  fit.close();
/************Needs work!!!!! Work in progress*****************************/
/*
//We need the partition function in order to calculate other thermo variables
  double partition_val;    //value for partition function at kT
  double partition[191];   //array for partition values at various 
  double index_Z[191];     //array for kT needed for fit
  int z=0;                 //counter to allow me to imput Z vals into array
   for (double i =20; i>=0.9;i-=0.1)
    {
       partition_val = calculate_partition_Z(reduced_g_array, counting_array, n, E_high, i);
       z++;
       Z << "  " << i << "     " << partition_val<< endl;
       partition[z-1] = log(partition_val);
     }
    
    for(int i=0; i<191;i++)
    {
       index_Z[i] = i;
    }

//we need to fit the partition function so we can access it derivative with
//respect to kT
    Spline my_cubic_spline_Z (index_Z, partition, 191, type);

   for(int i = 0; i<191;i++)
   {
      
      // Evaluate the spline
      double yval_prime = my_cubic_spline_Z.yp (double(i));
         
      fit_Z << " "<<i<< "      " << -yval_prime<<endl;
      
    }
*/    


return(0);   
}

/*******************Calculating Energy of 2D Ising Model**********************/
//************************* calculate_energy ******************************
//
//  Given the array of integers configuration[0...num_sites-1], which 
//   specifies the spin at each lattice point, find the energy of that 
//   configuration [eq.(13.6) in Session 13 notes].
//  Note that a freeboundary condition is specified here.
//
//*************************************************************************
int 
calculate_energy (int configuration[])
{
  int nearest = 0;
  double energy = 0.;

// go through the 2d lattice with free boundary conditions
    for (int i = 0; i < L - 1; i++)
    {
      for (int j = 0; j < L - 1; j++)
      {
        int id = i + j * L;

        // x direction
        nearest = id + 1;
        energy += - J_ising * double(configuration[id]*configuration[nearest]);

        // y direction;
        nearest = id + L;
        energy += - J_ising * double(configuration[id]*configuration[nearest]);
      }
    }


  return (energy);
}   

/**************************calculate_partition_Z*****************************/
double 
calculate_partition_Z(double g[], double energy[], int size,double E_high, double kT)
{
    double sum = 0; 
    for (int j =0; j<size;j++)
    {
       sum+= exp((g[j]-log(E_high))-(energy[j]+log(E_high))/(kT));

    }
  return (sum);
}








