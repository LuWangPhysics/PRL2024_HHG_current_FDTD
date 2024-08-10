#pragma once
#include <complex>
struct CONST{
      // ===================================================================
    //                     FUNDAMENTAL CONSTANTS
    // ===================================================================
        

    const double c     = 299792458;       // speed of light in free space in [m/s]
    const double M     = 9.10938215e-31; // mass of electron in [kg]
    const double pi   = 3.141592654;
    const double hbar  = 1.0545718e-34;
    const double kb    = 1.38064852e-23; // boltzmann constant
    const double eps0  = 8.85418782e-12; // permitivitty of free space
    const double mu0   = 1.25663706212e-6;
    const  std::complex<double> i{std::complex<double>(0.0,1.0)};
    const double e     = 1.60217649e-19;  // elementary charge in [C]
    const double T     = 300;     //room temperature
    const double g     =2; // degeneracy for Weyl nodes with one chirality
    double lambda;
    std::string name_tail;
    double omega;
    
    
    CONST(double& lambda_){
           lambda=lambda_;
        omega= 2*pi*c/lambda_;   
    
    }
    ~CONST(){};
}; 
