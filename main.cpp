#include <stdio.h>
#include <iostream>
#include <stdlib.h>                                                                 // include random number
#include <sstream>                                                                  // for precision of the number to string
#include <fstream>                                                                  // read and write file
#include <typeinfo>                                                                 // templete of cwise operation of vector
#include <string>	
#include <time.h>
#include <vector>
#include <thread>
#include <mutex>
#include <iomanip>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "omp.h"
#include "my_func/include.hpp"
#include "my_mesh/include.hpp"
#include "my_num_method/include.hpp"
#include <unistd.h>
#include <iostream>

int main(int argc, char* argv[]){ 

//---------------------------------------------------------
//define input pulse parameters
//---------------------------------------------------------
//if the cluster is dual-core multi-threading, then divide by two. 
int N_thread=omp_get_max_threads()/2;  
//int N_thread=omp_get_max_thread();

omp_set_num_threads(N_thread);
std::cout<<"total threads="<<N_thread<<std::endl;
// -----------------------------------------------------------
//the omp for nested loops
// -----------------------------------------------------------
//allow nested loops if flag set to zero, nested omp is disabled
int nest_omp_flag=1;
if (nest_omp_flag==0){}
else{omp_set_nested(1);}

// -----------------------------------------------------------
//for IR
// -----------------------------------------------------------
double lambda   = 1000e-9;
double tau = 10e-15;                                                        
double E0=2e8;            //2e8
int t_check_num=3000;      
// -----------------------------------------------------------                                                
//for THz
// -----------------------------------------------------------
// double lambda   = 10*1e-6;
// double tau = 100e-15;                                                        
// double E0=1.5e6;    //2e7
// int t_check_num=20000;


CONST C(lambda);
//---------------------------------------------------------
//construct saving name and time count
//---------------------------------------------------------
C.name_tail=(std::string) "test_";
mytime time_count;
time_count.start();

// ---------------------------------------------------------
// define input electric fields
// ---------------------------------------------------------
// for 3D identical mesh size dt=dz/sqrt(3)/c, for 1D is dz/c
double L_wsm=100e-9;//lambda/100;                                                     
double L_air=0*100e-9+(100e-9)*atof(argv[1]);
double dx,dt;
if ((L_air<51e-9)&&(L_air>0)){
   dx=2e-9;
}
else{
     dx=3e-9;
    }
//4dt cross one cell
dt=dx/(C.c*2);
double phi_loop=0;
double L_total=L_wsm*2+L_air;
E_field E(dt,tau,E0,C,phi_loop,L_total);

//---------------------------------------------------------
//define propagation direction  mesh
//---------------------------------------------------------          
MESH::spatial_z mesh_x(E,L_wsm,L_air,dx,tau,C);                                
// ---------------------------------------------------------
//define material response
//include numerical method
// ---------------------------------------------------------

My_material WSM(mesh_x,E,C,nest_omp_flag);
WSM.init(mesh_x,C);

// ---------------------------------------------------------
// initialize saving path 
// ---------------------------------------------------------
std::string save_path=make_string(E0,tau,lambda,L_wsm,L_air,C,dx);

// --------------------------------------------------------------
// save output data
// ---------------------------------------------------------
save_data my_data(save_path,E,mesh_x,N_thread);

int t_iter;
for (t_iter=0;t_iter<E.t_stop;t_iter++){

               WSM.EH_advance(t_iter,E,mesh_x,my_data,C);
               //manage the checkout part                 
               if((t_iter%500)==0){
                    my_print(t_iter); 
                    if(((t_iter)%t_check_num)==0){ my_data.check_point( WSM,t_iter);}
               }
     
                       
}

my_data.save_E(WSM,t_iter,mesh_x,C);

// -----------------------------------------------------------
// REMOVE all the input pointers if exist
//     VERY IMPORTANT
// -----------------------------------------------------------

time_count.end();


return 0;
}
