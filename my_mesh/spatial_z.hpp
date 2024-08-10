#pragma once
#include <math.h>

struct  spatial_z
{      

	 
        double dx;
        double L_wsm;
        double L_air;
        Eigen::ArrayXd x;
        std::vector<int> J_flag;  

        std::vector<int> Nx_save;
     	int Nx, Nj,N_front,N_wsm,N_air;
       

	
  template<typename T1,typename T2>
	spatial_z(T1& E,double L_wsm_, double L_air_,  double dx_  , double tau_,T2& C)
        {   L_wsm=L_wsm_;
            L_air=L_air_;
            dx=dx_;
           
           //make sure that the reflected pulse does not interference with the boundary
           //make sure it is exactly E.left_bc_on*E.dt*C.c/dx/2; 
           N_front=1*(E.left_bc_on*E.dt*C.c/dx/2);
           int N_end=200;
           N_wsm=L_wsm/dx;
           N_air=L_air/dx;
           std::vector<int>::const_iterator p_first,p_last;
          int N_start=0.4*tau_*C.c/dx;
           
           //if air gap =0
           if(L_air<1e-14){
              //remember counting from zero, the last index of arry is Nz-1
              Nx=N_end+N_front+2*N_wsm;
              J_flag=std::vector<int>(Nx,0);
              //this part has current J
              std::fill(J_flag.begin()+N_front,J_flag.begin() + N_front+2*N_wsm,1);

            
            }else{
             //remember counting from zero, the last index of arry is Nz-1
           
              Nx=N_end+N_front+2*N_wsm+N_air;
              J_flag=std::vector<int>(Nx,0);
              //this part has current J
              //std::fill(J_flag.begin()+1 ,J_flag.begin()+2,1);.begin()+1 is location= [1] .begin()+2 means ends (open bracket at [2]) 
              // thus the final assigend value is J[1]=1
              std::fill(J_flag.begin()+N_front ,J_flag.begin() + N_front+N_wsm,1);
              std::fill(J_flag.begin()+N_front+N_air+N_wsm ,J_flag.begin() + N_front+2*N_wsm+N_air,1);
                  
              //-----------------------------
              //define the saving location  
              //-----------------------------
              //starting position 5fs apart from left boundary
          
             
            };

         // if L_air =0, then the center save is inside the WSM
          int N_middle=N_front+N_wsm+(int) N_air/2;
           Nx_save={N_start,N_middle,Nx-(N_end-5)};
           x=Eigen::ArrayXd::Zero(Nx); 
           Nj=0;

           
           for (int z_i=0;z_i<Nx;z_i++){
                x(z_i)=dx*z_i;
                Nj+= J_flag[z_i];
       
           }
        
         
            
        };
	~spatial_z(){} ;
 	
};
