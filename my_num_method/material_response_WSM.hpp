#ifndef MATERIAL_RESPONSE_WSM_HPP
#define MATERIAL_RESPONSE_WSM_HPP
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <math.h> 
#include <numeric>
#include <cmath>
#include "my_num_method/include.hpp"


struct My_material{
    //get the material response for each position z
    Eigen::ArrayXXd ME,McB,A,bc_holdr,bc_holdl;
    Eigen::ArrayXd kz,dkx_m,dkx_p,kx_p,kx_m;
    Eigen::RowVectorXd ky;
    double dky,dkz;

    double a,b,c,g,q;
    Eigen::Array3d v,wsm_b;
    double ef,M,mu0,eps0,c_light,hbar;
    double dt,dx, mur_coe,tau;
    int Nx,N_k;
    int chi_shift;
    //omp iter number
    int id_sec;
    int flag_dim;
    int inner_threads;
    int thread_count;
   //---------------------------------------
   // get the current related variables
   //---------------------------------------
  

    Eigen::ArrayXd t;
    Eigen::ArrayXXd Epre, Eend;
    Eigen::ArrayXXd kyy;
    


    Eigen::ArrayXXd Rho,Rho0;
    Eigen::ArrayXXcd  Gamma, Energy_omega;
    
  //check the convergence

   int max_iter=300;
   double err_tolerate=0.01;
   const  std::complex<double> I{std::complex<double>(0.0,1.0)}; 
   const double PI  =3.141592653589793238463;

template<typename T1, typename T2,typename T3>
        My_material(T1 & mesh_x ,T2 &E, T3& C,int nest_opm_flag){
        dt=E.dt;
        flag_dim=E.flag_dim;
       //this tau is the decay time for gamma and rho
        tau=10e-15;//E.tau;
       // tau=E.tau;
        dx=mesh_x.dx;
        c_light=C.c;
        hbar=C.hbar;
       //---------------------------
       // Weyl nodes in TaAs
       //---------------------------
        a=3.4e-10;
        b=3.4e-10;
        c=11.6e-10;
        //fermi energy
        ef=100e-3*C.e;
        q=-C.e;

        v<< 1e6,1e6,1e6;
 
        wsm_b<<0.06*PI/a,0,0;
        g=2;
        mu0=C.mu0;
        eps0=C.eps0;
 
        mur_coe=(1-dx/(dt*C.c))/(1+dx/(C.c*dt));

      //---------------------------
       // omp set
       //---------------------------   
      id_sec=0;
      thread_count=0;
      //this must be an even number 
      if (nest_opm_flag==0){
       inner_threads=omp_get_max_threads(); }
      else{inner_threads=(omp_get_max_threads())/2; }
      
      //for not neseted for loop
         //   
      //---------------------------
       // K space mesh
       //---------------------------        

        double k_lim;
        if(wsm_b(0)!=0){
            k_lim=3.0*wsm_b(0);
        }
        else{
           //k_lim=1.5*C.omega/v(0);
           //the BZ
           k_lim=PI/a;
        }

        double k_margin=1e7;
        N_k=60;
        Nx=mesh_x.Nx;

        if(N_k*0.1<1){
               kx_p=Eigen::ArrayXd::LinSpaced(N_k,wsm_b(0)-k_lim,wsm_b(0)+k_lim);
        }
        else{

                kx_p=Eigen::ArrayXd::Zero(N_k);
                if (wsm_b(0)!=0){
                int kx_sec1= N_k*0.9;
                int kx_sec2=N_k*0.1;
            
                N_k=kx_sec1+kx_sec2;
                kx_p.head(kx_sec1)= Eigen::ArrayXd::LinSpaced(kx_sec1,wsm_b(0)-k_lim,wsm_b(0)+k_lim);
                kx_p.tail(kx_sec2)=Eigen::ArrayXd::LinSpaced(kx_sec2,wsm_b(0)+k_lim+k_margin,C.pi/a); 
                }
                else{
                    kx_p=Eigen::ArrayXd::LinSpaced(N_k,-k_lim, k_lim);}
        }

        kx_m=-kx_p.reverse();
       // kx=Eigen::ArrayXd::Zero(2*N_k);
        ky=Eigen::RowVectorXd::LinSpaced(N_k,-k_lim, k_lim); 
        kz=Eigen::ArrayXd::LinSpaced(N_k,-k_lim, k_lim);     


        dky=ky(1)-ky(0);
        dkz=kz(1)-kz(0);
        dkx_m=kx_m.tail(kx_m.size()-1)-kx_m.head(kx_m.size()-1);
        dkx_p=kx_p.tail(kx_p.size()-1)-kx_p.head(kx_p.size()-1);
        //kxx is vertically aligned the first index corresponds to x
       
        kyy= ky.replicate(kx_p.size(),1);
        
        


  
         
    }
    ~My_material(){}


    template<typename T1,typename T2>
    void init(T1& mesh_x,T2& C);
    Eigen::ArrayXXd  theta_phi_init( int & chi_iter,int& x_iter,int & kz_iter);
    template<typename B1,typename B2,typename B3,typename B4>
    void EH_advance(int & t_iter,B1& E,B2& mesh_x, B3& my_data,B4& C);
    template<typename T1,typename T2>
    std::tuple<double, double> J_wsm(int& J_flag,int& n_x,int& t_iter,T1& A_half,T2& my_data);
    template<typename T1>
    std::tuple<Eigen::ArrayXXd , Eigen::ArrayXXd ,Eigen::ArrayXXd ,Eigen::ArrayXXd ,Eigen::ArrayXXd, 
    Eigen::ArrayXXd ,Eigen::ArrayXXd > theta_phi( int & chi_iter,int& x_iter,int & kz_iter,const T1& A_loop);
    std::tuple<Eigen::ArrayXXcd,Eigen::ArrayXXd> gamma_rho(int& J_flag,Eigen::ArrayXXd& stheta,Eigen::ArrayXXd&ctheta,Eigen::ArrayXXd&sphi,
                               Eigen::ArrayXXd& cphi,Eigen::ArrayXXd&dtheta,Eigen::ArrayXXd& dphi,
                               Eigen::ArrayXXd& pi_scaler,int & chi_iter,int& kz_iter);
   std::tuple<double, double> J_calculate(int& J_flag,Eigen::ArrayXXd& stheta,Eigen::ArrayXXd&ctheta,Eigen::ArrayXXd&sphi,
                               Eigen::ArrayXXd& cphi,Eigen::ArrayXXcd& temp,Eigen::ArrayXXd& rho_temp,int & chi_iter,int& kz_iter);
    
    //calculate A
    template<typename T1>
    void A_calculate(Eigen::Ref<Eigen::ArrayXd> A, const T1& ME);
    template<typename T1,typename T2>
    void A_calculate(Eigen::Ref<Eigen::ArrayXd> A, const T1& E,T2& A_loop);

    template<typename T1>
    double J_kxky_sum(T1& j_kxky,int&chi_iter);
    template<typename T1>
    void   boundary_conditions(int& t_iter,T1& E);



    template<typename T1,typename T2,typename T3>
    void    half_converge_method(int& n_x,int& J_flag, int& t_iter,T1&E_guess,T1&E_loop,T2&my_data, T3 & converge);
    template<typename T1,typename T2>
    void source_converge(int& J_flag,int& t_iter,int& n_x,T1& C,T2& my_data);
    template<typename T1,typename T2>
    void source_direct(int& J_flag,int& t_iter,int& n_x,T1& C,T2& my_data);
  

};



template<typename T1,typename T2>
void My_material::init(T1& mesh_x,T2& C){
    
  //define variables will be modifined in main  
    
            //storeage for y, z dimension
        ME=Eigen::ArrayXXd::Zero(2, mesh_x.Nx);
        McB=Eigen::ArrayXXd::Zero(2, mesh_x.Nx);
        Epre=Eigen::ArrayXXd::Zero(2, 2);
        Eend=Eigen::ArrayXXd::Zero(2, 2);

    
       //---------------------------
       // Current calculation initialization
       //---------------------------     
       //only store the part with source 
       int Nj=mesh_x.Nj;
        //converge current guess iteration 
    

        A=Eigen::ArrayXXd::Zero(2,Nj);
     
        
        bc_holdr=Eigen::ArrayXXd::Zero(2,2);
        bc_holdl=Eigen::ArrayXXd::Zero(2,2);
        
      

        Rho=Eigen::ArrayXXd::Zero(N_k*Nj*2,ky.size()*kz.size());
        Gamma=Eigen::ArrayXXcd::Zero(N_k*Nj*2,ky.size()*kz.size());
        Energy_omega=Eigen::ArrayXXcd::Zero(N_k*Nj*2,ky.size()*kz.size());

         chi_shift=Nj*N_k;
        int pas=0;
         # pragma omp parallel for collapse(3)  
                    for (int chi_iter=0;chi_iter<2;chi_iter++){
                            for(int J_flag=0; J_flag<Nj ;J_flag++){
                                for(int kz_iter=0;kz_iter<kz.size();kz_iter++){ 
                                   
                                    auto pi_scaler= theta_phi_init( chi_iter, pas,kz_iter);
                              Rho.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(),N_k, ky.size())=1/(1+exp((pi_scaler-ef)/(C.kb*C.T)))-1/(1+exp((-pi_scaler-ef)/(C.kb*C.T)));                                      
                                }
                                
                            }

                    }
      Rho0=Rho;
};

template<typename B1,typename B2,typename B3,typename B4>
void My_material::EH_advance(int & t_iter,B1& E,B2& mesh_x, B3& my_data,B4& C){
       
 
//---------------------------------    
//update H, H in x mesh has one less element than E
//---------------------------------  
# pragma omp parallel for
   for ( int ii = 0; ii < Nx-1; ii++){
        McB(0,ii) +=(c_light*dt/dx)* (ME(1,ii+1) - ME(1,ii));
        McB(1,ii) += -(c_light*dt/dx)* (ME(0,ii+1) - ME(0,ii));   
    } 

   
//---------------------------------  
//storage E H field at 0, 1 position at previous time step
//---------------------------------  

     Epre=ME.bottomLeftCorner(2, 2);
     Eend=ME.bottomRightCorner(2, 2);

//---------------------------------  
  //update E value for the t+1 for all x
  //aprat from the boundary
//---------------------------------  

// the front part is only vacuum, parallel that first
# pragma omp parallel for
    for(int n_x=1;n_x<mesh_x.N_front-1;n_x++)
    {   
        
               ME(0,n_x)+=-(dt*c_light/dx)*(McB(1,n_x)-McB(1,n_x-1));
               ME(1,n_x)+=(dt*c_light/dx)*(McB(0,n_x)-McB(0,n_x-1));    
               my_data.save_location(t_iter,n_x, ME(0,n_x), ME(1,n_x));
    }
   


#pragma omp parallel num_threads(2) private(id_sec)
{
  id_sec=omp_get_thread_num();
  int iter_start,iter_end;
  iter_start=(1-id_sec)*(mesh_x.N_front-1)+id_sec*mesh_x.Nx_save[1];
  iter_end=(1-id_sec)*mesh_x.Nx_save[1]+id_sec*(Nx-1);
   

                            for(int n_x=iter_start;n_x<iter_end;n_x++)
                            {       
                                if (mesh_x.J_flag[n_x]&&(abs(ME.col(n_x)).sum()>0)){
                                    //count the t_iter when it starts the parallel
                                            if(thread_count==0){
                                                std::string message_save="t_iter starts to loop";
                                                test_line(t_iter);
                                                my_data.save_note_struc(t_iter,message_save);
                                                thread_count=1;}

                                    
                                            int J_flag = std::accumulate(mesh_x.J_flag.begin(), mesh_x.J_flag.begin()+n_x, 0);
                                        
                                            //  ---------------------------------  
                                            // calculate the J for all x array
                                            //  ---------------------------------   
                                            source_converge(J_flag,t_iter,n_x,C,my_data);
                                            //source_direct(J_flag,t_iter,n_x,C,my_data);
                                            
                                            
                                    }
                                else{
                                    ME(0,n_x)+=-(dt*c_light/dx)*(McB(1,n_x)-McB(1,n_x-1));
                                    ME(1,n_x)+=(dt*c_light/dx)*(McB(0,n_x)-McB(0,n_x-1));  

                                    my_data.save_location(t_iter,n_x, ME(0,n_x), ME(1,n_x));
                                    }
                                
                               
                            }
}

//  ---------------------------------      
 //calculate boundaries at two sides
//  ---------------------------------   

 boundary_conditions(t_iter,E);


};




template<typename T1>
void My_material::boundary_conditions(int& t_iter,T1& E){
//----------------------------------------
//----------------------------------------
    //mur boundary 
//----------------------------------------
//----------------------------------------

    
//---------------------------------  
  //define left boundary of E
//---------------------------------  
   if (t_iter<E.left_bc_on) {
        // assign Eyz  to the left boundary   
        ME(0,0)=E.E_vec(1,t_iter);
        ME(1,0)=E.E_vec(2,t_iter);

    } else {
     // switch let boundary to absorbing boundary to avoid reflection
      ME.col(0)=Epre.col(1)+mur_coe*(ME.col(1)-Epre.col(0));

    }
  

     
//---------------------------------     
 //define right boundary   
//---------------------------------  
   ME.col(Nx-1)=Eend.col(0)+mur_coe*(ME.col(Nx-2)-Eend.col(1));
  
};


#endif
