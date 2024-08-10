#ifndef MATERIAL_RESPONSE_WSM_J_FUNC_HPP
#define MATERIAL_RESPONSE_WSM_J_FUNC_HPP
#include "material_response_WSM.hpp"
#include "my_func/include.hpp"
#include <complex>



Eigen::ArrayXXd  My_material::theta_phi_init( int & chi_iter,int& x_iter,int & kz_iter){
          double  chi=pow(-1,chi_iter);
          Eigen::ArrayXXd pix,piy;   
          double piz;
          Eigen::ArrayXXd pi_scaler;
          Eigen::ArrayXXd kxx;
          if(chi_iter==0){kxx= kx_p.replicate(1,ky.size());}
          else{kxx= kx_m.replicate(1,ky.size());}
         
           pix=v(0)*hbar*(kxx-chi*wsm_b(0));
           piy=v(1)*(hbar*(kyy-chi*wsm_b(1)));  
           piz=v(2)*(hbar*(kz(kz_iter)-chi*wsm_b(2)));
           
                    
           pi_scaler=(pow(pix,2)+pow(piy,2)+pow(piz,2)).sqrt();
       
     
 return pi_scaler;
};

template<typename T1>
std::tuple<Eigen::ArrayXXd , Eigen::ArrayXXd ,Eigen::ArrayXXd ,Eigen::ArrayXXd ,Eigen::ArrayXXd, 
Eigen::ArrayXXd ,Eigen::ArrayXXd > My_material::theta_phi( int & chi_iter,int& x_iter,int & kz_iter,const T1& A_loop){
          double  chi=pow(-1,chi_iter);
          Eigen::ArrayXXd pi_xy,pix,piy;   
          double piz,py_dot,pz_dot;
          Eigen::ArrayXXd pi_scaler ,stheta,ctheta,sphi,cphi,dtheta,dphi;

          Eigen::ArrayXXd kxx;
          if(chi_iter==0){kxx= kx_p.replicate(1,ky.size());}
          else{kxx= kx_m.replicate(1,ky.size());}
         
           pix=v(0)*hbar*(kxx-chi*wsm_b(0));
           piy=v(1)*(hbar*(kyy-chi*wsm_b(1))-q*A_loop(0));  
           piz=v(2)*(hbar*(kz(kz_iter)-chi*wsm_b(2))-q*A_loop(1));
           
                    
           py_dot=q*v(1)*ME(0,x_iter);
           pz_dot=q*v(2)*ME(1,x_iter);

           pi_xy=(pow(pix,2)+pow(piy,2)).sqrt();
           pi_scaler=(pow(pix,2)+pow(piy,2)+pow(piz,2)).sqrt();
       

           stheta=pi_xy/pi_scaler;
           ctheta=piz/pi_scaler;
           sphi=piy/pi_xy;
           cphi=pix/pi_xy;

       
          dtheta=-(stheta*pz_dot-ctheta*(sphi*py_dot))/pi_scaler;
          dphi=(py_dot*cphi)/pi_xy;
    
    
                      
 return {stheta,ctheta,sphi,cphi,dtheta,dphi,pi_scaler};
};




std::tuple<Eigen::ArrayXXcd,Eigen::ArrayXXd> My_material::gamma_rho(int& J_flag,Eigen::ArrayXXd& stheta,Eigen::ArrayXXd&ctheta,Eigen::ArrayXXd&sphi,Eigen::ArrayXXd& cphi, Eigen::ArrayXXd&dtheta,Eigen::ArrayXXd& dphi,Eigen::ArrayXXd&pi_scaler,int & chi_iter,int& kz_iter){
  
          double  chi=pow(-1,chi_iter);
      
            Eigen::ArrayXXd rho_temp,rho0_temp,f1_temp;
            Eigen::ArrayXXcd  ew_temp,gamma_temp,temp,f2_temp;
            
            Energy_omega.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size())+=dt*pi_scaler/hbar;
            ew_temp=Energy_omega.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size());
            //get the related section to update
            rho_temp= Rho.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size());
            rho0_temp=Rho0.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size());
            gamma_temp=Gamma.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size());
            temp=gamma_temp*exp(-2.0*I*ew_temp);
            
   
            // f1_temp=-chi*2.0*dtheta*temp.real()+2.0*dphi*stheta*temp.imag()-(rho_temp-rho0_temp)/tau;
            // f1_temp=((rho_temp>1.0)&&(f1_temp>0)).select(0,f1_temp);
            // f1_temp=((rho_temp<-1.0)&&(f1_temp<0)).select(0,f1_temp);
            // rho_temp+=f1_temp*dt;



             rho_temp+=(-chi*2.0*dtheta*temp.real()+2.0*dphi*stheta*temp.imag()-(rho_temp-rho0_temp)/tau)*dt;
             gamma_temp+=(chi*0.5*rho_temp*exp(2.0*I*ew_temp)*(dtheta-chi*I*stheta*dphi)
                 +chi*I*ctheta*dphi*gamma_temp-gamma_temp/tau)*dt;
             //gamma_temp=(rho_temp.abs()>1).select(0,gamma_temp);
          
                
             
            temp=gamma_temp*exp(-2.0*I*ew_temp);
            Rho.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(), N_k, ky.size())=rho_temp;
            Gamma.block(J_flag*N_k+chi_iter*chi_shift, kz_iter*ky.size(),N_k, ky.size())=gamma_temp;

            return {temp, rho_temp};
}

std::tuple<double, double> My_material::J_calculate(int& J_flag,Eigen::ArrayXXd& stheta,Eigen::ArrayXXd&ctheta,Eigen::ArrayXXd&sphi,Eigen::ArrayXXd& cphi,Eigen::ArrayXXcd & temp,Eigen::ArrayXXd & rho_temp,int & chi_iter,int& kz_iter){
          double  chi=pow(-1,chi_iter);
            //JxJy 
            Eigen::ArrayXXd jkxky;
             double jy,jz;
            jkxky=(rho_temp+1)*stheta*sphi-2.0*(temp*(chi*sphi*ctheta+I*cphi)).real();
            jy=J_kxky_sum(jkxky,chi_iter)*(g*q*v(1)*dkz/pow(2*PI,3));
            
            jkxky=(rho_temp+1)*ctheta+2.0*(temp*chi*stheta).real();
            jz=J_kxky_sum(jkxky,chi_iter)*(g*q*v(2)*dkz/pow(2*PI,3));
            

            return {jy,jz};
}



template<typename T1,typename T2>
std::tuple<double, double>  My_material::J_wsm(int& J_flag,int& n_x,int& t_iter, T1& A_half,T2& my_data){
  
double Jy_sum=0;
double Jz_sum=0;


# pragma omp parallel for collapse(2) num_threads(inner_threads)  reduction(+:Jy_sum,Jz_sum)
//# pragma omp parallel for collapse(2) reduction(+:Jy_sum,Jz_sum)
 for (int chi_iter=0;chi_iter<2;chi_iter++){
            for (int kz_iter=0;kz_iter<kz.size();kz_iter++){
                           
                              
                           //theta phi analytical so everything on E (t)
                            auto [stheta,ctheta,sphi,cphi,dtheta,dphi,pi_scaler]=theta_phi(chi_iter,n_x,kz_iter,A.col(J_flag));
               
                            //gamma rho on (t+0.5dt)
                            auto [temp, rho_temp]= gamma_rho(J_flag,stheta,ctheta,sphi,cphi,dtheta,dphi,pi_scaler,chi_iter,kz_iter);
                  
                           //get the theta phi on (t+0.5dt)
                            auto [stheta1,ctheta1,sphi1,cphi1,dtheta1,dphi1,pi_scaler1]=theta_phi(chi_iter,n_x,kz_iter,A_half);
                            // J is analytical with all so it is  on (t+0.5dt)
                            auto [jy,jz]= J_calculate(J_flag,stheta1,ctheta1,sphi1,cphi1,temp,rho_temp,chi_iter,kz_iter);
                            Jy_sum+=jy;
                            Jz_sum+=jz;
                            //-----------------------
                            //check thread iterations;
                            //-----------------------
                            // int t_id = id_sec*inner_threads+omp_get_thread_num();
                            // my_data.thread_iteration[t_id]+=1;
                           
            } 
 }
    

return {Jy_sum,Jz_sum};
    
}


template<typename T1,typename T2>
void My_material::source_converge(int& J_flag,int& t_iter,int& n_x,T1& C,T2&  my_data){

         int iter_converge=0;
        // std::vector<int> converge={0,0};
        //  Eigen::Array<bool, 1, 2> converge(2);
        //  converge = Eigen::Array<bool, 1, 2>::Zero(2);
        Eigen::ArrayXd converge=Eigen::ArrayXd::Zero(2);
 
         Eigen::ArrayXXd f=Eigen::ArrayXXd::Zero(2,2);
         Eigen::ArrayXd E_guess=ME.col(n_x);
         Eigen::ArrayXd E_loop=ME.col(n_x);
       
         Eigen::ArrayXd A_loop,temp;
       

        //--------------------------------
         //first guess to start with
        //--------------------------------
         A_calculate(A.col(J_flag),E_loop, A_loop);
         auto [Jy_sum,Jz_sum]=J_wsm(J_flag,n_x,t_iter,0.5*(A.col(J_flag)+A_loop),my_data);

         E_guess(0)=ME(0,n_x)-dt*((c_light/dx)*(McB(1,n_x)-McB(1,n_x-1))+Jy_sum/eps0);
         E_guess(1)=ME(1,n_x)+dt*((c_light/dx)*(McB(0,n_x)-McB(0,n_x-1))-Jz_sum/eps0);
   
         
      
      //check tthe J and E(1/2) need to be consistent, but only the driving field converge is important, the other dimension depend all on E_drive so 
      //dont need to get convergence
       while((converge[0]*converge[1] != 1) && (iter_converge <max_iter)){
 
                //---------------------------------------
                //the current must integrate over the entire space!
                //---------------------------------------    
                half_converge_method(n_x,J_flag,t_iter,E_guess,E_loop,my_data,converge);
 
                //---------------------------------------
                // check convergence 
                //---------------------------------------
                for(int dim_loop=0;dim_loop<2;dim_loop++){
         
                                     if(abs(E_guess[dim_loop])<1){
                                         if(abs(E_guess[dim_loop]-E_loop[dim_loop])<err_tolerate){converge[dim_loop]=1;}
                                        
                                    }
                                    else{
                                         if(abs(E_guess[dim_loop]-E_loop(dim_loop))<abs(E_guess[dim_loop]*err_tolerate)){converge[dim_loop]=1;}
                                         
                                    }
                                       //if(abs(E_guess[dim_loop]-E_loop(dim_loop))<abs(E_guess[dim_loop]*err_tolerate)){converge[dim_loop]=1;}
                                    
                }
                
    
                my_data.test.block(id_sec*2,iter_converge,2,1)=E_guess;
                my_data.testconv.block(id_sec*2,iter_converge,2,1)=converge;
                iter_converge++;           
      }
      
              
      //case 1.if doesnt converge break
      //case 2. if doesnt converge, get average of all the convergence iteration
      //update E(x,t) 
      if(iter_converge==max_iter&&((converge[0]*converge[1])==0)){
              //----------------------
              //break the code if doesnt converge
              //----------------------
              // std::cout<<"does not converge, code break, t_iter=" << t_iter<<std::endl; 
              // my_data.save_check();
              // std::cout<<converge[0]<<","<<converge[1]<<std::endl;
              // exit(0);

              //----------------------
              //get the average
              //----------------------
              std::cout<<"does not converge, t_iter=" << t_iter<<", nx="<<n_x<<", "<<"id_parallel="<<id_sec<<", "<<converge[0]<<converge[1]<<std::endl; 
              //std::string message_save="does not converge "+std::to_string(converge[0])+std::to_string(converge[1]);                                 
              //my_data.save_note_struc(t_iter,message_save);
              //save the convergence trance every 200 points
              if (t_iter%20==0){
              my_data.save_check(t_iter);
              }
              //count time of the main one not converge 
             // my_data.E_z_converge_check[n_x]+=1-converge[flag_dim];
              E_guess(0)=my_data.test.row(id_sec*2).sum()/max_iter;
              E_guess(1)=my_data.test.row(id_sec*2+1).sum()/max_iter;
              
              ME.col(n_x)=E_guess;
              A_calculate(A.col(J_flag),E_guess);
         }
         else{ 
          ME.col(n_x)=E_guess;
          A_calculate(A.col(J_flag),E_guess);
        
        }

        
    
}










#endif
