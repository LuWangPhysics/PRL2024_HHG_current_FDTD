#ifndef MATERIAL_RESPONSE_WSM_SUPPORT_FUNC_HPP
#define MATERIAL_RESPONSE_WSM_SUPPORT_FUNC_HPP
#include "material_response_WSM.hpp"
#include "my_func/include.hpp"


template<typename T1>
void My_material::A_calculate(Eigen::Ref<Eigen::ArrayXd> A,const T1& E){

      A+=-dt*E;
}

template<typename T1,typename T2>
void My_material::A_calculate(Eigen::Ref<Eigen::ArrayXd> A, const T1& E,T2& A_loop){

     A_loop=A-dt*E;
}

 template<typename T1>
 double My_material::J_kxky_sum(T1& j_kxky,int& chi_iter){

    
        double j_xy_sum;
        Eigen::ArrayXd temp=j_kxky.rowwise().sum()*dky;
        if(chi_iter==0){   j_xy_sum=0.5*((temp.head(N_k-1)+temp.tail(N_k-1))*dkx_p).sum();}
        else{   j_xy_sum=0.5*((temp.head(N_k-1)+temp.tail(N_k-1))*dkx_m).sum();}
        return j_xy_sum;
 }
 




template<typename T1,typename T2,typename T3>
void My_material::half_converge_method(int &n_x,int &J_flag,int& t_iter,T1&E_guess,T1&E_loop,T2& my_data,T3 & converge ){

                        Eigen::ArrayXd A_half;
                        //A  in next 0.5dt time step
                        A_calculate(A.col(J_flag),(E_guess+ME.col(n_x))*0.5, A_half);
                        auto [Jy_sum,Jz_sum]=J_wsm(J_flag,n_x,t_iter,A_half,my_data);

                        // my_data.jy(J_flag,t_iter)=Jy_sum;
                        // my_data.jz(J_flag,t_iter)=Jz_sum;
                        E_loop=E_guess;
                       
                        //---------------------
                        //case 1 yz separate converge
                        //--------------------
                   
                          // if(converge[0]!=1){
                          //    E_guess[0]=ME(0,n_x)-dt*((c_light/dx)*(McB(1,n_x)-McB(1,n_x-1))+Jy_sum/eps0);
                          // }
                          // if(converge[1]!=1){
                          //    E_guess[1]=ME(1,n_x)+dt*((c_light/dx)*(McB(0,n_x)-McB(0,n_x-1))-Jz_sum/eps0);
                          // }
                        //---------------------
                        //case 2 yz together converge
                        //--------------------
                        E_guess[0]=ME(0,n_x)-dt*((c_light/dx)*(McB(1,n_x)-McB(1,n_x-1))+Jy_sum/eps0);
                        E_guess[1]=ME(1,n_x)+dt*((c_light/dx)*(McB(0,n_x)-McB(0,n_x-1))-Jz_sum/eps0);
                
                        //---------------------
                        //bisection
                        //--------------------
                         E_guess=(E_guess+E_loop)*0.5;

            
}

template<typename T1,typename T2>
void My_material::source_direct(int& J_flag,int& t_iter,int& n_x,T1& C,T2& my_data){
  
       Eigen::ArrayXd A_half=A.col(J_flag);
         //auto [Jy_sum,Jz_sum]=J_wsm(J_flag,n_x,t_iter,A_half,my_data);
         double Jy_sum=0;
         double Jz_sum=0;
         ME(0,n_x)+=-dt*((c_light/dx)*(McB(1,n_x)-McB(1,n_x-1))+Jy_sum/eps0);
         ME(1,n_x)+=dt*((c_light/dx)*(McB(0,n_x)-McB(0,n_x-1))-Jz_sum/eps0);

         A_calculate(A.col(J_flag), ME.col(n_x));
 }



#endif

