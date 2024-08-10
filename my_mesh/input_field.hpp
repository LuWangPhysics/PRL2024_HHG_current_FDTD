#ifndef INPUT_FIELD_H
#define INPUT_FIELD_H
#include <Eigen/Dense>
//#include <fftw3.h>
#include <math.h> 
//#include "my_func/include.hpp"
//CONST C;
struct E_field{
    Eigen::ArrayXXd E_vec;
    Eigen::ArrayXd t; 
    Eigen::ArrayXd t_arr,tempE; 

    double tau;
    double E0;
    double dt;
    double t_pulse_window,t_shift;
    double phi; //relative phase
    int Nt;
    int left_bc_on;
    int flag_dim;
    int t_stop;

template<typename T1>
        E_field( double & dt_,double & tau_ ,double &E0_, T1& C,double phi_loop,double L_total){
        dt=dt_;
        tau=tau_;
       
        E0=E0_;
        phi=phi_loop;
       
//define the input margin of the pulse
        double t_period=2*C.pi/C.omega;
        double aim_shift=3.5*tau;
        //distance from zero
        int all_period=aim_shift/t_period;
        int remain=(aim_shift-all_period*t_period)/(t_period/4);
        //make sure the reman is odd
        if((remain%2)==0){remain+=1;}
        t_shift=aim_shift;
        //contain the half pulse + /4 for the oscillation reach zero
        t_pulse_window=all_period*t_period+remain*(t_period/4);


//define the array of the pulse
        double t_alter=t_pulse_window+t_shift+L_total/3e8+0.5*tau;
        double t_bound;
        if(tau<15e-15){
             t_bound=7.4*tau;}
        else{
             t_bound=7.8*tau;
        }
        if (t_bound<t_alter){t_bound=t_alter;}
        t=Eigen::ArrayXd::LinSpaced((int) 2*t_bound/dt,-t_bound, t_bound); 
        Nt=t.size();
     
        E_vec=Eigen::ArrayXXd::Zero(3, Nt);
        tempE=Eigen::ArrayXd::Zero(Nt);
   
   
        dt=t[1]-t[0];



        t_arr=t+t[Nt-1]-t_shift;   //center of the pulse distance to t[0] is t_shift
       
        //my_print(t_stop);
        //for thz input make sure at 0 frequncy, spectrum goes to zero
        double tau_fwhm=tau;
        double s=(pow(C.omega*tau_fwhm,2)+(C.omega*tau_fwhm)*sqrt(pow(C.omega*tau_fwhm,2)-8*log(2)))/(8*log(2));

        auto field_define = [&] (const Eigen::ArrayXd & t, const double  phi,Eigen::ArrayXd & tempE, double & lambda) {
        if(C.lambda>5e-4){
            tempE= E0*exp(-pow(t,2)/pow(tau,2))*exp(C.i*(C.omega*t+phi)).real();
         }
        else{

            tempE=E0*(pow(1+C.i*C.omega*t/s,-s-1)*exp(C.i*phi)).real();

         }
 
        };

        // field_define(t_arr,phi,tempE,C.lambda);
        // my_print(s);
        // save_array(tempE, "Ein.txt" );


        //select the input driving dimension
        flag_dim=0;
        if (flag_dim==0){
        field_define(t_arr,phi,tempE,C.lambda);
        E_vec.row(1)=1*tempE;

        field_define(t_arr,0*phi,tempE,C.lambda);
        E_vec.row(2)=0*tempE;
        }
        else{
        field_define(t_arr,0.0,tempE,C.lambda);
        E_vec.row(1)=0*tempE;

        field_define(t_arr,0*phi,tempE,C.lambda);
        E_vec.row(2)=1*tempE;
        }


       //--------------------------------------
       //find the minimum value to get the switch of bc
       //--------------------------------------
        left_bc_on=(t_shift+t_pulse_window)/dt;
        int N_win=t_period/dt/2;
        //find minimum E value at left bc on
        double minE=abs(E_vec.block(1,left_bc_on-N_win,1,2*N_win)).minCoeff();
        int iter=left_bc_on-N_win;

        while(abs(E_vec(1,iter))!=minE){
           iter++;
        }
        left_bc_on=iter+1;
        //convert from x to t
        double x_front=left_bc_on*dt*C.c/2;
        t_stop=left_bc_on+(2*L_total/C.c+x_front/C.c+0.3*tau)/dt;

    };
    ~E_field(){};
};


#endif
