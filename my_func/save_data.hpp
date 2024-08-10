#ifndef SAVE_DATA_H
#define SAVE_DATA_H
#include"save_print.hpp"
#include <math.h> 
#include <numeric>
#include <cmath>
#include "my_func/include.hpp"


struct save_data{
    
      std::string save_path;
      Eigen::ArrayXXd E_store,test,testconv;
       Eigen::ArrayXXd  jz,jy;
      //Eigen::Array<bool,4,1000> testconv;
      Eigen::ArrayXd E_z_converge_check,t;
      std::vector<int> Nx_save,thread_iteration;

      int Nt,Nx,save_iter,N_p,N_front;
      
    
      template<typename T1,typename T2>
      save_data(std::string save_path_,T1& E,T2 &mesh_x,int& N_thread){
          save_path=save_path_;
          Nt=E.Nt;
          t=E.t;
          N_front=mesh_x.N_front;
          Nx=mesh_x.Nx;
          Nx_save=mesh_x.Nx_save;
          N_p=mesh_x.Nx_save.size();
          int N_flag=std::accumulate(mesh_x.J_flag.begin(), mesh_x.J_flag.begin()+(Nx-1), 0);
          E_store=Eigen::ArrayXXd::Zero(2*N_p,Nt);

          int assum_conv=400;
          test=Eigen::ArrayXXd::Zero(4,assum_conv);
          testconv=Eigen::ArrayXXd::Zero(4,assum_conv);
          jy=Eigen::ArrayXXd::Zero(N_flag,Nt);
          jz=Eigen::ArrayXXd::Zero(N_flag,Nt);
          //testconv= Eigen::Array<bool, 4,1000>::Zero(4,assum_conv);
          
          thread_iteration=std::vector<int>(N_thread,0);
          E_z_converge_check=Eigen::ArrayXd::Zero(Nx)+1.0;


        save_array(mesh_x.x,save_path+ "x.txt" );
        save_array(E.t,save_path+ "t.txt" );
        save_array(E.E_vec,save_path+ "Eyz_t.txt");
        //save_array(mesh_x.J_flag,save_path+ "J_flag.txt");
        save_array(Nx_save,save_path+"saveNx.txt");
        std::string save_str= "N_front="+std::to_string(mesh_x.N_front)+", N_wsm="+std::to_string(mesh_x.N_wsm)+", N_air="+std::to_string(mesh_x.N_air); 

        save_note(save_path+"save_note.txt",save_str);
        save_str="end numerical time iter="+std::to_string(E.t_stop); 
        save_note(save_path+"save_note.txt",save_str);
          
    };
      ~save_data(){}

    template<typename T1,typename T2>
    void save_location(int t_iter,int iter_x, T1& Ey, T2& Ez);
    template<typename T1,typename T2,typename T3>
    void  save_E(T1& WSM,int& t_iter,T2& mesh_x,T3& C);  


    int getIndex(std::vector<int>& v, int& K); 
    void save_check(int& t_iter); 
    template<typename T1>
    void check_point(T1& WSM,int& t_iter);
    void save_note_struc(int& a,std::string&b);
};

template<typename T1,typename T2>
void save_data::save_location(int t_iter,int n_x, T1& Ey, T2& Ez){
    
    int save_iter=getIndex(Nx_save, n_x) ;
    if (save_iter<N_p) {
    
          E_store(save_iter,t_iter)=Ey;
          E_store(save_iter+3,t_iter)=Ez;
       
     }
      
}


int save_data::getIndex(std::vector<int>& v, int& K) 
{ 
    auto it = find(v.begin(), v.end(), K); 
    int index;
    // If element was found 
    if (it != v.end())  
    { 
      
        index = it - v.begin(); 
     
    } 
    else { 
          index=N_p+1;
    } 
return index;
} 



template<typename T1,typename T2,typename T3>
void save_data::save_E(T1& WSM,int &t_iter,T2 &mesh_x,T3&C){

if (t_iter<(Nt-100)){

   //convert the Ex to Et in vacuum 
   Eigen::ArrayXd Ey_loss=WSM.ME.row(0).segment(Nx_save[0],N_front+1-Nx_save[0]);
   Eigen::ArrayXd Ez_loss=WSM.ME.row(1).segment(Nx_save[0],N_front+1-Nx_save[0]);
   Eigen::ArrayXd t_loss=mesh_x.x.head(N_front+1-Nx_save[0])/C.c;
   double Ey,Ez;
   double t_pre;
   double dt_loss=t_loss(1)-t_loss(0);
   double t_loss_end=t_loss(t_loss.size()-1);
   int t_x_loop=0;
    for(int t_loop=t_iter-1;t_loop<Nt;t_loop++){

        t_pre=t(t_loop)-t(t_iter-1);
        if (t_pre>t_loss_end) {break;}
        if(t_pre>t_loss(t_x_loop+1)){t_x_loop++;}

        Ey=(Ey_loss(t_x_loop)*(t_loss(t_x_loop+1)-t_pre)+Ey_loss(t_x_loop+1)*(t_pre-t_loss(t_x_loop)))/dt_loss; 
        Ez=(Ez_loss(t_x_loop)*(t_loss(t_x_loop+1)-t_pre)+Ez_loss(t_x_loop+1)*(t_pre-t_loss(t_x_loop)))/dt_loss;
        E_store(0,t_loop)=Ey;
        E_store(3,t_loop)=Ez;
      
    }

}
    std::cout<<"save all the files to "+save_path<<std::endl;
    save_array(E_store,save_path+ "Eyz_3location_t.txt");
    //at the end of t it should be zero since the entire pulse has passed
    save_array(WSM.ME,save_path+ "Eyz_x_dis.txt");
    //save_array(thread_iteration,save_path+ "thread_check.txt");

   
}
template<typename T1>
void save_data::check_point(T1& WSM,int& t_iter){
    std::cout<<"save check point iteration"<<t_iter<<  " to "+save_path<<std::endl;
    save_array(E_store,save_path+ "Eyz_3location_t"+std::to_string(t_iter)+".txt");
    //at the end of t it should be zero since the entire pulse has passed
    save_array(WSM.ME,save_path+ "Eyz_x_dis"+std::to_string(t_iter)+".txt");
    // save_array(jy,save_path+ "jyt"+std::to_string(t_iter)+".txt");
    // save_array(jz,save_path+ "jzt"+std::to_string(t_iter)+".txt");
  //  save_array(thread_iteration,save_path+ "thread_check.txt");
   // save_array(E_z_converge_check,save_path+"saveConv"+".txt");
   
}

void save_data::save_check(int &t_iter){
    std::cout<<"save all the files to "+save_path<<std::endl;

   // save_array(thread_iteration,save_path+ "thread_check.txt");
    save_array(test,save_path+ "test"+std::to_string(t_iter)+".txt");
      save_array(testconv,save_path+ "testconv"+std::to_string(t_iter)+".txt");


   
}
void save_data::save_note_struc(int& a,std::string& b){
    save_note(a,save_path+"save_note.txt",b);
}
#endif
