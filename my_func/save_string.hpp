#include <sys/stat.h>


template<typename T1>
std::string make_string(double & E_02,double& tau,double & lambda,double & L_wsm,double &L_air,T1 & C,double& dx){

std::stringstream E_text,s_lambda;
std::string s_efield2;
if(E_02>1e5){
E_text << std::fixed << std::setprecision(2) <<E_02*1e-5 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text.str()+"E5V_per_m";}else{
E_text << std::fixed << std::setprecision(2) <<E_02*1e-3 ;                                  // 2 digit significant number of the E_01
s_efield2 = E_text.str()+"E3V_per_m";
}




s_lambda<<std::fixed << std::setprecision(2) <<lambda*1e6;
std::string extra_name="E"+s_efield2+"tau"+std::to_string(int (tau*1e15) )+"fs_"+s_lambda.str()+"um_Lwsm";
extra_name+=std::to_string(int(L_wsm*1e9))+"nm_Lair"+std::to_string(int (L_air*1e9))+"nm_dx"+std::to_string(int (dx*1e9) )+"nm";

std::string save_path="my_output/"+C.name_tail+extra_name+"/";
if (mkdir(save_path.c_str(),0777)!= 0){mkdir(save_path.c_str(),0777);}


    return save_path;
}
