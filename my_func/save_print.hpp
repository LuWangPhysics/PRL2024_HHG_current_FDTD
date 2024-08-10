//all my print functions
#ifndef SAVE_PRINT_HPP
#define SAVE_PRINT_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <sstream>
#include <vector>

//print which thread with how many iterations 
void my_print(const std::vector<int>& iterations)
{
    int index = 0;
    for (const auto& i : iterations)
    {
        std::cout << "Thread " << index++ << " -> " << i << " iterations" << std::endl;
    }
}
void my_print(std::vector<double>& iterations)
{
    int index = 0;
    for (const auto& i : iterations)
    {
        std::cout <<  i << std::endl;
        index++ ;
    }
}

void my_print(const Eigen::ArrayXd & a)
{
          int N_f=a.size();
       
         
            for (int i=0;i<N_f;i++){
              std::cout<<a(i)<<"\n "<<std::endl; }
}
template<typename T>
void my_print(const T& x)
{std::cout<<x<<std::endl;}

template<typename T>
void test_line(const T& n){
    std::cout<<"****************************************************************"<<std::endl;
     std::cout<<n<<std::endl;
     std::cout<<"****************************************************************"<<std::endl;
}

void test_line(std::string n){
    std::cout<<"****************************************************************"<<std::endl;
    std::cout<<n<<std::endl;
    std::cout<<"****************************************************************"<<std::endl;
}



//---------------------------------------------------------------------------------------------------------
//my save functions
//---------------------------------------------------------------------------------------------------------

      
        void save_note(const Eigen::ArrayXd& a,const std::string & title,const std::string& message  )
        {         std::ofstream myfile;

                myfile.open(title, std::ios_base::app); // append instead of overwrite
                myfile <<message; 
               
                int N_f=a.size();

                for (int i=0;i<N_f;i++){
                    myfile <<a[i]<<"\n "; }
                myfile.close();

        }
      
        void save_note(const std::string & title,const  std::string & message  )
        {         std::ofstream myfile;

                myfile.open(title, std::ios_base::app); // append instead of overwrite
                myfile <<message<<"\n "; 

                myfile.close();

        }
       template<typename T1>
         void save_note(T1& i_iter,const std::string & title, const std::string& message  )
        {         std::ofstream myfile;

                myfile.open(title, std::ios_base::app); // append instead of overwrite
                myfile <<message<<"\n "; 
                myfile <<i_iter<<"\n ";
                myfile.close();

        }


        void save_array(Eigen::ArrayXcd & a, std::string & b  )
        {
            int N_f=a.size();
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::showpoint);
	        myfile.setf(std::ios::scientific);
            for (int i=0;i<N_f;i++){
                myfile <<a[i].real()<<","<< a[i].imag()<<"\n "; }
            myfile.close(); 
      

        }

template<typename T1>
  void save_array(const T1 & a, std::string b  )
  {
      int N_f=a.size();
      std::ofstream myfile(b);
      myfile.precision(15);
      myfile.setf(std::ios::showpoint);
      myfile.setf(std::ios::scientific);
      for (int i=0;i<N_f;i++){
          myfile <<a[i]<<"\n "; }
      myfile.close();


  }
          
  void save_array(const Eigen::ArrayXd & a, std::string b  )
  {
      int N_f=a.size();
      std::ofstream myfile(b);
      myfile.precision(15);
      myfile.setf(std::ios::showpoint);
      myfile.setf(std::ios::scientific);
      for (int i=0;i<N_f;i++){
          myfile <<a[i]<<"\n "; }
      myfile.close();


  }
  
    void save_array(std::vector<int> & a, std::string b  )
  {
      int N_f=a.size();
      std::ofstream myfile(b);

      for (int i=0;i<N_f;i++){
          myfile <<a[i]<<"\n "; }
      myfile.close();


  }

        void save_array(int N_f,double* a,std::string b )
        {
         

            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);
            for (int i=0;i<N_f;i++){
                myfile <<a[i]<<"\n ";
            }
            myfile.close();
            
        }
        void save_array( Eigen::MatrixXcd& a, std::string b )
        {   int N_f=a.rows();
            int N_x=a.cols();
        
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);

                for (int i=0;i<N_f;i++)
                {
                    for (int j=0;j<N_x;j++)
                    {
                        myfile<< a(i,j).real()<<","<<a(i,j).imag() << ","; 
                    }
                  myfile << "\n";
                }
            myfile.close();
            
        }
        
        
     void save_array( Eigen::ArrayXXcd& a, std::string b )
        {   int N_f=a.rows();
            int N_x=a.cols();
        
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);

                for (int i=0;i<N_f;i++)
                {
                    for (int j=0;j<N_x;j++)
                    {
                        myfile<< a(i,j).real()<<","<<a(i,j).imag() << ","; 
                    }
                  myfile << "\n";
                }
            myfile.close();
            
        }

        
         void save_array( Eigen::ArrayXXd& a, std::string b )
        {   int N_f=a.rows();
            int N_x=a.cols();
        
            std::ofstream myfile(b);
            myfile.precision(15);
            myfile.setf(std::ios::scientific);
            myfile.setf(std::ios::showpoint);

                for (int i=0;i<N_f;i++)
                {
                    for (int j=0;j<N_x;j++)
                    {
                        myfile<< a(i,j)<< ","; 
                    }
                  myfile << "\n";
                }
            myfile.close();
            
        }
        
#endif
