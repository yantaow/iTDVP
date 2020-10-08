#ifndef IMPS_HPP_INCLUEDED
#define IMPS_HPP_INCLUEDED 1

#include"itensor/all.h"
#include<vector>
#include<iostream>
#include<random>
#include "util.h"

using namespace itensor;
using namespace std; 
typedef ITensor T; 
/*------------------------------------------------------------*/
class iMPS{
  public: 
    iMPS(); 
    iMPS(T A, Index il, Index ir, Index is, T W, Index wl, Index wr, string s); 

    int m_; 
    T A_, AL_, AR_, C_, AC_, lw_, rw_, L_, R_;   
    T XR_, XL_; //guess vectors when computing lw and rw
    Index il_, ir_; 
    Index is_;
    T W_; 
    Index wl_, wr_; 
    double entropy_, energy_; 
    ofstream out_fs_; 

    void init_A(); 
    void set_canonical(); 
    void set_MPO(T W, Index wl, Index wr); 
    void mixed_gauge(); 
    void iTDVP(Cplx dt); 
    double get_energy(bool is_print); 
    double get_entropy(bool is_print); 
    double check_canonical(bool is_print); 
    T left_normal(double err_goal); //compute AL from A 
    T right_normal(double err_goal); //compute AR from A  
    T get_l();
    T get_r(); 
    T get_lw(); 
    T get_rw(); 
    T split_AL(); 
    T split_AR();   

    void to_screen(int i, double t); 
    void to_file(int i, double t); 
}; 

#endif



