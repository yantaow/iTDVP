#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED 1

#include "itensor/all.h"
//  #define ARMA_DONT_USE_CXX11
#include <armadillo>
#include <vector>

using namespace itensor; 
using namespace arma; 
using namespace std;
typedef ITensor T; 

/*-------------------------------------------------------*/
ITensor subTensor(ITensor T,std::vector<Index> indices,std::vector<int> values); 
vector<vector<ITensor> > extract_Tw(T& A, T& W, Index a, Index b);  
void print_MPO(T W);  
/*-------------------------------------------------------*/
class BigM_L{
  public: 
    ITensor E_; 
    ITensor r_; 
    Index il_; 
    Index ir_; 
    ITensor Il_; 
    bool isR_; 
    double lambda_; 

    BigM_L(ITensor E, Index i_l, Index i_r, ITensor r, double lambda, bool isR)
      :E_(E), il_(i_l), ir_(i_r), r_(r), lambda_(lambda), isR_(isR)
    {
      Il_ = ITensor(il_, prime(il_)); 
      for(int i = 1; i <= dim(il_); i++)
      {
        Il_.set(il_(i), prime(il_)(i), 1); 
      }
    }

    int size() const
    {
//        return dim(il_) * dim(il_) * dim(il_) * dim(il_); 
      return dim(il_) * dim(il_); 
    }

    void product(ITensor& x, ITensor& Ax) const
    {
      ITensor L1, L2, L3; 
      L1 = x;  
      L2 = (x*lambda_*E_) * delta(il_, ir_) * delta(prime(il_), prime(ir_)); 
      if(!isR_)
      {
        Ax = L1 - L2; 
      }
      else
      {
        L3 = x * delta(il_, ir_) * delta(prime(il_), prime(ir_)) * r_ * Il_; 
        Ax = L1 - L2 + L3; 
      }

      return; 
    }
  };
/*-------------------------------------------------------*/
class BigM_R{
  public: 
    ITensor E_; 
    ITensor l_; 
    ITensor Ir_; 
    Index il_; 
    Index ir_; 
    bool isL_; 
    double lambda_; 

    BigM_R(ITensor E, Index i_l, Index i_r, ITensor l, double lambda, bool isL)
      :E_(E), il_(i_l), ir_(i_r), l_(l), lambda_(lambda), isL_(isL)
    {
      Ir_ = ITensor(ir_, prime(ir_)); 
      for(int i = 1; i <= dim(ir_); i++)
      {
        Ir_.set(ir_(i), prime(ir_)(i), 1); 
      }
    }

    int size() const
    {
//        return dim(ir_) * dim(ir_) * dim(ir_) * dim(ir_); 
      return dim(ir_) * dim(ir_); 
    }

    void product(ITensor& x, ITensor& Ax) const
    {
      ITensor L1, L2, L3; 
      L1 = x;  
      L2 = (lambda_*E_*x) * delta(il_, ir_) * delta(prime(il_), prime(ir_)); 
      if(!isL_)
      {
        Ax = L1 - L2; 
      }
      else 
      {
        L3 = l_ * delta(il_, ir_) * delta(prime(il_), prime(ir_)) * x * Ir_; 
        Ax = L1 - L2 + L3; 
      }

      return; 
    }
  };
/*-------------------------------------------------------*/
class K_op{
  public: 
    ITensor lw_; 
    ITensor rw_; 
    Index il_, ir_; 

    K_op(ITensor lw, ITensor rw, Index il, Index ir)
      :lw_(lw), rw_(rw), il_(il), ir_(ir)
    {
    }

    int size() const
    {
      return dim(il_) * dim(ir_); 
    }

    void product(ITensor& C, ITensor& K_C) const
    {
      K_C = C * lw_ * rw_; 
      K_C.noPrime(); 
      return; 
    }
  };
/*-------------------------------------------------------*/
class BigMatrix{
  public: 

    Index i_open_; 
    Index i_contact_; 
    T E_; 

    BigMatrix(T& E, Index i_open, Index i_contact)
      :E_(E), i_open_(i_open), i_contact_(i_contact)
    {
      /*nothing*/
    }

    int size() const
    {
      return dim(i_open_); 
    }

    void product(ITensor& x, ITensor& Ex) const
      //i_contact-x
    {
      T x_tmp = x * delta(i_contact_, prime(i_contact_)); 
      Ex = (E_ * delta(i_contact_, prime(i_contact_))) * x_tmp; //iopen-Ex 
      Ex = Ex * delta(i_open_, i_contact_); 
      return; 
    }
  };
/*------------------------------------------------------------*/
#endif

