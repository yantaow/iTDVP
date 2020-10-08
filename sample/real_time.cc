#include "itensor/all.h"
#include <chrono>
#include <random>
#include <omp.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "util.h"
#include "imps.h"

using namespace itensor;
typedef ITensor T; 
Cplx _i = Cplx(0, 1); 
/*----------------------------------*/
void init_MPS(SiteSet& sites, T& A, Index& il, Index& ir, Index& is, int m);   
MPS expand_dim(MPS psi, int m);  
void init_MPO(SiteSet& sites, T& W, Index& wl, Index& wr, double gx);   
void schur(T& W);  
/*----------------------------------*/
void take_input(int argc, char *argv[]); 
Cplx dt0, dt1; 
int m; //bond dimension  
int n_sw0, n_sw1;    
double J, g0, g1; 
/*----------------------------------*/
int main(int argc, char *argv[])
{
  take_input(argc, argv); 
  //initialize an open chain with 3 sites. 
  //take the middle one as our system
  SiteSet sites = SpinHalf(3, {"ConserveQNs=", false}); 
  /*-------initialize the system-------*/
  T A; 
  Index il, ir, is;
  init_MPS(sites, A, il, ir, is, m);   
  T W; 
  Index wl, wr, wl_Z, wr_Z; 
  init_MPO(sites, W, wl, wr, g0);   
  iMPS iA = iMPS(A, il, ir, is, W, wl, wr, "A"); 
  /*--------------*/
  Cplx t = 0; 
  for(int i = 0; i < n_sw0; i++, t+=abs(dt0))
  {
    bool is_screen = n_sw0 < 10 ? true : i%(n_sw0/10) == 0; 
    if(is_screen)   
    {
      cout << "---------search for ground state-----------" << endl; 
      cout << "i: " << i << "\t" << "t: " << abs(t) << endl; 
      iA.to_screen(i, abs(t)); 
    }
    iA.iTDVP(-dt0); 
  }

  if(abs(g0-1.5) < 10E-5)
  {
    cout << "Exact energy for gx 1.5:" << -1.671926221536197 << endl;  
  }
  else if(abs(g0-1.0) < 10E-5)
  {
    cout << "Exact energy for gx 1.0:" << -1.273239544735164 << endl;  
  }
  cout << "---------------------------------------------" << endl; 

  init_MPO(sites, W, wl, wr, g1);   
  iA.set_MPO(W, wl, wr); 
  t = 0; 
  for(int i = 0; i < n_sw1; i++, t+=abs(dt1))
  {
    bool is_screen = n_sw1 < 10 ? true : i%(n_sw1/10) == 0; 
    if(is_screen)   
    {
      cout << "---------post-quench-----------" << endl; 
      cout << "i: " << i << "\t" << "t: " << abs(t) << endl; 
      iA.to_screen(i, abs(t)); 
    }
    iA.to_file(i, abs(t));  
    iA.iTDVP(-dt1*_i); 
  }

  return 0;
}
/*----------------------------------------------*/
void init_MPS(SiteSet& sites, T& A, Index& il, Index& ir, Index& is, int m)   
  //initialize an open chain with 3 sites. 
  //take the middle one as our system
{
  auto state = InitState(sites);
  state.set(1,"Up");
  state.set(2,"Up");
  state.set(3,"Up"); 

  auto psi = MPS(state);
  psi = expand_dim(psi, m); 
  A = psi(2); 
  il = commonIndex(A, psi(1));  
  ir = commonIndex(A, psi(3));  
  is = findIndex(A, "Site"); 

  return; 
}
/*----------------------------------------------*/
void init_MPO(SiteSet& sites, T& W, Index& wl, Index& wr, double gx)   
{
  AutoMPO ampo = AutoMPO(sites); 
  //initialize an open chain with 3 sites. 
  //take the middle one as our system
  for(int j = 1; j < 3; ++j)
  {
    ampo += 4*J,"Sz",j,"Sz",j+1;
    ampo += 2*gx,"Sx",j;
  }
  ampo += 2*gx,"Sx",3;
  
  auto H = toMPO(ampo); 
  W = H.ref(2); 
  schur(W); 
  wl = commonIndex(W, H(1)); 
  wr = commonIndex(W, H(3));
  /*------long range TFIM-------*/
  print_MPO(W); 

  return; 
}
/*----------------------------------------------*/
MPS expand_dim(MPS psi, int m) 
//expand the bond dimension of psi to a new MPS while keep all the tags
{
  int l = length(psi); 
  vector<Index> links_old(l+1); 
  links_old[0] = findIndex(uniqueInds(psi.ref(1).inds(), psi.ref(2).inds()), "Link");  
  links_old[l] = findIndex(uniqueInds(psi.ref(l).inds(), psi.ref(l-1).inds()), "Link");  

  for(int j = 1; j < l; j++)
  {
    links_old[j] = commonIndex(psi.ref(j), psi.ref(j+1)); 
  }
  
  vector<Index> links_new(l+1); 
  for(int j = 0; j <= l; j++)
  {
    auto i = Index(m); 
    i.setTags(links_old[j].tags()); 
    links_new[j] = i; 
  }

  for(int j = 1; j <= l; j++)
  {
    auto is = findIndex(psi.ref(j).inds(), "Site");
    T A; 
//        A = randomITensorC(links_new[j-1], is, links_new[j]); 
    A = T(links_new[j-1], is, links_new[j]); 
    A.set(1, 1, 1, 1.0);      
    psi.ref(j) = A; 
  }

  return psi; 
}
/*----------------------------------------------*/
void schur(T& W) 
//AutoMPO doesn't produce MPO in Schur form
//we have to manualy put W in Schur form....
//the code is written for TFIM
//The Index id does not change
{
  int i1 = 2; 
  int i2 = 3; 
  Index a = W.inds().index(1);  
  Index b = W.inds().index(2); 
  /*---------------------------------*/
  T swap_a = T(a, prime(a)); 
  swap_a.set(a(i1), prime(a)(i2), 1); 
  swap_a.set(a(i2), prime(a)(i1), 1); 
  for(int i = 1; i <= dim(a); i++)
  {
    if(i != i1 && i != i2)
    {
      swap_a.set(a(i), prime(a)(i), 1); 
    }
  }
  /*---------------------------------*/
  T swap_b = T(b, prime(b)); 
  swap_b.set(b(i1), prime(b)(i2), 1); 
  swap_b.set(b(i2), prime(b)(i1), 1); 
  for(int i = 1; i <= dim(b); i++)
  {
    if(i != i1 && i != i2)
    {
      swap_b.set(b(i), prime(b)(i), 1); 
    }
  }
  /*---------------------------------*/
  W = (W * swap_a) * swap_b; 
  W.noPrime(prime(a)); 
  W.noPrime(prime(b)); 

  return; 
}
/*----------------------------------------------*/
void take_input(int argc, char *argv[])
{
  int i = 1;
  m      =  stod(argv[i++]); 
  dt0    =  stod(argv[i++]); 
  dt1    =  stod(argv[i++]); 
  n_sw0  =  atoi(argv[i++]); 
  n_sw1  =  atoi(argv[i++]); 
  J      =  stod(argv[i++]); 
  g0     =  stod(argv[i++]); 
  g1     =  stod(argv[i++]); 


  return; 
}

