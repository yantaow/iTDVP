#include "itensor/all.h"
#include <vector>
#include "util.h"

using namespace itensor; 
using namespace arma; 
using namespace std;
typedef ITensor T; 

extern complex<double> _i; 
/*----------------------------------------------*/
ITensor subTensor(ITensor T, vector<Index> indices, vector<int> values)
{
  ITensor retT;
  for(unsigned l=0; l<indices.size();l++)
  {
    if(l==0)
      retT = T*setElt(indices[l]=values[l]);
    else
      retT*=setElt(indices[l]=values[l]);
  }
  return retT;
}
/*----------------------------------------------*/
vector<vector<ITensor> > extract_Tw(T& A, T& W, Index a, Index b) 
{
  Index is = findIndex(A, "Site"); 
  int dw = dim(a); 
  vector<vector<ITensor> > Tw(dw, vector<ITensor>(dw));  
  for(int i = 0; i < dw; i++)
  {
    for(int j = 0; j < dw; j++)
    {
      T Wab = subTensor(W, {a, b}, {i+1, j+1}); 
      Tw[i][j] = Wab * A * prime(dag(A)); 
    }
  }

  return Tw; 
}
/*----------------------------------------------*/
void print_MPO(T W) 
{
  auto ia = findInds(W, "Link")(1);  
  auto ib = findInds(W, "Link")(2);  
  auto is = findInds(W, "Site")(1); 
  PrintData(ia);
  PrintData(ib); 
  PrintData(is); 
  cout << "The MPO tensor W is: " << endl; 
  int width = 12; 
  for(int i = 1; i <= dim(ia); i++)
  {
    for(int s = 1; s <= dim(is); s++)
    {
      for(int j = 1; j <= dim(ib); j++)
      {
        T Wab = subTensor(W, {ia, ib}, {i, j}); 
        for(int ss = 1; ss <= dim(is); ss++)
        {
          cout.width(width); 
          std::ostringstream out;  
          out.precision(2); 
          if(abs(imag((eltC(Wab, is(s), prime(is)(ss))))) < 1E-5)
          {
            double num = real(eltC(Wab, is(s), prime(is)(ss))); 
            if(abs(num-int(num)) < 1E-5)
            {
              out.precision(0); 
              out << std::fixed << num; 
              out.precision(2); 
            }
            else 
            {
              out << std::fixed << num; 
            }
          }
          else 
          {
            Cplx num = (eltC(Wab, is(s), prime(is)(ss))); 
            out << std::fixed << real(num) << std::showpos << imag(num) << "i"; 
          }
          
          string str; 
          if(ss == 1)
            str = "|" + out.str(); 
          else 
            str = out.str(); 

          cout << left << str;   
        }
      }
      cout << endl; 
    }
    for(int k = 1; k <= dim(ib)*dim(is)*width; k++)
    {
      cout << "-"; 
    }
    cout << endl; 
  }

  return; 
}

