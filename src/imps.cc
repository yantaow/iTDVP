//The code is based on PHYSICAL REVIEW B 97, 045145 (2018)
//and SciPost Phys. Lect. Notes 7 (2019), and follows their
//notation. One should cite these papers when using this code. 

# include<vector>
# include<iostream>
# include<iomanip>
# include<string>
# include "imps.h"

using namespace std; 
double gmres_err = 10E-15, gauge_err = 10E-15; 
/*------------------------------------------------------------*/
T identity(Index a, Index b)
{
  T id(a, b); 
  for(int i = 1; i < dim(a); i++)
  {
    id.set(a(i), b(i), 1); 
  }
  return id; 
}
/*----------------------------------------------*/
iMPS::iMPS()
{
  cout << "Hello World!" << endl; 
}
/*----------------------------------------------*/
iMPS::iMPS(T A, Index il, Index ir, Index is, T W, Index wl, Index wr, string s)
  :m_(dim(il)), A_(A), il_(il), ir_(ir), is_(is), W_(W), wl_(wl), wr_(wr)
{
  A_ = T(il_, is_, ir_); 
  A_.set(1, 1, 1, 1.0);      
//    A = randomITensorC_(il_, is_, ir_); 
  XR_ = randomITensor(ir_, prime(ir_)); 
  XL_ = randomITensor(il_, prime(il_)); //to be used in gmres 
  set_canonical(); 
  lw_ = get_lw(); 
  rw_ = get_rw(); 
  get_energy(false); 
  get_entropy(false); 

  out_fs_.open("imps_"+s+".out"); 
}
/*----------------------------------------------*/
void iMPS::set_MPO(T W, Index wl, Index wr)
{
  W_ = W; 
  wl_ = wl; 
  wr_ = wr; 
  
  lw_ = get_lw(); 
  rw_ = get_rw(); 
  
  return; 
}
/*----------------------------------------------*/
void iMPS::set_canonical()
{
  mixed_gauge(); 
  //il-AL-il', il'-C_-ir', ir'-AR-ir
  /*----------------------------*/
  AC_ = AL_ * delta(ir_, prime(ir_)) * (C_ * delta(il_, prime(ir_))); // il-AC-ir 
  return; 
}
/*----------------------------------------------*/
void iMPS::mixed_gauge()
  //Algorithm 2 in SciPost Phys. Lect. Notes 7 (2019)
  //input: il-A-ir; 
  //return: il-AL-ir, il-C-ir, il-AR-ir, il'-L-il, ir-R-ir'
{
  L_ = identity(il_, prime(il_)); 
  AL_ =  left_normal(gauge_err);  
  R_ = identity(ir_, prime(ir_)); 
  AR_ = right_normal(gauge_err); 

  C_ = L_ * R_ * delta(il_, ir_);  
  C_ = C_ * delta(il_, prime(il_)) * delta(ir_, prime(ir_)); 
  C_ = C_/norm(C_); 
  return; 
}
/*----------------------------------------------*/
void iMPS::iTDVP(Cplx dt) 
  //"Integrating the TDVP equations" on page 35 of SciPost Phys. Lect. Notes 7 (2019)
  //input: il, ir, il-AC-ir, il-C-ir, il-AL-ir, il-AR-ir
  //output: il-AC'-ir, il-C'-ir, il-AL'-ir, il-AR'-ir
{
  LocalOp Heff(W_, lw_, rw_, {"numCenter=",1}); 
  applyExp(Heff,AC_,dt);
  
  lw_ *= delta(wl_, wr_); 
  K_op Keff(lw_, rw_, il_, ir_); //K_op is defined in util.h
  applyExp(Keff,C_,dt); 

  AL_ = split_AL(); //il_-AL-ir_
  AR_ = split_AR(); //il_-AR-ir_ 
  C_ /= norm(C_); 
  AC_ /= norm(AC_); 

  lw_ = get_lw(); 
  rw_ = get_rw(); 

  return; 
}
/*----------------------------------------------*/
double iMPS::get_energy(bool is_print) 
{
  T Ew = W_ * AL_ * prime(dag(AL_), il_, ir_, is_); 
  T lw_1 = lw_ * Ew * delta(il_, ir_) * delta(prime(il_), prime(ir_)) * delta(wl_, wr_); 
  
  Ew = W_ * AR_ * prime(dag(AR_), il_, ir_, is_); 
  T rw_1 = rw_ * Ew * delta(il_, ir_) * delta(prime(il_), prime(ir_)) * delta(wl_, wr_); 
  if(is_print)
  {
    PrintData(lw_1 - lw_); 
    PrintData(rw_1 - rw_); 
  }
  int dw = dim(findInds(W_, "Link")(1)); 
  double energy_l = real(eltC(lw_1 - lw_, il_(1), prime(il_)(1), wl_(1)));  
  double energy_r = real(eltC(rw_1 - rw_, ir_(1), prime(ir_)(1), wr_(dw)));  

  energy_ = energy_l; 
  return energy_l; 
}
/*----------------------------------------------*/
double iMPS::get_entropy(bool is_print)
{
  T U(il_),S,V;        
  auto spec = svd(C_, U, S, V); 
  if(is_print)
  {
    cout << "itdvp spectrum: " << endl; 
    PrintData(S); 
  }
  auto eigs = spec.eigsKept(); 
  double entropy = 0.0; 
  for(auto& p: eigs)
  {
    if(p > 1E-13)
    {
      entropy += -p * log(p); 
    }
  }
  entropy_ = entropy; 

  return entropy; 
}
/*----------------------------------------------*/
double iMPS::check_canonical(bool is_print) 
{
  Index ii(dim(il_)); 
  T AL_C = AL_ * delta(ir_, ii) * (C_ * delta(il_, ii)); 
  T C_AR = C_ * delta(ir_, ii) * (AR_ * delta(il_, ii)); 
  if(is_print)
  {
    cout << "=============check_canonical===========" << endl; 
    PrintData(AL_ * prime(dag(AL_), ir_));   
    PrintData(AR_ * prime(dag(AR_), il_));   
    PrintData(AC_ - AL_C);  
    PrintData(AC_ - C_AR); 
    cout << "=============done checking  ===========" << endl; 
  }
  return norm(AC_ - AL_C); 
}
/*----------------------------------------------*/
T iMPS::left_normal(double err_goal) 
  //Algorithm 1 in SciPost Phys. Lect. Notes 7 (2019)
  //input: il_'-L_guess-il_, il_-A-ir_  
  //return: il_-AL-ir_, il_'-L-il_
{
  L_ /= norm(L_); 
  T L_old = L_; 

  T L_A = (L_ * delta(il_, prime(ir_))) * (A_ * delta(il_, prime(ir_))); //il_'-L_A-ir_ 
  T AL = T(is_, prime(il_)); 

  auto [AL_temp, L_temp] = qr(L_A, AL.inds()); 
  AL = AL_temp;
  L_ = L_temp; 

  Index iQR = findIndex(L_, "QR"); 
  AL = AL * delta(prime(il_), il_) * delta(iQR, ir_); 
  L_ = L_ * delta(iQR, prime(il_)) * delta(ir_, il_); //il_'-L-il_ 
  L_ /= norm(L_);
  Real err = norm(L_-L_old); 
  int i = 1;  
  while(err > err_goal)
  {
    T E = A_ * prime(dag(AL), il_, ir_);  

    auto [CC, ic] = combiner(prime(il_), il_);//contact 
    auto [O, io] = combiner(prime(ir_), ir_);//open 
    auto E_comb = E * CC * O; 

    auto E_big = BigMatrix(E_comb, io, ic);  
    T guess = L_ * CC; //guess-ic 
    arnoldi(E_big, guess, {"ErrGoal=", err/10.});  
    L_ = guess * CC; 
    L_ /= norm(L_); 
    L_old = L_; 

    L_A = (L_ * delta(il_, prime(ir_))) * (A_ * delta(il_, prime(ir_))); //il_'-L_A-ir_ 
    AL = T(is_, prime(il_)); 


    auto [AL_temp, L_temp] = qr(L_A, AL.inds()); 
    AL = AL_temp;
    L_ = L_temp; 

    iQR = findIndex(L_, "QR"); 
    AL = AL * delta(prime(il_), il_) * delta(iQR, ir_); 
    L_ = L_ * delta(iQR, prime(il_)) * delta(ir_, il_); //il_'-L-il_ 
    L_ /= norm(L_);
    err = norm(L_-L_old); 
    i++; 
  }

  return AL; 
}
/*----------------------------------------------*/
T iMPS::right_normal(double err_goal)
  //right version of Algorithm 1 in SciPost Phys. Lect. Notes 7 (2019)
  //input: il_-A-ir_, ir_-R_guess-ir_' 
  //return: il_-AR-ir_, ir_-R-ir_'
{
  R_ /= norm(R_); 
  T R_old = R_; 

  T A_R = (A_ * delta(ir_, prime(il_))) * (R_ * delta(prime(il_), ir_)); //il_-A_R-ir_' 
  T AR = T(is_, prime(ir_)); 

  auto [AR_temp, R_temp] = qr(A_R, AR.inds()); 
  AR = AR_temp;
  R_ = R_temp; 

  Index iQR = findIndex(R_, "QR"); 
  AR = AR * delta(prime(ir_), ir_) * delta(iQR, il_); 
  R_ = R_ * delta(iQR, prime(ir_)) * delta(il_, ir_); //ir_-R-ir_' 
  R_ /= norm(R_);
  Real err = norm(R_-R_old); 
  while(err > err_goal)
  {
    T E = A_ * prime(dag(AR), ir_, il_);  

    auto [CC, ic] = combiner(prime(ir_), ir_);//contact 
    auto [O, io] = combiner(prime(il_), il_);//open 
    auto E_comb = E * CC * O; 

    auto E_big = BigMatrix(E_comb, io, ic);  
    T guess = R_ * CC; //guess-ic 
    arnoldi(E_big, guess, {"ErrGoal=", err/10.});  
    R_ = guess * CC; 
    R_ /= norm(R_); 
    R_old = R_; 

    A_R = (A_ * delta(ir_, prime(il_))) * (R_ * delta(prime(il_), ir_)); //il_-A_R-ir_' 
    AR = T(is_, prime(ir_)); 
    auto [AR_temp, R_temp] = qr(A_R, AR.inds()); 
    AR = AR_temp;
    R_ = R_temp; 

    iQR = findIndex(R_, "QR"); 
    AR = AR * delta(prime(ir_), ir_) * delta(iQR, il_); 
    R_ = R_ * delta(iQR, prime(ir_)) * delta(il_, ir_); //ir_-R-ir_' 
    R_ /= norm(R_);
    err = norm(R_-R_old); 
  }

  return AR; 
}
/*----------------------------------------------*/
T iMPS::get_l()
  //return the leading left eigenvector of AR_*dag(AR_) 
{
  T r = identity(ir_, prime(ir_)); 
  auto l_AR = C_ * delta(il_, prime(il_)) * (dag(C_) * delta(il_, prime(il_)) * delta(ir_, prime(ir_)));//ir_-l-ir_' 
  l_AR = l_AR * delta(ir_, il_) * delta(prime(ir_), prime(il_)); 
  Cplx tr = eltC((r * l_AR * delta(ir_, il_)) * delta(prime(il_), prime(ir_)));  
  return l_AR/tr; 
}
/*----------------------------------------------*/
T iMPS::get_r() 
  //return the leading right eigenvector of AL_*dag(AL_) 
{
  T l = identity(il_, prime(il_)); 
  auto r_AL = C_ * delta(ir_, prime(ir_)) * (dag(C_) * delta(ir_, prime(ir_)) * delta(il_, prime(il_))); //il_-r-il_' 
  r_AL = r_AL * delta(il_, ir_) * delta(prime(il_), prime(ir_)); 
  Cplx tr = eltC((l * r_AL * delta(ir_, il_)) * delta(prime(il_), prime(ir_)));  
  return r_AL/tr; 
}
/*----------------------------------------------*/
T iMPS::get_lw()
  //Algorithm 6 in PRB 97, 045145 (2018)
  //return the left leading generalized eigenvector of EW_AL
{
  T E = AL_ * prime(dag(AL_), il_, ir_);  
  Index is_ = commonIndex(W_, AL_); 
  T l = ITensor(il_, prime(il_)); 
  for(int i = 1; i <= dim(il_); i++)
  {
    l.set(il_(i), prime(il_)(i), 1); 
  }
  T r_AL = get_r(); 
  /*----------------------*/
  int dw = dim(wl_); 
	vector<vector<ITensor> > Tw = extract_Tw(AL_, W_, wl_, wr_);   

  vector<ITensor> Lw(dw, ITensor(il_, prime(il_)));  
	for(int a = dw; a > 0; a--)
  {
    if(a == dw) 
    {
      Lw[a-1] = l; 
    }
    else
    {
      ITensor Ya(ir_, prime(ir_)); 

      for(int b = dw; b > a; b--)
      {
        Ya += Lw[b-1] * Tw[b-1][a-1]; 
      }
      Ya *= delta(il_, ir_);
      Ya *= delta(prime(il_), prime(ir_)); 

      if(a > 1)
      {
        //lambda = 0 for nn-TFIM
        double lambda = 0; 

        //solve (L1|(1- lambda*E) =  (Ya| 
        T guess = Ya; 
        BigM_L M(E, il_, ir_, r_AL, lambda, /*isR*/false); 
				gmres(M, Ya, guess, {"ErrGoal=", gmres_err});	

        Lw[a-1] = guess;  
      }
      else if(a == 1)
      {
        //solve (L1|(1- E +|R)(1|) =  (Ya| - (Ya|R)(1|
        T Y1 = Ya; 
        T Y2 = (Ya * r_AL * delta(ir_, il_)) * delta(prime(il_), prime(ir_)) * l;  
        T Y = Y1-Y2; 
        
				//solve M*X = Y for X
        //XL_ is_ the guess before, and the solution after
        BigM_L M(E, il_, ir_, r_AL, /*lambda*/1.0, /*isR*/true); 
				gmres(M, Y, XL_, {"ErrGoal=", gmres_err});	

        Lw[a-1] = XL_;   
      }
      else 
      {
        cerr << "wrong" << endl; 
      }
    }
  }

  T lw(il_, wl_, prime(il_)); 
  for(int i = 1; i <= dim(il_); i++)
  {
    for(int j = 1; j <= dim(prime(il_)); j++)
    {
      for(int a = 1; a <= dim(wl_); a++)
      {
        lw.set(il_(i), prime(il_)(j), wl_(a), eltC(Lw[a-1], il_(i), prime(il_)(j))); 
      }
    }
  }

  return lw; 
}
/*----------------------------------------------*/
T iMPS::get_rw() 
  //return the right leading generalized eigenvector of EW_AL
{
  T E = AR_ * prime(dag(AR_), il_, ir_);  
  T r = ITensor(ir_, prime(ir_)); 
  for(int i = 1; i <= dim(ir_); i++)
  {
    r.set(ir_(i), prime(ir_)(i), 1); 
  }
  T l_AR = get_l(); 
  /*-------------------------------*/
  int dw = dim(wr_); 
	vector<vector<ITensor> > Tw = extract_Tw(AR_, W_, wl_, wr_);  
  vector<ITensor> Rw(dw, ITensor(ir_, prime(ir_)));  
	for(int a = 1; a <= dw; a++)
  {
    if(a == 1) 
    {
      Rw[a-1] = r; 
    }
    else
    {
      ITensor Ya(il_, prime(il_)); 

      for(int b = 1; b < a; b++)
      {
        Ya += Tw[a-1][b-1] * Rw[b-1]; 
      }
      Ya *= delta(il_, ir_);
      Ya *= delta(prime(il_), prime(ir_)); 

      if(a < dw)
      {
        //lambda = 0 for nn-TFIM
        double lambda = 0; 

        T guess = Ya; 
        BigM_R M(E, il_, ir_, l_AR, lambda, /*isL*/false); 
				gmres(M, Ya, guess, {"ErrGoal=", gmres_err});	
        Rw[a-1] = guess;  
      }
      else if(a == dw)
      {
        //solve (L1|(1- TL_ +|R)(1|) =  (Ya| - (Ya|R)(1|
        T Y1 = Ya; 
        T Y2 = (Ya * l_AR * delta(ir_, il_)) * delta(prime(il_), prime(ir_)) * r;  
 
        T Y = Y1-Y2; 
        
				//solve A X = Y for x
        BigM_R M(E, il_, ir_, l_AR, /*lambda*/1.0, /*isL*/true); 
				gmres(M, Y, XR_, {"ErrGoal=", gmres_err});	
        Rw[a-1] = XR_;   
      }
      else 
      {
        cerr << "wrong" << endl; 
      }
    }
  }
  T rw(ir_, wr_, prime(ir_)); 
  for(int i = 1; i <= dim(ir_); i++)
  {
    for(int j = 1; j <= dim(prime(ir_)); j++)
    {
      for(int a = 1; a <= dim(wr_); a++)
      {
        rw.set(ir_(i), prime(ir_)(j), wr_(a), eltC(Rw[a-1], ir_(i), prime(ir_)(j))); 
      }
    }
  }

  return rw; 
}
/*----------------------------------------------*/
T iMPS::split_AL() 
  //Eq. 139-142 in SciPost Phys. Lect. Notes 7 (2019) 
  //input: il_-AC_-ir_, il_-C_-ir_; 
  //return: il_-AL-ir_ 
{
  T AC_Cdag = AC_ * dag(C_ * delta(il_, prime(il_))); 
  Index is = findIndex(AC_Cdag, "Site"); 
  T U(il_, is), S, V; 
  svd(AC_Cdag, U, S, V); //AC_Cdag = U*S*V 
  Index iU = commonIndex(U, S); 
  Index iV = commonIndex(S, V); 
  T AL = U * V * delta(iU, iV); 
  AL *= delta(prime(il_), ir_); 

  return AL; 
}
/*----------------------------------------------*/
T iMPS::split_AR()   
  //Eq. 139-142 in SciPost Phys. Lect. Notes 7 (2019) 
  //input: il_-AC_-ir_, il_-C_-ir_; 
  //return: il_-AR-ir_ 
{
  T Cdag_AC = dag(C_ * delta(ir_, prime(ir_))) * AC_;  
  T U(prime(ir_)), S, V;  
  svd(Cdag_AC, U, S, V); 
  Index iU = commonIndex(U, S); 
  Index iV = commonIndex(S, V); 
  T AR = U * V * delta(iU, iV);
  AR *= delta(prime(ir_), il_); 

  return AR; 
}
/*----------------------------------------------*/
void iMPS::to_screen(int i, double t)
{
  double entropy = get_entropy(true); 
  double energy = get_energy(false); 
  cout << setprecision(16); 
  cout << "itdvp entropy: " << entropy << endl; 
  cout << "itdvp energy_l: " << energy << endl;

  return; 
}
/*----------------------------------------------*/
void iMPS::to_file(int i, double t)
{
  T Sx = T(is_, prime(is_));  
  Sx.set(is_(1), prime(is_)(2), 1); 
  Sx.set(is_(2), prime(is_)(1), 1); 

  Cplx mx = eltC(delta(il_, prime(il_)) * AC_ * prime(dag(AC_)) * Sx * delta(ir_, prime(ir_)));  

  entropy_ = get_entropy(false);  
  energy_ = get_energy(false); 
  out_fs_ << setprecision(16); 
  out_fs_ << i << "\t" << abs(t) << "\t" << entropy_ << "\t" << energy_ << "\t" << real(mx) << endl; 
  return; 
}
