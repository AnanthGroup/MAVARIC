/* Author:Elliot Eklund

   Theta contains all data and functions necessary for computations relating to theta of
   the MV-RPMD hamiltonian. This includes calculating all the derivatives of theta needed
   for force calculations.

 */

#ifndef THETA_H
#define THETA_H

#include "C_Matrix.h"
#include "M_Matrix.h"
#include "functions.h"

#include <iostream>
#include <valarray>

struct Theta{

  int num_beads; //number of ring polymer beads
  int num_states; //number of electronic states
  double theta; //theta
  double dTheta_dBeta; //derivative of theta w.r.t to beta
  std::complex<double> trGam; //used in update_dTheta_dBeta

  std::vector<std::vector<std::complex<double> > > theta_prod2; //used in update_theta
  std::vector<std::vector<std::complex<double> > > theta_prod1; //used in update_theta

  /* forward and backward chain are a collection of pre-multiplied matricies (a chain of matrix
     multiplication) that optimize force calculations  */
  std::vector<std::vector<std::vector<std::complex<double> > > > forward_chain_dC; 
  std::vector<std::vector<std::vector<std::complex<double> > > > backward_chain_dC;
  std::vector<std::vector<std::vector<std::complex<double> > > > forward_chain_dM;
  std::vector<std::vector<std::vector<std::complex<double> > > > backward_chain_dM;

  std::valarray<double> dTheta_dx; //derivate of theta w.r.t x
  std::valarray<double> dTheta_dp; //derivate of theta w.r.t p
  std::valarray<double> dTheta_dQ; //derivative of theta w.r.t Q

  std::vector<std::vector<std::complex<double> > > identity; // num_states x num_states identity 
  std::vector<std::vector<std::complex<double> > > dummy; //junk matrix used in update_dTheta_dBeta
  std::vector<std::vector<std::complex<double> > > product; //junk matrix used in update_dTheta_dBeta
  std::vector<std::vector<std::complex<double> > > gamma; //gamm matrix used in update_dTheta_dBeta

  C_Matrix C; //C_Matrix object
  M_Matrix M; //M_Matrix object

  /* Initialize all data above. */
  void initialize(int num_beadsIN, int num_statesIN, double betaIN, double massIN){
    num_beads = num_beadsIN;
    num_states = num_statesIN;
    
    C.initialize(num_beads,num_states);
    M.initialize(num_beads,num_states,betaIN,massIN);

    construct_tensor(theta_prod1,num_states,num_states);
    construct_tensor(theta_prod2,num_states,num_states);

    construct_tensor(identity,num_states,num_states);
    construct_tensor(product,num_states,num_states);
    construct_tensor(dummy,num_states,num_states);
    construct_tensor(gamma,num_states,num_states);

    dTheta_dQ.resize(num_beads);
    dTheta_dx.resize(num_states*num_beads);
    dTheta_dp.resize(num_states*num_beads);

    construct_tensor(forward_chain_dM,num_beads,num_states,num_states);
    construct_tensor(backward_chain_dM,num_beads,num_states,num_states);
    construct_tensor(forward_chain_dC,num_beads,num_states,num_states);
    construct_tensor(backward_chain_dC,num_beads,num_states,num_states);


    for(int i=0; i<num_states; i++){
      identity[i][i] = 1.0;
    }
  
    theta = 0;
  
  }

  /* Update theta data member given Q, x, and p. */
  void update_Theta(const std::valarray<double> &Q,const std::valarray<double> &x,
    const std::valarray<double> &p){

    C.update_C(x,p);
    M.update_M(Q);
    theta = compute_Theta(C.C,M.M);

  }

  /* Update theta data member and all forward and backward chains. This call is used during a 
     dynamics calculation. */
  void update_Theta_dyn(const std::valarray<double> &Q,const std::valarray<double> &x,
    const std::valarray<double> &p){

    C.update_C(x,p);
    M.update_M(Q);
    theta = compute_Theta(C.C,M.M);
    update_forward_chain();
    update_backward_chain();
    
  }

  /* Give, C and M, compute theta. That is, multiply C1 * M1 ... Cnum_beads * Mnum_beads, and take
     the real part of the trace. */
  double compute_Theta(const std::vector<std::vector<std::vector<std::complex<double> > > > &C,
    const std::vector<std::vector<std::vector<std::complex<double> > > > &M){

    matrix_mult(C[0],M[0],theta_prod1);

    for(int bead=1; bead<num_beads; bead++){
    
      matrix_mult(theta_prod1,C[bead],theta_prod2);
      matrix_mult(theta_prod2,M[bead],theta_prod1);
    }

    return trace(theta_prod1,num_states).real(); 
  }

  /* Update data member dTheta_dx given x and p.  */
  void update_dTheta_dx(const std::valarray<double> &x,const std::valarray<double> &p){
  
    C.update_dC_dx(x,p);
    int s = 0;
 
    for(int bead=0; bead<num_beads; bead++){
      for(int state=0; state<num_states; state++){
        s = bead*num_states + state;
  
        matrix_mult(forward_chain_dC[bead],C.dC_dx[s],product);
        matrix_mult(product,backward_chain_dC[bead],gamma);
  
        dTheta_dx[s] = trace(gamma,num_states).real();
      }
    }

  }

  /* Update data member dTheta_dp given x and p.  */
  void update_dTheta_dp(const std::valarray<double> &x,const std::valarray<double> &p){
  
    C.update_dC_dp(x,p);
  
    int s = 0;
  
    for(int bead=0; bead<num_beads; bead++){
      for(int state=0; state<num_states; state++){
        s = bead*num_states + state;
  
        matrix_mult(forward_chain_dC[bead],C.dC_dp[s],product);
        matrix_mult(product,backward_chain_dC[bead],gamma);
  
        dTheta_dp[s] = trace(gamma,num_states).real();
      }
    }

  }

  /* Update data member dTheta_dQ given Q.  */
  void update_dTheta_dQ(const valarray<double> &Q){
  
    M.update_dM_dQ(Q);
 
    for(int bead=0; bead<num_beads; bead++){
  
      matrix_mult(forward_chain_dM[bead],M.dM_dQ[bead],product);
      matrix_mult(product,backward_chain_dM[bead],gamma);
  
      dTheta_dQ[bead] = trace(gamma,num_states).real();
    }
  }

  /* Update dTheta_dBeta given Q vector. This functions implicitly calls update_dM_dBeta*/
  void update_dTheta_dBeta(const std::valarray<double> &Q){
 
    trGam = 0;
    M.update_dM_dBeta(Q);
  
    for(int bead=0; bead<num_beads; bead++){
      gamma = identity;
  
      for(int i=0; i<num_beads; i++){
  
        if(bead != i){
          matrix_mult(gamma,C.C[i],dummy);
          matrix_mult(dummy,M.M[i],gamma);

        }
        else{
          matrix_mult(gamma,C.C[i],dummy);
          matrix_mult(dummy,M.dM_dBeta[bead],gamma);
        }
      }

     trGam = trGam + trace(gamma,num_states);

    }

    dTheta_dBeta = trGam.real();
  }

  /* Update both forward_chain_dM and forward_chain_dC. */
  void update_forward_chain(){
  
    forward_chain_dM[0] = C.C[0];
    forward_chain_dC[0] = identity;
  
    for(int i=1; i<num_beads; i++){
      matrix_mult(forward_chain_dM[i-1],M.M[i-1],product);
      forward_chain_dC[i] = product;
      matrix_mult(product,C.C[i],forward_chain_dM[i]);
    }
  }

  /* Update both backward_chain_dM and backward_chain_dC. */
  void update_backward_chain(){
  
    backward_chain_dM[num_beads - 1] = identity;
    backward_chain_dC[num_beads - 1] = M.M[num_beads - 1];
  
    for(int i=num_beads - 2; i>-1; i--){
      matrix_mult(C.C[i+1],backward_chain_dC[i+1],product);
      backward_chain_dM[i] = product;
      matrix_mult(M.M[i],product,backward_chain_dC[i]);
    }
  }

};

#endif
