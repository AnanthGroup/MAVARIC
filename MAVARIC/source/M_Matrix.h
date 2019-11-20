/* Author:Elliot Eklund

 M_Matrix contains all data and functions necessary for computations relating to the
   "M Matrix". */

#ifndef M_MATRIX_H
#define M_MATRIX_H

#include "Potentials.h"
#include "functions.h"

#include <complex>
#include <valarray>

struct M_Matrix{

  int num_beads; //number of ring polymer beads
  int num_states; //number of electronic states
  double beta; //1.0/temp
  double beta_n; //beta/num_beads

  std::vector<std::vector<std::vector<std::complex<double> > > > M; //collection of num_bead M matricies

  /* Collection of num_beads derivative w.r.t. beta M matricies. The derivative of M w.r.t to beta
     corresponding to the alpha-th bead is sotred in element alpha.*/
  std::vector<std::vector<std::vector<std::complex<double> > > > dM_dBeta;

  /* Collection of num_beads derivative w.r.t. Q M matricies. The derivative of M w.r.t to Q
     corresponding to the alpha-th bead is sotred in element alpha.*/
  std::vector<std::vector<std::vector<std::complex<double> > > > dM_dQ;

  std::valarray<double> Vstate; //used in update_M
  std::valarray<double> dVstate; //used in update_dM_dQ
  std::valarray<double> exp_Vstate; //used in update_M
  std::complex<double> M_off_diag; //= 0; //used in update_M
  std::complex<double> M_diag;// = 0; //used in update_dM_dBeta
  std::complex<double> dM_diag;// = 0; //used in update_dM_dBeta
  std::complex<double> trGam; // = 0; //used in update_dTheta_dBeta
  double ONE_n; //1.0/num_beads

  Potentials V; //potential energy structure

  /* Initialize all data above. */
  void initialize(int num_beadsIN, int num_statesIN, double betaIN, double massIN){
    num_beads = num_beadsIN;
    num_states = num_statesIN;
    beta = betaIN;
    beta_n = betaIN/num_beadsIN;

    construct_tensor(M,num_beads,num_states,num_states);
    construct_tensor(dM_dQ,num_beads,num_states,num_states);
    construct_tensor(dM_dBeta,num_beads,num_states,num_states);
 
    Vstate.resize(num_beads*num_states);
    dVstate.resize(num_beads*num_states);
    exp_Vstate.resize(num_beads*num_states);
    ONE_n = 1.0/num_beads;

    M_off_diag = 0;
    M_diag = 0;
    dM_diag = 0;
    trGam = 0;

    V.initialize(num_beads,num_states,massIN);
  }

  /* Update M given Q vector*/
  void update_M(const std::valarray<double> &Q){

    V.Velec(Q,Vstate);
    exp_Vstate = exp(-beta_n * Vstate);

    /* Fill diagonal elements */
    for(int bead=0; bead<num_beads; bead++){
    
      for(int row=0; row<num_states; row++){
        M[bead][row][row] = exp_Vstate[bead + row*num_beads];
      }

    /* Fill upper right triangle */
      for(int row=0; row<num_states-1; row++){
        for(int col=1+row; col<num_states; col++){
          M_off_diag = -beta_n * M[bead][row][row] * V.V_couple(Q[bead],row,col);
          M[bead][row][col] = M_off_diag;
        }
      }
  
    /* Fill lower left triangle */
      for(int row=1; row<num_states; row++){
        for(int col=0; col<row; col++){
          M_off_diag = -beta_n * M[bead][row][row] * V.V_couple(Q[bead],row,col);
          M[bead][row][col] = M_off_diag;
        }
      }
    }

  }

  /* Update dM_dBeta given Q vector */
  void update_dM_dBeta(const std::valarray<double> &Q){

    V.Velec(Q,Vstate);

    /* Fill Diagonal elements */
    for(int bead=0; bead<num_beads; bead++){
      for(int state=0; state<num_states; state++){

          dM_dBeta[bead][state][state] = -ONE_n * Vstate[bead + state*num_beads] * M[bead][state][state];
      }
    }

    /* Fill upper right triangle */
    for(int bead=0; bead<num_beads; bead++){
      for(int row=0; row<num_states-1; row++){
        M_diag = M[bead][row][row];
        dM_diag = dM_dBeta[bead][row][row];

        for(int col=1+row; col<num_states; col++){
          dM_dBeta[bead][row][col] = -ONE_n * V.V_couple(Q[bead],row,col)*(M_diag + beta*dM_diag);
        }
      }
    }

    /* Fill lower left triangle */
    for(int bead=0; bead<num_beads; bead++){
      for(int row=1; row<num_states; row++){
        M_diag = M[bead][row][row];
        dM_diag = dM_dBeta[bead][row][row];

        for(int col=0; col<row; col++){
          dM_dBeta[bead][row][col] = -ONE_n * V.V_couple(Q[bead],row,col)*(M_diag + beta*dM_diag);
        }
      }
    }
  }

  /* Update dM_dQ given Q vector */
  void update_dM_dQ(const valarray<double> &Q) {
  
    V.dVelec(Q,dVstate);

    /* Fill diagonal elements */

    for(int bead=0; bead<num_beads; bead++){
 
      for(int row=0; row<num_states; row++){
        dM_dQ[bead][row][row] = -beta_n * dVstate[bead + row*num_beads]*M[bead][row][row];
      }

    /* Fill upper right triangle */
      for(int row=0; row<num_states-1; row++){
        for(int col=1+row; col<num_states; col++){
	
	  dM_dQ[bead][row][col] = -beta_n * (V.dV_couple(Q[bead],row,col) * M[bead][row][row] +
	     V.V_couple(Q[bead],row,col) * dM_dQ[bead][row][row]);
        }
      }
  
    /* Fill lower left triangle */
      for(int row=1; row<num_states; row++){
        for(int col=0; col<row; col++){

	  dM_dQ[bead][row][col] = -beta_n * (V.dV_couple(Q[bead],row,col) * M[bead][row][row] +
	     V.V_couple(Q[bead],row,col) * dM_dQ[bead][row][row]);
        }
      }
    }

  }

};

#endif
