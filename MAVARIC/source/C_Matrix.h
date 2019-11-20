/* Author:Elliot Eklund

   C_Matrix contains all data and functions necessary for computations relating to the
   "C Matrix". */

#ifndef C_MATRIX_H
#define C_MATRIX_H

#include "functions.h"

#include <complex>
#include <valarray>

struct C_Matrix{

  int num_beads; //number of ring polymer beads
  int num_states; //number of electronic states
  double real, imag; //real and imaginary part of C

  std::valarray<double> x_bead; //used in update_C
  std::valarray<double> p_bead; //used in update_C

  std::vector<std::vector<std::vector<std::complex<double> > > > C; //collection of num_bead C matricies

  /* Collection of num_beads*num_states derivative w.r.t. x C matricies. The derivative of C w.r.t x
     of bead alpha and state k is stored in element alpha*num_beads + k*/
  std::vector<std::vector<std::vector<std::complex<double> > > > dC_dx;

  /* Collection of num_beads*num_states derivative w.r.t. p C matricies. The derivative of C w.r.t p
     of bead alpha and state p is stored in element alpha*num_beads + k*/
  std::vector<std::vector<std::vector<std::complex<double> > > > dC_dp;

  /* Pre-multiplication of M matricies with C matricies used for optimized performance. */
  std::vector<std::vector<std::vector<std::complex<double> > > > forward_chain_dC;

  /* Pre-multiplication of M matricies with C matricies used for optimized performance. */
  std::vector<std::vector<std::vector<std::complex<double> > > > backward_chain_dC;


  /* Initialize all data types above. */
  void initialize(int num_beadsIN,int num_statesIN){
    num_beads = num_beadsIN;
    num_states = num_statesIN;
    real = 0;
    imag = 0;

    construct_tensor(C,num_beads,num_states,num_states);
 
    construct_tensor(dC_dx,num_beads*num_states,num_states,num_states);
    construct_tensor(dC_dp,num_beads*num_states,num_states,num_states);

    x_bead.resize(num_beads);
    p_bead.resize(num_beads);

  }

  /* Update C given x and p vectors */
  void update_C(const std::valarray<double> &x, const std::valarray<double> &p){
  
    for(int bead=0; bead<num_beads; bead++){

      x_bead = x[slice(bead*num_states,num_states,1)];
      p_bead = p[slice(bead*num_states,num_states,1)];

      for(int row=0; row<num_states; row++){
        for(int col=0; col<num_states; col++){
  
          real = x_bead[row]*x_bead[col] + p_bead[row]*p_bead[col];
          imag = p_bead[row]*x_bead[col] - p_bead[col]*x_bead[row];

          C[bead][row][col] = complex<double>(real,imag);
        }
      }

      for(int state=0; state<num_states; state++){
        C[bead][state][state] = C[bead][state][state] - 0.5;
      }
    }

  }

  /* Update dC_dp given x and p vectors */
  void update_dC_dp(const std::valarray<double> &x, const std::valarray<double> &p){
  
    int s = 0;
    int s_row = 0;
    int s_col = 0;
  
    for(int bead=0; bead<num_beads; bead++){
      for(int state=0; state<num_states; state++){
        s = bead*num_states + state;
  
        for(int row=0; row<num_states; row++){
          s_row = bead*num_states + row;
  
          for(int col=0; col<num_states; col++){
            s_col = bead*num_states + col;
  
            real = p[s_col]*kronecker_delta(row,state) + p[s_row]*kronecker_delta(col,state);
            imag = x[s_col]*kronecker_delta(row,state) - x[s_row]*kronecker_delta(col,state);
  
            dC_dp[s][row][col] = complex<double>(real,imag);
  
          }
        }
  
      }
    }
  }

  /* Update dC_dx given x and p vectors */
  void update_dC_dx(const std::valarray<double> &x, const std::valarray<double> &p){
  
    int s = 0;
    int s_row = 0;
    int s_col = 0;
  
    for(int bead=0; bead<num_beads; bead++){
      for(int state=0; state<num_states; state++){
        s = bead*num_states + state;
  
        for(int row=0; row<num_states; row++){
          s_row = bead*num_states + row;
  
          for(int col=0; col<num_states; col++){
            s_col = bead*num_states + col;
  
            real = x[s_col]*kronecker_delta(row,state) + x[s_row]*kronecker_delta(col,state);
            imag = p[s_row]*kronecker_delta(col,state) - p[s_col]*kronecker_delta(row,state);
  
            dC_dx[s][row][col] = complex<double>(real,imag);
  
          }
        }
  
      }
    }

  }

};

#endif
