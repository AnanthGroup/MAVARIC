/* Author: Elliot Eklund
 *
 * Potential energy functions and their derivatives. */

#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <iostream>
#include <valarray>

#define DELTA 0.1
#define A 1.0
#define D 0.0


struct Potentials{

  /*Parameters passed to initialize. */
  int num_beads;
  int num_states;
  int elec_size;
  double mass;

  void initialize(double num_beadsIN,double num_statesIN,double massIN){
    num_beads = num_beadsIN;
    num_states = num_statesIN;
    elec_size = num_beads*num_states;
    mass = massIN;
  }

  /* State Independent Potential */
  inline double V0(const std::valarray<double> &Q){
    return 0.5*(Q*Q).sum();
  }

  /* Derivative of State Independent Potential*/
  inline std::valarray<double> dV0_dQ(const std::valarray<double> &Q){
    return Q;
  }
 
  /* State Specific Potential Energy Surfaces */
  inline void Velec(const std::valarray<double> &Q,std::valarray<double> &Vout){
  
    //state 1
    Vout[std::slice(0,num_beads,1)] = A*Q + D;
  
    //state2
    Vout[std::slice(num_beads,num_beads,1)] = -A*Q;
  }
  
  /* Derivative of State Specific Potential Energy Surfaces */
  inline void dVelec(const std::valarray<double> &Q,std::valarray<double> &Vout){
  
    //state 1
    Vout[std::slice(0,num_beads,1)] = A;
  
    //state2
    Vout[std::slice(num_beads,num_beads,1)] = -A;
  }
  
  /* Off-diagonal Coupling */
  inline double V_couple(const double &Q,int state1, int state2){
    return DELTA;
  }
  
  /* Derivative of off-diagonal Coupling */
  inline double dV_couple(const double &Q, int state1, int state2){
    return 0;
  }

};

#endif
