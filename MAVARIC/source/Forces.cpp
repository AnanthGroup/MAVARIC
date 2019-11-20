/* Author:Elliot Eklund*/ 
 
#include "Forces.h"

using namespace std;

Forces::Forces(int num_beads,int num_states,double beta,double mass)
	:num_beads(num_beads),num_states(num_states),beta(beta),mass(mass),
	 ONE_beta_n(1.0/(beta/num_beads)),TWO_beta_n(2.0/(beta/num_beads)),
	 ONE_mass(1.0/mass),spring_coeff(num_beads*num_beads*mass/(beta*beta)),
	 beta_n(beta/num_beads),ONE_TWO_m(1.0/(2.0*mass))
{
  myTheta.initialize(num_beads,num_states,beta,mass);
  V.initialize(num_beads,num_states,mass);
}

Forces::Forces(const Forces &F)
	:num_beads(F.num_beads),num_states(F.num_states),beta(F.beta),mass(F.mass),
	 ONE_beta_n(F.ONE_beta_n),TWO_beta_n(F.TWO_beta_n),
	 ONE_mass(F.ONE_mass),spring_coeff(F.spring_coeff),
	 beta_n(F.beta_n),ONE_TWO_m(F.ONE_TWO_m)
{
  myTheta.initialize(F.num_beads,F.num_states,F.beta,mass);
  V.initialize(F.num_beads,F.num_states,F.mass);
}

void Forces::update_forces(const valarray<double> &Q, const valarray<double> &x,
  const valarray<double> &p){

  myTheta.update_Theta_dyn(Q,x,p);
  myTheta.update_dTheta_dx(x,p);
  myTheta.update_dTheta_dp(x,p);
  myTheta.update_dTheta_dQ(Q);
  sgn_theta = sign(myTheta.theta);
}

valarray<double> Forces::dH_dQ(const valarray<double> &Q, const valarray<double> &x,
  const valarray<double> &p){

  return dspring_dQ(Q) + V.dV0_dQ(Q) - ONE_beta_n * myTheta.dTheta_dQ / myTheta.theta;
}

valarray<double> Forces::dH_dP(const valarray<double> &P){
  return ONE_mass * P;
} 

valarray<double> Forces::dH_dp(const valarray<double> &Q, const valarray<double> &x,
  const valarray<double> &p){
  return TWO_beta_n * p - ONE_beta_n * myTheta.dTheta_dp/myTheta.theta;
}

valarray<double> Forces::dH_dx(const valarray<double> &Q, const valarray<double> &x,
  const valarray<double> &p) {

  return  TWO_beta_n * x - ONE_beta_n * myTheta.dTheta_dx/myTheta.theta;
}

valarray<double> Forces::dspring_dQ(const valarray<double> &Q){

  valarray<double> force (num_beads);

  int i_up = 0; //index one bead up
  int i_down = 0; //index one bead down

  for(int bead=0; bead<num_beads; bead++){

    i_up = (bead + 1) % num_beads;
    i_down = (bead + num_beads - 1) % num_beads;

    force[bead] = 2*Q[bead] - Q[i_up] - Q[i_down];
  }

  return spring_coeff * force;
}

double Forces::get_sgn_theta() {return sgn_theta;}

inline double Forces::G(const valarray<double> &x, const valarray<double> &p){
  return ((x*x).sum() + (p*p).sum()) / beta_n;
}

double Forces::spring(const valarray<double> &Q){

  double sum = 0; //total sum of spring term
  int b_up = 0; //index one bead up
  double diff = 0; //difference between Q[i] and Q[i+1]

  for(int bead=0; bead<num_beads; bead++){

    b_up = (bead + 1) % num_beads;
    diff = Q[bead] - Q[b_up];
    sum += diff*diff;
  }

  return 0.5 * spring_coeff * sum;
}

double Forces::get_energy(const valarray<double> &Q, const valarray<double> &P,
  const valarray<double> &x, const valarray<double> &p){

  return 0.5 * (P*P).sum() + spring(Q) + V.V0(Q) + G(x,p) - ONE_beta_n * log(abs(myTheta.theta));
}
