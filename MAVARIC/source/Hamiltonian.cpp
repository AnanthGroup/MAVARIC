/* Author:Elliot Eklund */

#include "Hamiltonian.h"

using namespace std;

Hamiltonian::Hamiltonian(double mass,int num_beads,int num_states,double beta)
	:mass(mass),num_beads(num_beads),num_states(num_states),beta(beta),
	 beta_n(beta/num_beads),spring_coeff(mass / (beta_n*beta_n)),ONE_beta_n(1.0/beta_n),
	 ONE_n(1.0/num_beads),ONE_mass(1.0/mass),TWO_beta_n(2.0/beta_n)
{
  myTheta.initialize(num_beads,num_states,beta,mass);
  V.initialize(num_beads,num_states,mass);
}

Hamiltonian::Hamiltonian(const Hamiltonian & h)
	:mass(h.mass),num_beads(h.num_beads),num_states(h.num_states),beta(h.beta),
	 beta_n(h.beta_n),spring_coeff(h.spring_coeff),ONE_beta_n(h.ONE_beta_n),
	 ONE_n(h.ONE_n),ONE_mass(h.ONE_mass),TWO_beta_n(h.TWO_beta_n)
{
  myTheta.initialize(h.num_beads,h.num_states,h.beta,mass);
  V.initialize(h.num_beads,h.num_states,h.mass);
}

double Hamiltonian::get_energy(const valarray<double> &Q, const valarray<double> &x, 
  const valarray<double> &p){

  myTheta.update_Theta(Q,x,p);
  return spring(Q) + V.V0(Q) + G(x,p) - ONE_beta_n * log(abs(myTheta.theta));
}

void Hamiltonian::update_energy(const valarray<double> &Q, const valarray<double> &x, 
  const valarray<double> &p){

  myTheta.update_Theta(Q,x,p);
  myTheta.update_dTheta_dBeta(Q);

  double sys_spring = spring(Q);
  double v0 = V.V0(Q);
  double theta = myTheta.theta;
  sgn_theta = sign(theta);

  energy = sys_spring + v0 + G(x,p) - ONE_beta_n * log(abs(theta));
  estimator = (0.5*ONE_beta_n + ONE_n*(v0 - sys_spring) - myTheta.dTheta_dBeta/theta)*sgn_theta;

}

double Hamiltonian::spring(const valarray<double> &Q){

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

inline double Hamiltonian::G(const valarray<double> &x, const valarray<double> &p){
  return ONE_beta_n*((x*x).sum() + (p*p).sum());
}

double Hamiltonian::get_energy(){return energy;}

double Hamiltonian::get_estimator(){return estimator;}

int Hamiltonian::get_sgn(){return sgn_theta;}
