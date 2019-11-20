#include "RK4_MV_RPMD.h"

RK4_MV_RPMD::RK4_MV_RPMD(Forces F,int num_beads,int num_states,double dt)
	:F(F),
	 num_beads(num_beads),num_states(num_states),
	 dt(dt), dt_half(0.5*dt),
	 coeff1(dt/6.0),coeff2(2.0*dt/6.0),coeff3(2.0*dt/6.0),coeff4(dt/6.0)
{
  zero_ks();
}

void RK4_MV_RPMD::zero_ks(){

  k1Q.resize(num_beads);
  k2Q.resize(num_beads);
  k3Q.resize(num_beads);
  k4Q.resize(num_beads);

  k1P.resize(num_beads);
  k2P.resize(num_beads);
  k3P.resize(num_beads);
  k4P.resize(num_beads);

  k1x.resize(num_states*num_beads);
  k2x.resize(num_states*num_beads);
  k3x.resize(num_states*num_beads);
  k4x.resize(num_states*num_beads);

  k1p.resize(num_states*num_beads);
  k2p.resize(num_states*num_beads);
  k3p.resize(num_states*num_beads);
  k4p.resize(num_states*num_beads);
}

void RK4_MV_RPMD::take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);
  update_k1(Q,P,x,p);

  F.update_forces(Q + dt_half*k1Q, x + dt_half*k1x, p + dt_half*k1p);
  update_k2(Q + dt_half*k1Q, P + dt_half*k1P, x + dt_half*k1x, p + dt_half*k1p);

  F.update_forces(Q + dt_half*k2Q, x + dt_half*k2x, p + dt_half*k2p);
  update_k3(Q + dt_half*k2Q, P + dt_half*k2P, x + dt_half*k2x, p + dt_half*k2p);

  F.update_forces(Q + dt*k3Q, x + dt*k3x , p + dt*k3p );
  update_k4(Q + dt*k3Q, P + dt*k3P, x + dt*k3x, p + dt*k3p);

  update_final(Q,P,x,p);
}

void RK4_MV_RPMD::update_k1(const valarray<double> &Q, const valarray<double> &P,
  const valarray<double> &x, const valarray<double> &p){

  k1Q = F.dH_dP(P);
  k1P = - F.dH_dQ(Q,x,p);
  k1x =  F.dH_dp(Q,x,p);
  k1p = - F.dH_dx(Q,x,p);
}

void RK4_MV_RPMD::update_k2(const valarray<double> &Q, const valarray<double> &P,
  const valarray<double> &x, const valarray<double> &p){

  k2Q = F.dH_dP(P);
  k2P = - F.dH_dQ(Q,x,p);
  k2x =  F.dH_dp(Q,x,p);
  k2p = - F.dH_dx(Q,x,p);
}

void RK4_MV_RPMD::update_k3(const valarray<double> &Q, const valarray<double> &P,
  const valarray<double> &x, const valarray<double> &p){

  k3Q = F.dH_dP(P);
  k3P = - F.dH_dQ(Q,x,p);
  k3x =  F.dH_dp(Q,x,p);
  k3p = - F.dH_dx(Q,x,p);
}

void RK4_MV_RPMD::update_k4(const valarray<double> &Q, const valarray<double> &P,
  const valarray<double> &x, const valarray<double> &p){

  k4Q = F.dH_dP(P);
  k4P = - F.dH_dQ(Q,x,p);
  k4x =  F.dH_dp(Q,x,p);
  k4p = - F.dH_dx(Q,x,p);
}

void RK4_MV_RPMD::update_final(valarray<double> &Q, valarray<double> &P,
  valarray<double> &x,valarray<double> &p){

  Q = Q + (coeff1 * k1Q) + (coeff2 * k2Q) + (coeff3 * k3Q) + (coeff4 *k4Q);
  P = P + (coeff1 * k1P) + (coeff2 * k2P) + (coeff3 * k3P) + (coeff4 *k4P);
  x  = x  + (coeff1 * k1x)  + (coeff2 * k2x)  + (coeff3 * k3x)  + (coeff4 *k4x);
  p  = p  + (coeff1 * k1p)  + (coeff2 * k2p)  + (coeff3 * k3p)  + (coeff4 *k4p);
}

