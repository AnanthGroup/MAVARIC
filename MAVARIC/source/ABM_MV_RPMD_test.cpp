/* Author:Elliot Eklund */

#include "ABM_MV_RPMD_test.h"

ABM_MV_RPMD_test::ABM_MV_RPMD_test(Forces F,int num_beads,int num_states, double dt)
	:F(F), my_rk4(F,num_beads,num_states,-dt),
	 num_beads(num_beads),num_states(num_states),dt(dt),
	 h1_p(dt*55.0/24.0),h2_p(dt*59.0/24.0),h3_p(dt*37.0/24.0),h4_p(dt*9.0/24.0),
  	 h1_c(dt*9.0/24.0),h2_c(dt*19.0/24.0),h3_c(dt*5.0/24.0),h4_c(dt*1.0/24.0)
{
  initialize();
}

ABM_MV_RPMD_test::ABM_MV_RPMD_test(const ABM_MV_RPMD_test &ABM)
	:F(ABM.F), my_rk4(ABM.F,ABM.num_beads,ABM.num_states,ABM.dt),
	 num_beads(ABM.num_beads),num_states(ABM.num_states),dt(ABM.dt),
	 h1_p(ABM.h1_p),h2_p(ABM.h2_p),h3_p(ABM.h3_p),h4_p(ABM.h4_p),
	 h1_c(ABM.h1_c),h2_c(ABM.h2_c),h3_c(ABM.h3_c),h4_c(ABM.h4_c)
{
  initialize();
}

void ABM_MV_RPMD_test::initialize(){

  /* Variables relevent to RK4 integrator */
  f_Q_alpha.resize(num_beads);
  f_Q_beta.resize(num_beads);
  f_Q_gamma.resize(num_beads);
  f_Q_delta.resize(num_beads);
  f_Q_p1.resize(num_beads);

  f_P_alpha.resize(num_beads);
  f_P_beta.resize(num_beads);
  f_P_gamma.resize(num_beads);
  f_P_delta.resize(num_beads);
  f_P_p1.resize(num_beads);

  f_x_alpha.resize(num_states*num_beads);
  f_x_beta.resize(num_states*num_beads);
  f_x_gamma.resize(num_states*num_beads);
  f_x_delta.resize(num_states*num_beads);
  f_x_p1.resize(num_states*num_beads);

  f_p_alpha.resize(num_states*num_beads);
  f_p_beta.resize(num_states*num_beads);
  f_p_gamma.resize(num_states*num_beads);
  f_p_delta.resize(num_states*num_beads);
  f_p_p1.resize(num_states*num_beads);

  pred_Q.resize(num_beads);
  pred_P.resize(num_beads);
  pred_x.resize(num_states*num_beads);
  pred_p.resize(num_states*num_beads);

  initialize_pointers();
}

void ABM_MV_RPMD_test::initialize_pointers(){

  p_Q_free = &f_Q_alpha; 
  p_Q_1 =    &f_Q_beta; 
  p_Q_2 =    &f_Q_gamma; 
  p_Q_3 =    &f_Q_delta; 

  p_P_free = &f_P_alpha; 
  p_P_1 =    &f_P_beta; 
  p_P_2 =    &f_P_gamma; 
  p_P_3 =    &f_P_delta; 

  p_x_free = &f_x_alpha; 
  p_x_1 =    &f_x_beta; 
  p_x_2 =    &f_x_gamma; 
  p_x_3 =    &f_x_delta; 

  p_p_free = &f_p_alpha; 
  p_p_1 =    &f_p_beta; 
  p_p_2 =    &f_p_gamma; 
  p_p_3 =    &f_p_delta; 
}

void ABM_MV_RPMD_test::take_first_steps(valarray<double> &Q,valarray<double> &P,
  valarray<double> &x,valarray<double> &p){

  valarray<double> Q_copy = Q;
  valarray<double> P_copy = P;
  valarray<double> x_copy = x;
  valarray<double> p_copy = p;

  /* Use RK4 to take steps backwards. These steps will seed ABM. */
  RK4_MV_RPMD my_rk4(F,num_beads,num_states,-dt);

  /* Take step backward to t_(-1)*/
  my_rk4.take_step(Q,P,x,p);
  update_f_beta(Q,P,x,p);

  /* Take step backward to t_(-1)*/
  my_rk4.take_step(Q,P,x,p);

  update_f_gamma(Q,P,x,p);
  
  /* Take step backward to t_(-2)*/
  my_rk4.take_step(Q,P,x,p);

  update_f_delta(Q,P,x,p);

  /* Copy original PSV at t_0 back to outgoing variables. */
  Q = Q_copy; 
  P = P_copy; 
  x = x_copy;
  p = p_copy;

}

void ABM_MV_RPMD_test::take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
}

void ABM_MV_RPMD_test::take_sgn_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p,double &sgn_theta){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
  sgn_theta = F.get_sgn_theta();
}

void ABM_MV_RPMD_test::take_energy_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p,double &energy){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
  energy = F.get_energy(Q,P,x,p);
}

void ABM_MV_RPMD_test::predict(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  update_p_free(Q,P,x,p);

  pred_Q = Q + (h1_p *  *p_Q_free) - (h2_p *  *p_Q_1) + (h3_p *  *p_Q_2) - (h4_p *  *p_Q_3);
  pred_P = P + (h1_p *  *p_P_free) - (h2_p *  *p_P_1) + (h3_p *  *p_P_2) - (h4_p *  *p_P_3);
  pred_x = x + (h1_p *  *p_x_free) - (h2_p *  *p_x_1) + (h3_p *  *p_x_2) - (h4_p *  *p_x_3);
  pred_p = p + (h1_p *  *p_p_free) - (h2_p *  *p_p_1) + (h3_p *  *p_p_2) - (h4_p *  *p_p_3);

}

void ABM_MV_RPMD_test::correct(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  update_f_p1(pred_Q,pred_P,pred_x,pred_p);

  Q = Q + (h1_c * f_Q_p1) + (h2_c * *p_Q_free) - (h3_c * *p_Q_1) + (h4_c * *p_Q_2);
  P = P + (h1_c * f_P_p1) + (h2_c * *p_P_free) - (h3_c * *p_P_1) + (h4_c * *p_P_2);
  x = x + (h1_c * f_x_p1) + (h2_c * *p_x_free) - (h3_c * *p_x_1) + (h4_c * *p_x_2);
  p = p + (h1_c * f_p_p1) + (h2_c * *p_p_free) - (h3_c * *p_p_1) + (h4_c * *p_p_2);

  swap();

}

void ABM_MV_RPMD_test::update_f_delta(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_delta =   F.dH_dP(P);
  f_P_delta = - F.dH_dQ(Q,x,p);
  f_x_delta  =   F.dH_dp(Q,x,p); 
  f_p_delta  = - F.dH_dx(Q,x,p);
}
  
void ABM_MV_RPMD_test::update_f_gamma(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_gamma =   F.dH_dP(P);
  f_P_gamma = - F.dH_dQ(Q,x,p);
  f_x_gamma  =   F.dH_dp(Q,x,p); 
  f_p_gamma  = - F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD_test::update_f_beta(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_beta =   F.dH_dP(P);
  f_P_beta = - F.dH_dQ(Q,x,p);
  f_x_beta  =   F.dH_dp(Q,x,p); 
  f_p_beta  = - F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD_test::update_p_free(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);

  *p_Q_free =  F.dH_dP(P);
  *p_P_free = -F.dH_dQ(Q,x,p);
  *p_x_free =  F.dH_dp(Q,x,p); 
  *p_p_free = -F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD_test::update_f_p1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);

  f_Q_p1 =  F.dH_dP(P);
  f_P_p1 = -F.dH_dQ(Q,x,p);
  f_x_p1 =  F.dH_dp(Q,x,p); 
  f_p_p1 = -F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD_test::swap(){

  p_Q_3    = p_Q_2;
  p_Q_2    = p_Q_1;
  p_Q_1    = p_Q_free;
  p_Q_free = p_Q_3;

  p_P_3    = p_P_2;
  p_P_2    = p_P_1;
  p_P_1    = p_P_free;
  p_P_free = p_P_3;

  p_x_3    = p_x_2;
  p_x_2    = p_x_1;
  p_x_1    = p_x_free;
  p_x_free = p_x_3;

  p_p_3    = p_p_2;
  p_p_2    = p_p_1;
  p_p_1    = p_p_free;
  p_p_free = p_p_3;
}
