/* Author:Elliot Eklund */

#include "ABM_MV_RPMD.h"

ABM_MV_RPMD::ABM_MV_RPMD(Forces F,int num_beads,int num_states, double dt)
	:F(F), my_rk4(F,num_beads,num_states,-dt),
	 num_beads(num_beads),num_states(num_states),dt(dt),
	 h1_p(dt*55.0/24.0),h2_p(dt*59.0/24.0),h3_p(dt*37.0/24.0),h4_p(dt*9.0/24.0),
  	 h1_c(dt*9.0/24.0),h2_c(dt*19.0/24.0),h3_c(dt*5.0/24.0),h4_c(dt*1.0/24.0)
{
  initialize();
}

ABM_MV_RPMD::ABM_MV_RPMD(const ABM_MV_RPMD &ABM)
	:F(ABM.F), my_rk4(ABM.F,ABM.num_beads,ABM.num_states,ABM.dt),
	 num_beads(ABM.num_beads),num_states(ABM.num_states),dt(ABM.dt),
	 h1_p(ABM.h1_p),h2_p(ABM.h2_p),h3_p(ABM.h3_p),h4_p(ABM.h4_p),
	 h1_c(ABM.h1_c),h2_c(ABM.h2_c),h3_c(ABM.h3_c),h4_c(ABM.h4_c)
{
  initialize();
}

void ABM_MV_RPMD::initialize(){

  /* Variables relevent to RK4 integrator */
  f_Q_3.resize(num_beads);
  f_Q_2.resize(num_beads);
  f_Q_1.resize(num_beads);
  f_Q_0.resize(num_beads);
  f_Q_p1.resize(num_beads);

  f_P_3.resize(num_beads);
  f_P_2.resize(num_beads);
  f_P_1.resize(num_beads);
  f_P_0.resize(num_beads);
  f_P_p1.resize(num_beads);

  f_x_3.resize(num_states*num_beads);
  f_x_2.resize(num_states*num_beads);
  f_x_1.resize(num_states*num_beads);
  f_x_0.resize(num_states*num_beads);
  f_x_p1.resize(num_states*num_beads);

  f_p_3.resize(num_states*num_beads);
  f_p_2.resize(num_states*num_beads);
  f_p_1.resize(num_states*num_beads);
  f_p_0.resize(num_states*num_beads);
  f_p_p1.resize(num_states*num_beads);

  pred_Q.resize(num_beads);
  pred_P.resize(num_beads);
  pred_x.resize(num_states*num_beads);
  pred_p.resize(num_states*num_beads);

}

void ABM_MV_RPMD::take_first_steps(valarray<double> &Q,valarray<double> &P,
  valarray<double> &x,valarray<double> &p){

  valarray<double> Q_copy = Q;
  valarray<double> P_copy = P;
  valarray<double> x_copy = x;
  valarray<double> p_copy = p;

  /* Use RK4 to take steps backwards. These steps will seed ABM. */
  RK4_MV_RPMD my_rk4(F,num_beads,num_states,-dt);


  my_rk4.take_step(Q,P,x,p);
  update_f_1(Q,P,x,p);

  /* Take step backward to t_(-1)*/
  my_rk4.take_step(Q,P,x,p);

  update_f_2(Q,P,x,p);

  /* Take step backward to t_(-2)*/
  my_rk4.take_step(Q,P,x,p);

  update_f_3(Q,P,x,p);

  /* Copy original PSV at t_0 back to outgoing variables. */
  Q = Q_copy; 
  P = P_copy; 
  x = x_copy;
  p = p_copy;

}

void ABM_MV_RPMD::take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
}

void ABM_MV_RPMD::take_sgn_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p,double &sgn_theta){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
  sgn_theta = F.get_sgn_theta();
}

void ABM_MV_RPMD::take_energy_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
    valarray<double> &p,double &energy){

  predict(Q,P,x,p);
  correct(Q,P,x,p);
  energy = F.get_energy(Q,P,x,p);
}

void ABM_MV_RPMD::predict(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  update_f_0(Q,P,x,p);

  pred_Q = Q + (h1_p * f_Q_0) - (h2_p * f_Q_1) + (h3_p * f_Q_2) - (h4_p * f_Q_3);
  pred_P = P + (h1_p * f_P_0) - (h2_p * f_P_1) + (h3_p * f_P_2) - (h4_p * f_P_3);
  pred_x = x + (h1_p * f_x_0) - (h2_p * f_x_1) + (h3_p * f_x_2) - (h4_p * f_x_3);
  pred_p = p + (h1_p * f_p_0) - (h2_p * f_p_1) + (h3_p * f_p_2) - (h4_p * f_p_3);

}

void ABM_MV_RPMD::correct(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  update_f_p1(pred_Q,pred_P,pred_x,pred_p);

  Q = Q + (h1_c * f_Q_p1) + (h2_c * f_Q_0) - (h3_c * f_Q_1) + (h4_c * f_Q_2);
  P = P + (h1_c * f_P_p1) + (h2_c * f_P_0) - (h3_c * f_P_1) + (h4_c * f_P_2);
  x  = x  + (h1_c * f_x_p1)  + (h2_c * f_x_0)  - (h3_c * f_x_1)  + (h4_c * f_x_2);
  p  = p  + (h1_c * f_p_p1)  + (h2_c * f_p_0)  - (h3_c * f_p_1)  + (h4_c * f_p_2);

  f_Q_3 = f_Q_2;
  f_Q_2 = f_Q_1;
  f_Q_1 = f_Q_0;

  f_P_3 = f_P_2;
  f_P_2 = f_P_1;
  f_P_1 = f_P_0;

  f_x_3 = f_x_2;
  f_x_2 = f_x_1;
  f_x_1 = f_x_0;

  f_p_3 = f_p_2;
  f_p_2 = f_p_1;
  f_p_1 = f_p_0;

}

void ABM_MV_RPMD::update_f_3(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_3 =   F.dH_dP(P);
  f_P_3 = - F.dH_dQ(Q,x,p);
  f_x_3  =   F.dH_dp(Q,x,p); 
  f_p_3  = - F.dH_dx(Q,x,p);
}
  
void ABM_MV_RPMD::update_f_2(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_2 =   F.dH_dP(P);
  f_P_2 = - F.dH_dQ(Q,x,p);
  f_x_2  =   F.dH_dp(Q,x,p); 
  f_p_2  = - F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD::update_f_1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);  

  f_Q_1 =   F.dH_dP(P);
  f_P_1 = - F.dH_dQ(Q,x,p);
  f_x_1  =   F.dH_dp(Q,x,p); 
  f_p_1  = - F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD::update_f_0(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);

  f_Q_0 =  F.dH_dP(P);
  f_P_0 = -F.dH_dQ(Q,x,p);
  f_x_0 =  F.dH_dp(Q,x,p); 
  f_p_0 = -F.dH_dx(Q,x,p);
}

void ABM_MV_RPMD::update_f_p1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
  valarray<double> &p){

  F.update_forces(Q,x,p);

  f_Q_p1 =  F.dH_dP(P);
  f_P_p1 = -F.dH_dQ(Q,x,p);
  f_x_p1 =  F.dH_dp(Q,x,p); 
  f_p_p1 = -F.dH_dx(Q,x,p);
}
