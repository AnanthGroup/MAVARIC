/* Author: Elliot Eklund */

#include "MonteCarlo.h"

using namespace std;

MonteCarlo::MonteCarlo(Hamiltonian H,vector<double>MC_params,vector<double> sys_params,
  vector<double> elec_params)
	:H(H), num_beads(sys_params[1]), num_states(elec_params[0]), 
	 beta_n((1.0/sys_params[2])/num_beads),
	 MC_steps(MC_params[1]),
	 esti_rate(MC_params[2]),
	 PSV_W(MC_params[3]), data_W(MC_params[4]), PSV_R(MC_params[5]), data_R(MC_params[6]),
	 elec_size(num_states*num_beads),mt(myRand()),sys_dist(-sys_params[3],sys_params[3]),
	 rand_bead(0,sys_params[1]-1),unif_1(0.0,1.0),elec_dist(-elec_params[1],elec_params[1])
{
  Q.resize(num_beads);
  x.resize(elec_size);
  p.resize(elec_size);
  estimator_t.resize(MC_steps/esti_rate);
  if(PSV_R){read_PSV(num_beads,num_states,Q,x,p);}
}

MonteCarlo::~MonteCarlo(){

}

void MonteCarlo::set_PSV(valarray<double> Qin, valarray<double> xin,valarray<double> pin){
  Q = Qin;
  x = xin;
  p = pin;
}

void MonteCarlo::runMC(){

  Q_prop.resize(num_beads);
  x_prop.resize(elec_size);
  p_prop.resize(elec_size);

  Q_prop = Q;
  x_prop = x;
  p_prop = p;

  sys_steps_accpt = 0;
  sys_steps = 0;

  elec_steps_accpt = 0;
  elec_steps = 0;

  H.update_energy(Q,x,p);
  energy = H.get_energy();

  int esti_samples = 0;
  double estimator_total = 0;
  double sgn_total = 0;

  if(data_R){
    read_MC_data(sgn_total,estimator_total);
  }
  else{
   estimator = H.get_estimator();
  }

  for(int step=0; step<MC_steps; step++){
    //cout << "step: " << step << endl;


    if(unif_1(mt) > 0.5){
      sample_system();
    }

    else {
      sample_elec();
    }
 
    estimator_total += estimator;//*sgn_theta;
    sgn_total += sgn_theta;
    if(step % esti_rate == 0){
      estimator_t[esti_samples] = estimator_total/(sgn_total + 1);
      esti_samples += 1;
    }
  }

  print_sys_accpt(sys_steps,sys_steps_accpt);
  print_elec_accpt(elec_steps,elec_steps_accpt);
  print_avg_energy(estimator_total,sgn_total);
  write_estimator(estimator_t,esti_rate);
  if(PSV_W){write_PSV(num_beads,num_states,Q,x,p);}
  if(data_W){write_MC_data(sgn_total,estimator_total);}
}

void MonteCarlo::sample_system(){

  int mcMove = 0;

  /* Propose new system moves.*/
  for (int i=0; i<num_beads; i++) {
    mcMove = rand_bead(mt);
    Q_prop[mcMove] = Q[mcMove] + sys_dist(mt);
  }

  H.update_energy(Q_prop,x,p);
  double energ_prop = H.get_energy();
  double esti_prop = H.get_estimator();
  int sgn_prop = H.get_sgn();

  /* Accept new system moves if energ_prop < energy*/
  if(energ_prop < energy){
    Q = Q_prop;
    energy = energ_prop;
    estimator = esti_prop;
    sgn_theta = sgn_prop;
    sys_steps_accpt += 1;
  }

  /* Accept new system moves if inequality is met*/
  else if (unif_1(mt) <= exp(-beta_n * (energ_prop - energy))){
    Q = Q_prop;
    energy = energ_prop;
    estimator = esti_prop;
    sgn_theta = sgn_prop;
    sys_steps_accpt += 1;
  }
  else{Q_prop = Q;}

  sys_steps += 1;
}

void MonteCarlo::sample_elec(){

  int mcMove = 0;

  /* Propose new electronic moves.*/
  for (int i=0; i<num_beads; i++) {
    mcMove = num_states * rand_bead(mt);

   for(int state=0; state<num_states; state++){
     mcMove = mcMove + state; 

     x_prop[mcMove] = x[mcMove] + elec_dist(mt);
     p_prop[mcMove] = p[mcMove] + elec_dist(mt);
    }
  }
  
  H.update_energy(Q,x_prop,p_prop);
  double energ_prop = H.get_energy();
  double esti_prop = H.get_estimator();
  int sgn_prop = H.get_sgn();

  /* Accept new electronic move if energ_prop < energy */
  if(energ_prop < energy){
    x = x_prop;
    p = p_prop;
    energy = energ_prop;
    estimator = esti_prop;
    sgn_theta = sgn_prop;
    elec_steps_accpt += 1;
  }

  /* Accept new system moves if inequality is met*/
  else if (unif_1(mt) <= exp(-beta_n * (energ_prop - energy))){
    x = x_prop;
    p = p_prop;
    energy = energ_prop;
    estimator = esti_prop;
    sgn_theta = sgn_prop;
    elec_steps_accpt += 1;
  }

  else{
    x_prop = x;
    p_prop = p;
  }

  elec_steps += 1;
}
 
valarray<double> MonteCarlo::get_Q(){return Q;}

valarray<double> MonteCarlo::get_x(){return x;}

valarray<double> MonteCarlo::get_p(){return p;}
