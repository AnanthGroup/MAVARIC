/* Author:Elliot Eklund */

#include "Dynamics.h"

using namespace std;

Dynamics::Dynamics(ABM_MV_RPMD abm,vector<double> Dyn_params,vector<double> Samp_params,
  vector<double> sys_params,vector<double> elec_params)
	:abm(abm),
	 num_trajs(Samp_params[1]),num_states(elec_params[0]),num_beads(sys_params[1]),
	 run_time(Dyn_params[1]), dt(Dyn_params[2]),
	 elec_size(num_states*num_beads),
	 check_energy(Dyn_params[4]),tol(Dyn_params[5])
{
}

void Dynamics::get_PSV(valarray<double> QIN, valarray<double> PIN,
  valarray<double> xIN, valarray<double>pIN){

  Q.resize(num_trajs*num_beads);
  P.resize(num_trajs*num_beads);
  x.resize(num_trajs*elec_size);
  p.resize(num_trajs*elec_size);

  Q = QIN;
  P = PIN;
  x = xIN;
  p = pIN;
}

void Dynamics::run_dynamics(){

  /* _temp vectors store the current trajectory being propogated.*/
  valarray<double> Q_temp(num_beads);
  valarray<double> P_temp(num_beads);
  valarray<double> x_temp(elec_size);
  valarray<double> p_temp(elec_size);

  int num_steps = run_time/dt;

  double sgn_theta = 0;
  double centroid_0 = 0;
  double energy_0 = 0;
  
  sgn_theta_total = 0;
  CQQ.resize(num_steps);

  for(int traj=0; traj<num_trajs; traj++){

    load_temp(Q_temp,P_temp,x_temp,p_temp,traj);
    centroid_0 = Q_temp.sum()/num_beads;
    abm.take_first_steps(Q_temp,P_temp,x_temp,p_temp);

    abm.take_sgn_step(Q_temp,P_temp,x_temp,p_temp,sgn_theta);
    update_CQQ(0,centroid_0,sgn_theta,centroid_0);
    sgn_theta_total += sgn_theta;

    for(int step=1; step<num_steps; step++){

      abm.take_step(Q_temp,P_temp,x_temp,p_temp);
      update_CQQ(step,Q_temp.sum()/num_beads,sgn_theta,centroid_0);
    }
  }

  write_CQQ();
}

void Dynamics::check_energ_conserv(){

  /* _temp vectors store the current trajectory being propogated.*/
  valarray<double> Q_temp(num_beads);
  valarray<double> P_temp(num_beads);
  valarray<double> x_temp(elec_size);
  valarray<double> p_temp(elec_size);

  int num_steps = run_time/dt;
  double energy_0 = 0;
  double energy_t = 0;
  bool broken;  
  int total_broken = 0;

  for(int traj=0; traj<num_trajs; traj++){
    broken = false;
    load_temp(Q_temp,P_temp,x_temp,p_temp,traj);
    abm.take_first_steps(Q_temp,P_temp,x_temp,p_temp);
    
    abm.take_energy_step(Q_temp,P_temp,x_temp,p_temp,energy_0);

    for(int step=1; step<num_steps; step++){
      abm.take_energy_step(Q_temp,P_temp,x_temp,p_temp,energy_t);

      if ((100 * abs(energy_0 - energy_t)/energy_0) >= tol){
        total_broken += 1;
        break;
      }

    }
  }

  cout << endl;
  cout << "Percentage of Trajectories Broken: " << 100 * total_broken/num_trajs << endl;
}

void Dynamics::load_temp(valarray<double> &Q_temp, valarray<double> &P_temp,
  valarray<double> &x_temp, valarray<double> &p_temp, int traj){

  Q_temp = Q[slice(traj*num_beads,num_beads,1)];
  P_temp = P[slice(traj*num_beads,num_beads,1)];
  x_temp = x[slice(traj*elec_size,elec_size,1)];
  p_temp = p[slice(traj*elec_size,elec_size,1)];
}

void Dynamics::update_CQQ(int step, double centroid, double sgn_theta,double centroid_0){

  CQQ[step] += centroid * centroid_0 * sgn_theta;
}

void Dynamics::write_CQQ(){

  ofstream myFile;
  myFile.open("./Results/pos_auto_corr");

  for(int i=0; i<CQQ.size(); i++){

    myFile << i*dt << " " << CQQ[i]/sgn_theta_total << endl;

  }

  myFile.close();

  cout << "Successfully wrote pos_auto_corr to Results." << endl;
}

void Dynamics::read_Q(){

  Q.resize(num_trajs*num_beads);

  ifstream myFile;
  myFile.open("./Results/Trajectories/Q");

  for(int i=0; i<num_trajs*num_beads; i++){
    myFile >> Q[i];    
  }

  myFile.close();
}

void Dynamics::read_P(){

  P.resize(num_trajs*num_beads);

  ifstream myFile;
  myFile.open("./Results/Trajectories/P");

  for(int i=0; i<num_trajs*num_beads; i++){
    myFile >> P[i];    
  }

  myFile.close();
}

void Dynamics::read_x(){

  x.resize(num_trajs*num_beads*num_states);

  ifstream myFile;
  myFile.open("./Results/Trajectories/xelec");

  for(int i=0; i<num_trajs*num_beads*num_states; i++){
    myFile >> x[i];    
  }

  myFile.close();
}

void Dynamics::read_p(){

  p.resize(num_trajs*num_beads*num_states);

  ifstream myFile;
  myFile.open("./Results/Trajectories/pelec");

  for(int i=0; i<num_trajs*num_beads*num_states; i++){
    myFile >> p[i];    
  }

  myFile.close();
}

void Dynamics::read_trajectories(){

  read_Q();
  read_P();
  read_x();
  read_p();

  cout << "Successfully read in all trajectories from Trajectories folder." << endl;
}
