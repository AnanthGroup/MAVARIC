/* Author:Elliot Eklund */

#include "Sampling.h"

using namespace std;


Sampling::Sampling(Hamiltonian H,vector<double>Samp_params,vector<double> sys_params,
  vector<double> elec_params)
	:H(H),mass(sys_params[0]),num_beads(sys_params[1]),beta(1.0/sys_params[2]),
	 beta_n((1.0/sys_params[2])/num_beads),
	 num_states(elec_params[0]),elec_size(num_beads*num_states),
	 num_trajs(Samp_params[1]),decorr_len(Samp_params[2]),
	 mt(myRand()),unif_1(0.0,1.0),
     sys_dist(-sys_params[3],sys_params[3]),
     elec_dist(-elec_params[1],elec_params[1]),
	 rand_bead(0,sys_params[1]-1),
     get_hist_pos(Samp_params[4]),num_bins(Samp_params[5]),
	 save_samp(Samp_params[3]),
     read_PSV(Samp_params[6])
{
  Q.resize(num_beads);
  x.resize(elec_size);
  p.resize(elec_size);
}

Sampling::~Sampling(){
}

void Sampling::set_PSV(valarray<double> Qin, valarray<double> xin,valarray<double> pin){
  Q = Qin;
  x = xin;
  p = pin;
}

void Sampling::runSampling(){

  if(read_PSV){
    read_PSV_fxn();    
  }

  Q_prop.resize(num_beads);
  x_prop.resize(elec_size);
  p_prop.resize(elec_size);

  Q_samp.resize(num_trajs*num_beads);
  P_samp.resize(num_trajs*num_beads);
  x_samp.resize(num_trajs*elec_size);
  p_samp.resize(num_trajs*elec_size);

  Q_prop = Q;
  x_prop = x;
  p_prop = p;

  energy = H.get_energy(Q,x,p);

  for(int traj=0; traj<num_trajs; traj++){
    for(int i=0; i<decorr_len; i++){    

      if(unif_1(mt) > 0.5){
        sample_system();
      }

      else {
        sample_elec();
      }
    }

    Q_samp[slice(traj*num_beads,num_beads,1)] = Q;
    x_samp[slice(traj*elec_size,elec_size,1)] = x;
    p_samp[slice(traj*elec_size,elec_size,1)] = p;

  }

  sample_boltzman();

  if(get_hist_pos){
    histogram_pos(); 
  }

  if(save_samp){
    save_trajectories();
  }

}

void Sampling::sample_system(){

  int mcMove = 0;

  /* Propose new system moves.*/
  for (int i=0; i<num_beads; i++) {
    mcMove = rand_bead(mt);
    Q_prop[mcMove] = Q[mcMove] + sys_dist(mt);
  }

  double energ_prop = H.get_energy(Q_prop,x,p);

  /* Accept new system moves if energ_prop < energy*/
  if(energ_prop < energy){
    Q = Q_prop;
    energy = energ_prop;
  }

  /* Accept new system moves if inequality is met*/
  else if (unif_1(mt) <= exp(-beta_n * (energ_prop - energy))){
    Q = Q_prop;
    energy = energ_prop;
  }

  else{Q_prop = Q;}
}

void Sampling::sample_elec(){

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
  
  double energ_prop = H.get_energy(Q,x_prop,p_prop);

  /* Accept new system moves if energ_prop < energy*/
  if(energ_prop < energy){
    x = x_prop;
    p = p_prop;
    energy = energ_prop;
  }

  /* Accept new system moves if inequality is met*/
  else if (unif_1(mt) <= exp(-beta_n * (energ_prop - energy))){
    x = x_prop;
    p = p_prop;
    energy = energ_prop;
  }

  else{
    x_prop = x;
    p_prop = p;
  }

}

void Sampling::sample_boltzman(){

  double stdev = sqrt(num_beads*mass/beta);
  std::normal_distribution<double> boltz_dist(0,stdev);

  for(int i=0; i<num_trajs*num_beads; i++){
    P_samp[i] = boltz_dist(mt);
  }
}

void Sampling::histogram_pos(){

  valarray<double> Q_centroid (num_trajs); 
  valarray<double> Q_temp (num_beads);


  for(int i=0; i<num_trajs; i++){
    Q_temp = Q_samp[slice(i*num_beads,num_beads,1)]; 
    Q_centroid[i] = Q_temp.sum()/num_beads;
  }

  string name = "./Results/pos_histogram";

  histogram(name,num_bins,Q_centroid); 

}

void Sampling::save_trajectories(){
  save_Q();
  save_P();
  save_x();
  save_p();
    
  cout << "\t Successfully saved sampled trajectories to Results/Trajectories." << endl;
}

void Sampling::save_Q(){

  ofstream myFile;
  myFile.open("./Results/Trajectories/Q");

  for(int i=0; i<num_trajs*num_beads; i++){
    myFile << Q_samp[i] << endl;
  }

  myFile.close();
}

void Sampling::save_P(){

  ofstream myFile;
  myFile.open("./Results/Trajectories/P");

  for(int i=0; i<num_trajs*num_beads; i++){
    myFile << P_samp[i] << endl;
  }

  myFile.close();
}

void Sampling::save_x(){

  ofstream myFile;
  myFile.open("./Results/Trajectories/xelec");

  for(int i=0; i<num_trajs*num_beads*num_states; i++){
    myFile << x_samp[i] << endl;
  }

  myFile.close();
}

void Sampling::save_p(){

  ofstream myFile;
  myFile.open("./Results/Trajectories/pelec");

  for(int i=0; i<num_trajs*num_beads*num_states; i++){
    myFile << p_samp[i] << endl;
  }

  myFile.close();
}

valarray<double> Sampling::get_Q(){return Q_samp;}

valarray<double> Sampling::get_P(){return P_samp;}

valarray<double> Sampling::get_x(){return x_samp;}

valarray<double> Sampling::get_p(){return p_samp;}

void Sampling::read_PSV_fxn(){

  ifstream myFile;
  myFile.open("Results/PSV");

  int num_beads_temp;
  int num_states_temp;

  myFile >> num_beads_temp;
  myFile >> num_states_temp;

  for(int i=0; i<num_beads; i++){
    myFile >> Q[i];
  }

  for(int i=0; i<num_beads*num_states; i++){
    myFile >> x[i];
  }

  for(int i=0; i<num_beads*num_states; i++){
    myFile >> p[i];
  }

  try{
    if(num_beads_temp != num_beads){
      throw "WARNING: Previous calculation used a different num_beads.";
    }
  }

  catch (const char *c){cout << c;}

  myFile.close();
  cout << "\t Successfully read PSV from Results." << endl;
}
