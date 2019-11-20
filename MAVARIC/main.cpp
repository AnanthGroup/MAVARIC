/* Author: Elliot Eklund*/

#include "MainHelper.h"
#include "Hamiltonian.h"
#include "MonteCarlo.h"
#include "Sampling.h"
#include "Forces.h"
#include "Dynamics.h"
#include "ABM_MV_RPMD.h"

#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <valarray>

using namespace std;

int main(int argc, char ** argv){

  /* Vectors used to store parameters from InputFiles */
  vector<double> sys_parameters;
  vector<double> elec_parameters;
  vector<double> MC_parameters;
  vector<double> Samp_parameters;
  vector<double> Dyn_parameters;

  int abort = input_file_handler(sys_parameters,elec_parameters,MC_parameters,Samp_parameters,Dyn_parameters);
 
  if (abort == -1){
    return -1;
  }

  /* From MonteCarlo */
  const bool runMC = MC_parameters[0];
  
  /* From Sampling */
  const bool runSamp = Samp_parameters[0];
  const int num_trajs = Samp_parameters[1];
  
  /* From Dynamics */
  const bool runDyn = Dyn_parameters[0];
  const double dt = Dyn_parameters[2];
  const bool read_trajs = Dyn_parameters[3];
  const bool check_energy = Dyn_parameters[4];

  /* From SystemParameters */
  const double mass = sys_parameters[0];
  const int num_beads = sys_parameters[1];
  const double temp = sys_parameters[2];
  const double beta = 1.0/temp;

  /* From ElecParameters */
  const int num_states = elec_parameters[0];

  /* Vectors containing phase space variables PSV for MonteCarlo and Sampling. */
  valarray<double> Q (num_beads);
  valarray<double> P (num_beads);
  valarray<double> x (num_beads*num_states);
  valarray<double> p (num_beads*num_states);

  /* Vectors containing initial trajectories used for Dynamics. */
  valarray<double> Q_samp;
  valarray<double> P_samp;
  valarray<double> x_samp;
  valarray<double> p_samp;

  /* Random Number Generation Process */
  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<double> unif(-1.0,1.0);

  /* Randomly distributed PSV from [-1,1] */
  for(int i=0; i<num_beads; i++){
    Q[i] = unif(mt);
    P[i] = unif(mt);
  }

  for(int i=0; i<num_beads*num_states; i++){
    x[i] = unif(mt); 
    p[i] = unif(mt);
  }

  cout << endl << endl << endl;
  cout << "Calculating..." << endl << endl;

  /* Run MonteCarlo simulation if runMC is true. */
  if(runMC){

    cout << "Running Monte Carlo calculations ..." << endl;
    cout << endl;   

    cout << endl;
    cout << "MonteCarlo Results:" << endl;
    Hamiltonian H(mass,num_beads,num_states,beta);
    MonteCarlo myMC(H,MC_parameters,sys_parameters,elec_parameters);
    myMC.set_PSV(Q,x,p);

    auto start = std::chrono::high_resolution_clock::now();

    myMC.runMC();

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    cout << endl;  
    cout << "Monte Carlo Simulation Run Time: " << elapsed.count() << endl;

    Q = myMC.get_Q();
    x = myMC.get_x();
    p = myMC.get_p();
  }

  /* Run Sampling if runSamp is true. */
  if(runSamp){

    cout << endl << endl;
    cout << "Running Sampling ..." << endl;

    Hamiltonian H(mass,num_beads,num_states,beta);
    Sampling mySamp(H,Samp_parameters,sys_parameters,elec_parameters); 
    mySamp.set_PSV(Q,x,p);

    auto start = std::chrono::high_resolution_clock::now();

    mySamp.runSampling();

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
 
    cout << endl; 
    cout << "Sampling Run Time: " << elapsed.count() << endl;

    Q_samp.resize(num_trajs*num_beads);
    P_samp.resize(num_trajs*num_beads);
    x_samp.resize(num_trajs*num_beads*num_states);
    p_samp.resize(num_trajs*num_beads*num_states);

    Q_samp = mySamp.get_Q();
    P_samp = mySamp.get_P();
    x_samp = mySamp.get_x();
    p_samp = mySamp.get_p();
  }


  /* Run Sampling if runSamp is true. */
  if(runDyn){
    cout << endl << endl;
    cout << "Running Dynamics ..." << endl;

    Forces myForces(num_beads,num_states,beta,mass);
    ABM_MV_RPMD my_abm(myForces,num_beads,num_states,dt);
    Dynamics myDyn(my_abm,Dyn_parameters,Samp_parameters,sys_parameters,elec_parameters);

   if(read_trajs){
     myDyn.read_trajectories();
   }
   else{
     myDyn.get_PSV(Q_samp,P_samp,x_samp,p_samp);
   }

   auto start = std::chrono::high_resolution_clock::now();

   if(!check_energy){
     myDyn.run_dynamics();
   }
   else{
     myDyn.check_energ_conserv();
   }

   auto stop = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = stop - start;
  
   cout << "Dynamics Run Time: " << elapsed.count() << endl;
  }

  cout << endl;
  cout << "Finished Calculations" << endl;

  return 0;
}
