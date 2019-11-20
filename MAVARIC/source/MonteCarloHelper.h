/* Author:Elliot Eklund
   MonteCarloHelper are helper functions for MonteCarlo. They are defined here
   to prevent MonteCarlo from becoming too clutered. */

#ifndef MONTECARLOHELPER
#define MONTECARLOHELPER

#include <fstream>
#include <string>
#include <valarray>
#include <vector>

inline void write_MC_data(double sgn_total,double estimator_total){

  ofstream myFile;
  myFile.open("Results/mc_data");

  myFile << sgn_total << endl;
  myFile << estimator_total << endl;

  myFile.close();

  cout << "\t Successfully saved MC data to Results." << endl;

}

inline void read_MC_data(double &sgn_totalGlobal, double &estimator_total){

  ifstream myFile;
  myFile.open("Results/mc_data");

  myFile >> sgn_totalGlobal;
  myFile >> estimator_total;

  myFile.close();

  cout << "\t Successfully read MC datat from Results." << endl;
}

/* Print system Monte Carlo acceptance ratios after calling runMC.  */
inline void print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt){
  std::cout << "\t System Acceptance Ratio: " << 100*(double)sys_steps_accpt/sys_steps << std::endl;
}

/* Print electronic Monte Carlo acceptance ratios after calling runMC.  */
inline void print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt){
  std::cout << "\t Electronic Acceptance Ratio: " << 100*(double)elec_steps_accpt/elec_steps << std::endl;
}

/* Print the average energy after calling runMC */
inline void print_avg_energy(double estimator_total, double sgn_total){
  std::cout << "\t Average Energy: " << estimator_total/sgn_total << endl;
}

/* Write estimator data to file "estimator" after calling runMC.  */
inline void write_estimator(std::valarray<double> estimator,int interval){

  int size = estimator.size();

  ofstream myFile;
  myFile.open("./Results/energy_estimator");

  for(int i=0; i<size; i++){
    myFile << i*interval << " " <<  estimator[i] << endl;
  }

  myFile.close();

  cout << "\t Successfully wrote energy_estimator file to Results." << endl;
}

/* Read in PSV.*/
inline void read_PSV(int num_beads,int num_states, std::valarray<double> Q, std::valarray<double> x,
  std::valarray<double> p){

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

  cout << "\t Successfully read PSV file from Results." << endl;
}

inline void write_PSV(int num_beads,int num_states, std::valarray<double> Q, std::valarray<double> x,
  std::valarray<double> p){

  ofstream myFile;
  myFile.open("./Results/PSV");

  myFile << num_beads << endl;
  myFile << num_states << endl;

  for(int i=0; i<num_beads; i++){
    myFile << Q[i] << endl;
  }

  for(int i=0; i<num_beads*num_states; i++){
    myFile << x[i] << endl;
  }

  for(int i=0; i<num_beads*num_states; i++){
    myFile << p[i] << endl;
  }

  cout << "\t Successfully saved PSV to Results." << endl;
  myFile.close();
}

#endif 
