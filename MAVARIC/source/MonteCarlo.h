/* Author:Elliot Eklund

  MonteCarlo performs a Monte Carlo(MC) simulation using the Metropolis Hastings algorithm.
  From this object, acceptance ratios for system and electronic steps can be calculated.
  This object also allows the user to calculate energy estimators.

*/

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "functions.h"
#include "Hamiltonian.h"
#include "MonteCarloHelper.h"

#include <fstream>
#include <random>
#include <valarray>

class MonteCarlo{

  private:

    Hamiltonian H;

    std::random_device myRand;
    std::mt19937 mt;
    std::uniform_real_distribution<double> sys_dist; 
    std::uniform_real_distribution<double> elec_dist;
    std::uniform_real_distribution<double> unif_1;
    std::uniform_int_distribution<int> rand_bead;

    /* Physical parameters */
    const int num_beads; //number of ring polymer beads
    const int num_states; //number of electronic states
    const double beta_n; // beta/num_beads
    const int elec_size; //num_states * num_beads

    /* Monte Carlo parameters */
    const unsigned long long MC_steps; //number of MC steps for current simulation
    const int esti_rate; //rate at which energy estimator is collected
    const bool PSV_W; //save PSV if true
    const bool data_W; //save Monte Carlo data if true
    const bool PSV_R; //read in PSV if true
    const bool data_R; //read in Monte Carlo if true

    unsigned long long sys_steps_accpt; //system steps accepted
    unsigned long long sys_steps; //systems steps tried

    unsigned long long elec_steps_accpt; //electronic steps accepted
    unsigned long long elec_steps; //electronic steps tried

    double energy; //current energy
    int sgn_theta; //current sign of theta;
    double estimator; //current energy estimator
    std::valarray<double> estimator_t; //estimator at each step

    /* PSV Variables */
    std::valarray<double> Q; //system positions
    std::valarray<double> x; //electronic x coordinate
    std::valarray<double> p; //electronic p coordinate

    std::valarray<double> Q_prop; //Q proposed by MC algorithm
    std::valarray<double> x_prop; //x proposed by MC algorithm
    std::valarray<double> p_prop; //p proposed by MC algorithm


    /** Private Functions **/

    /* Propose a system move; accept/reject based on Metropolis Hastings criteria. */
    void sample_system();

    /* Propose an electronic move; accept/reject based on Metropolis Hastings criteria. */
    void sample_elec();

  public:

    /* MonteCarlo constructor. */
    MonteCarlo(Hamiltonian H,std::vector<double> MC_parameters,std::vector<double> sys_parameters,
      std::vector<double> elec_parameters);

    /* MonteCarlo destructor */
    ~MonteCarlo();

    /* Set PSV (Q,...) equal incoming variables (QIN,...). */
    void set_PSV(std::valarray<double> QIN, std::valarray<double> xIN, std::valarray<double> pIN);

    /* Run Monte Carlo calculation. Acccpetance ratios and energy estimator will be calculated.*/
    void runMC();

    /* Return current values stored in Q. */
    std::valarray<double> get_Q();

    /* Return current values stored in x. */
    std::valarray<double> get_x();

    /* Return current values stored in p. */
    std::valarray<double> get_p();


};

#endif
