/* Author:Elliot Eklund 

  Sampling is used to generate trajectories to be used during the Dynamics calculation.
  It should be used after bead convergence and equilibrium has been established with a
  MonteCarlo object. From Sampling, the user can histogram distributions.
*/

#ifndef SAMPLING_H
#define SAMPLING_H

#include "functions.h"
#include "Hamiltonian.h"

#include <ctime>
#include <fstream>
#include <random>
#include <valarray>

class Sampling{

  private:

    Hamiltonian H;

    std::random_device myRand;
    std::mt19937 mt;
    std::uniform_real_distribution<double> sys_dist; /* sample Unif(-sys_SS,sys_SS), SS = step size */
    std::uniform_real_distribution<double> elec_dist; /* sample Unif(-elec_SS,elec_SS))*/
    std::uniform_real_distribution<double> unif_1; /* sample Unif(0,1) */
    std::uniform_int_distribution<int> rand_bead; /* generate random bead between [0,num_beads-1] */

    /* Physical parameters */
    const int num_beads; //number of ring polymer beads
    const int num_states; //number of electronic states
    const double beta; //1.0/beta
    const double beta_n; // beta/num_beads
    const double mass; //system mass
    const int elec_size; //num_states * num_beads

    /* Sampling parameters */
    const int num_trajs; //number of trajectories that will be sampled
    const int decorr_len; //decorrelation between two adjacent trajectories
    double energy; //current energy

    bool get_hist_pos; //histogram_pos will be called if this is true
    bool save_samp; //save sampled trajectories if true
    bool read_PSV; //read in PSV if true
    int num_bins; //number of bins for histogram

    /* PSV Variables */
    std::valarray<double> Q; //system positions
    std::valarray<double> x; //electronic x coordinate
    std::valarray<double> p; //electronic p coordinate

    std::valarray<double> Q_prop; //Q proposed by MC algorithm
    std::valarray<double> x_prop; //x proposed by MC algorithm
    std::valarray<double> p_prop; //p proposed by MC algorithm

    std::valarray<double> Q_samp; //collection of final system positions
    std::valarray<double> P_samp; //collection of final system momenta
    std::valarray<double> x_samp; //collection of final electronic x coordinates
    std::valarray<double> p_samp; //collection of final electronic p coordinates
    

    /** Private Functions **/

    /* Propose a system move; accept/reject based on Metropolis Hastings criteria. */
    void sample_system();

    /* Propose an electronic move; accept/reject based on Metropolis Hastings criteria. */
    void sample_elec();

    /* Sample the boltzman distribution and generate a new momenta.*/
    void sample_boltzman();

    void histogram_pos();

  public:

    /* Sampling constructor*/
    Sampling(Hamiltonian H,std::vector<double> MC_parameters,std::vector<double> sys_parameters,
      std::vector<double> elec_parameters);

    /* Sampling destructor */
    ~Sampling();


    /* Set PSV (Q,...) equal incoming variables (QIN,...). */
    void set_PSV(std::valarray<double> QIN, std::valarray<double> xIN,std::valarray<double> pIN);

    /* Run sampling procedure. After runSampling is called, Q_samp, P_samp, x_samp, p_samp will filled
       with num_trajs initial conditions. These can be used for a Dynamics calculation.*/
    void runSampling();

    /* Save sampled trajectories to Results/Trajectories/  */
    void save_trajectories();

    /* Save trajectory positions.  */
    void save_Q();

    /* Save trajectory momenta. */
    void save_P();

    /* Save trajectory electronic positions */
    void save_x();

    /* Save trajectory electronic momenta */
    void save_p();

    /* Return the positions of the sampled trajectories. */
    std::valarray<double> get_Q();

    /* Return the momenta of the sampled trajectories. */
    std::valarray<double> get_P();

    /* Return the electronic positions of the sampled trajectories. */
    std::valarray<double> get_x();

    /* Return the electronic momenta of the sampled trajectories. */
    std::valarray<double> get_p();

    /* Read in PSV varibes stored in Results/ and assign them to 
     * Q, P, x, and p. */
    void read_PSV_fxn();

};

#endif
