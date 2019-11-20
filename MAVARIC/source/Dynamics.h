/* Author:Elliot Eklund */

#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "ABM_MV_RPMD.h"

#include <valarray>
#include <vector>

class Dynamics{

  private:
    ABM_MV_RPMD abm;

    /* Parameters relating to Dynamics calculation. */
    const int num_trajs; //number of trajectories
    const int num_states; //number of electronic states
    const int num_beads; //number of ring polymer beads
    const int elec_size; // num_states x num_beads
    const double dt; //time step
    const double run_time; //simulation run time

    std::valarray<double> Q; //vector bath modes size: num_trajs x num_modes x num_beads
    std::valarray<double> P; //vector bath modes size: num_trajs x num_modes x num_beads
    std::valarray<double> x; //vector elec x's size: num_trajs x num_states x num_beads
    std::valarray<double> p; //vector elec p's size: num_trajs x num_states x num_beads

    std::valarray<double> CQQ; //position auto-correlation function
    double sgn_theta_total;//sum of all sgn theta 
    bool check_energy; //check energy conservation if true
    double tol; //energy conservation tolerance

    /** Private Functions **/

    /* Update the position auto-correlation function at time "step" given centroid,
     * sgn_theta, and centroid_0.*/
    void update_CQQ(int step, double centroid, double sgn_theta, double centroid_0);

    /* Load a trajectory from (Q,P,x,p) into (Q_temp,P_temp,x_temp,p_temp). */
    void load_temp(std::valarray<double> &Q_temp,std::valarray<double> &P_temp,
      std::valarray<double> &x_temp, std::valarray<double> &p_temp,int traj);

    /* Write position auto-correlation function to file /Results/pos_auto_corr*/
    void write_CQQ();

    /* Read in trajectory positions from /Results/Trajectories/Q */
    void read_Q();
   
    /* Read in trajectory momenta from /Results/Trajectories/P */
    void read_P();

    /* Read in trajectory electronic positions from /Results/Trajectories/x */
    void read_x();

    /* Read in trajectory electronic momenta from /Results/Trajectories/p */
    void read_p();


  public:
    /* Dynamics constructor */
    Dynamics(ABM_MV_RPMD abm, std::vector<double>Dyn_params,std::vector<double> Samp_params,
      std::vector<double> sys_params, std::vector<double> elec_params);

    /* Pass PSV to Dynamics. */
    void get_PSV(std::valarray<double> QIN, std::valarray<double> PIN,
      std::valarray<double> xIN, std::valarray<double>pIN);

    /* Run algorithm for checking energy conservation. */
    void check_energ_conserv();

    /* Run Dynamics simulation and calculate requested quantities */
    void run_dynamics();

    /* Read in saved trajectories from /Results/Trajectories/ */
    void read_trajectories();

};

#endif
