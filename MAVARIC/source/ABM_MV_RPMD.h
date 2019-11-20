/* Author:Elliot Eklund 

   ABM_MV_RPMD implements an Adams-Bashforth-Moulton fourth order predictor-corrector integration scheme
   by using the MV-RPMD Hamiltonian to provide forces.

*/

#ifndef ABM_MV_RPMD_H
#define ABM_MV_RPMD_H

#include "functions.h"
#include "RK4_MV_RPMD.h"
#include "Forces.h"

#include <iostream>
#include <valarray>

using namespace std;

class ABM_MV_RPMD{

  private:
    Forces F;
    RK4_MV_RPMD my_rk4;
   
    /* Variables relevent to the physical system.*/
    const int num_states; // number of electronic states
    const int num_beads; // number of beads

    /* Variables relevent to RK4 integrator. */
    const double dt; // time step

    const double h1_p; //predictor coefficient 1
    const double h2_p; //predictor coefficient 2
    const double h3_p; //predictor coefficient 3
    const double h4_p; //predictor coefficient 4

    const double h1_c; //corrector coefficient 1
    const double h2_c; //corrector coefficient 2
    const double h3_c; //corrector coefficient 3
    const double h4_c; //corrector coefficient 4

    valarray<double> f_Q_3; //derivative of Hamiltonian wrt Qb, evaluated at time t-3
    valarray<double> f_Q_2; //derivative of Hamiltonian wrt Qb, evaluated at time t-2
    valarray<double> f_Q_1; //derivative of Hamiltonian wrt Qb, evaluated at time t-1
    valarray<double> f_Q_0; //derivative of Hamiltonian wrt Qb, evaluated at time t
    valarray<double> f_Q_p1; //derivative of Hamiltonian wrt Qb, evaluated at time t+1

    valarray<double> f_P_3; //derivative of Hamiltonian wrt Pb, evaluated at time t-3
    valarray<double> f_P_2; //derivative of Hamiltonian wrt Pb, evaluated at time t-2
    valarray<double> f_P_1; //derivative of Hamiltonian wrt Pb, evaluated at time t-1
    valarray<double> f_P_0; //derivative of Hamiltonian wrt Pb, evaluated at time t
    valarray<double> f_P_p1; //derivative of Hamiltonian wrt Pb, evaluated at time t+1

    valarray<double> f_x_3; //derivative of Hamiltonian wrt x, evaluated at time t-3
    valarray<double> f_x_2; //derivative of Hamiltonian wrt x, evaluated at time t-2
    valarray<double> f_x_1; //derivative of Hamiltonian wrt x, evaluated at time t-1
    valarray<double> f_x_0; //derivative of Hamiltonian wrt x, evaluated at time t
    valarray<double> f_x_p1; //derivative of Hamiltonian wrt x, evaluated at time t+1

    valarray<double> f_p_3; //derivative of Hamiltonian wrt p, evaluated at time t-3
    valarray<double> f_p_2; //derivative of Hamiltonian wrt p, evaluated at time t-2
    valarray<double> f_p_1; //derivative of Hamiltonian wrt p, evaluated at time t-1
    valarray<double> f_p_0; //derivative of Hamiltonian wrt p, evaluated at time t
    valarray<double> f_p_p1; //derivative of Hamiltonian wrt p, evaluated at time t+1

    valarray<double> pred_Q; //Q given by predictor phase
    valarray<double> pred_P; //P given by predictor phase
    valarray<double> pred_x; //x given by predictor phase
    valarray<double> pred_p; //p given by predictor phase


    /** Private Functions **/

    /* All data is initialized after calling this function.*/
    void initialize();

    /* Perform prediction phase of ABM method.*/
    void predict(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Perform correcton phase of ABM method. */
    void correct(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_0.*/
    void update_f_0(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-1).*/
    void update_f_1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-2).*/
    void update_f_2(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-3).*/
    void update_f_3(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force of predicted step */
    void update_f_p1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);


  public:

    /* ABM_MV_RPMD constructor */
    ABM_MV_RPMD(Forces F, int num_beads,int num_states,double dt);

    /* Copy constructor */
    ABM_MV_RPMD(const ABM_MV_RPMD &ABM);

    /* RK4_MV_RPMD is used to take two steps backward and seed the ABM_MV_RPMD integrator.
       f_v_1, f_v_2, and f_v_3 are all set after calling this function. */
    void take_first_steps(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Advance the current trajectory one time step forward using ABM_MV_RPMD integrator. */
    void take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p);

    /* Advance the current trajectory one time step forward using ABM_MV_RPMD integrator. Return the sign of
       theta in sgn_theta that is passed */
    void take_sgn_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p,
      double &sgn_theta);

    /* Advnace the current trajectory one time step forward using ABM_MV_RPMD integrator. Return the energy of
     * of the trajectory in energy that is passed */
    void take_energy_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p,
      double &energy);
};

#endif
