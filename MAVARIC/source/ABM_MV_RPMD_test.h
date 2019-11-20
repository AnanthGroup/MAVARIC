/* Author:Elliot Eklund 

   ABM_MV_RPMD_test implements an Adams-Bashforth-Moulton fourth order predictor-corrector integration scheme
   by using the MV-RPMD Hamiltonian to provide forces.

*/

#ifndef ABM_MV_RPMD_TEST_H
#define ABM_MV_RPMD_TEST_H

#include "functions.h"
#include "RK4_MV_RPMD.h"
#include "Forces.h"

#include <iostream>
#include <valarray>

using namespace std;

class ABM_MV_RPMD_test{

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

    valarray<double> f_Q_alpha; //derivative of Hamiltonian wrt Qb, evaluated at time t-3
    valarray<double> f_Q_beta; //derivative of Hamiltonian wrt Qb, evaluated at time t-2
    valarray<double> f_Q_gamma; //derivative of Hamiltonian wrt Qb, evaluated at time t-1
    valarray<double> f_Q_delta; //derivative of Hamiltonian wrt Qb, evaluated at time t
    valarray<double> f_Q_p1; //derivative of Hamiltonian wrt Qb, evaluated at time t+1

    valarray<double> f_P_alpha; //derivative of Hamiltonian wrt Pb, evaluated at time t-3
    valarray<double> f_P_beta; //derivative of Hamiltonian wrt Pb, evaluated at time t-2
    valarray<double> f_P_gamma; //derivative of Hamiltonian wrt Pb, evaluated at time t-1
    valarray<double> f_P_delta; //derivative of Hamiltonian wrt Pb, evaluated at time t
    valarray<double> f_P_p1; //derivative of Hamiltonian wrt Pb, evaluated at time t+1

    valarray<double> f_x_alpha; //derivative of Hamiltonian wrt x, evaluated at time t-3
    valarray<double> f_x_beta; //derivative of Hamiltonian wrt x, evaluated at time t-2
    valarray<double> f_x_gamma; //derivative of Hamiltonian wrt x, evaluated at time t-1
    valarray<double> f_x_delta; //derivative of Hamiltonian wrt x, evaluated at time t
    valarray<double> f_x_p1; //derivative of Hamiltonian wrt x, evaluated at time t+1

    valarray<double> f_p_alpha; //derivative of Hamiltonian wrt p, evaluated at time t-3
    valarray<double> f_p_beta; //derivative of Hamiltonian wrt p, evaluated at time t-2
    valarray<double> f_p_gamma; //derivative of Hamiltonian wrt p, evaluated at time t-1
    valarray<double> f_p_delta; //derivative of Hamiltonian wrt p, evaluated at time t
    valarray<double> f_p_p1; //derivative of Hamiltonian wrt p, evaluated at time t+1

    valarray<double> pred_Q; //Q given by predictor phase
    valarray<double> pred_P; //P given by predictor phase
    valarray<double> pred_x; //x given by predictor phase
    valarray<double> pred_p; //p given by predictor phase

    /* Pointers */
    valarray<double> *p_Q_free;
    valarray<double> *p_Q_1;
    valarray<double> *p_Q_2;
    valarray<double> *p_Q_3;

    valarray<double> *p_P_free;
    valarray<double> *p_P_1;
    valarray<double> *p_P_2;
    valarray<double> *p_P_3;

    valarray<double> *p_x_free;
    valarray<double> *p_x_1;
    valarray<double> *p_x_2;
    valarray<double> *p_x_3;

    valarray<double> *p_p_free;
    valarray<double> *p_p_1;
    valarray<double> *p_p_2;
    valarray<double> *p_p_3;

    /** Private Functions **/

    /* All data is initialized after calling this function.*/
    void initialize();

    void initialize_pointers();

    /* Perform prediction phase of ABM method.*/
    void predict(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Perform correcton phase of ABM method. */
    void correct(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_0.*/
    void update_p_free(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-1).*/
    void update_f_beta(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-2).*/
    void update_f_gamma(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force at t_(-3).*/
    void update_f_delta(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    /* Update force of predicted step */
    void update_f_p1(valarray<double> &Q, valarray<double> &P, valarray<double> &x,
      valarray<double> &p);

    void swap();

  public:

    /* ABM_MV_RPMD_test constructor */
    ABM_MV_RPMD_test(Forces F, int num_beads,int num_states,double dt);

    /* Copy constructor */
    ABM_MV_RPMD_test(const ABM_MV_RPMD_test &ABM);

    /* RK4_MV_RPMD is used to take two steps backward and seed the ABM_MV_RPMD_test integrator.
       f_v_1, f_v_2, and f_v_3 are all set after calling this function. */
    void take_first_steps(valarray<double> &Q,valarray<double> &P,valarray<double> &x,
      valarray<double> &p);

    /* Advance the current trajectory one time step forward using ABM_MV_RPMD_test integrator. */
    void take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p);

    /* Advance the current trajectory one time step forward using ABM_MV_RPMD_test integrator. Return the sign of
       theta in sgn_theta that is passed */
    void take_sgn_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p,
      double &sgn_theta);

    void take_energy_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p,
      double &energy);
};

#endif
