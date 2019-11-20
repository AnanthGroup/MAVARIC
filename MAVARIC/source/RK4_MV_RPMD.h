/* Author:Elliot Eklund

   RK4_MV_RPMD implements a 4th order Runge Kutta integration scheme by using the  MV-RPMD
   Hamiltonian to provide forces. */

#ifndef RK4_MV_RPMD_H
#define RK4_MV_RPMD_H

#include "functions.h"
#include "Forces.h"

#include <iostream>
#include <valarray>

using namespace std;

class RK4_MV_RPMD{

  private:
    /* Hamiltonian used by RK4. */
    Forces F;

    /* Variables relevent to the physical system.*/
    const int num_states; // number of electronic states
    const int num_beads; // number of beads
    
    /* Variables relevent to RK4 integrator. */
    const double dt; // time step
    const double dt_half; // 0.5 * dt
    const double coeff1, coeff2, coeff3, coeff4; // coefficients for final step of RK4

    /* Four k-parameters for each degree of freedom in the system. */
    valarray<double> k1Q, k2Q, k3Q, k4Q;
    valarray<double> k1P, k2P, k3P, k4P;
    valarray<double> k1x, k2x, k3x, k4x;
    valarray<double> k1p, k2p, k3p, k4p;


    /** Private Functions **/

    /* All k parameters are filled with the appropriate number of zeros after zero_ks
       is called. */
    void zero_ks();

    /* k1 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k1(const valarray<double> &Q,const valarray<double> &P,const valarray<double> &x,
      const valarray<double> &p);

    /* k2 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k2(const valarray<double> &Q,const valarray<double> &P,const valarray<double> &x,
      const valarray<double> &p);

    /* k3 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k3(const valarray<double> &Q,const valarray<double> &P,const valarray<double> &x,
      const valarray<double> &p);

    /* k4 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k4(const valarray<double> &Q,const valarray<double> &P,const valarray<double> &x,
      const valarray<double> &p);

    /* All k parameters must be updated before calling update_final.
       This function ultimately advances all degrees of freedom using the k parameters. */
    void update_final( valarray<double> &Q, valarray<double> &P, valarray<double> &x,
       valarray<double> &p);

  public:
    /* Constructor. All variables defined above are initialized after the constructor is called. */
    RK4_MV_RPMD(Forces F,int num_beads,int num_states,double dt);

    /* Advance the current trajectory one time step forward using RK4-MV-RPMD-SB integrator. */
    void take_step(valarray<double> &Q, valarray<double> &P, valarray<double> &x, valarray<double> &p);

};

#endif
