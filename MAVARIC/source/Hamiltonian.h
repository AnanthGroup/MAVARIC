/* Author:Elliot Eklund

   Hamiltonian contains all functions necessary to compute MV-RPMD energy and energy estimators. 
   This object works in conjunction with MonteCarlo and Sampling. */

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "Potentials.h"
#include "Theta.h"
#include "functions.h"

#include <cmath>
#include <complex>
#include <valarray>
#include <vector>

class Hamiltonian{

  private:

    Theta myTheta; 
    Potentials V;

    /* Model Specific Parameters */
    const double mass; //system mass [a.u.]
    const int num_beads; //number of ring polymer beads
    const int num_states; //number of electronic states
    const double beta; // 1.0 / temp [a.u.^-1]

    /* Parameters defined for optimized performance*/
    const double beta_n; // beta/num_beads
    const double spring_coeff; // mass/beta_n^2
    const double ONE_beta_n; // 1.0/beta_n
    const double TWO_beta_n; // 2.0/beta_n
    const double ONE_n; //1.0/num_beads
    const double ONE_mass; //1.0/mass

    double energy; //energy of the current state of the system
    double estimator; //energy estimator of the current state of the system
    int sgn_theta; //sign of theta of the current state of the system


    /** Private functions **/

    /* Return ring polymer spring energy given Q.*/
    double spring(const std::valarray<double> &Q);

    /* Return G-term given x and p. */
    inline double G(const std::valarray<double> &x, const std::valarray<double> &p);

    /* Return derivative of spring term w.r.t Q. */
    std::valarray<double> dspring_dQ(const valarray<double> &Q);


  public:

    /* Hamiltonian Constructor  */
    Hamiltonian(double mass,int num_beads,int num_states,double beta);

    /* Hamiltonian copy-constructor */
    Hamiltonian(const Hamiltonian & h);

    /* Optimized version of calculating MV-RPMD energy. Returns MV-RPMD energy. */
    double get_energy(const std::valarray<double> &Q, const std::valarray<double> &x, 
      const std::valarray<double> &p);

    /* Both energy and estimator are updated after this function is called. This is an optimization
       for when the estimator needs to be calculated. It reduces redundant calculations that would
       occur if energy and estimator where to be calculated separately.*/
    void update_energy(const std::valarray<double> &Q, const std::valarray<double> &x, 
      const std::valarray<double> &p);

    /* Return energy variable.*/
    double get_energy();

    /* Return estimator variable. */
    double get_estimator();

    /* Return sgn_theta */
    int get_sgn();
};

#endif
