/* Author:Elliot Eklund

   Forces contains everything needed to calculate forces corresponding to an MV-RPMD
   Hamiltonian. Forces also contains functions need to compute the energy at a given
   trajectory.

*/

#ifndef FORCES_H
#define FORCES_H

#include "functions.h"
#include "Potentials.h"
#include "Theta.h"

#include <valarray>
#include <vector>

class Forces{

  private:
    Theta myTheta;
    Potentials V;

    /* Model specific variables */
    const int num_beads; //number of ring polymer beads
    const int num_states; //number of electronic states
    const double beta; //1.0/temp
    const double beta_n; // beta/num_beads
    const double mass; //system mass
    double sgn_theta; // sign of theta

    /* Parameters defined for efficiency */
    const double spring_coeff; //num_beads^2 * mass /(beta^2)
    const double TWO_beta_n; //2.0/beta_n
    const double ONE_beta_n; //1.0/beta_n
    const double ONE_mass; //1.0/mass
    const double ONE_TWO_m; //1.0/(2.0 * mass)

    /* Return the derivative of the spring term w.r.t Q*/
    std::valarray<double> dspring_dQ(const valarray<double> &Q);

    /* Return energy of spring term */
    double spring(const valarray<double> &Q);

    /* Return energy of the G term. */
    double G(const valarray<double> &x, const valarray<double> &p);

  public:

    /* Forces constructor */
    Forces(int num_beads, int num_states, double beta, double mass);

    /* Forces copy constructor */
    Forces(const Forces &F);

    /* Updates all data necessary for dH_dQ, dH_dP, dH_dp, and dH_dx to be called. Note: This 
       function should always be called before calling the functions previously listed! */
    void update_forces(const std::valarray<double> &Q, const std::valarray<double> &x,
      const std::valarray<double> &p);

    /* Return the derivative of the MV-RPMD hamiltonian w.r.t Q.*/
    std::valarray<double> dH_dQ(const std::valarray<double> &Q, const std::valarray<double> &x,
      const std::valarray<double> &p);
  
    /* Return the derivative of the MV-RPMD hamiltonian w.r.t P.*/
    std::valarray<double> dH_dP(const std::valarray<double> &P);
  
    /* Return the derivative of the MV-RPMD hamiltonian w.r.t p.*/
    std::valarray<double> dH_dp(const std::valarray<double> &Q, const std::valarray<double> &x,
      const std::valarray<double> &p);
  
    /* Return the derivative of the MV-RPMD hamiltonian w.r.t x.*/
    std::valarray<double> dH_dx(const std::valarray<double> &Q, const std::valarray<double> &x,
      const std::valarray<double> &p);

    /* Return sgn_theta data member. */
    double get_sgn_theta();

    /* Return the energy of the MV-RPMD Hamiltonian given a trajectory */
    double get_energy(const std::valarray<double> &Q, const std::valarray<double> &P,
      const std::valarray<double> &x, const std::valarray<double> &p);

};

#endif
