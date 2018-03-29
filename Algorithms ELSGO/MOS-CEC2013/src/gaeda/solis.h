/**
* Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
*
* This file is part of software Realea
*
* Realea is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Realea is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _SOLIS_H
#define _SOLIS_H 1

#include "solisrand.h"
#include "genomes/GA1DArrayGenome.h"

#define tChromosomeReal GA1DArrayAlleleGenome<double>
#define tFitness long double
#define tGen double
#define tReal double

class SolisParams {
public:
  SolisParams(unsigned dim) : bias(dim) {}

  unsigned dim(void) {
    return bias.size();
  }

  double delta;
  vector<double> bias;
  unsigned numFailed;
  unsigned numSuccess;

  virtual ~SolisParams(void) {}

  virtual void store(tGen **aparams, unsigned *psize) {
    unsigned size = 3+bias.size();
    tGen *params = new tGen[size];
    params[0] = delta;
    params[1] = numFailed;
    params[2] = numSuccess;
    copy(bias.begin(), bias.end(), params+3);
    *aparams = params;
    *psize = size;
  }

  virtual void recover(tGen *params, unsigned size) {
    assert(size > 3);
    delta = params[0];
    numFailed = (unsigned) ceil(params[1]);
    numSuccess = (unsigned) ceil(params[2]);
    copy(params+3, params+size, bias.begin());
  }
};

class SolisWets {
public:

  SolisWets(void);
  void setDelta(double maxdelta);
  void setDelta(double mindelta, double maxdelta);

  SolisParams *getInitOptions(tChromosomeReal &sol, GAPopulation* m_pop);
  unsigned apply(SolisParams *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned& evals, double& fit_inc_acum, unsigned& improvements, unsigned maxeval);

private:
  SolisParams *recoverOptions(tGen *params, unsigned size);

  /**
   * Store in a real vector the values of a parameters. 
   *
   * @param params parameter to store
   * @param paparams reference to sequencial parameters (can be NULL)
   * @param size reference to size (always is right, even if params == NULL)
   *
   * @see recoverOptions
   */
  void storeOptions(SolisParams *params, tGen **paparams, unsigned *psize);

private:
  /**
   * Obtain a neighbour solution from a solution given.
   *
   * @param p Solis Wets parameters
   * @param sol current solution
   * @param dif  calculate ddifferences vector, output
   * @param newsol new solution, output
   *
   * @return fitness of newsol
   */
  tFitness getNeighbour(SolisParams *param, tChromosomeReal &sol, tChromosomeReal &dif, tChromosomeReal &newsol);

  void clip(tChromosomeReal &sol);

  double normal(double desv);

  double distreal(const tChromosomeReal &x, const tChromosomeReal &y);
  double distanceMin(const tChromosomeReal &x, GAPopulation *pop, unsigned *posmin);

  double m_maxdelta;
  double m_mindelta;

  SRandom m_rand;
};

#endif
