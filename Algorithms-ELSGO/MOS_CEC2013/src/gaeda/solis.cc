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

#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "solis.h"
#include "genomes/MOSGenome.h"
#include <algorithm>
#include <cassert>
#include <cmath>

using namespace std;

SolisWets::SolisWets(void) : m_rand(12345679) {
  m_maxdelta = m_mindelta = -1;
}

void SolisWets::setDelta(double maxdelta) {
  assert(maxdelta > 0);
  m_maxdelta = maxdelta;
}

void SolisWets::setDelta(double mindelta, double maxdelta) {
  assert(mindelta > 0 && mindelta <= maxdelta);
  m_mindelta = mindelta;
  m_maxdelta = maxdelta;
}

SolisParams *SolisWets::recoverOptions(tGen *params, unsigned size) {
  SolisParams *option;
  unsigned dim = GAEDAConfig::handle()->getProblemSize();

  option = new SolisParams(dim);
  option->recover(params, size);
  return ((SolisParams *) option);
}

void SolisWets::storeOptions(SolisParams *params, tGen **paparams, unsigned *psize) {
  unsigned dim = GAEDAConfig::handle()->getProblemSize();
  unsigned size = dim+3;

  if (params != NULL) {
    SolisParams *p = (SolisParams *) params;
    p->store(paparams, psize);
    assert(size == *psize);
  }
  else {
    *paparams = NULL;
  }

  *psize = size;
}

void SolisWets::clip(tChromosomeReal &sol) {
  for (int i = 0; i < sol.length (); i++) {
    if (sol.gene(i) > sol.alleleset(i).upper())
      sol.gene(i, sol.alleleset(i).upper());
    else if (sol.gene(i) < sol.alleleset(i).lower())
      sol.gene(i, sol.alleleset(i).lower());
  }
  return;
}

double SolisWets::normal(double desv) {
  double u1, u2;
  double result;

  do {
    u1=m_rand.rand();
  } while (u1 == 0);

  u2=m_rand.rand();
  result = desv * sqrt (-2*log(u1)) * sin (2*PI*u2);

  return result;
}

double SolisWets::distreal(const tChromosomeReal &x, const tChromosomeReal &y) {
  double dist = 0;
  unsigned size = x.size();
  assert(x.size() == y.size());

  for (unsigned i = 0; i < size; i++)
    dist += (x[i]-y[i])*(x[i]-y[i]);

  return sqrt(dist)/size;
}

double SolisWets::distanceMin(const tChromosomeReal &x, GAPopulation *pop, unsigned *posmin) {
  double dist, lowest;
  unsigned i;

  if (pop->size() == 0) throw new string("dist:Error, popsize is zero");

  *posmin = 0;
  lowest = 0;

  for (i = 0; i < pop->size(); i++) {
    MOSGenome& gen = dynamic_cast<MOSGenome&>(pop->individual(i));
    tChromosomeReal& ind = dynamic_cast<tChromosomeReal&>(*gen.getGenome(GAID::RealEncoding));
    dist = distreal(x, ind);

    if (dist > 0 && (lowest == 0 || dist < lowest) ) {
      *posmin = i;
      lowest = dist;
    }
  }

  return lowest;
}

SolisParams *SolisWets::getInitOptions(tChromosomeReal &sol, GAPopulation* m_pop) {
  SolisParams *option;
  unsigned dim = GAEDAConfig::handle()->getProblemSize();
  option = new SolisParams(dim);
  option->numFailed = option->numSuccess = 0;
  unsigned nearest;

  // if (m_pop == NULL) {
  //   assert(m_maxdelta > 0);
  //   option->delta = m_maxdelta;
  // }
  // else {
  //   /**
  //    * Calculo sigma como la mitad de la distancia al mas cercano
  //    */
  //   double dist = distanceMin(sol, m_pop, &nearest);
  //   option->delta = dist/2.0;

  //   if      (m_maxdelta > 0 && option->delta > m_maxdelta) option->delta = m_maxdelta;
  //   else if (m_mindelta > 0 && option->delta < m_mindelta) option->delta = m_mindelta;
  // }

  option->delta = GAEDAConfig::handle()->getSWDelta();

  fill(option->bias.begin(), option->bias.end(), 0.0);

  return ((SolisParams *) option);
}

tFitness SolisWets::getNeighbour(SolisParams *param, tChromosomeReal &actual, tChromosomeReal &dif, tChromosomeReal &newsol) {
  unsigned i;
  SolisParams *p = (SolisParams *) param;
  unsigned ndim = actual.size();
  //DomainRealPtr domain = m_problem->getDomain();

  for (i = 0; i < ndim; i++) {
    //if (domain->canBeChanged(i)) {
    dif.gene(i, normal(p->delta));
    newsol.gene(i, actual.gene(i) + p->bias[i] + dif.gene(i));
    //}
    //else 
    //  newsol[i] = actual[i];
  }

  clip(newsol);
  return newsol.evaluate();
}

static double increm_bias(const double &bias, const double &dif) {
  return 0.2*bias+0.4*(dif+bias);
}

static double dec_bias(const double &bias, const double &dif) {
  return bias-0.4*(dif+bias);
}

unsigned SolisWets::apply(SolisParams *params, tChromosomeReal &sol, tFitness &sol_perf, unsigned& evals, double& fit_inc_acum, unsigned& improvements, unsigned maxeval) {

  static const double maxSuccess    = GAEDAConfig::handle()->getSWMaxSuccess   ();
  static const double maxFailed     = GAEDAConfig::handle()->getSWMaxFailed    ();
  static const double adjustSuccess = GAEDAConfig::handle()->getSWAdjustSuccess();
  static const double adjustFailed  = GAEDAConfig::handle()->getSWAdjustFailed ();

  SolisParams *p = (SolisParams *) params;
  unsigned ndim = sol.size();
  //tChromosomeReal dif(ndim), newsol(ndim);
  tChromosomeReal dif(sol), newsol(sol);
  //DomainRealPtr domain = m_problem->getDomain();
  tFitness newsol_perf;
  unsigned gen;
  unsigned numEval = 0;

  fit_inc_acum = 0;

  //for (numEval = 0; numEval < maxeval && !m_running->isFinish(); ) {
  for (numEval = 0; numEval < maxeval;) {
    newsol_perf = getNeighbour(p, sol, dif, newsol);
    numEval++;
    fit_inc_acum += newsol.computeFitnessIncrement(sol_perf);

    // Si lo mejoro lo copio
    if (GAGenome::compareScores(newsol_perf, sol_perf) == GAGenome::BETTER) {
      //copy(newsol.begin(), newsol.end(), sol.begin());
      sol.copy(newsol);
      sol_perf = newsol_perf;
      improvements++;

      // Adapto bias
      //transform(p->bias.begin(), p->bias.end(), dif.begin(), p->bias.begin(), ptr_fun(increm_bias));
      for (unsigned i = 0; i < p->bias.size(); i++)
        p->bias[i] = increm_bias(p->bias[i], dif.gene(i));

      p->numSuccess++;
      p->numFailed = 0;
    }
    //else if (numEval < maxeval && !m_problem->isBetter(newsol_perf, sol_perf) && !m_running->isFinish()) {
    else if (numEval < maxeval && GAGenome::compareScores(newsol_perf, sol_perf) != GAGenome::BETTER) {
      // Invierto la direccion
      for (gen = 0; gen < ndim; gen++)
        newsol.gene(gen, sol.gene(gen) - p->bias[gen] - dif.gene(gen));

      clip(newsol);
      newsol_perf = newsol.evaluate();
      numEval++;
      fit_inc_acum += newsol.computeFitnessIncrement(sol_perf);

      if (GAGenome::compareScores(newsol_perf, sol_perf) == GAGenome::BETTER) {
        //copy(newsol.begin(), newsol.end(), sol.begin());
        sol.copy(newsol);
        sol_perf = newsol_perf;
        improvements++;

        // Disminuyo bias
        //transform(p->bias.begin(), p->bias.end(), dif.begin(), p->bias.begin(), ptr_fun(dec_bias));
        for (unsigned i = 0; i < p->bias.size(); i++)
          p->bias[i] = dec_bias(p->bias[i], dif.gene(i));

        p->numSuccess++;
        p->numFailed = 0;
      }
      else {
        // Se reduce bias a la mitad
        transform(p->bias.begin(), p->bias.end(), p->bias.begin(), bind2nd(multiplies<double>(), 0.5));
        p->numFailed++;
        p->numSuccess = 0;
      }

    }

    if (p->numSuccess >= maxSuccess) {
      p->numSuccess = 0;
      p->delta *= adjustSuccess;
    }
    else if (p->numFailed >= maxFailed) {
      p->numFailed = 0;
      p->delta *= adjustFailed;
    }

  } // De comprobar las iteraciones

  // For compatibility with other LSs
  evals = numEval;

  return numEval;
}
