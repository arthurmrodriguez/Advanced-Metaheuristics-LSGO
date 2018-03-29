/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
 * GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
 * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
 * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/


/* GASimpleEDA.cc: Pure EDA algorithm (implementation) */

/* INCLUDES */
#include "GASimpleEDA.h"

#include "garandom.h"
#include "GAGraphModel.h"
#include "GABayesianNetwork.h"
#include "GAGaussianNetwork.h"

/* FUNCTIONS */

/**
 * This is an static function. It provides the configuration variables along
 * with its default values
 */
GAParameterList &GASimpleEDA::registerDefaultParameters( GAParameterList &p ) {

	// The function at the base class does most of the job
	GAGeneticAlgorithm::registerDefaultParameters( p );


	// Add elitism to the list
	p.add( gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &DefaultParameters::elitism );

	// Add the selection percentage (percentage of the population used to
	// train the EDA) to the list
	p.add( gaNselectionPercentage, gaSNselectionPercentage, GAParameter::DOUBLE, &DefaultParameters::EDASelectionPerc );

	// Learning method for discrete EDAs
	p.add( gaNdiscreteLearningMethod, gaSNdiscreteLearningMethod,
			GAParameter::INT, &DefaultParameters::bayesianLearningMethod );

	// Score method for the EBNA_LOCAL learning method in
	// discrete EDAs
	p.add( gaNdiscreteEBNAScoring, gaSNdiscreteEBNAScoring,
			GAParameter::INT, &DefaultParameters::bayesianEBNAScoring );

	// Simulation method for discrete EDAs
	p.add( gaNdiscreteSimulationMethod, gaSNdiscreteSimulationMethod,
			GAParameter::INT, &DefaultParameters::bayesianSimulationMethod );

	// Learning methods for continuous EDAs
	p.add( gaNcontinuousLearningMethod, gaSNcontinuousLearningMethod,
			GAParameter::INT, &DefaultParameters::gaussianLearningMethod );

	// Scoring method for continuous EDAs
	p.add( gaNcontinuousScoreMethod, gaSNcontinuousScoreMethod,
			GAParameter::INT, &DefaultParameters::gaussianScoringMethod );

	return p;
}



/**
 * Constructs the class given an individual. The population is formed cloning
 * the given genome
 *
 * @param[in] c Genome used to build the population
 *
 */
GASimpleEDA::GASimpleEDA( GAGenome &c ) : GAGeneticAlgorithm( c ) {

	// Bayesian or gaussian network?
	if(GAGraphModel::GuessModelType( c ) == GAGraphModel::BAYESIAN_NETWORK){
		network = (GAGraphModel *) new GABayesianNetwork( c, DefaultParameters::bayesianLearningMethod,
			   DefaultParameters::bayesianEBNAScoring, DefaultParameters::bayesianSimulationMethod );

	} else if(GAGraphModel::GuessModelType( c ) == GAGraphModel::GAUSSIAN_NETWORK){
	        network= (GAGraphModel *) new GAGaussianNetwork( c, DefaultParameters::gaussianLearningMethod,
			   DefaultParameters::gaussianScoringMethod	);
	}

	// The constructor for the GA base class (GAGeneticAlgorithm())
	// has already created a population (pop) by cloning the genome c
	// We initialize the holder for the old population by cloning
	// the actual population (as created by the mentioned constructor)
	mOldPop = pop->clone();

	// are we elitists?. Note that originally this was set directly to
	// true, but I think we should use the default variable
	// (after all it is there for something, right?)
	this->el = DefaultParameters::elitism;
	params.add(gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &el);


	// Set the selection percentage
	this->mSelPerc = DefaultParameters::EDASelectionPerc;
	params.add( gaNselectionPercentage, gaSNselectionPercentage,
			GAParameter::DOUBLE, &mSelPerc);


	// Set the parameters for discrete EDAs
	params.add( gaNdiscreteLearningMethod, gaSNdiscreteLearningMethod,
			GAParameter::INT, &DefaultParameters::bayesianLearningMethod );
	params.add( gaNdiscreteEBNAScoring, gaSNdiscreteEBNAScoring,
			GAParameter::INT, &DefaultParameters::bayesianEBNAScoring );
	params.add( gaNdiscreteSimulationMethod, gaSNdiscreteSimulationMethod,
			GAParameter::INT, &DefaultParameters::bayesianSimulationMethod );

	// Set the parameters for continuous EDAs
	params.add( gaNcontinuousLearningMethod, gaSNcontinuousLearningMethod,
			GAParameter::INT, &DefaultParameters::gaussianLearningMethod );
	params.add( gaNcontinuousScoreMethod, gaSNcontinuousScoreMethod,
			GAParameter::INT, &DefaultParameters::gaussianScoringMethod );

}

/**
 * Constructs the class given an existing population.
 *
 * @param[in] p Initial population
 *
 */

GASimpleEDA::GASimpleEDA( const GAPopulation &p ) : GAGeneticAlgorithm( p ) {
	// We get an individual from the given population to check whether
	// we are using gaussian or bayesian networks
	GAGenome &c = p.individual(0);

	// Bayesian or gaussian network?
	if(GAGraphModel::GuessModelType( c ) == GAGraphModel::BAYESIAN_NETWORK) {
		network = (GAGraphModel *) new GABayesianNetwork( c, DefaultParameters::bayesianLearningMethod,
			   DefaultParameters::bayesianEBNAScoring, DefaultParameters::bayesianSimulationMethod);
	} else if(GAGraphModel::GuessModelType( c ) == GAGraphModel::GAUSSIAN_NETWORK) {
		network= (GAGraphModel *) new GAGaussianNetwork( c, DefaultParameters::gaussianLearningMethod,
			   DefaultParameters::gaussianScoringMethod	);
	}

	// We can learn from the given population
	network->learn(p, mSelPerc);

	// The constructor for the GA base class (GAGeneticAlgorithm())
	// has already created a population (pop) by cloning the genome c
	// We initialize the holder for the old population by cloning
	// the actual population (as created by the mentioned constructor)
	mOldPop = pop->clone();

	// Are we elitists?
	el = DefaultParameters::elitism;
	params.add(gaNelitism, gaSNelitism, GAParameter::BOOLEAN, &el);


	// Set the selection percentage
	this->mSelPerc = DefaultParameters::EDASelectionPerc;
	params.add( gaNselectionPercentage, gaSNselectionPercentage,
			GAParameter::DOUBLE, &mSelPerc);


	// Set the parameters for discrete EDAs
	params.add( gaNdiscreteLearningMethod, gaSNdiscreteLearningMethod,
			GAParameter::INT, &DefaultParameters::bayesianLearningMethod );
	params.add( gaNdiscreteEBNAScoring, gaSNdiscreteEBNAScoring,
			GAParameter::INT, &DefaultParameters::bayesianEBNAScoring );
	params.add( gaNdiscreteSimulationMethod, gaSNdiscreteSimulationMethod,
			GAParameter::INT, &DefaultParameters::bayesianSimulationMethod );

	// Set the parameters for continuous EDAs
	params.add( gaNcontinuousLearningMethod, gaSNcontinuousLearningMethod,
			GAParameter::INT, &DefaultParameters::gaussianLearningMethod );
	params.add( gaNcontinuousScoreMethod, gaSNcontinuousScoreMethod,
			GAParameter::INT, &DefaultParameters::gaussianScoringMethod );


}


/**
 * Class destructor
 */
GASimpleEDA::~GASimpleEDA(){

	// Note: For some reason we don't delete the current population (?)
	// Maybe it is because we need it for the results (wild guess)
	delete mOldPop;

	delete network;
}


/**
 * Implementation of the = operator
 */
GASimpleEDA &GASimpleEDA::operator= ( const GASimpleEDA& ga ){
	if(&ga != this)
		this->copy(ga);
	return *this;
}



// Regular member functions


// Copy the current algorithm to the provided parameter
void GASimpleEDA::copy ( const GAGeneticAlgorithm &ga ) {

	// The base class does most of the job
	GAGeneticAlgorithm::copy(ga);

	// Obtain an appropiate reference....
	const GASimpleEDA& gaeda = DYN_CAST(const GASimpleEDA&,ga);

	// ... so we can fill the EDA-specific values
	el = gaeda.el;

	// Manage the old population
	if(mOldPop) {
		mOldPop->copy(*(gaeda.mOldPop));
	} else  {
		mOldPop = gaeda.mOldPop->clone();
	}
	mOldPop->geneticAlgorithm(*this);
}





// Function to set configuration values.
int GASimpleEDA::setptr( const char* name, const void* value ){

	// We don't know yet if this can be handled in the
	// base class. We try it first
	int status = GAGeneticAlgorithm::setptr(name, value);

	// Check if the request can be handled here
	if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0 ) {
		el = (*((int*)value) != 0 ? gaTrue : gaFalse);
		status = 0;
	} else if (strcmp(name, gaNselectionPercentage) == 0 || strcmp(name, gaSNselectionPercentage) == 0 ) {
		mSelPerc = (*((double*)value));
		status = 0;

	} else if ( strcmp(name, gaNdiscreteLearningMethod) == 0 || strcmp(name, gaSNdiscreteLearningMethod) == 0 ||
		strcmp(name, gaNdiscreteEBNAScoring) == 0 || strcmp(name, gaSNdiscreteEBNAScoring) == 0 ||
		strcmp(name, gaNdiscreteSimulationMethod) == 0 || strcmp(name, gaSNdiscreteSimulationMethod) == 0 ||
		strcmp(name, gaNcontinuousLearningMethod) == 0 || strcmp(name, gaSNcontinuousLearningMethod) == 0 ||
		strcmp(name, gaNcontinuousScoreMethod) == 0 || strcmp(name, gaSNcontinuousScoreMethod) == 0 ) {

		status = network->setptr( name, value );
	}
	return status;
}

// Get a configuration value.
int GASimpleEDA::get( const char* name, void* value ) const {

	// As with setptr, we try with the base class first
	int status = GAGeneticAlgorithm::get(name, value);

	// Check for elitism
	if(strcmp(name, gaNelitism) == 0 || strcmp(name, gaSNelitism) == 0){
		*((int*)value) = (el == gaTrue ? 1 : 0);
		status = 0;

	} else if (strcmp(name, gaNselectionPercentage) == 0 || strcmp(name, gaSNselectionPercentage) == 0 ) {
		(*((double*)value)) = mSelPerc;
		status = 0;

	} else if ( strcmp(name, gaNdiscreteLearningMethod) == 0 || strcmp(name, gaSNdiscreteLearningMethod) == 0 ||
		strcmp(name, gaNdiscreteEBNAScoring) == 0 || strcmp(name, gaSNdiscreteEBNAScoring) == 0 ||
		strcmp(name, gaNdiscreteSimulationMethod) == 0 || strcmp(name, gaSNdiscreteSimulationMethod) == 0 ||
		strcmp(name, gaNcontinuousLearningMethod) == 0 || strcmp(name, gaSNcontinuousLearningMethod) == 0 ||
		strcmp(name, gaNcontinuousScoreMethod) == 0 || strcmp(name, gaSNcontinuousScoreMethod) == 0 ) {

		status = network->get( name, value );
	}

	return status;
}

// Set the objetive function
void GASimpleEDA::objectiveFunction(GAGenome::Evaluator f){
	// As usual, the base class does most of the job
	GAGeneticAlgorithm::objectiveFunction(f);

	// The old pop is created on this class, so
	// we make the changes here
	for(unsigned int i=0; i<pop->size(); i++)
		mOldPop->individual(i).evaluator(f);
}


// Set the objective data member on all individuals used by the genetic
// algorithm. This can be changed during the course of an evolution.
void GASimpleEDA::objectiveData( const GAEvalData &v ){
	// As usual, the base class does most of the job
	GAGeneticAlgorithm::objectiveData(v);

	for(unsigned int i=0; i<pop->size(); i++)
		// BUG??. Changed to oldpop
		//pop->individual(i).evalData(v);
		mOldPop->individual(i).evalData(v);
}

// Set/Get the population. Returns a reference to the current population.
const GAPopulation& GASimpleEDA::population(const GAPopulation &p) {

	// Stupidity test
	if(p.size() < 1) {
		GAErr(GA_LOC, className(), "population", gaErrNoIndividuals);
		return *pop;
	}

	// Let the base class to the work...
	GAGeneticAlgorithm::population( p );

	// ...and we fix the mOldPop part
	//mOldPop->copy(*pop->clone());

	GAPopulation* tmpPop = pop->clone();
	mOldPop->copy (*tmpPop);
	delete tmpPop;




	mOldPop->geneticAlgorithm(*this);

	return *pop;
}

// Set/Get the population size. This can be changed during the course of an evolution.
unsigned int GASimpleEDA::populationSize( unsigned int value ) {

	// Base class
	GAGeneticAlgorithm::populationSize(value);

	// And old population
	mOldPop->size(value);
	return value;
}


// Initialize the population,  do a few stupidity
// checks, reset the stats.  We must initialize the old pop because there is no
// guarantee that each individual will get initialized during the course of our
// operator++ operations.  We do not evaluate the old pop because that will
// happen as-needed later on.
void GASimpleEDA::initialize() {
	pop->initialize();
	pop->evaluate(gaTrue);	// the old pop will get it when the pops switch

	stats.reset(*pop);

	if(!scross)
		GAErr(GA_LOC, className(), "initialize", gaErrNoSexualMating);

	printStats("Initial Stats");

}

// Create offsprings from the current population
void GASimpleEDA::offspring(GAPopulation* offpop) {
	//TODO Aprendemos de TODA la antigua poblaci�. XXX: ??
	network->learn(*pop, mSelPerc );

	// Generate the individuals in the temporary population from individuals in
	// the main population.

	for(unsigned int i=0; i<offpop->size(); i++)
	  network->simulate(offpop->individual(i));

	stats.numeval += offpop->size();
	stats.numrep  += offpop->size();
	stats.numcro  += offpop->size();
	stats.numsel  += offpop->size();
}


//   Evolve a new generation of genomes.  When we start this routine, pop
// contains the current generation.  When we finish, pop contains the new
// generation and mOldPop contains the (no longer) current generation.  The
// previous old generation is lost.  We don't deallocate any memory, we just
// reset the contents of the genomes.
//
//   The selection routine must return a pointer to a genome from the old
// population.
void GASimpleEDA::step() {

	unsigned int i;

	offspring(mOldPop);

	for(i=0; i<mOldPop->size(); i++)
		pop->add(&mOldPop->individual(i));

	pop->evaluate();	// get info about current pop for next time
	pop->scale();		// remind the population to do its scaling

	for(i=0; i<mOldPop->size(); i++)
		mOldPop->replace(pop->remove(GAPopulation::WORST, GAPopulation::RAW), i);

	stats.update(*pop);		// update the statistics by one generation

	printStats("End of Step Stats");

}
