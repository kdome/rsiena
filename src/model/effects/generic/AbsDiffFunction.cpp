/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: AbsDiffFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * AbsDiffFunction.
 *****************************************************************************/

#include <cstdlib>

#include "AbsDiffFunction.h"


namespace siena
{

/**
 * Creates a new function as a difference between the value of the two
 * give functions.
 */
AbsDiffFunction::AbsDiffFunction(AlterFunction * pFirstFunction,
	AlterFunction * pSecondFunction)
{
	this->lpFirstFunction = pFirstFunction;
	this->lpSecondFunction = pSecondFunction;
}


/**
 * Deallocates this function.
 */
AbsDiffFunction::~AbsDiffFunction()
{
	delete this->lpFirstFunction;
	delete this->lpSecondFunction;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void AbsDiffFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	AlterFunction::initialize(pData, pState, period, pCache);
	this->lpFirstFunction->initialize(pData, pState, period, pCache);
	this->lpSecondFunction->initialize(pData, pState, period, pCache);
}


/**
 * Does the necessary preprocessing work for calculating the alter
 * function for a specific ego. This method must be invoked before
 * calling AbsDiffFunction::value(...).
 */
void AbsDiffFunction::preprocessEgo(int ego)
{
	AlterFunction::preprocessEgo(ego);
	this->lpFirstFunction->preprocessEgo(ego);
	this->lpSecondFunction->preprocessEgo(ego);
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double AbsDiffFunction::value(int alter)
{
	return std::abs(this->lpFirstFunction->value(alter) -
		this->lpSecondFunction->value(alter));
}

}
