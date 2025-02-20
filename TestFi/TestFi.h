/*--------------------------------------------------------------------------*/
/*---------------------------- File TestFi.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * TestFi is a simple example of a "concrete" class which implements the
 * interface defined by the abstract base class FiOracle; it is mainly meant
 * to provide a template where a user can cut & paste its code to produce a
 * first simple but working FiOracle.
 *
 * For the sake of the example, the function implemented by TestFi is the
 * (convex, differentiable) quadratic function
 *
 *   Fi( Lambda ) = 1/2 \sum{i=0}^{NumVar-1} || Lambda[ i ] - L0 ||_2^2 + 1
 *
 * with (unique sub)gradient
 *
 *   Gi( Lambda )[ i ] = Lambda[ i ] - L0.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __TestFi
 #define __TestFi  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS TestFi -------------------------------*/
/*--------------------------------------------------------------------------*/

class TestFi : public FiOracle
               // the class may also derive from other class
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** NV is the number of variables, and L0 the "center" of the function;
 * these parameters are only needed for the sake of! this example: remove
 * them and replace with your own. */

 TestFi( cIndex NV , cLMNum L0 );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

 virtual HpNum Fi( cIndex wFi = Inf< Index >() ) override final;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 virtual Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
		      cIndex Name = Inf< Index >() , cIndex strt = 0 ,
		      Index stp = Inf< Index >() ) override final;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

 virtual HpNum GetLowerBound( cIndex wFi = Inf< Index >() ) override final
 {
  return( wFi > 1 ? 0 : - Inf< HpNum >() );  // the optimal value is 1
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/// destructor of the class: it must be virtual

 virtual ~TestFi();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

 LMNum Cntr;  /**< the minimum of the function (L0): this protected data
	       * structure is only needed for the sake of this example,
	       * remove it and replace with your own. */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 };  // end( class TestFi )

/*--------------------------------------------------------------------------*/

};  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* TestFi.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File TestFi.h ------------------------------*/
/*--------------------------------------------------------------------------*/
