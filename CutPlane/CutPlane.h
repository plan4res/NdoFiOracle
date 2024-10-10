/*--------------------------------------------------------------------------*/
/*---------------------------- File CutPlane.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the CutPlane class, which implements the NDOSolver interface
 * for NonDifferentiable Optimization Solvers, as described in NDOSlver.h,
 * using Kelley's "Cutting Plane" algorithm.
 *
 * The class requires that the function to be minimized be available under
 * the form of a FiOracle object, as described in FiOracle.h.
 *
 * \author Vitor Barbosa \n
 *         Algoritmi Research Unit \n
 *         Universidade do Minho \n
 *
 * \author Filipe Alvelos \n
 *         Grupo de Optimizacao e Investigacao Operacional \n
 *         Departamento de producao e Sistemas \n
 *         Universidade do Minho \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Vitor Barbosa, Filipe Alvelos, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __CutPlane
 #define __CutPlane  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NDOSlver.h"

#include "OsiSolverInterface.hpp"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS CutPlane -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** This class implements Kelley's Cutting Plane method for minimization of
 *  a convex function Fi() (possibly, the sum of k convex functions) over box
 *  constraints. At each step, a Restricted Master Problem (RMP) is
 *  constructed that provides a valid lower bound on the minimum value of the
 *  function using all (or part of) the first-prder information found so far,
 *  the RMP is solved and the function is evaluated in the minimum, until the
 *  gap betwen the RMP solution and the Fi solution is less than the required
 *  precision.
 *
 *  This implementation *very partially* follows the NDOSolver interface for
 *  "general" NonDifferentiable Optimization solvers, in the following sense:
 *
 *  1) Does not implement a "Phase 0", meaning that it requires the RMP to
 *     be bounded below, either by box constraints, or by the subgradients
 *     fetched at the first call, or at least by any crude but finite lower
 *     bound on the minimum of Fi().
 *  2) Does not support more complicated constraints set than boxes (although
 *     it would be very simple).
 *  3) Does not correctly handle adding/removing variables etc.
 *  4) Does not correctly handle functions Fi() that are only "approximately"
 *     computed. */

class CutPlane : public NDOSolver
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
/** @name Constructor
    @{ */

  CutPlane( std::istream *iStrm = 0 );

/**< Constructor of the class.

   The parameter `iStrm', if provided, is taken as a pointer to a istream from
   which the algorithmic parameters for the Cutting plane algorithm are
   sequentially read in the following order. Each parameter must be placed at
   the beginning of a separate line, max 255 carachters long, with all the
   rest of the line up to the first newline carachter '\n' (apart from a
   separating whitespace) being available for comments. Any line whose first
   carachter is '#' and any blank line is ignored. If 0 is passed, the file
   ends before reaching a given parameter, or some parameter is in the wrong
   format, each non-specified parameter is given a default value, shown in []
   below.

   Note that `iStrm' is passed to the constructor of NDOSolver [see
   NDOSolver.h], which reads the general algorithmic parameters out of it;
   since the constructor of the CutPlane class is executed after the one of
   NDOSolver, the following parameters specific for the CutPlane have to be
   found in the stream *after* those of the base class.

  - Index NumNames    [100] number of items (re)allocated each time it
                      is needed.

  - HpNum PurgeRows   [-1] Purge (delete) rows options:
                       < 0  never purge rows;
                       = 0  purge rows with slack greater than gap;
                       > 0  purge rows with slack greater than gap * <this>.

  - Index PurgeInvl   [5] controls how often rows are purged.

  - Index HeurSubGU   [1] number of subgradients asked to the FiOracle,
                      at every iteration and for every Fi-component,
		      when the RMP is "artificially" bounded.

  - Index HeurSubGF   [1] number of subgradients asked to the FiOracle,
                      every HeurInvl [see below] and for every
		      Fi-component, when the RMP is "naturally" bounded.

  - Index HeurInvl    [30] how often "extra" subgradients are asked to
                      the FiOracle, for every Fi-component, when the RMP
		      is "naturally" bounded.

  - Index KpPrimals   [0] if == 0, do not care for subgradient "names".

  - Index PrintInvl   [50] how often log information is printed. */

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   void SetOsiSolver( OsiSolverInterface* osi = 0 );

///< Sets the OsiSolver to be used by CutPlane to solve the RMP

/*--------------------------------------------------------------------------*/

   void SetRadius( LMNum rad = Inf< LMNum >() );

/**< Allows a very primitive form of stabilization by forcing all iterates
   from now on to lie in a ball (in the INF-Norm) of radius 2 * rad around
   the current point Lambda. This can be called before Solve(), possibly
   having already set a "smart" Lambda with SetLambda() [see below]. By
   calling SetRadius( Inf< LMNum >() ), the original box constraints (if any)
   on the variables are restored. */

/*--------------------------------------------------------------------------*/

   void SetFiOracle( FiOracle *Fi = 0 );

/*--------------------------------------------------------------------------*/

   void SetLambda( cLMRow tLambda = 0 );

/*--------------------------------------------------------------------------*/

   void KeepBestLambda( const bool KBL = true );

/*--------------------------------------------------------------------------*/

   void SetNDOLog( ostream *outs = 0 , const char lvl = 0 );

/**< Set the log file and the level of log to be used. lvl controls the
   "level of verbosity" of the code. 

   - 0  =>  no log except error messages;

   - 1  =>  just final results are written;

   - 2  =>  succint results are written every PrintInvl iterations;

   - 3  =>  results are written every iteration and are *much* more verbose.
    */

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   NDOStatus Solve( void );

/*--------------------------------------------------------------------------*/

//!!   void ReSetAlg( unsigned char RstLvl = 0 ) {}

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   cLMRow ReadSol( cIndex_Set &I , Index &D );

/*--------------------------------------------------------------------------*/

   cLMRow ReadBestSol( cIndex_Set &I , Index &D );

/*--------------------------------------------------------------------------*/
      
   HpNum ReadFiVal( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   HpNum ReadBestFiVal( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/
   
   bool IsOptimal( HpNum eps = 0 ) const;
  
/*--------------------------------------------------------------------------*/

  cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf< Index >() );
  
/*--------------------------------------------------------------------------*/

   HpNum ReadLBMult( cIndex wFi = Inf< Index >() );

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

   void AddVariables( Index NNwVrs , cLMRow IVs = 0 );

/*--------------------------------------------------------------------------*/

   void RemoveVariables( cIndex_Set whch = 0 , Index hwmny = 0 );

/*--------------------------------------------------------------------------*/

   void ChgFiV( cIndex wFi = Inf< Index >() )
   {
    throw NDOException( "CutPlane::ChgFiV not implemented yet" );
    }

/*--------------------------------------------------------------------------*/

   void ChgSbG( cIndex strt = 0 , Index stp = Inf< Index >() ,
		cIndex wFi = Inf< Index >() )
   {
    throw NDOException( "CutPlane::ChgSbG not implemented yet" );
    }

/** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~CutPlane();

///< Destructor. It *must* be virtual.

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/


 private:

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void FiAndGi( void );

/* Solve Fi subproblem and retrieve corresponding subgradients.
   The Fi problem is solved with the last Lambda found and if some
   solution is atractive, corresponding sugradients is retrieved and
   a new line is added in the RMP. The number of lines added depends
   of the number of subproblems in the Oracle. */

/*--------------------------------------------------------------------------*/

   inline Index DelRMPRows( HpNum lim );

// Delete RMP rows with slacks greater than lim.

/*--------------------------------------------------------------------------*/

   inline void OutRsts( void );

// Function to output results when LOG >=1.

/*--------------------------------------------------------------------------*/

   inline void MemDealloc( void );

// Function to free the allocated arrays. 

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   // algorithmic parameters- - - - - - - - - - - - - - - - - - - - - - - - -

   Index NumNames;     // number of items (re)allocated each time

   HpNum PurgeRows;    // purge (delete) rows options
   Index PurgeInvl;    // how often the purge rows is done.

   Index HeurSubGU;    /* number of subgradients asked to the FiOracle
			  (every iteration) when the restricted master
			  problem is artificially feasible */
   Index HeurSubGF;    /* number of subgradients asked to the FiOracle
			  (every HeurInvl - next parameter) when the
			  restricted master problem is feasible */
   Index HeurInvl;     /* how often "extra" subgradients are asked by the
			  CutPlane (when the restricted master problem is
			  feasible) */

   Index KpPrimals;    // if == 0, do not care for primal information

   Index PrintInvl;    // how often do we print log information

   // other stuff - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Index MaxNumVar;
   LMRow Lambda;           // The current point
   LMRow BestLambda;       // The best point found so far
   Index_Set Base;         // Vector of the indices of nonzeroes in Lambda
   Index BaseSize;         // Size of Base
   bool KpBstL;            // If LmbdBst has to be kept
   HpRow FiLambda;         // Fi( Lambda ) of each component of Fi
   HpRow BestFiLambda;     // Fi( BestLambda ) of each component of Fi
   HpNum AFiLambda;        // Value of Fi( Lambda )
   HpNum ABestFiLambda;    // Best value of Fi found so far.
   SgRow Gi;               // SubGradient returned by the FiOracle

   OsiSolverInterface *OsiS;  // Instance of OsiSolverInterface.

   HpNum AValueRMP;        // Value of the approximation of Fi in the RMP
   HpRow ValueRMP;         // Optimal values of the RMP, for each component
   Index RowsAdded;        // Rows added to the RMP in the current iteration
   Index MaxNRows;         // Maximum number of rows of the RMP.
   HpNum Gap;              // Gap between FiLambda and ValueRMP
   HpNum LowerBound;       // Lower Bound over Fi
   bool TrueLB;            // true if LowerBound is a "true" lower bound
			   // rather than just the "minus infinity"
   bool loptimal;          // true if the solution found is optimal to the
			   // current master problem. When adding of
			   // variables is being used, note that the solution
			   // may not be optimal to the full problem.

   // oracle solutions management related - - - - - - - - - - - - - - - - - -
   // CutPlane is responsible to give a name to each of the primal solutions
   // on which the subgradients are calculated. it is for him to "know" which
   // names are being used and which are not. the following structures are
   // relative to that management of names. also, each name that is being
   // used has a corresponding row in the RMP. the primal solution
   // associated with row i is kept in FullNames[ i ]

   Index_Set EmptyNames;   // Vector of Empty Gi names.
   Index_Set FullNames;    // Vector of Full Gi names.
   Index_Set SubPNames;    // Vector with Fi-components of subgradients

   SIndex LastFull;        // Pointer to the last Full name.
   Index FirstEmpty;       // Pointer to the first Empty name.
   Index AllocNumber;      // Number of allocations to RowNames.

/*--------------------------------------------------------------------------*/

 };  // end( class CutPlane )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

};  //end( namespace( NDO_di_unipi_it ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* CutPlane.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File CutPlane.h -----------------------------*/
/*--------------------------------------------------------------------------*/
