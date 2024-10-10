/*--------------------------------------------------------------------------*/
/*--------------------------- File MPSolver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * MPSolver is an abstract base class which defines the interface for
 * algorithms solving the (Stabilized) Master Problem of "Generalized
 * Bundle" algorithms for NonDifferentiable Optimization. Objects of derived
 * classes, actually implementing the solvers, are thought to be used by an
 * algorithm conforming to the NDOSolver/FiOracle abstract interface [see
 * NDOSolver.h and FiOracle.h], such as the Bundle class [see Bundle.h].
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

#ifndef __MPSolver
 #define __MPSolver  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS MPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class MPSolver provides a standard interface between "Generalized
 *  Bundle" algorithms for NonDifferentiable Optimization and the (complex)
 *  algorithms that solve their Master Problems.
 *
 *  The user is assumed to be familiar with generalized Bundle algorithms:
 *  refer to
 *
 *    A. Frangioni "Generalized Bundle Methods"
 *    SIAM Journal on Optimization 13(1), p. 117-156, 2002
 *
 *  available at
 *
 *   \link ftp://ftp.di.unipi.it/pub/project/orgroup/frangio/GenBundle.pdf
 *   \endlink
 *
 *  As a quick description, the data of the Master Problem is a set B of 
 *  vectors z[ i ] that are alfa_{Lambda}[ i ]-subgradients of some (convex)
 *  function Fi in some current point Lambda. This set of vectors (bundle)
 *  is used to construct a "model" Fi_{B,Lambda} of the "translated function"
 *
 *    Fi_{Lambda}( d ) = Fi( Lambda + d ) - Fi( Lambda )
 *
 *  such that Fi_{Lambda}( 0 ) = 0. The "classical" model is the (translated)
 *  Cutting Plane model of Fi in the current point Lambda
 *
 *  \hat{Fi}_{B,Lambda}( d ) = max{ z[ i ] * d - alfa_{Lambda}[ i ], i \in B }
 *
 *  that is a polyhedral lower approximation of the "true" (translated)
 *  function Fi_{Lambda}( d ). alfa_{Lambda}[ i ] (>= 0) is called the
 *  linearization error of the subgradient z[ i ] in the point Lambda, and
 *  it depends linearly on Lambda.
 *
 *  Fi_{B,Lambda} is used to find a tentative descent direction for Fi at
 *  Lambda by solving the (Stabilized) Primal Master Problem
 *
 *    P_{B,Lambda,t}:   inf{ Fi_{B,Lambda}( d ) + D_t( d ) }
 *
 *  where D_t() is a family of "stabilizing functions", enjoying some
 *  properties [see the reference], which is indiced over the "proximal
 *  parameter" t, a positive real number. P_{B,Lambda,t} has a (Fenchel)
 *  dual, the (Stabilized) Dual Master Problem
 *
 *    D_{B,Lambda,t}:   inf{ Fi_{B,Lambda}*( z ) + D_t*( -z ) }
 *
 *  where "*" indicates the Fenchel's conjugate operator.
 *
 *  A number of different Master Problems may be constructed by properly
 *  choosing the actual model Fi_{B,Lambda} and the "stabilizing term" D_t:
 *  this base class is designed to define a common interface for all them.
 *
 *  The above presentation is only schematic: the interface of MPSolver offers
 *  support for advanced features like constraints on the primal space,
 *  "disaggregated bundles" for decomposable functions, on-line variables
 *  generation and active-set variable reduction strategies, sparse
 *  subgradients and others. For instance, if box constraints
 *
 *    LB <= Lambda <= UB
 *
 *  are imposed on (some of the) variables Lambda [see GetUC() and GetUB() in
 *  FiOracle.h], then the Stabilized Primal Master Problem becomes
 *
 *    P_{B,Lambda,t}:   inf{ Fi_{B,Lambda}( d ) + D_t( d ) ,
 *                           LB - Lambda <= d <= UB - Lambda }
 *
 *  since it must guarantee that LB <= Lambda + d <= UB, being its variable d
 *  a displacement from Lambda. For further details, see the references and
 *  the interface below, as well as the interface of Bundle, NDOSolver and
 *  FiOracle [in Bundle.h, NDOSolver.h, and FiOracle.h, respectively]. */

class MPSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

/** Public enum, meant to be used in SetPar() [see below] for selecting the
    algorithmic parameter to be set. kLastMPParam is provided for allowing
    derived classes to extend the set of parameters they can handle while
    keeping the same signature for SetPar(): by setting
    kDerivedClassFirstParam = kLastMPParam, one is guaranteed that the two
    sets of parameters will not collide. */

   enum MPParam { kMaxTme = 0 ,
		  kOptEps , kFsbEps , // kZero ,
		  kLastMPParam
                  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Public enum describing the possible return values of SolveMP() [see
    below]. */

   enum MPStatus { kOK = 0 ,
		   kUnbndd , kUnfsbl ,
		   kStppd , kError
                   };

/** @} ---------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class: it should only have few initializations to do,
 * since the actual data structures can only be constructed after that the
 * "size" of the problem is known [see SetDim()]. */

 MPSolver( void ) : MaxBSize( 0 ) , MaxSGLen( 0 ) , CrrSGLen( 0 ) ,
    NrFi( 1 ) , MPLog( nullptr ) , MPLLvl( 0 ) , MPt( nullptr ) {}

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

/** Provides the MPSolver with basic information about the "size" of the
 * Master problem:
 *
 * - MxBSz is the maximum number of different items (subgradients and
 *   constraints) that can be managed by the class, i.e. the maximum size of
 *   the bundle: the protected field MaxBSize is offered by  the base class
 *   to store this information;
 *
 * - Oracle is a pointer to the FiOracle object that computes the function
 *   Fi to be minimized: the other data describing the size of the problem
 *   can be obtained by calls to methods of the FiOracle class, and three
 *   protected fields are offered by the base class to store some of this
 *   information:
 *
 *   = MaxSGLen the maximum length of any item, i.e., the maximum number of
 *     variables in the Primal Master Problem (PMP);
 *
 *   = CrrSGLen (<= MaxSGLen) is the current length of any item: the
 *     interface supports on-line creation and destruction of primal
 *     variables;
 *
 *    = NrFi the function Fi to be minimized may be "decomposable", i.e.,
 *      the sum of k functions; in this case, its model Fi_{B,Lambda} is
 *      usually decomposable as well (for instance, the Cutting Plane model
 *      is), For decomposable functions, it is possible (although not
 *      necessary) to keep k "separate bundles", each one describing a model
 *      Fi[ k ]_{B,Lambda} of each component Fi[ k ] of Fi, which may be
 *      beneficial to improve the "exactness" of the model. Special support
 *      is offered for the case when one of the components is a linear
 *      function, by considering it as the "0-th component" of Fi [see
 *      ReadFiBLambda() and SetItem() and GetItem()].
 *
 * - UsAvSt is true if the Active Set technique will be exploited by the
 *   solver, and false otherwise. In the first case, the initial "active set"
 *   of variables is *empty*, and a nonempty set of active variables must be
 *   set with SetActvSt(). In the second case, SetActvSt() *must not* be
 *   called, and all the variables are considered to be always "active".
 *
 * This information may be necessary for allocating the data structures of
 * the class, so this method *must* be called prior to any other method.
 * Other information from the FiOracle is likely to be useful for the
 * MPSolver, such as the max expected density of the subgradients [see
 * FiOracle::GetMaxNZ()], which variables are constrained in sign or has
 * upper bounds [see FiOracle::Get[UC/UB]()], and which of the components of
 * Fi(), if any, is "easy" [see FiOracle::GetBNC() ...].
 *
 * \note "easy" components of Fi() have basically to be dealt with by the
 *       MPSolver, by inserting their description in the Master Problem;
 *	 if this is not possible, the Master Problem has to signal it (e.g.
 *	 by throwing an exception) because the Bundle algorithm relies on
 *	 this and the FiOracle is not going to provide "ordinary" black-box
 *	 information (function values, subgradients, ...) for these
 *	 components.
 *
 * \note *Important*: "easy" components of Fi() in the Master Problem are
 *       *not* translated *by value* (but they are by argument), meaning that
 *	 the function that is included in the MP is
 *
 *	   Fi_{Lambda}[ k ]( d ) = Fi[ k ]( Lambda + d )
 *
 *	 and *not*
 *
 *	   Fi_{Lambda}[ k ]( d ) = Fi[ k ]( Lambda + d ) - Fi[ k ]( Lambda )
 *
 *	 as one would expect by analogy with the ordinary "difficult"
 *	 components. This is done to spare the Bundle algorithm with the
 *	 need to compute Fi[ k ]( Lambda ) for all "easy" components k and
 *	 each current point Lambda, which may be problematic especially for
 *	 the *initial* current point---consider the case where Lambda *is
 *	 not feasible*, so Fi[ k ]( Lambda ) = +INF! However, the value of
 *	 Fi[ k ] is actually computed by the MPSolver at Lambda + d* (d*
 *	 being the optimal solution of the Primal MP), so the value of
 *	 Fi[ k ] is known for the current point at least after every Serious
 *	 Step with step 1 along d*. Yet, this choice has some impact on the
 *	 "output" methods [see ReadFiBLambda() and ReadSigma() below].
 *
 * This method can be called more than once to modify the settings, but expect
 * the implementation to be quite expensive in time and/or memory. There are
 *  essentially three different ways for calling SetDim():
 *
 * - SetDim( 0 , ... ) makes the MPSolver to deallocate all its memory and to
 *   quietly wait for new instructions.
 *
 * - SetDim( n , 0 , ASV ) with n != 0 sets the max bundle size to n and
 *   activate/deactivate the Active Set Mechanism without changing anything
 *   else. The existing items in the bundle (if any) are all kept if n is
 *   larger than the previous setting, but a smaller value will force deletion
 *   of all the items with "name" [see [Get/Set]Item() and RmvItem() below]
 *   greater than or equal to n. Also, if the active set technique is
 *   initialized (ASV == true while it was false previously), the set of
 *   "active" variables is set to be *empty*.
 *
 * - SetDim( n , Orcl , ASV ) with n != 0 and Orcl != 0 discards all the
 *   previous settings and re-allocates everything; all the existing items in
 *   the bundle are lost. Note that such a call also resets every parameter
 *   of the algorithm, such as the starting point (which is set to 0).
 *
 * In general, calling this method with a non-empty bundle could be costly, so
 * if the items are to be discarded anyway, this should be done *before* the
 * call to SetDim() [see RmvItems() below]. */

   virtual void SetDim( cIndex MxBSz = 0 , FiOracle *Oracle = nullptr ,
			const bool UsAvSt = false )
   {
    MaxBSize = MxBSz;

    if( Oracle ) {
     MaxSGLen = Oracle->GetMaxNumVar();
     CrrSGLen = Oracle->GetNumVar();
     NrFi = Oracle->GetNrFi();
     }
    }

/*--------------------------------------------------------------------------*/
/** Sets the value of the proximal parameter, t, to tt. Prior to the first
 * call to Sett(), t is taken to be equal to 1. */

   virtual void Sett( cHpNum tt = 1 ) = 0;

/*--------------------------------------------------------------------------*/
/** @defgroup Prcsn Setting algorithmic parameters for the MP solution

    The solution of the Master Problem may be a costly task; hence, especially
    in the first iterations it may just not be smart to solve it exactly. This
    requires properly setting algorithmic parameters of the MPSolver. At least
    the algorithmic parameters defined in the base class need to be handled by
    the derived ones, although they can possibly add their own. @{ */

/** Change "float" algorithmic parameters of the MPSolver. This method is
   pure virtual, i.e., it has to be defined in derived classes, which may
   extend this it to allow setting the algorithmic-specific parameters.

   The enum MPParam is used for selecting the parameter to be set, with the
   following intended meaning:

   - kMaxTme: maximum running time (in seconds) for each call to SolveMP()
              [see below]. A 0 or negative value means "no time limit",
	      which should be the default. If more time required, SolveMP()
	      should stop returning kStppd. Note that this check can only be
	      performed if timing of the code is activated [see SetMPTime()
	      below].

   - kOptEps: tells the solver to solve the MP only with *at least* the
              specified relative precision (>= 0), the exact meaning of this
	      being solver-dependent. The solver can always decide to use
	      higher precision than that, if so it prefers. By default a
	      solver-specific "minimal" precision is used.

   - kFsbEps: tells the solver to tolerate the specified relative violation 
              of in the (primal) constraints; the solver can always decide
	      to use higher precision than that, if so it prefers. By default
	      a solver-specific "minimal" precision is used. 

   Note that the type of `wp' is int rather than MPParam, which is not a
   problem as enums can always be used where intsare required (i.e.,
   SetPar( MPSolver::kMaxTime , ... ) is a valid call). This is done in order
   to allow derived classes to extend the set of parameters they can handle
   while keeping the same signature for SetPar(). */

   virtual void SetPar( const int wp , cHpNum value ) = 0;

/*--------------------------------------------------------------------------*/
/** Set the maximum number of threads that the MPSolver is allowed to use.
 * This method is given a default empty implementation (doing nothing) for
 * those MPSolver that do no have any parallel processing capability. */

   virtual void SetThreads( int nthreads ) {}
   
/** @} ---------------------------------------------------------------------*/
/** Sets a lower bound on the value of the (translated) model Fi_{B,Lambda};
   since Fi_{B,Lambda}( 0 ) = 0, LwBnd must clearly be a nonpositive value.

   \note The translated model is Fi_{B,Lambda} is only translated w.r.t. the
         function value in Lambda *of the non-easy components*, while the
	 easy components retain their original value. This means that LwBnd
	 must be the original value of the lower bound minus the value
         \sum_{ k is not an easy component } Fi[ k ]( Lambda )

   \note In the above formula, the 0-th linear component is considered to be
         a non-easy component, i.e., its value has to be subtracted from the
	 original value of the lower bound. There can hardly be an easier
	 component than the 0-th one, but the variable in the master problem
	 is d = Lambda1 - Lambda and the term in the objective function is
	 RHS * d, therefore the 0-th component is naturally translated
	 w.r.t. Fi[ 0 ]( Lambda ) like the non-easy ones.

   When Fi is a decomposable function [see SetDim() above], each of its
   components may have a separate model Fi[ k ]_{B,Lambda}, thus this method
   has different meanings according to the value of wFi. For 1 <= wFi <= NrFi,
   the lower bound is on the *individual* model Fi[ wFi ]_{B,Lambda}, while
   for wFi > NrFi the lower bound is on the *global* model
   Fi_{B,Lambda} = sum_{k = 0 ... NrFi} Fi[ k ]_{B,Lambda}. Note that the
   0-th component is taken into account for the global model, while, being
   a linear function, it clearly cannot have any "individual" lower bound.

   \note Individual lower bounds are meaningless for "easy" components of
         Fi() [see FiOracle::GetBNC() ...], and therefore they are not
	 allowed.

   These lower bounds are automatically updated when the current point Lambda
   is changed [see ChangeCurrPoint() below], and they are equivalent to
   inserting in the bundle of the proper component an all-0 subgradient with
   alfa_{Lambda}[ i ] = - LwBnd (which is nonnegative). The same holds for
   the "global" lower bound, except that it is an "aggregated" all-zero
   subgradient.

   Passing - Inf< HpNum >() as the lower bound means that no such bound is
   available; this is the default. */

   virtual void SetLowerBound( cHpNum LwBnd = - Inf< HpNum >() ,
			       cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** When a new item is inserted in the bundle, checks can be done to
   ensure that it is not identical to other items already present, i.e., an
   useless duplicate [see CheckItem() below]. Since these checks can be
   costly, they are done only if CheckIdentical( true ) is called; by
   default, or if CheckIdentical( false ) is called, they are not. */

   virtual void CheckIdentical( const bool Chk = true ) {
    // the derived class may ignore the setting of CheckIdentical(), that's
    // why a default implementation doing nothing is provided
    }

/*--------------------------------------------------------------------------*/
/** The class outputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetMPLog( std::ostream *outs = 0 , const char lvl = 0 ) {
    if( ( MPLog = outs ) )
     MPLLvl = lvl;
    else
     MPLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** SetMPTime() allocates an OPTtimers object [see OPTtypes.h] that should be
   used for timing the calls to relevant methods of the class. The time can be
   read with MPTime() [see below]. By default, or if SetMPTime( false ) is
   called, no timing is done. Note that, since all the relevant methods of the
   class are pure virtual, MPSolver can only manage the OPTtimers object, but
   it is due to derived classes to actually implement the timing.

   Note that time accumulates over the calls: calling SetMPTime(), however,
   resets the counters, allowing to time specific groups of calls. */

   virtual void SetMPTime( const bool TimeIt = true ) {
    if( TimeIt )
     if( MPt )
      MPt->ReSet();
     else
      MPt = new OPTtimers();
    else
     delete MPt;
    }

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

/** Attempts to solve the current Master Problem.

   \note The bundle for some components of Fi may be empty; if this happens,
         those components must just be temporarily ignored (like they did not
         exist) in the solution of the MP. If all the "sub-bundles" for each
         of the components are empty, what is solved is just a "pure
         feasibility" MP which returns just any point in the polyhedron
         defined by the constraints, if any. If there are neither subgradient
         nor constraints, then the optimal (both primal and dual) solution
         must be all-0.

   \note The last sentence is not valid if there there are "easy" components
         of Fi() [see FiOracle::GetBNC() ...], since in this case the primal
	 solution d* of the MP is typically non-0 even if there is no
	 subgradient for any of the "difficult" components; what happens is
	 that the bundle is not empty since "those of the easy components
	 always contain all the possible subgradients".

   The method must return:

   kOK       if the solver found an optimal solution;

   kUnbndd   if the MP is primal unbounded: this can only happen with some
             "not-really-stabilizing" terms D_t, and should always be
	     solvable by properly increasing the proximal parameter t;

   kUnfsbl   if the MP is primal unfeasible (dual unbounded or empty, see
             ReadMult() below);

   kStppd    the algorithm has been stopped because the time/resources
             allotted to solve the MP have been depleted; this is not a
	     "serious" condition, in that it can be "solved" by providing
	     the MPSolver with more time/resources (usually, just calling
	     it again, possibly multiple times);

   kError    in case of a failure in the MP solver.

   If the "active set technique" is being used [see SetActvSt() below], the
   problem that is solved is actually a *restriction* of the "real" MP, where
   (possibly many) primal variables d[ i ] are forced to 0; this "restricted"
   MP may be *empty* (so that the solver might have to return kUnfsbl) while
   the full MP is nonempty. Hence, after a kUnfsbl return a "price in" [see
   Readd() below] should be done to find if new variables can to be added to
   the "restricted" MP in order to try to make it nonempty. Thus, in the
   unfeasible case ReadMult() [see below] must return the multipliers that
   make a *feasible ascent extreme ray* in the dual, so that the corresponding
   direction returned by Readd() [see below] is nonzero on variable `i'
   (strictly positive if `i' is constrained) if and only if that variable
   has to be added. */

   virtual MPStatus SolveMP( void ) = 0;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading results
    @{ */

/** Returns the value of the model Fi_{B,Lambda} in d*, the optimal primal
   solution of the MP; this value is often called "v*". This is the value of
   the *translated* model

      v* = Fi_{B,Lambda}( d* ) = Fi_B( Lambda + d* ) - Fi( Lambda ) ,

   hence it is a negative value that indicates the improvement in the
   Fi-value that the model Fi_B() predicts for a step to Lambda + d*.

   When Fi is a decomposable function [see SetDim() above], each of its
   components has a separate model Fi[ k ]_{B,Lambda}; ReadFiBx() must
   return:

   - Fi[ wFi ]_{B,Lambda}( d* ), if 0 <= wFi <= NrFi;

   - \sum_{k = 1}^{NrFi} Fi[ k ]_{B,Lambda}( d* ) if Inf< Index >() > wFi >
     NrFi, i.e., the value of all the models *excluding* the linear one;

   - the value of the aggregated model if wFi == Inf< Index >().

   Note that if there are no subgradients (of a given component) in B, the
   problem that is solved is different from the standard one, as there is
   no such thing as a model for that component. Thus the model is temporarily
   removed from the master problem (like adding a constraint v[ wFi ] >= 0
   for the usual representation where v[ wFi ] is the value of the model).
   In this case, ReadFiBLambda( wFi ) must return Inf< HpNum >().

   The above changes if Fi has "easy" components [see FiOracle::GetBNC() ...],
   because the model of each "easy" component is *not* translated by value.
   Hence, if wFi is the index of an easy component then the method returns

     Fi[ wFi ]( Lambda + d* )

   rather than

     Fi[ wFi ]( Lambda + d* ) - Fi[ wFi ]( Lambda ) .

   Note that in both cases this is the *actual value* of the (corresponding
   component of the) function, because for easy components *the model is the
   actual function* (that is, Fi[ wFi ] == Fi_B[ wFi ]). This information is
   easily extracted out of the solution of the (Dual) Master Problem: if
   x[ wFi ]* is the optimal solution of the Master Problem in terms of the
   "original variables" of the wFi-th component, the required value is just
   ( c[ wFi ] - (Lambda + d* ) * A[ wFi ] ) x[ wFi ]*).

   Note that the 0-th component is an "easy" one but it is dealt with as a
   "difficult" one, in that the translated model (the only possible one,
   the true linear function) is used. */

   virtual HpNum ReadFiBLambda( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns the value of the primal stabilizing term D() in the primal
   optimal solution d* if tt is taken as the proximal parameter. */

   virtual HpNum ReadDt( cHpNum tt = 1 ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns the value of the conjugate Fi_{B,Lambda}* of the model
   Fi_{B,Lambda} in z*, the optimal dual solution of the MP; this value is
   often called "sigma*".

   When Fi is a decomposable function [see SetDim() above], each component
   (of the model) has a separate conjugate function Fi[ k ]_{B,Lambda}* and
   a separate dual solution z[ k ]*: the dual stabilizing term is then

     D_t*( - \sum_{k = 0}^{NrFi} z[ k ] ).

   If primal constraints Lambda \in L are added to the primal, this is viewed
   as adding the indicator function I_L of the feasible set to the stabilized
   Primal Master Problem; this corresponds to having the support function
   \sigma_L appearing in the dual objective function. Thus, in this case the
   dual objective function is

     \sum_{k = 1}^{NrFi} Fi[ k ]_{B,Lambda}*( z[ k ] ) + \sigma_L( w ) +
                         D_t*( - \sum_{k = 0}^{NrFi} z[ k ] )

   i.e., it has a term for each component, the stabilizing term and an
   extra term for constraints (this term is always zero in the primal
   because the indicator function I_L is zero in L).

   Thus, ReadSigma() must return:

   - Fi[ wFi ]_{B,Lambda}*( z[ i ]* ) if 1 <= wFi <= NrFi: note that the
     value for the 0-th component of Fi is always 0, since the conjugate of
     a linear function is 0 iif its argument is equal to the gradient and
     +Infinity otherwise---in other words, z[ 0 ] is always equal to the
     gradient of the 0-th component;

   - \sum_{k = 1}^{NrFi} Fi[ k ]_{B,Lambda}*( z[ k ] ) + \sigma_L( w ) if
     wFi > NrFi, i.e., the dual objective function except the (dual)
     stabilizing term.

  \note For "easy" components of Fi() [see FiOracle::GetBNC() ...] the value
        that should be returned by this method is the *actual* value of the
	conjugate of the (corresponding component of) Fi(). Although it may
	be possible to extract that information out of the solution of the
	MP, a Bundle algorithm typically will never use it, so the MPSolver
	is allowed to assume that ReadSigma( wFi ) will never be called with
	wFi the index of an "easy" component. Furthermore, as in the case of
	ReadFiBLambda() [see above], in an "aggregated" call (that is, with
	wFi > NrFi) *only the difficult components are counted*, that is, the
	value of the conjugate of the model of an "easy" component *never*
	has to be computed. */

   virtual HpNum ReadSigma( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns the value of the dual stabilizing term D*() in the dual optimal
   solution z* if tt is taken as the proximal parameter. */

   virtual HpNum ReadDStart( cHpNum tt = 1 ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns a read-only pointer d to a vector containing the primal optimal
   solution d*; the entry of d* relative to variable `i', if required (see
   below), can be found in d[ i ] for i in [ 0 .. CrrSGLen - 1 ] (see
   [Add/Remove/Subst]Var() below for a discussion on variable "names").

   If Fulld == true, the solution d[] returned by Readd() is a "full" one,
   even if an "active set" of variables has been defined [see SetActvSt()
   below]; that is, even the entries of d* relative to variables that are not
   "active" are computed in response to a call to Readd( true ). This may be
   costly, and it is not always necessary: thus, a call to Readd( false )
   computes only the entries relative to "active" variables, the contents
   of d[ i ] for an "inactive" variable `i' being unpredictable.

   Note that the returned pointer may point to a temporary data structure
   rather than to a "persistent" one: the pointer is no more valid (it may
   even point to no-more-allocated memory) after any call to any other method
   of the class. */

   virtual cLMRow Readd( bool Fulld = false ) = 0;

/*--------------------------------------------------------------------------*/
/** Writes in tz[] the dual optimal solution z* in its "natural" format,
   i.e., a single vector z*, with as many components as d*. The format is
   "dense" if I == 0 and "sparse" otherwise: if I != 0, then only the
   first D entries of tz[] contain relevant information, the i-th entry being
   relative to the variable `I[ i ]', all other entries being zero. If
   I == 0, a full CrrSGLen-vector is returned. Note that all the entries
   of z are "significant", even the ones corresponding to variables which
   are not in the "active set" [see SetActvSt() below].

   If Fi is separable, then there are as many dual optimal solutions z[ i ]*
   as components [see ReadSigma() above], the optimal solution z being their
   sum. Thus, ReadZ() must return:

   - z[ 0 ]* == the gradient of the 0-th component if wFi == 0;

   - z[ wFi ]* if 1 <= wFi <= NrFi;

   - the sum of all the things that would be reported by calls with wFi going
     from 1 to NrFi if NrFi < wFi < Inf< Index >();

   - the sum of all the things that would be reported by calls with wFi going
     from 0 to NrFi if wFi == INF (this is the same as summing what is
     separately obtained with wFi == 0 and wFi == NrFi + 1).

  \note In the constrained case, "extra" dual variables w (w[ i ])
        corresponds to constraints, and the optimal solution is actually
	z* + w* (z[ i ]* + w[ i ]*); ReadZ() only returns the "z* part" of
	the dual solution, leaving the "w* part" out.

  \note The vector pointed by I[] may be that of a temporary data structure;
        hence, it is no more valid (it may even be dangling) after any call to
	any other method of the class. There is one notable exception: I[] can
	be passed to SetItemBase(), and be used as set of indices af a sparse
	subgradient to be inserted in the Master Problem. In fact, this is the
	typical reason why a z[ i ]* is computed: to perform "aggregation".

  \note For "easy" components of Fi() [see FiOracle::GetBNC() ...] the vector
        that should be returned by this method is a subgradient of Fi[ wFi ]
	in Lambda. Although it is possible to extract that information out of
	the solution of the Master Problem (if x[ wFi ]* is the optimal
	solution of the Master Problem in terms of the "original variables" of
	the wFi-th component, this is just - A[ wFi ] x[ wFi ]*), a Bundle
	algorithm typically will never use it except in its "aggregated" form
	(that is, in a call with wFi > NrFi), so the MPSolver is allowed to
	assume that ReadZ( wFi ) will never be called with wFi the index of an
	"easy" component. */

   virtual void ReadZ( LMRow tz , cIndex_Set &I , Index &D ,
		       cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** The dual optimal solution z* lies in the domain of Fi_{B,Lambda}*; at
   least with the "classical" Cutting Plane model \hat{Fi}_{B,Lambda}, this is
   the convex hull of all the subgradients in B. Hence, there is another
   possible forms in which the dual optimal solution can be represented, i.e.,
   as a set of (convex, for \hat{Fi}) multipliers Theta[] (>= 0 and summing
   to 1) such that

      z* = Sum{ i \in B } Theta[ i ] * z[ i ]

   where z[ i ] is the item with "name" `i' [see [Set/Get]Item() below].
   ReadMult() returns these multipliers. The format is "dense" if I == 0
   and "sparse" otherwise: if I != 0, then only the first D entries of the
   returned vector contain relevant information, the i-th entry being relative
   to item with name `I[ i ]', all other entries being zero. If I == 0, a
   full MaxBSize-vector is returned.

   If Fi is separable, then there are as many dual optimal solutions z[ i ]*
   as components [see ReadSigma() above], the optimal solution z being their
   sum; thus, ReadMult() must return:

   - the multipliers making z[ wFi ]* (those relative to the subgradients
     of the i-th component) if 1 <= wFi <= NrFi (wFi = 0 makes no sense since
     z[ 0 ]* is a known constant vector, i.e., its multiplier is always 1);

   - the union of all the things that would be reported by calls with wFi
     going from 1 to NrFi otherwise.

   If the parameter IncldCnst is left to true, then the returned multipliers
   are both of subgradients and constraints (corresponding to component wFi).
   In the constrained case, "extra" dual variables w (w[ i ]) corresponds to
   constraints, and the optimal solution is actually z* + w*
   (z[ i ]* + w[ i ]*): the multipliers corresponding to constraints are just
   nonnegative---they are not constrained to sum to 1---as they represent
   "extreme rays" in the dual feasible set. If IncldCnst == false, only the
   multipliers corresponding to subgradients are returned.

   \note The pointers returned by ReadMult() may point to a temporary data
         structure; the pointers are no more valid (they may even be dangling)
	 after any call to any other method of the class.

   \note No subgradient of "easy" components of Fi() [see FiOracle::GetBNC()
         ...] are ever produced, so a call to ReadMult( wFi ) with wFi the
	 index of an "easy" component may appear to make no sense. However,
	 something else naturally "takes the place" of the multipliers: the
	 optimal solution x[ wFi ]* of the Master Problem in terms of the
	 "original variables" of the wFi-th component. Thus, calls to
	 ReadMult( wFi ) with an "easy" wFi are indeed allowed,	 with the
	 following different meaning: the vector returned by the method (that
	 may be "dense" or "sparse" as in the standard case) represents
	 x[ wFi ]*, and it is therefore a FiOracle::GetBNC( wFi )-vector. This
	 information, however, can *only* be accessed when calling
	 ReadMult( wFi ) for wFi <= NrFi: in "global" calls (wFi > NrFi) only
	 the multipliers corresponding to non-easy components are returned. */

   virtual cHpRow ReadMult( cIndex_Set &I , Index &D ,
			    cIndex wFi = Inf< Index >() ,
			    const bool IncldCnst = true ) = 0;

/*--------------------------------------------------------------------------*/
/** Must return the information required by NDOSolver::ReadLBMult(), see
   NDOSlver.h. */

   virtual HpNum ReadLBMult( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Must return the scalar product between the item with name `Nm' [see
   SetItem() below] and the primal optimal solution d*; this information is
   typically computed anyway, so it may be made available at little cost.
   When Nm >= MaxBSize, then the scalar product to be returned is that with
   the constant subgradient of the "linear part" of Fi, the 0-th component. */

   virtual HpNum ReadGid( cIndex Nm = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Writes in (the array pointed by) Lmbd1 the vector
   Lmbd + ( Tau / t ) * d*, where d* is the optimal primal solution of the MP.

   The "format" of Lmbd and Lmbd1 depends on whether or not the "active set"
   is being used [see SetActvSt() below]; if not, then Lmbd is (and Lmbd1 is
   expected to be) just a "dense" CrrSGLen-vectors. Otherwise, both Lmbd and
   Labd1 are intended to be "restricted" to the set of "active" variables
   that are defined by SetActvSt(); that is, only the first k entries of
   Lmbd/Lmbd1 are significant, where k is the size of the active set as
   resulting from all the previous calls to *ActvSt() [see below], and
   Lmbd/Lmbd1[ i ] contain the value of the i-th "active" variable (all the
   other variables being fixed to zero).

   Note that, if an "active set" is defined, this method can compute only the
   entries of d* that are strictly necessary for the purpose (like Readd()
   does, see above), thereby possibly saving some time. */

   virtual void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau ) = 0;

/*--------------------------------------------------------------------------*/
/**< Gives a linear lower approximation on the optimal value of v (the
   predicted decrease from the Cutting Plane model) as a function of the
   parameter t in the subproblem. That is, returns two numbers lp and cp such
   that

     v( t ) = Fi_{B,x}( x + t * d ) >= t * lp + cp. */

   virtual void SensitAnals( HpNum &lp , HpNum &cp ) = 0;

/*--------------------------------------------------------------------------*/
/** If this method is called within any of the methods of the class that are
    "actively timed" (this depends on the subclasses), it returns the user
    and sistem time (in seconds) since the start of that method. If methods
    that are actively timed call other methods that are actively timed, these
    methods return the (...) time since the beginning of the *outer* actively
    timed method. If these methods are called outside of any actively timed
    method, they return the (...) time spent in all the previous executions
    of all the actively timed methods of the class.

    Implementing the proper calls to MPt->Start() and MPt->Stop() is due to
    derived classes; these should at least be placed at the beginning and at
    the end, respectively, of SolveMP() that is, at least SolveMP() should be
    "actively timed". */

   void MPTime( double &t_us , double &t_ss ) {
    t_us = t_ss = 0;
    if( MPt ) MPt->Read( t_us , t_ss ); 
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** If this method is called within any of the methods of the class that are
    "actively timed" (this depends on the subclasses), it returns the total
    time (in seconds) since the start of that method. If methods that are
    actively timed call other methods that are actively timed, these methods
    return the (...) time since the beginning of the *outer* actively timed
    method. If these methods are called outside of any actively timed method,
    they return the (...) time spent in all the previous executions of all
    the actively timed methods of the class.

    Implementing the proper calls to MPt->Start() and MPt->Stop() is due to
    derived classes; these should at least be placed at the beginning and at
    the end, respectively, of SolveMP() that is, at least SolveMP() should be
    "actively timed". */

   double MPTime( void ) { return( MPt ? MPt->Read() : 0 );  }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   /// returns the maximum bundle size

   Index GetMaxBSize( void ) { return( MaxBSize ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// returns the maximum subgradient length

   Index GetMaxSGLen( void ) { return( MaxSGLen ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// returns the current subgradient length

   Index GetCrrSGLen( void ) { return( CrrSGLen ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /// returns the number of components of Fi()

   Index GetNrFi( void ) { return( NrFi ); }

/*--------------------------------------------------------------------------*/
  /// returns the current number of items (of either type) in the bundle.

   virtual Index BSize( cIndex wFi = Inf< Index >() ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /// returns the current number of constraints in the bundle.

   virtual Index BCSize( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns 1 + the maximum "name" of an item in the bundle if the bundle is
   nonempty, and 0 otherwise. wFi "restricts" the count to the items
   corresponding to that particular component, with the usual convention. */

   virtual Index MaxName( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns the name of the component that the item `i' is relative to; a
   return value of Inf< Index >() means that no item with name `i' is in the
   bundle. Note that a return value of 0 is possible only if the item is a
   constraint, since the 0-th component may have constraints but no
   subgradients: its only (sub)gradient is dealt with in a special way and it
   has no name. */

   virtual Index WComponent( cIndex i ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns true if item `i' is in the Bundle and it is a subgradient, and
   false otherwise; hence, an item `i' for which wFi == WComponent( i ) > 0
   [see above] and IsSubG( i ) == false is a constraint relative to the wFi-th
   component of Fi (again, wFi here can be zero). */

   virtual bool IsSubG( cIndex i ) = 0;

/*--------------------------------------------------------------------------*/
/** Return the number of variables that are constrained in sign; this must of
    course be <= CrrSGLen. The standard implementation is fine for MP Solvers
    not accepting constrained problems. */

   virtual Index NumNNVars( void ) { return( 0 ); }

/*--------------------------------------------------------------------------*/
/** Return the number of variables (<= CrrSGLen) that have either a lower
    bound constraint or an upper bound constraint, or both. Hence,
    NumBxdVars() >= NumNNVars(), and the difference gives the number of
    variables with u.b. but not l.b.. The standard implementation is fine for
    MP Solvers not accepting constrained problems. */

   virtual Index NumBxdVars( void ) { return( 0 ); }

/*--------------------------------------------------------------------------*/
/** This method shall return true if the variable `i' is constrained to be
    NonNegative, and false otherwise; the standard implementation is fine
    for MP Solvers not accepting constrained problems. */

   virtual bool IsNN( cIndex i ) { return( false ); }

/*--------------------------------------------------------------------------*/
/** Must return a read-only pointer to a vector containing in position i the
   current value of the linearization error/right hand side of item `i'. */

   virtual cHpRow ReadLinErr( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Reads back the current value of the Lower Bound; see SetLowerBound() for
 * the meaning of wFi. Note that the value is relative to the function values
 * in the current point, and therefore changes when these (say, the current
 * point itself) do. The updating need be done by the MPSolver, which is why
 * it makes sense to be able to read the current value. */

   virtual HpNum ReadLowerBound( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/* Writes in Gi[] the description of the item with "name" Nm. It is an error
   if there is no item with that name in the MP. The format of the returned
   vector depends on what is returned in GB. If GB == 0, then the vector in
   Gi[] is "dense". Otherwise, the vector in Gi[] is "sparse", GB[] contains
   the read-only vector of indices of the nonzero elements (ordered in
   increasing sense and Inf< Index >()-terminated), and the number of these is
   returned by the method. That is, if GB != 0 and NNZGi is the value returned
   by the method, then the "true" Gi is

     Gi[ j ] = / Gi[ i ]   if GB[ i ] == j for some 0 <= i < NNZGi
               \   0       otherwise. */

//   virtual Index ReadGi( cIndex Nm , SgRow Gi , cIndex_Set &GB ) = 0;

/*--------------------------------------------------------------------------*/
/** Must return the threshold used by the solver to decide when a variable is
   zero; Lambda[ i ] <= EpsilonD() is considered as == 0. */

   virtual HpNum EpsilonD( void ) = 0;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

/** GetItem( wFi ) must return a pointer SG to a vector of SgNum's where a
   new item has to be written. wFi tells to which of the components of Fi()
   [see SetDim() above] the new item belongs, with the usual format:

   - when wFi == 0, the new item is relative to the linear 0-th component of
     Fi, i.e., it is either its constant (sub)gradient, or a constraint which
     defines its domain Feas0;

   - when 1 <= wFi <= NrFi, the new item is relative to the wFi-th component
     of Fi;

   - when Inf< Index >() > wFi > NrFi, the new item is relative to all the
     components *except* the 0-th;

   - when wFi == Inf< Index >(), the new item is relative to all the components,
     *comprised* the 0-th.

   This information allows the MPSolver to properly dimension the vector to
   be returned [see GetMaxNZ() in FiOracle.h]. */

   virtual SgRow GetItem( cIndex wFi = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** When the correct values have been written in the vector SG returned by
   GetItem() [see above], SetItemBse() can be used to tell the MPSolver
   the "format" of the information. If SGBse == 0, then SG is in "dense"
   format, i.e., SG[ i ] contains the entry of the subgradient relative to
   the variable with name `i', 0 <= i < CrrSGLen - 1. If SGBse != 0 instead,
   then SG is in "sparse" format, i.e., SG[ i ] is the entry of the
   subgradient relative  to the variable with name `SGBse[ i ]', 0 <= i <
   SGBDm - 1. SGBse is assumed to be ordered in increasing sense and
   Inf< Index >()-terminated, i.e., SGBse[ SGBDm ] == Inf< Index >(); if
   SGBse == 0, then the value of SGBDm is ignored. If SetItemBse() is *not*
   called, then the "dense" format is assumed.

   Note that all the CrrSGLen entries of the item must be given (although
   in the sparse format those that are 0 can be given only implicitly);
   that is, the entry relative to the variable `i' must be given even if
   `i' is currently out of the "active set" [see SetActvSt() below].

   The MPSolver is allowed *not to copy* the memory pointed by SGBse (if any)
   and to keep using it throughout all the sequence of Check[SubG/Const](),
   ChangesMPSol() and SetItem() [see below] that is needed to decide the
   "fate" of the item. As soon as SetItem() returns, or GetItem() [see above]
   is called again, the memory pointed by SGBse must not be accessed again. */


   virtual void SetItemBse( cIndex_Set SGBse = 0 , cIndex SGBDm = 0 ) = 0;

/*--------------------------------------------------------------------------*/
/** After that SetItemBse() has been called, CheckSubG() (if the item is
   a subgradient, otherwise see CheckCnst() below] must be called to:
 
   - verify that the new item is not already in the bundle, and

   - compute some information about the item, i.e., the linearization
     error of the item and the scalar product with d*.

   The return value must be the "name" of an item in the bundle that is
   identical (element-wise) to GS, if such an item exists, or Inf< Index >()
   otherwise. This allows to check for useless copies (but this check is
   costly and can be deactivated, see CheckIdentical() above), which is not
   an unlikely possibility e.g. when Fi is the Lagrangian function of a
   combinatorial optimization problem (and therefore the set of possible
   subgradients is finite). Note that subgradients that are element-wise
   identical but *belong to different components* are *not* identical and
   should *not* be substituted for one another; in other words, the check
   *only has to be performed with subgradients of the same component as SG*
   (which, besides being correct, is also likely to be less costly). Also,
   *subgradients should not be checked against constraints*, as they are
   two entirely different kinds of items.

   Also, the scalar product between the primal optimal solution d* of the MP
   and SG must be calculated and returned in ScPri.

   DFi must contain the valueof Fi[ wFi ]( Lambda1 ) - Fi[ wFi ]( Lambda ),
   where Lambda1 = Lambda + ( Tau / t ) * d*, and must Ai contain its
   linearization error w.r.t. *Lambda1*, i.e., SG is an Ai-subgradient of
   Fi[ wFi ]() in Lambda1. Ai must be updated to the linearization error of
   SG w.r.t. *the current point Lambda*, that is

     Ai = Ai - ( DFi - ( Tau / t ) * SG * d*  ).

   A special case, that has a separate treatment, is that of Tau <= 0. This is
   taken to mean that SG is an Ai-subgradient *in Lambda* (typically, DFi
   then is also 0). This happens when adding subgradients for whom the
   linearization error in Lambda is already known, an important case being
   that of aggregation, i.e., adding to the bundle the  dual optimal solution
   z* as provided by ReadZ() and ReadSigma(). In such a case, Ai is *not*
   modified, and the original one is directly used.

   A further specific difference is whether CheckSubG() with Tau <= 0 is
   called inside the normal loop "compute function value and subgradient at
   next iterate, put them in the master problem, re-solve" or outside it;
   in other words, if it is possible that ChangeCurrPoint() is called
   before the next master problem solution. If not, then the computation of
   the scalar product (that anyway is not used to update Ai)  can be
   avoided, saving some cost: this is done if Tau == 0. If, instead,
   ChangeCurrPoint() may be called, then the scalar product need be computed
   anyway, even if it is not used to update Ai: this is done if Tau < 0.

   The Ai computed by CheckSubG() are typically used in the MP. Note that
   CheckSubG() *must* be called before SetItem() for all items *except* the
   constant subgradient of the "linear part" of Fi (wFi == 0 in GetItem()),
   where CheckSubG() *must not* be called. */

   virtual Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			    HpNum &ScPri ) = 0;

/*--------------------------------------------------------------------------*/
/** After that SetItemBse() has been called, CheckCnst() [if the item is a
   constraint, otherwise see CheckSubG() above] must be called to:
 
   - verify that the new item is not already in the bundle, and

   - compute some information about the item, i.e., the Right Hand Side of
     the item and the scalar product with d*.

   The return value must be the "name" of a constraint in the bundle that is
   identical (element-wise) to GS, if such an item exists, or Inf< Index >()
   otherwise. This allows to check for useless copies (but this check is
   costly and can be deactivated, see CheckIdentical() above). Note, however,
   that *the component to which the constraint belongs is now irrelevant*,
   and therefore that the checks *have to be performed across all the
   components* (comprised the 0-th). The reason is that while a constraint
   nominally belongs to one component (in Lagrangian terms, the corresponding
   primal unbounded ray is for one particular subproblem) this has no impact
   on the master problem: unlike subgradients, two identical constraints
   (columns) corresponding to two different components actually are the same
   constraint (column) in the master problem, and therefore one is redundant.
   Thus one can (and therefore should) safely declare a constraint a copy of
   another even if they belong to different components. Note that this
   corresponds, from the primal side, to choose between two different primal
   rays belonging to different subproblems, and therefore it does have an
   impact on the generated primal (convexified) solution. However, if one
   would keep the two constraints then it would be the Master Problem solver
   which would arbitrarily "split the corresponding primal variables among
   them", so a decision would still be somewhat taken in this respect that is
   entirely out of control; better, then, to at least save some running time
   and memory. If this turns out to be a problem for some weird reason, the
   check can be deactivated with CheckIdentical() and the problem solved.

   As for CheckSubG(), the scalar product between the primal optimal solution
   d* of the MP and SG must be calculated and returned in ScPri.

   Ai must contain the Right Hand Side, i.e., the constraint on the original
   variables is

     SG * Lambda <= Ai.

   However, the MP is defined over the d variables, that is the constraint
   in the MP must be

     SG * ( CrrPnt + d ) <= Ai

   where CrrPnt is the current point of the Bundle algorithm. The format
   of CrrPnt depends on whether or not the "active set" strategy has been
   activated with SetDim() [see above]. If the active set is used, then
   CrrPnt[] is in "sparse" format, i.e., `CrrPnt[ i ]' is the value for the
   variable with name `AVrs[ i ]' for i = 0 ..  AVDm - 1, all the other
   values being zero. Otherwise, CrrPnt[] is in "dense" format. Upon return,
   Ai must be updated to the "scaled" Right Hand Side w.r.t. the current
   point, that is

     Ai = Ai - SG * CurrPnt.

   The Ai computed by CheckCnst() are typically used in the MP. Note that
   CheckCnst() *must* be called before SetItem() for all constraints,
   *comprised* those of the 0-th component (wFi == 0 in GetItem()). */

   virtual Index CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt ) = 0;

/*--------------------------------------------------------------------------*/
/** After that Check[SubG/Cnst]() [see above] has been called, this method
   can be used to check whether or not the introduction of the new item
   changes the solution of the MP. If not, the Bundle algorithm may decide
   not to insert it in the bundle, and however knows if the newly obtained
   information is enough for obtaining a different d* at the next iteration.
   This allows to cope with "non-exact functions", i.e., with a FiOracle that
   is not capable of - or is instructed not to - computing the Fi-values
   and/or the subgradients exactly. Note that ChangesMPSol() *must not* be
   called for the constant subgradient of the linear 0-th component of Fi. */

   virtual bool ChangesMPSol( void ) = 0;

/*--------------------------------------------------------------------------*/
/** If -- after all the above checks -- the Bundle code decides to actually
   insert the item in the bundle, it has just to call SetItem( Nm ). `Nm' is
   the "name" of the new item: the names of items are the integers between 0
   and MaxBSize - 1. `Nm' must be a currently unused name, that is, either a
   name that has not yet been used or a name that has been used in RmvItem()
   [see below] and not yet used in GetItem() ever since. However, if wFi == 0
   has been passed to GetItem(), then Nm == Inf< Index >() has to be passed, as
   the only subgradient corresponding to the (linear) 0-th component of Fi has
   no name. */

   virtual void SetItem( cIndex Nm = Inf< Index >() ) = 0;

/*--------------------------------------------------------------------------*/
/** Although having multiple copies of the same item in the bundle is legal,
   it is also useless and it should be avoided. If CheckItem() returned a
   "valid" name (< Inf< Index >()), then SubstItem( Nm ) can be used to replace
   the current item in the bundle with name `Nm' with SG. It is assumed that
   `Nm' is the name of an item element-wise identical to SG (but possibly with
   different Ai), typically the name returned by CheckItem(). Typically, this
   operation makes sense when the Ai (Right Hand Side/linearization error) of
   the new item is *smaller* than that of the existing copy in the bundle.
   Thus, SubstItem( Inf< Index >() ) - that would operate on the (only)
   (sub)gradient of the 0-th component of Fi - *cannot* be called. */

   virtual void SubstItem( cIndex Nm ) = 0;

/*--------------------------------------------------------------------------*/
/** Removes the item with name `i' from the bundle. Items removed are unused,
   and their names can be used again in GetItem() [see above]. Initially, all
   the items are unused. */

   virtual void RmvItem( cIndex i ) = 0;

/*--------------------------------------------------------------------------*/
/** Removes *all* items from the bundle. Note that RmvItems() does *not*
   remove the constant subgradient of the linear 0-th component of Fi
   (RmvItem() [see above] also cannot be used for the purpose as that
   subgradient has no name). In principle, any Fi has a linear part, possibly
   all-0; thus, to remove that item one has just to use the
   [Get/Check/Set]Item() sequence passing an all-0 vector. */

   virtual void RmvItems( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Solving the MP can be costly; a possible way of decreasing the cost is to
   resort to an "active set" strategy over the primal variables d[ i ]. That
   is, one may guess that a (large) set of the variables will have 0 as the
   optimal value, and only let the others (the "active set") to vary; if it
   turns out that the guess has been incorrect, the active set must be
   revised by inserting or deleting variables [see [Add/Rmv]ActvSt() below].
   This way, a smaller MP is solved at each step.

   Whether or not an "active set" is used is decided in SetDim() [see above];
   if it is, then SetActvSt() sets it to the variables whose names are in the
   first AVDm entries of the vector AVrs[] (which must be ordered in
   increasing sense and Inf< Index >()-terminated, i.e. AVrs[ AVDm ] ==
   Inf< Index >()). It is an error to call SetActvSt() if the active set option
   has not previously been selected with SetDim().

   Note that changing the active set has *no impact* on the linearization
   errors of the items. Setting the active set dictates which among the
   primal variables d[ i ] of the MP can be different from zero, hence,
   which primal variables Lambda[ i ] of the original NDO problem may change
   their value w.r.t. the value that they have in the current point. The
   linearization errors of the items depend on the current point only,
   and therefore they change if and only if the current point changes [see
   ChangeCurrPoint() and [Add/Rmv]Vars() below]. Since setting the active
   set affects only the variables d[ i ] of the MP and does not affect the
   value of the primal variables Lambda[ i ] of the Bundle algorithm, setting
   (or changing) the active set has no impact on the linearization errors.

   Note that **the pointer passed to SetActvSt() will be retained and used**
   until it is changed, either by a new call to SetActvSt() - which just
   substitutes any previous active set with the new one - or by calls to
   [Add/Rmv]ActvSt() [see below]. Thus, the caller should not change it while
   it is "in use" by the MPSolver object; the vector can only be changed
   *after* the return of the call to SetActvSt() or [Add/Rmv]ActvSt() which
   makes it unused. This implies that any call to SetActvSt() which sets an
   entirely new active set must pass a pointer that is different from the one
   which contained the old active set [see also [Add/Rmv]ActvSt() below].
   Exception is of course AVrs == 0, that sets an *empty* active set: that
   is, the active set option is used (only SetDim() can change this), and it
   currently contains no variables. This is the default if SetActvSt() has
   not been called yet.

   Important note: after a call to SetActvSt(), all the solution information
   corresponding to the last call to SolveMP() is *lost*, and the results
   returned by all the query methods are unpredictable. This is true even for
   CheckSubG(), CheckCnst() and ChangesMPSol(), so adding items to the bundle
   after a call to this method is not advised. */

   virtual void SetActvSt( cIndex_Set AVrs = 0 , cIndex AVDm = 0 ) = 0;

/*--------------------------------------------------------------------------*/
/** Adds `AdDm' new variables to the current active set, whose names are
   contained in the first AdDm positions of the vector Addd (ordered in
   increasing sense and Inf< Index >()-terminated); those names must *not*
   already belong to the active set.

   AVrs must point to a vector containing the *new* active set; a pointer to
   that vector is kept and used from now on as the new active set, replacing
   the old vector set with the latest call to one of SetActvSt(), AddActvSt()
   or RmvActvSt(). Thus, the vector pointed by the "old" pointer can be
   freely modified after the return of AddActvSt(), while the one pointed by
   AVrs must not be modified until it is in turn substituted by a new one.

   Thus, AVrs must usually be different from the "old" pointer; one exception
   is permitted, however, in the case when all the variables to be added are
   to be found *after* all the variables that were present previously. In
   this case, the caller is authorized to use the old pointer as the new
   pointer (and checking this is exactly how the MPSolver detects the fact),
   overwriting the "Inf< Index >()" which previously terminated the active set
   with the names of the variables to be added *before* the call to
   AddActvSt() (i.e., Addd = AVrs +  <n. of variables in the old active set>).

   This method can be used only if an active set has been previously defined
   with SetActvSt( <non-0> ) and not yet removed with SetActvSt( 0 ).

   Important note: after a call to AddActvSt(), all the solution information
   corresponding to the last call to SolveMP() is *lost*, and the results
   returned by all the query methods are unpredictable. This is true even for
   CheckSubG(), CheckCnst() and ChangesMPSol(), so adding items to the bundle
   after a call to this method is not advised. */

   virtual void AddActvSt( cIndex_Set Addd , cIndex AdDm ,
                           cIndex_Set AVrs ) = 0;

/*--------------------------------------------------------------------------*/
/** Removes `RmDm' variables from the current active set, whose names are
   contained in the first RmDm positions of the vector Rmvd (ordered in
   increasing sense and Inf< Index >()-terminated); those names *must* be in
   the current active set.

   AVrs must point to a vector containing the *new* active set; a pointer to
   that vector is kept and used from now on as the new active set, replacing
   the old vector set with the latest call to one of SetActvSt(), AddActvSt()
   or RmvActvSt(). Thus, the vector pointed by the "old" pointer can be
   freely modified after the return of RmvActvSt(), while the one pointed by
   AVrs must not be modified until it is in turn substituted by a new one.

   This method can be used only if an active set has been previously defined
   with SetActvSt( <non-0> ) and not yet removed with SetActvSt( 0 ).

   Important note: after a call to RmvActvSt(), all the solution information
   corresponding to the last call to SolveMP() is *lost*, and the results
   returned by all the query methods are unpredictable. This is true even for
   CheckSubG(), CheckCnst() and ChangesMPSol(), so adding items to the bundle
   after a call to this method is not advised. */

   virtual void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm ,
                           cIndex_Set AVrs ) = 0;

/*--------------------------------------------------------------------------*/
/** This method adds `NNwVrs' new variables to the MP, in response to a call
   to NDOSolver::AddVariables( NNwVrs , IVs ) [see NDOSlvr.h]. The derived
   classes are required to use GetUC(), GetUB(), GetGi() and (possibly)
   GetADesc() [see FiOracle.h] to acquire directly from the FiOracle all the
   information about the newly created variables they need.

   The method is used when NDOSolver::AddVariables() is called, i.e., when
   new variables Lambda[ i ] of the original NDO problem are created. Note
   that adding variables Lambda[ i ] changes, in principle, the current point,
   and therefore the linearization errors of the items. For instance, in the 
   Lagrangian case, it is easy to check that the linearization error
   alfa_{Lambda}[ i ] of the subgradient named `i' with respect to the
   current point Lambda is simply

    alfa_{Lambda}[ i ] = Fi( Lambda )
                       - c( x[ i ] ) + Lambda * ( b - A( x[ i ] ) )

   where x[ i ] is the solution of the Lagrangian problem corresponding to
   the subgradient. Thus, adding a new variable Lambda_j changes the
   linearization error of Lambda_j * ( b_j - A_j( x[ i ] ) ). Hence, the
   linearization errors of all items change *except if Lambda_j == 0*. In
   order to avoid the recomputation of the linearization errors, it is
   assumed that all the variables d[ i ] of the Primal Master Problem added
   with AddVars() correspond to variables Lambda[ i ] that are *set to 0*;
   in other words, it is assumed that adding the variables does not change
   the linearization errors (changing them is always possible with
   ChangeCurrPoint()).

   Note that, if an "active set" has been defined [see SetDim() above], then
   the newly added variables are *not* added to the current active set; this
   *must* be done explicitly, if needed, with [Set/Add]ActvSt() [see above].
   If the active set option has not been selected, then the variables are
   automatically added to the MP.

   Important note: after a call to AddVars(), all the solution information
   corresponding to the last call to SolveMP() is *lost*, and the results
   returned by all the query methods are unpredictable. This is true even for
   CheckSubG(), CheckCnst() and ChangesMPSol(), so adding items to the bundle
   after a call to this method is not advised. */

   virtual void AddVars( cIndex NNwVrs ) = 0;

/*--------------------------------------------------------------------------*/
/** This method removes a set of variables from the MP, in response to a call
   to NDOSolver::RemoveVariables( which ) [see NDOSlvr.h]; the remaining
   variables in the MS have to be *renamed* as described in the comments to
   NDOSolver::RemoveVariables(). As for NDOSolver::RemoveVariables(),
   whch == 0 means "remove all variables" (hwmny is ignored in this case).

   As opposed to RmvActvSt() [see above], the effect of this removal is
   permanent; the variables are eliminated from the whole MP, rather than
   just "temporarily set aside" outside of the active set. Indeed, this
   "permanent" removal can only happen *after* the removal from the active
   set: no variable in the active set can be removed, i.e., RmvVars( 0 )
   can only be called if the active set is *empty* (if, of course, the
   active set has been at all defined).

   As for AddVars() [see above], it is assumed that the removal of the
   variables does not change the linearization errors of the items, i.e.,
   that all the variables d[ i ] of the Primal Master Problem that are
   removed correspond to variables Lambda[ i ] of the NDO problem that
   have value zero in the current point; it is responsibility of the
   caller to ensure this, using ChangeCurrPoint() if required.

   Important note: after a call to RmvVars(), all the solution information
   corresponding to the last call to SolveMP() is *lost*, and the results
   returned by all the query methods are unpredictable. This is true even for
   CheckSubG(), CheckCnst() and ChangesMPSol(), so adding items to the bundle
   after a call to this method is not advised. */

   virtual void RmvVars( cIndex_Set whch = 0 , Index hwmny = 0 ) = 0;

/*--------------------------------------------------------------------------*/
/** Changes the vector of the linearization errors of the subgradients.

   DeltaAlfa[] is a (NrFi + 1)-vector; the linearization error of each
   subgradient corresponding to a component wFi is *increased* by
   DeltaAlfa[ wFi ]. DeltaAlfa[ 0 ] is ignored, and the right hand sides of
   the constraints are unchanged.

   This is the typical correction that ought to be done when "nonexact"
   functions have to be handled. */

   virtual void ChgAlfa( cHpRow DeltaAlfa ) = 0;

/*--------------------------------------------------------------------------*/
/** Changes the vector of the linearization errors of the subgradients/right
   hand sides of the constraints.

   For every item i in the bundle, the new value for its linearization error/
   right hand side is taken from NewAlfa[ i ]. 0 <= wFi <= NrFi means that
   only the items corresponding to the wFi-th component of Fi are affected by
   the change. In particular, if wFi == 0 then what changes is the constant
   term `b0' in the linear 0-th component of Fi; the new value for b0 is to
   be found in NewAlfa[ MaxBSize ]. If Inf< Index >() > wFi > NrFi, then all the
   components except the 0-th have changed, while if wFi == Inf< Index >() then
   all the components, comprised the 0-th, have changed. */

   virtual void ChgAlfa( cHpRow NewAlfa , cIndex wFi ) = 0;

/*--------------------------------------------------------------------------*/
/** Changes the linearization error of the item with name i in the bundle to
   the new value Ai. */

   virtual void ChgAlfa( cIndex i , cHpNum Ai ) = 0;

/*--------------------------------------------------------------------------*/
/** Updates the data of the Master Problem in response to a change of the
   "current point" Lambda:

    DLambda = newLambda - Lambda,

   i.e., DLambda is the direction that is followed (with step 1) to go from
   the old current point Lambda to the new current point newLambda. DLambda
   is in "dense" format [see MakeLambda1() above], i.e., `DLambda[ i ]'
   contains the displacement relative to variable `i' for
   i = 0 .. CrrSGLen - 1. Note that the set of variables with nonzero
   displacement may *not* be a subset of the "active set" [see
   [Set/Add/Rmv]ActvSt() above], i.e., that the entry of the current point
   for variable `i' may vary even if i is currently forced to be zero by
   being out of the active set, In order for the Master problem to become
   "aware" of the change for these variables, the active set has to be
   changed.

   DFi[ k ] must be the change in the value of the function between the two
   points, i.e.,

    DFi[ k ] = Fi[ k ]( newLambda ) - Fi[ k ]( Lambda )

   for 1 <= k <= NrFi, while

    DFi[ 0 ] = Fi( newLambda ) - Fi( Lambda )

   (i.e., the value of the whole function). Note that the values for "easy"
   components are provided as well, altough the MPSolver in principle knows
   them already.

   The formulae for updating the data are the following. If z[ i ] is a
   alfa[ i ]-subgradient to the k-th component of Fi in Lambda, then it is
   a newalfa[ i ]-subgradient (to the k-th component of Fi) in newLambda,
   where

     newalfa[ i ] = alfa[ i ] + DFi[ k ] - DLambda * z[ i ] .

   If z[ i ] is a constraint

     z[ i ] * ( Lambda + d ) <= RHSi

   its right-hand side in the d-variables is RHSi - z[ i ] * Lambda. Hence,
   for newLambda = Lambda + DLambda the constraint has to become

    z[ i ] * ( newLambda + d ) <= RHSi ,

   that is, its right-hand side in the d-variables is

    RHSi - z[ i ] * newLambda =
           ( RHSi - z[ i ] * Lambda ) - z[ i ] * Dlambda .

  In particular, for box constraints on the Lambda variables [see GetUC()
  and GetUB() in FiOracle.h], the Master Problem has to contain box
  constraints on the d variables (assuming for simplicity every variable
  is box constrained)

    LB - Lambda <= d <= UB - Lambda .
 
  Clearly, when the current point is changed to newLambda = Lambda + DLambda
  these constraints have to become

    LB - newLambda <= d <= UB - newLambda ,

  that is, ( LB - Lambda ) - Dlambda <= d <= ( UB - Lambda ) - Dlambda.  */

   virtual void ChangeCurrPoint( cLMRow DLambda , cHpRow DFi ) = 0;

/*--------------------------------------------------------------------------*/
/** Equivalent to ChangeCurrPoint( ( Tau / t ) * d* , DFi ), where d* is the
   optimal solution of the direction finding subproblem (Tau is a *relative
   step* along d, relative w.r.t. t).  Here, however, no problem with the
   active set can arise, since d* is nonzero only for variables of the active
   set. */

   virtual void ChangeCurrPoint( cHpNum Tau , cHpRow DFi ) = 0;

/*--------------------------------------------------------------------------*/
/** This method allows to change the entries corresponding to variables with
   names between strt (included) and stp (excluded) for all the items in the
   bundle corresponding to the wFi-th component of Fi, in response to a call
   to NDOSolver::ChgSbG( strt , stp, wFi ) [see NDOSlvr.h]. The method has to
   use GetGi() and (possibly) GetADesc() [see FiOracle.h] to acquire directly
   from the FiOracle the new information about the changed variables. If 0 <=
   wFi <= NrFi, the NDO solver is told that only the wFi-th component of Fi
   has changed; if Inf< Index >() > wFi > NrFi then all the components except
   the 0-th have changed, while if wFi == Inf< Index >() then all the
   components, comprised the 0-th, have changed.

   Note that Fi-values, and therefore the Alfas, are assumed not to be
   involved in the changes. */

   virtual void ChgSubG( cIndex strt , Index stp , cIndex wFi ) = 0;

/*--------------------------------------------------------------------------*/
/** This method allows to change the costs of the "easy" component wFi. It is
 * illegal if it is called with wFi not an "easy" component. The method is
 * supposed to call GetBDesc() (with all parameters 0 save that of costs) to
 * retrieve the new costs and change the master problem accordingly.
 *
 * IMPORTANT NOTE: the number of columns in the "easy" component must *not*
 *                 have changed.
 *
 * An issue is that the costs that need be stored in the objective of the
 * master problem are not the original ones, but rather the reduced ones
 * w.r.t. the stability centre Lambda. Since the MPSolver typically has no
 * direct access to it, the second parameter is provided; it must point to
 * a "dense" representation of the stability centre.
 *
 * Since not all MPSolver may implement this advanced feature, a default
 * implementation is given throwing exception. */

   virtual void ChgCosts( Index wFi , cLMRow Lambda ) {
    throw( std::logic_error( "ChgCosts not implemented" ) );
    }

/*--------------------------------------------------------------------------*/
/** This method allows to change the LHS/RHS of the constraints of the "easy"
 * component wFi. It is illegal if it is called with wFi not an "easy"
 * component. The method is supposed to call GetBDesc() (with all parameters 0
 * save these of the LHS/RHS) to retrieve the new LHS/RHS and change the master
 * problem accordingly.
 *
 * IMPORTANT NOTE: the number of rows in the "easy" component must *not*
 *                 have changed.
 *
 * Since not all MPSolver may implement this advanced feature, a default
 * implementation is given throwing exception. */

   virtual void ChgRLHS( Index wFi ) {
    throw( std::logic_error( "ChgCosts not implemented" ) );
    }

/*--------------------------------------------------------------------------*/
/** This method allows to change the lower and upper bounds of the variables
 * of the "easy" component wFi. It is illegal if it is called with wFi not an
 * "easy" component. The method is supposed to call GetBDesc() (with all
 * parameters 0 save these of lower and upper bounds) to retrieve the new
 * lower and upper bounds and change the master problem accordingly.
 *
 * IMPORTANT NOTE: the number of columns in the "easy" component must *not*
 *                 have changed.
 *
 * Since not all MPSolver may implement this advanced feature, a default
 * implementation is given throwing exception. */

   virtual void ChgLUBD( Index wFi ) {
    throw( std::logic_error( "ChgCosts not implemented" ) );
    }

/** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   /// destructor: deletes the timer. it is virtual, as it must be

   virtual ~MPSolver() { delete MPt; }
    
/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  Index MaxBSize;     // max Bundle size
  Index MaxSGLen;     // max items length
  Index CrrSGLen;     // current items length
  Index NrFi;         // number of components of Fi

  std::ostream *MPLog;     // the output stream object
  char MPLLvl;        // the "level of verbosity"

  OPTtimers *MPt;     // total time

/*--------------------------------------------------------------------------*/

 };  // end( class MPSolver )

/*--------------------------------------------------------------------------*/

 }  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MPSolver.h -----------------------------*/
/*--------------------------------------------------------------------------*/
