/*--------------------------------------------------------------------------*/
/*-------------------------- File MinQuad.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Definition of the class MinQuad, implementing the TT (Tall and Thin)
 * algorithm for solving the Quadratic Problems arising (among others) as
 * tentative descent direction finding subproblem within Bundle algorithms
 * for the (linearly constrained) minimization of NonDifferentiable convex
 * functions.
 * The user is assumed to be familiar with the kind of problems that are
 * solved by this code; for a description of the algorithm and the basic
 * notations, refer to
 *
 *  A. Frangioni "Solving Semidefinite Quadratic Problems Within Nonsmooth
 *  Optimization Algorithms" Computers & O.R. 23(1), p. 1099-1118, 1996
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#COR96
 * \endlink
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

#ifndef __MinQuad
 #define __MinQuad  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MinQuad_MACROS Compile-time switches in MinQuad.h
 *  These macros control some important details of the class behavior.
 *  Although using macros for activating features of the interface is not
 *  very C++, switching off some unused features may allow some
 *  implementation to be more efficient in running time or memory.
 *  @{ */

/*------------------------------- LOG_MQ -----------------------------------*/
/** This macro controls how (if any) MinQuad produces a log of its activities
 * on a ostream object set with the method SetMQLog():
 *
 * 0  =>  no log at all (SetMQLog() is even removed from the interface);
 *
 * 1  =>  "basic" log: only the ERRORs are reported;
 *
 * 2  =>  as 1, plus the following indices are kept updated and *succinctly*
 *        (in a row, tab-separated) reported when the object is destroyed:
 *        - the number of calls
 *        - the average dimension of the base (in each step)
 *        - the average dimension of the bundle (in each call)
 *        - the total time spent by SolveQP() (if TIMERS_MQ > 0);
 *
 * 3  =>  as 2, but performances reports are more verbose and the "Faults" are
 *        also reported: Faults can happen in the normal run of the algorithm,
 *        but might also depend on erroneous settings of the parameters;
 *
 * 4  =>  a detailed step-by-step log the algorithm is displayed;
 *
 * 5  =>  as 4, plus some further debug informations are printed. */

#define LOG_MQ 0

/*------------------------------ VARCOEFF ----------------------------------*/
/** The typical (QP) that has to be solved has the linear equality constraint
 *
 *     Sum{ i } Mult[ i ] == 1.
 *
 * However, by setting VARCOEFF > 0 it is turned into the more general
 *
 *     Sum{ i } cf[ i ] * Mult[ i ] == 1
 *
 * where cf[ i ] are user-provided (nonnegative) coefficients, that can be
 * changed between two subsequent iterations. */

#define VARCOEFF 0

/*--------------------------- WASTE_MEM ------------------------------------*/
/** This macro controls the handling of memory within an object of class
 * MinQuad. There are three possible values:
 *
 * 0 = the minimum amount of memory is allocated, and no calls to new() and
 *     delete() are done within the normal execution of the algorithm: some
 *     tasks require considerably more computational effort;
 *
 * 1 = almost the minimum amount of memory is allocated, but some calls to
 *   new() and delete() are done (i.e. the amount of memory used can grow
 *   a little bit during a call to a method); this speeds up calculations;
 *
 * 2 = more memory is used, but calculations are even faster; this option can
 *     be convenient for small problems and/or if the C++ RTS is slow in
 *     allocating and deallocating memory. */

#define WASTE_MEM 2

/*------------------------------ LAZY_Q ------------------------------------*/
/** If LAZY_Q > 0, Q[ i , j ] is calculated, with a call to GiTGj( i , j )
 * [see], only the first time - if any - that it is required. This may save
 * the calculation of some scalar products, but slows down the access to
 * Q[]. If LAZY_Q == 0 instead, all the scalar products corresponding to item
 * i are required in one blow with a call to GiTG( i ... ) [see below] when i
 * is inserted in the bundle [see Add[SubGrad/Constr]() below]. */

#define LAZY_Q 0

/*------------------------------- EXACT ------------------------------------*/
/** If EXACT > 0, the scalar products x * y for x and y in {z1, z2, Lh, Lj}
 * are calculated from scratch at each iteration rather than kept updated in
 * O( 1 ). This is costly, but it avoids accumulation of rounding errors that
 * may lead to faults (and therefore to costly refactorizations) or even to
 * errors. */

#define EXACT 1

/**@}  end( group( MinQuad_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

#if( LOG_MQ )
 #include <iostream>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace MinQuad_di_unipi_it
{
 /** @namespace MinQuad_di_unipi_it
     The namespace MinQuad_di_unipi_it is defined to hold the MinQuad class
     for specially-structured QP optimization as well as its derived classes
     (which are necessary since MinQuad is abstract). It comprises the
     namespace OPTtypes_di_unipi_it. */

 using namespace OPTtypes_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS MinQuad -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/
/** MinQuad solves the following pair of dual (QP)s:
 *
 * (D)  min   f( x ) = (1/2) * ti * x^T * Q * x + Alfa * x
 *            s.t.   Sum{ i = 0 .. n - 1 } e_i x_i = 1 , x >= 0
 *
 * (P)  min   v + 1/(2 ti) || d ||_2^2
 *            s.t.   v * e_i >= G[ i ] * d - Alfa[ i ] , i = 0 .. n - 1
 *
 * where the n-times-n symmetric positive semidefinite matrix Q is given by
 * Q = G^T * G for some n vectors G[ i ] in some space. Each of these 
 * vectors is called an "item", and the set of all items is called "the
 * bundle". The item 'i' is called a "subgradient" if e_i is 1 and a
 * "constraint" if e_i is 0.
 *
 * Any pair of optimal solutions x^*, (v^*, d^*) for (D) and (P) has the
 * following properties:
 *
 *  d^* = - ti * z^*  where  z^* = Sum{ i = 0 .. n - 1 } G[ i ] x^*_i
 *  v^* = - ti * ( || d^* ||_2^2 + Sigma^* )
 *          where    Sigma^* = Sum{ i = 0 .. n - 1 } Alfa[ i ] x^*_i
 *
 * The code works with "virtual" items in the bundle (subgradients or
 * constraints); it never deals directly with any implementation of items
 * (such as the "naive" one, a vector of SgNum's), but it allows the user to
 * choose it. In fact, the code does not even need to know the dimension of
 * the space, that can be dynamically changed [see [Add/Cut]SGSpaceDim() 
 * below]. All the operations on the bundle are performed by referring to
 * "names" [see Add[SubGrad/Constr]() and ResteBundle() below].
 *
 * This is why MinQuad is an *abstract class*, i.e., it has *pure virtual
 * methods* [see the "pure virtual methods" part of the protected interface
 * below]; so that the user can provide it with the required information
 * about the items it cannot directly see. As a consequence, however,
 * instances of MinQuad cannot be built, and a derived class has to be
 * defined where these methods are actually implemented. */

class MinQuad {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

/** Public enum describing the possible return values of SolveQP(). */

 enum MQError { kOK =  0 ,         ///< everything went fine

		kNegativeDelta ,   ///< one kind of numerical issue
		kIncreasingStep ,  ///< another kind of numerical issue
		kNegVarInBase ,    ///< another kind of numerical issue
		kInvalidV ,        ///< another kind of numerical issue
		kLoop ,            ///< another kind of numerical issue

		kQPPrimUnbndd ,    ///< the problem is proven to be unbounded

		kStpTime ,         ///< stopped due to time limit
		kFatal             ///< a fatal numerical error
                };

/** @} ---------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor: initialize the object, but does not allocate memory. This is
 * done is SetMaxDim(), because only there inside the size of the various
 * data structures is known. */

 MinQuad( void );

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

/** This method must be called prior to any other method in the class (but the
 * constructor) to tell the MinQuad object the size of the problem to be
 * solved. The size is specified by the three numbers
 *
 * -  m    = max. dimension of the bundle;
 *
 * -  n    = max. ACTIVE dimension of the bundle (n <= m);
 *
 * - SDim = dim. (number of entries) of the items [optional].
 *
 * This method can be called more than once, but each call causes a complete
 * reset of the internal state of the object, and requires quite a lot of
 * memory allocation/deallocation. Calling the method with m == 0 makes all
 * the memory to be deallocated, and the object to sit down quietly in its
 * corner, waiting for either a new call to SetMaxDim() or its destruction.
 *
 * Note on memory requirements: m, n and SDim controls the way in which memory
 * is allocated by the object.
 *
 * - if SDim is != 0, it is taken as the dimension of the items: no more than
 *   SDim + 3 elements of matrix L are allocated, and therefore their max.
 *   length is SDim + 3; SDim has no impact if SDim + 3 >= m, in fact we
 *   consider SDim = min( SDim + 3 , m );
 *
 * - m cannot be changed: this is the actual max. dimension of the bundle,
 *   and 4 vectors of dimension m (+ 5 of dim. SDim) are allocated until the
 *   object is destroyed;
 *
 * - n can be changed (but must always be <= m), see ChangeCrrBDim(): at any
 *   time, no more than (the current value of) n items can be present, since
 *   only the memory for those items is allocated.
 *
 * If WASTE_MEM > 1, then O( l * n + n * n / 2 ) memory is allocated, 
 * otherwise O( n * n ) memory is allocated, where l = min( SDim + 3 , m ). */

 virtual void SetMaxDim( Index m , Index n , Index SDim = 0 );

/*--------------------------------------------------------------------------*/
/** Sets the max. *active* dimension of the bundle, that must always be <=
 * the max. *absolute* dimension. The max. active dimension can be both
 * increased or decreased: if it is decreased and there are items with "name"
 * >= Newn in the bundle, they are simply deleted. Since this method has no
 * other cost than the allocation/deallocation of the memory, the caller is
 * responsible for a proper memory management in case of memory shortage. */

 void SetCrrBDim( cIndex Newn );

/*--------------------------------------------------------------------------*/
/** Sets a lower bound on the value of the variable v in the QP (the predicted
 * value of the Cutting Plane model). This lower bound is automatically
 * updated when the current point is changed [see MoveAlongD() below]. Passing
 * -Inf< HpNum >()() as the lower bound means that no such bound is available.
 *
 * A finite LB is equivalent to inserting in the bundle an all-0 subgradient
 * with linearization error equals to - LB (LB must be negative), but it is
 * directly handled by the class without requiring to waste space. When a LB
 * is "in", the multipliers corresponding to subgradients returned by
 * ReadMult() may sum to a quantity *strictly smaller* than 1;
 * 1 - Sum{ multipliers } is the optimal multiplier that the all-0 subgradient
 * would get if it would be actually inserted in the bundle, and it is
 * returned by ReadLBMult(). */

 void SetLB( HpNum LB = -Inf< HpNum >() )
 {
  LwrBnd = -LB;
  if( LBIsInB ) {
   Lin = Inf< HpNum >();
   if( LwrBnd == Inf< HpNum >() )
    LBIsInB = 0;
   }
  }

/*--------------------------------------------------------------------------*/
/** In the code, a number eR is used as a measure of the relative precision
 * required from the algorithm: a number x is considered to be zero <=>
 *
 *     ABS( x ) <= eR * Lmu * size( x ).
 *
 * where Lmu is the condition number of the LL^T factorization of the
 * "nonsingular part" of the matrix Q restricted to the current base. In some
 * cases, slightly different formulae are better and they are used instead.
 *
 * eR is automatically changed (either increased or decreased) by SolveQP()
 * to try to face badly conditioned problems; this method allows to "manually"
 * set its current value. Passing 0 (clearly, an impossible value) resets the
 * value of eR to some "default" value. */

 void SetEpsilonR( HpNum NeweR = 0 );

/*--------------------------------------------------------------------------*/
/** At each "outer" iteration, i.e. when a restricted optima is reached, the
 * non-basic items are checked for insertion in the base. The MinQuad class
 * implements a "Steepest - Reverse Entering Order with Early Break" pricing,
 * that is:
 *
 * - in principle, all non-basic items are checked and the item with the *most
 *   negative* directional derivative (the mostly violated dual constraint),
 *   if any, is put in the base; (unlikely) ties are break since the items are
 *   (essentially) checked in Reverse Entering Order, i.e. items that have
 *   entered first are checked last;
 *
 * - however, pricing is interrupted as soon as an item with derivative <
 *   - Brk * BDim * | f( < current point > ) | is found (the BDim factor is
 *   because it is assumed that no multiplier will be larger than 1 / BDim,
 *   hence Brk * BDim is an upper bound on the expected decrease in f()).
 *
 * Setting Brk == 0 (the default if SetPricing() is not called) means that
 * the pricing is in fact a "pure" Reverse Entering Order; the most recently
 * entered item with negative directional derivative is selected. Note that
 * this avoids to calculate the directional derivative for all the items
 * entered before that one, hence it is less expensive. On the other hand,
 * setting Brk == Inf< HpNum >() means using a "pure" Steepest pricing, where
 * just the item with the most negative directional derivative is selected;
 * this *may* help in reducing the overall iterations count, but it is more
 * expensive, and therefore an intermediate value of Brk may help.
 *
 * The value of Brk can be changed at any time. */

 void SetPricing( cHpNum Brk = Inf< HpNum >() ) { EpsRo = Brk; }

/*----------------------------------------------------------------------------
   In typical conditions, obtaining an x such that f( x ) ~= 0 means
   terminating; hence, the extra stopping condition "f() <= MinFVal" can help
   to avoid useless operations and the numerical problems that can occur when
   f() is too near to the machine precision. This method allow to set MinFVal
   to an user-provided value; if it is never called,
   MinFVal == - Inf< HpNum >().

 void SetMinFVal( cHpNum NewMinFVal ) { MinFVal = NewMinFVal; }

----------------------------------------------------------------------------*/
/** Sets the maximum time (in seconds) that each call to SolveQP() is allowed
 * to last. If this method is never called, 0 is assumed, which means "no time
 * limit". Note that this only works if time is actually recorded, see
 * SetMQTime() below. */

 void SetMaxTime( cHpNum NewMaxTime )
 {
  MaxTime = NewMaxTime > 0 ? NewMaxTime : 0;
  }

/*--------------------------------------------------------------------------*/

#if( LOG_MQ )
/** The output of the code is directed onto the ostream object pointed by
 * log: if this method is never called, 'clog' (the standard error
 * bufferized) is used as default. Under Unix-like environments, it can be
 * redirected to a file by using ">&" from the command shell. */

 void SetMQLog( ostream *log )
 {
  MQLog = log;

  #if( LOG_MQ > 2 )
   *MQLog << std::endl << "MinQuad: #bundle = " << MaxBDim << " ~ active = "
	  << CrrBDim << "." << std::endl << std::endl; 
  #endif
  }

#endif

/*--------------------------------------------------------------------------*/
/** SetMQTime() decides whether or not an OPTtimers object is used to
 * record the solution time of the TT algorithm.
 *
 * Note that time accumulates over the calls; calling SetMQTime(), however,
 * resets the counters, allowing to time specific groups of calls. */

 virtual void SetMQTime( const bool TimeIt = true )
 {
  if( TimeIt )
   if( MQt )
    MQt->ReSet();
   else
    MQt = new OPTtimers();
  else {
   delete MQt;
   MQt = 0;
   }
  }

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

/** Add to the current bundle one more subgradient whose "name" in 'n'.
 * Note that nothing is known about the subgradient (item), except that
 * somewhere "outside" there is a subgradient whose name is 'n' and whose
 * linearization error is alfan. */

 virtual void AddSubGrad( cIndex n , cHpNum alfan );

/*--------------------------------------------------------------------------*/
/** Add to the current bundle one more constraint whose "name" in 'n'.
 * Note that nothing is known about the constraint (item), except that
 * somewhere "outside" there is a constraint whose name is 'n' and whose
 * Right Hand Side is alfan.
 *
 * IMPORTANT ADVICE: in constrained problem, a wrong scaling of constraints
 * may have an adverse influence on the behaviour of the algorithm. In fact,
 * even though the constraints C * d <= e and [ k * C ] * d <= k * e (for k a
 * positive scalar) are equivalent in theory, using a vector C whose norm is
 * much different from the "average" norm of the subgradients in the bundle
 * may result in a badly conditioned problem. Since error-recovering procedures
 * are costly, and may cause significant losses of time/precious information,
 * the caller should scale the constraints appropriately. */

 virtual void AddConstr( cIndex n , cHpNum alfan );

/*--------------------------------------------------------------------------*/
///  Returns true if there is an item with name 'n'.

 bool IsThere( cIndex n ) { return( G[ n ] != 0 ); }

/*--------------------------------------------------------------------------*/
///  Returns true if the item with "name" 'n' is a subgradient.

 bool IsASubG( cIndex n ) { return( ( G[ n ] & 1 ) == 1 ); }

/*--------------------------------------------------------------------------*/
///  Returns true if the item with "name" 'n' is a constraint.

 bool IsAConst( cIndex n ) { return( ( G[ n ] & 2 ) == 2 ); }

/*--------------------------------------------------------------------------*/
///  Returns true if a Lower Bound has been set

 bool IsLB( void ) { return( LwrBnd < Inf< HpNum >() ); }

/*--------------------------------------------------------------------------*/
/** There are no restrictions on the use of names, except that they must be
 * 0 <= n < CrrBDim and unique. MaxItemN() returns 1 + the maximum name of
 * an item currently in the bundle if the bundle is nonempty, and 0 if the
 * bundle is empty. That is, MaxItemN() returns the first name for which
 * calling IsThere() [see above] is useless because it will return false. */

 Index MaxItemN( void ) { return( NxtBIdx ); }

/*--------------------------------------------------------------------------*/
/** Eliminates the item named 'which' from the bundle: if which is in Base,
 * then all the data structures (L, z1, z2...) are properly updated. */

 void ResetBundle( cIndex which );

/*--------------------------------------------------------------------------*/
/// Complete reset of the bundle: all items are eliminated.

 void ResetBundle( void );

/*--------------------------------------------------------------------------*/
/** Changes the Alfa (linearization error of a subgradient or RHS of a
 * constraint) of the item with name 'i', adding it DeltaAlfai; that is, at
 * the end of the call
 *
 *    Alfa[ i ] = < old Alfa[ i ] > + DeltaAlfai. */

 virtual void ChangeAlfa( cIndex i , cHpNum DeltaAlfai )
 {
  if( G[ i ] && DeltaAlfai ) {
   Alfa[ i ] += DeltaAlfai;

   if( G[ i ] & 4 )       // IsInBse == 4
    z1Tz2 = Lin = Inf< HpNum >();  // Alfa (=> z2) is changed
   }
  }

/*--------------------------------------------------------------------------*/
/** Performs a "relative" change of all the vector Alfa[], i.e.
 *
 *                  / DeltaAlfa   'i' is a subgradient
 *    Alfa[ i ] += |
 *                 \     0       'i' is a constraint
 *
 * This is the typical correction that ought to be done in bundle algorithms
 * when "nonexact" functions have to be handled, and can be managed faster
 * than the "general" case above. */

 virtual void ChangeAlfa( cHpNum DeltaAlfa );

/*--------------------------------------------------------------------------*/
/** Has the same effect of ChangeAlfa() when the next current point is
 *
 *     next current point = current point - Tau * z^*
 *
 * i.e., - Tau is the absolute step taken along the optimal convex combination
 * of the items (in other words, Tau / ti is the relative step taken along
 * the optimal direction d^*). The method checks for the case Tau == ti (of
 * the latest call to SolveQP()) and in case exploits this to save some
 * computations.
 *
 * Note that a step Tau > ti leads to a point violating the "active"
 * constraints, if any, i.e. the constraints in the optimal "base" [see
 * ReadBase()]; hence, the new point is feasible if and only if the optimal
 * base does not contain constraints.
 *
 * The OBTAINED increase Fi( new current point ) - Fi( old current point )
 * has to be provided in DeltaFi; MoveAlongD() efficiently computes the new
 * Alfa's (O(1) for each of them), but it is faster in recalculating some of
 * the internal data structures.
 *
 * MoveAlongD() uses the well-known Alfa updating formula
 *
 *     Alfa[ i ] = Alfa[ i ] + Tau * z^* * G_i + DeltaFi
 *
 * (since Tau * z^* * G_i = - ( Tau / ti ) * d^* * G_i), hence it needs the
 * scalar product between z^* and the item 'i' to perform the task. Such a
 * scalar product is available for free as a byproduct of the latest
 * optimization if i was already in the bundle at that time, but it is
 * clearly *not* available if i has been inserted after the latest call to
 * SolveQP().
 *
 * Therefore, in order to make MoveAlongD() work properly, one needs
 *
 * - either that no new subgradients [constraints] be added between the end of
 *   SolveQP() and the call to MoveAlongD();
 *
 * - or that SetGTz( i , .. ) is invoked after Add[Subgrad/Constr]( i , .. )
 *   to provide the required information;
 *
 * - that Alfa has not been changed in the meantime, i.e. that no calls to
 *   *any version of* ChangeAlfa() have occurred;
 *
 * - that a previous direction d exists, i.e. SolveQP() has already been
 *   called and was successful. */

 virtual void MoveAlongD( cHpNum Tau , cHpNum DeltaFi );

/*--------------------------------------------------------------------------*/
/** Writes the current vector of Alfas, i.e., Alfa[ i ] is the item with
 * "name" i, in the vector NewAlfa. */

 virtual void ReadAlfa( HpRow NewAlfa );

/*--------------------------------------------------------------------------*/
/** Returns a read-only pointer to the current vector of Alfas, i.e.,
 * Alfa[ i ] is the item with "name" i. */

 virtual cHpRow ReadAlfa( void ) { return( Alfa ); }

/*--------------------------------------------------------------------------*/
/// returns the current value of the Lower Bound

 cHpNum ReadLB( void ) { return( - LwrBnd ); }

/*--------------------------------------------------------------------------*/
/** Provides the scalar product between z^* and the newely inserted item 'i'
 * prior to a call to MoveAlongD() [see above]. */

 virtual void SetGTz( cIndex i , cHpNum GTzi ) { GTz[ i ] = GTzi; }

/*--------------------------------------------------------------------------*/
/** Increases the dimension of the space of items, implicitly "adding" one
 * entry to all the items (i.e., a dual variable of the problem).
 *
 * NewDim[ i ] is supposed to be G[ i ][ h ], i.e. the "missing" component
 * of the item whose "name" is 'i'; the matrix Q[][] is updated for all the
 * already calculated entries (everywhere if ! LAZY_Q);
 *
 * If UpdtB == true, the current base is properly updated; this can at most
 * decrease (by one) Dependent, hence it cannot create problems. If UpdtB ==
 * false, instead, the current base is *not* updated, so after a(ny number of
 * consecutive) calls to AddSGSpaceDim( ... , false ), one among UpdateB()
 * and ResetB() *must* be called; failure to do that will generate incorrect
 * results. This feature is provided for the case when the size of the items
 * changes a lot, so that recomputing the base from scratch is more convenient
 * than keeping it updated. Since a base update costs ~ ReadBDim()^2 / 2 while
 * rebuilding it costs ~ ReadBDim()^3 / 3, rebuilding is convenient if the
 * dimension of the space changes more than ~ 2 ReadBDim() / 3. */

 void AddSGSpaceDim( cSgRow NewDim , const bool UpdtB = true );

/*--------------------------------------------------------------------------*/
/** Performs the same task of AddSGSpaceDim( cSgRow , const bool ) for the
 * case where MinQuad is used as the "core solver" for the box-constrained
 * case. It handles the case when a box-constrained active variable have to
 * be *deleted* from the Base, where lbh is its lower bound (or - its upper
 * bound). */

 void AddSGSpaceDim( cSgRow NewDim , cHpNum lbh , const bool UpdtB = true );

/*--------------------------------------------------------------------------*/
/** Decreases the dimension of the space of items, "cutting away" one of the
 * entries of the items (i.e., a dual variable of the problem).
 *
 * OldDim[ i ] is supposed to be G[ i ][ h ], i.e. the "undesidered"
 * component of the item whose "name" is 'i': the matrix Q[][] is updated
 * for all the already calculated entries (everywhere if ! LAZY_Q);
 *
 * If UpdtB == true, the current base is properly updated; this can increase
 * (by one) Dependent, so that if Dependent is == 2 it may become == 3,
 * causing the last item in base to be discarded. If UpdtB == false, instead,
 * the current Base is *not* updated, so after a(ny number of consecutive)
 * calls to CutSGSpaceDim( ... , false ), one among UpdateB() and ResetB()
 * *must* be called; failure to do that will generate incorrect results. This
 * feature is provided for the case when the size of the items changes a lot,
 * so that recomputing the base from scratch is more convenient than keeping
 * it updated. Since a base update costs ~ ReadBDim()^2 / 2 while rebuilding
 * it costs ~ ReadBDim()^3 / 3, rebuilding is convenient if the dimension of
 * the space changes more than ~ 2 ReadBDim() / 3. */

 void CutSGSpaceDim( cSgRow OldDim , const bool UpdtB = true );

/*--------------------------------------------------------------------------*/
/** Performs the same task of CutSGSpaceDim( cSgRow , const bool ) for the
 * case where MinQuad is used as the "core solver" for the box-constrained
 *case. It handles the case when a box-constrained active variable have to be
 * *deleted* from the Base, where lbh is its lower bound (or - its upper
 * bound). */

 void CutSGSpaceDim( cSgRow OldDim , cHpNum lbh , const bool UpdtB = true );

/*--------------------------------------------------------------------------*/
/** The current matrix Q is substituted with an entirely new Q, whose elements
 * are requested with calls to GiTG[j](). This can be a cheaper alternative to
 * the use of [Add/Cut]SGSpaceDim() in the case of many changes to the Q
 * matrix, or be simply used to insert an entirely uncorrelated problem.
 *
 * After a call to ChangeQ(), one among UpdateB() and ResetB() *must* be
 * called; failure to do that will generate incorrect results.
 *
 * This method can also be used to recover from failures of SolveQP() due to
 * accumulation of numerical errors in *the matrix Q*. These errors may
 * happen if QuNum is a floating point type and [Add/Cut]SGSpaceDim() are
 * used to update Q many times. MinQuad implements a sophisticated machinery
 * for handling this problem in the lower trapezoidal factor L of Q, but Q is
 * assumed to be "exact": one possibility for reacting to a failure of
 * SolveQP() is therefore to call ChangeQ(), for having Q "refreshed". */

 void ChangeQ( void );

/*--------------------------------------------------------------------------*/
/** After that data of the problem has been changed with any sequence of calls
 * to [Add/Cut]GSpaceDim( ... , false ) or ChangeQ(), either this method or
 * ResetB()  must be called to ensure that the data structures describing the
 * current base are updated.
 *
 * This method rebuilds (as much as possible) the data structures for the
 * latest optimal base, allowing a "warm start" of the method. This makes
 * sense if the changes were "few", so that one can expect the latest
 * optimal base to be a "good" one. */

 void UpdateB( void );

/*--------------------------------------------------------------------------*/
/** After that data of the problem has been changed with any sequence of calls
 * to [Add/Cut]GSpaceDim( ... , false ) or ChangeQ(), either this method or
 * UpdateB() must be called to ensure that the data structures describing the
 * current base are updated.
 *
 * This method starts from a clean new (empty) base, which may be better than
 * using UpdateB() if the changes were "many", so that one cannot expect the
 * latest optimal base to be a "good" one. */

 void ResetB( void );

/*--------------------------------------------------------------------------*/

#if( VARCOEFF )

/** If VARCOEFF == 1, the primal linear equality constraint is 
 *
 *     Sum{ i } cf[ i ] * Mult[ i ] == 1
 *
 * where the cf[ i ] are user-provided (nonnegative) coefficients; they are
 * defaulted to 1, but they can be changed with this method.
 *
 * The coefficient of the primal variable relative to the item whose "name"
 * is 'i' is set to ci: note that ci == 0 ==> 'i' is a constraint, and it is
 * actually possible to change a subgradient into a constraint and vice-versa
 * by setting ci == 0 (> 0) with this method. */

 void ChangeCoefficient( cIndex i , cHpNum ci );

#endif

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

/** Solves the problem with the current data, i.e. the current bundle, Alfa
 * and ti. The second level of dynamic tolerances adjustment is implemented
 * here inside; the relative precision parameter eR [see EpsilonR() above]
 * can be changed to accomplish with increase/decrease "requests".
 *
 * The method returns:
 *
 * - kOK for normal termination: a finite optimal solution has been found.
 *
 * - kQPPrimUnbndd if the problem is proven to be primal unbounded. This can
 *   only happen if there are constraints with *negative* right hand side in
 *   the bundle. In this case, the method ReadInfDir() can be used together
 *   with the solution returned by ReadMult() to construct an unbounded ray
 *   of the problem.
 *
 * - kStpTime if the maximum time limit set with SetMaxTime() has elapsed
 *   before the algorithm can find an optimal solution (or prove the problem
 *   to be unbounded). This is not a "serious" condition, in that it can be
 *   "cured" by just calling again (possibly repeatedly) the method to give it
 *   more time (since warm starts are done, this has a limited overhead w.r.t.
 *   directly calling the method once with a longer time limit).
 *
 * - kFatal if even adjusting eR was not sufficient to handle the numerical
 *   problems, i.e., either eR is increased/decreased out of the bounds or
 *   two contrasting requests (increase eR and then decrease it or 
 *   vice-versa) have been detected. */

 virtual MQError SolveQP( HpNum ti );

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading results
    @{ */

/** After (correct) termination of SolveQP(), ReadzNorm() returns
 * || z^* ||^2. Together with ReadSigma() this allows to obtain
 * the optimal value of the (dual and primal) objective function. */

 virtual HpNum ReadzNorm( void ) { return( Quad ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** After (correct) termination of SolveQP(), ReadSigma() returns Sigma^*.
 * Together with ReadzNorm() this allows to obtain the optimal
 * value of the (dual and primal) objective function.
 *
 * If IncldCnst == false, ReadSigma() returs the value of Sigma restricted to
 * only those items that are actually subgradients, excluding the
 * constraints.*/

 virtual HpNum ReadSigma( const bool IncldCnst = true )
 {
  if( IncldCnst || ( ! CBDim ) )
   return( Lin );
  else {
   HpNum tS = 0;
   cHpRow tM = Mult;
   cIndex_Set tB = Base;
   for( Index h ; ( h = *(tB++) ) < Inf< Index >() ; tM++ )
    if( ! IsAConst( h ) )
     tS += (*tM) * Alfa[ h ];

   return( tS );
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Returns v^*. 
 *
 * Note that if there are no subgradients in the bundle, the problem that is
 * solved is different from the standard one, as the simplex constraint cannot
 * be satisfied. Thus the constraint is temporarily removed (in dual terms, a
 * constraint v >= 0 is added, which renders that a redundant <= constraint).
 * In this case v^* is not significant (it is the dual optimal value of a
 * non-existent constraint); this is signalled by having Readv() returning
 * Inf< HpNum >().
 *
 * However, note that a finite lower bound counts as a subgradient. */

 HpNum Readv( void )
 {
  HpNum tv = - LastRo;

  if( tv < Inf< HpNum >() )
   tv *= PrvsTi;

  return( tv );
  }

/*--------------------------------------------------------------------------*/
/** After (correct) termination of SolveQP(), returns a pointer to a
 * read-only vector containing the non-zero entries of the optimal solution
 * x^*. The length of the vector and the items to which each entry belongs
 * are returned respectively by ReadBDim() and ReadBase().
 *
 * It is required that ResetBundle( i ) is *not* called for any i in the
 * optimal base between the latest call to SolveQP() and this call to
 * ReadMult(). */

 cHpRow ReadMult( void ) { return( Mult ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Returns the number of non-zero entries of x^*, i.e., the length of the
 * vector returned by ReadMult(). */

 Index ReadBDim( void ) { return( BDim ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Returns the optimal base. That is, ReadMult()[ i ] contains the multiplier
 * relative to theitem whose "name" is ReadBase()[ i ] for
 * i = 0 .. ReadBDim() - 1. The vector pointed by ReadBase() is
 * Inf< Index >()-terminated, i.e., ReadBase()[ ReadBDim() ] == Inf< Index >().
 */

 cIndex_Set ReadBase( void ) { return( Base ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// Returns true if the "name" 'n' is in the optimal base.

 bool IsInBase( cIndex n ) { return( ( G[ n ] & 4 ) == 4 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Returns the number of items in the optimal base that are constraints.
 * Note that ReadCBDim() <= ReadBDim(), and the inequality holds strictly
 * (i.e., at least one subgradient must be in the optimal base) whenever
 * there are subgradients in the bundle and the solution is finite. */

 Index ReadCBDim( void ) { return( CBDim ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** If SolveQP() returns kQPPrimUnbndd [see above], then this method returns
 * an unbounded descent direction v^* to the problem, meaning that
 *
 *      x^* + \beta v^*
 *
 * is a feasible solution to the problem for each non-negative value of \beta,
 * and that the objective function value is strictly decreasing as \beta
 * increases. Note that v^* is in the same "sparse" format as x^*, i.e.,
 * only the nonzero entries are returned, the corresponding indices being
 * obtainable by ReadBase(). Clearly, v^* must be zero on all the entries
 * relative to subgradients; however, there will in general be subgradients
 * in ReadBase()[], hence the corresponding entries of v^* will be == 0 (up
 * to the precision set by SetEpsilonR(), so be sure to actually zero-out
 * those entries manually if you are really going to use v^* for "large"
 * values of \beta). Note that v^* is contained into a temporary, which may
 * be lost after any call to any method of the class (but ReadinfDir(),
 * obviously); so, be sure to save it if you need it. */

 cHpRow ReadInfDir( void ) { return( tmpv ); }

/*--------------------------------------------------------------------------*/
/** Returns the optimal multiplier that the all-0 subgradient corresponding
 * to the Lower Bound [see SetLB() above] would get if it would be inserted
 * in the bundle. */

 HpNum ReadLBMult( void ) { return( LBIsInB ? LBMult : 0 ); }

/*--------------------------------------------------------------------------*/
/** After (correct) termination of SolveQP(), returns
 *
 *   z^* * G[ i ]       where, as usual, G[ i ] = item whose "name" is i. */

 virtual HpNum ReadGTz( cIndex i ) { return( GTz[ i ] ); }

/*--------------------------------------------------------------------------*/
/** Let x( t ) be the optimal solution of (D) as a function of the parameter
 * t, and let v( t ), z( t ) be the corresponding solutions of (P) (actually,
 * d( t ) is a solution of (P), but d( t ) = - t * z( t )); there exists an
 * interval [tMin, tMax] [see SensitAnals3()] in which the currently optimal
 * base stays optimal. This function reports three numbers v1, v2 and v3 such
 * that within that interval
 *
 *     v( t )               =    t * v1 + v2
 *
 *     || z( t ) ||^2       =  - v1 + ( 1 / t )^2 * v3
 *
 *     Alfa * x( t )        =  - v2 - ( 1 / t ) * v3
 *
 * Outside the "stability interval", the above expressions give *lower
 * approximations* of the real figures. */

 virtual void SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Given v1 and v2 reported by SensitAnals1(), calculates the two vectors
 * x1 and x2 such that x( t ) = x1 + ( 1 / t ) * x2 is the primal optimal
 * solution of the problem as a function of t, within the "stability interval"
   [tMin, tMax] [see SensitAnals3()]. */

 void SensitAnals2( HpNum v1 , HpNum v2 , HpRow x1 , HpRow x2 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/** Given v1 and v2 reported by SensitAnals1() and the two vectors x1 and x2
 * calculated by SensitAnals2(), calculates the "stability interval"
 * [tMin, tMax] in which the currently optimal base stays optimal, so that
 * v( t ), x( t ) etc. as reported by SensitAnals[1/2]() are "exact". */

 void SensitAnals3( cHpNum v1 , cHpNum v2 , HpRow x1 , HpRow x2 ,
		    HpNum &tMin , HpNum &tMax );

/*--------------------------------------------------------------------------*/
/// Return the user and sistem time (in seconds) spent so far by SolveQP().

 void MQTime( double &t_us , double &t_ss )
 {
  t_us = t_ss = 0;
  if( MQt )
   MQt->Read( t_us , t_ss ); 
  }

/*--------------------------------------------------------------------------*/
/// Return the total time (in seconds) spent so far by SolveQP().

 double MQTime( void ) { return( MQt ? MQt->Read() : 0 ); }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/** Returns the current value of eR (actually, the returns eR * Lmu * BDim),
 * see SetEpsilonR() above. */

 HpNum EpsilonR( void ) { return( eR * Lmu * ( BDim ? BDim : 1 ) ); }

/*--------------------------------------------------------------------------*/
/** Returns Q[ i ][ j ]. With the "lazy" calculation method, the first time
 * (if any) that Q[ i ][ j ] is requested the appropriate scalar product is
 * fetched by invoking GiTG[j]() [see below]. Note that invoking LowQ() with
 * either 'i' or 'j' being names of items currently not in the bundle is
 * wrong. */

 QuNum LowQ( cIndex i , cIndex j )
 {
  #if( LAZY_Q )
   if( i >= j ) {
    if( Q[ i ][ j ] == -QuINF )
     Q[ i ][ j ] = GiTGj( j , i );

    return( Q[ i ][ j ] );
    }
   else {
    if( Q[ j ][ i ] == -QuINF )
     Q[ j ][ i ] = GiTGj( i , j );

    return( Q[ j ][ i ] );
    }
  #else
   return( i >= j ? Q[ i ][ j ] : Q[ j ][ i ] );
  #endif
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// LowQ( i ) = LowQ( i , i )

 QuNum LowQ( cIndex i )
 {
  #if( LAZY_Q )
   if( Q[ i ][ i ] == -QuINF )
    Q[ i ][ i ] = GiTGj( i , i );
  #endif

  return( Q[ i ][ i ] );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

/// Memory deallocation. Statistics (if any) are printed.

 virtual ~MinQuad();

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PURE VIRTUAL METHODS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for providing the scalar products.
 *
 * These methods are "pure virtual", hence they *must* be implemented in
 * any derived class.
 * @{ */

/** If Q is calculated "lazyly", calculation of its elements is performed
 * "asynchronously" when (and only if) it is needed by a call to
 *
 *     GiTGj( i , j ) = Sum{ k = 0 .. n - 1 } G[ i ][ k ] * G[ j ][ k ].
 *
 * In the simplest case (there exist a "naive" form of items as n-vectors of
 * SgNum's, and the item whose "name" is 'i' is just the i-th row of the
 * matrix G) the above scalar product can be calculated with [see OPTvect.h]
 *
 *     ScalarProduct( G[ i ] , G[ j ] , n )
 *
 * However, the objects of class MinQuad doesn't need to know anything about
 * the "real" implementation of the items; therefore, the implementation of
 * GiTG[j]() is *not* provided, and it's due to the derived class. The
 * rationale is that the items may have some structure (e.g. sparsity) that
 * can be exploited to fasten the calculation of the scalar product.
 *
 * To ease the implementation, it is guaranteed that GiTGj( i , j ) will
 * always be called with i <= j.
 *
 * Conversely, if Q is not calculated "lazily", all the elements of row Q[ i ]
 * are computed together when the item of "name" 'i' is inserted in the bundle
 * (with Add[SubGrad/Constr]()); hence, in this case the derived class is
 * rather required to implement a method that gives them all together - this
 * may fasten the implementation, e.g. in the case that the matrix G is stored
 *  by rows (variables-wise) instead of columns (items-wise).
 *
 * The required effect of GiTG( i , Qi , iMax ) is described by the following
 * pseudo-code:
 *
 *   for( j = 0 ; j < iMax ; j++ )
 *    if( IsThere( j ) )
 *     Qi[ j ] = GiTGj( i , j );
 *
 * where GiTGj( i , j ) here means "whatever GiTGj() would have returned
 * if called with i , j". This method is called each time a new item is
 * inserted in the bundle, and in this case iMax will be just the largest
 * "name" among all items currently in the bundle (comprised i). However,
 * the method can also be called during a "from-scratch" recomputation of
 * the matrix Q [see ChangeQ()], and in this case iMax == i + 1. */

#if( LAZY_Q )

 virtual QuNum GiTGj( cIndex i , cIndex j ) = 0;

#else

 virtual void GiTG( cIndex i , QuRow Qi , cIndex iMax ) = 0;

#endif

/** @} ---------------------------------------------------------------------*/
/*------------------ "REAL" PROTECTED PART OF THE CLASS --------------------*/
/*--------------------------------------------------------------------------*/
/** @name "Normal" protected part of the class
 *
 * The standard user should not care about the following part: users who
 * need to extend the code by deriving a new class may use these methods
 * and data structures. It is *dangerous* to *modify* the data structures,
 * while it safe to read them.
 * @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** For efficiency, the Alfa vector is protected and it can be modified by
 * derived classes: however, when the Alfa[] of "basic" items are changed,
 * AlfaChanged() must be called to warn the base class. */

 void AlfaChanged( void ) { z1Tz2 = Lin = Inf< HpNum >(); }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

 Index MaxBDim;        ///< Max. dimension of the bundle
 Index ActBDim;        ///< Number of items currently in the bundle
 Index ActCNum;        /**< <= ActBDim: how many of the items currently in
			  the bundle are constraints */
 Index NxtBIdx;        /**< Next available name for a new item == max name
			  of an item in the bundle + 1 */
 MQError QPStatus;     ///< termination code of SolveQP()

 HpRow Mult;           ///< Primal variables
 Index_Set Base;       ///< Base (set of strictly positive variables)
 Index BDim;           ///< n. of items in Base
 HpNum Quad;           ///< Mult^T * Q * Mult
 HpNum Lin;            ///< Alfa^T * Mult

 HpRow GTz;            /**< Scalar product between the (optimal) direction
			  and the items in the bundle */
 HpRow Alfa;           ///< Linearization errors

 HpNum PrvsTi;         ///< Value of ti in the latest call of CalcOptDir()
 HpNum eR;             ///< Relative error for "== 0" tests
 // HpNum MinFVal;     < Stop if f() <= MinFVal

 HpRow tmpa;           /**< Temporary vectors, at least as long as the max
			  bundle dimension: can be freely used by derived
			  classes, since its content need not to be
			  conserved between two calls to any method */
 HpRow tmpv;           ///< like tmpa

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void GjTGBk( cIndex j , cIndex_Set Bk , Index Dim , HpRow v );

/*--------------------------------------------------------------------------*/

 void RebuildBase( void );

/*--------------------------------------------------------------------------*/

 void CalcOptDir( cHpNum ti );

/*--------------------------------------------------------------------------*/

 bool Feasible( void );

/*--------------------------------------------------------------------------*/

 HpNum vPerAlfa( cHpRow v );

/*--------------------------------------------------------------------------*/

 HpNum CalculateDelta( cIndex i );

/*--------------------------------------------------------------------------*/

 void AddSubGradToBase( Index BD , cHpNum Delta );

/*--------------------------------------------------------------------------*/

 void ChkNwLinInd( void );

/*--------------------------------------------------------------------------*/

 void CutSubGradFromBase( cIndex s );

/*--------------------------------------------------------------------------*/

 void MoveSubgradToLastPos( cIndex s );

/*--------------------------------------------------------------------------*/

 void UpdateLmu( cIndex s );

/*--------------------------------------------------------------------------*/

 void ComputeQuad( void );

/*--------------------------------------------------------------------------*/

 void MemDealloc( void );

/*--------------------------------------------------------------------------*/

 void CheckB( void );

 void CheckQ( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 Index_Set Order;      // Order[ 0 .. ActBDim - 1 ] contains the "names" of
                       // all the items currently in the bundle
 char *G;              // Which position contains what in the bundle

 QuMat Q;              // Q = G^T * G (the lower triangular part)

 HpMat L;              // L_B * L_B^T = Q_B (L may be trapezoidal)

 HpRow z1;             // L_{B'}^{-1} * e_{B'}    \ B' the "linearly
 HpRow z2;             // L_{B'}^{-1} * Alfa_{B'} / independent" part of B
 #if( VARCOEFF )
  HpRow cf;            // Coefficients of the linear equality constraint
 #endif

 HpNum f;              // f( x ) = (1/2) * ti * Quad + Lin
 HpNum Minf;           // Min. value of f() got in this run
 HpNum LastRo;         // Latest value of ro = - v / ti
 HpNum z1Tz2;          // z1 * z2
 HpNum z1Norm;         // || z1 ||^2
 HpNum LhTz1;          // L[ BDim - 1 ] * z1  \ used by SensytAnalsX()
 HpNum LhTz2;          // L[ BDim - 1 ] * z2  / if Dependent == 1
 #if( ! EXACT )
  HpNum LhNorm;        // || Lh ||^2  \ set by CalculateDelta, on the base
  HpNum LjNorm;        // || Lj ||^2  / of the current value of Dependent
 #endif
 HpNum LwrBnd;         // the Lower Bound on v
 HpNum LBMult;         // the multiplier corresponding to LwrBnd

 Index Dependent;      // How many items in Base are lin. dep. (0, 1 or 2)
 Index LBIsInB;        // 1 <=> if the Lower Bound is "in Base"
 Index SgSpDim;        // Max. dimension of a base
 Index CrrBDim;        // Max. *active* dimension of the bundle

 Index Ename;          // The name of the item that has entered last
 Index Ewher;          // The position in Base of the item Ename
 Index Eordr;          // The position in Order of the item Ename
 Index TLDim;          // The number of "taboo" items
 Index LiMax;          // Index of the max ...
 Index LiMin;          // ... and of the min. diagonal element of L
 HpNum Lmu;            // Lmu = L[ iMax ][ iMax ] / L[ iMin ][ iMin ], i.e.
                       // the condition number of L
 HpNum EpsRo;          // The "steepest" pricing is interrupted if
                       // df < - EpsRo * BDim * | f |
 Index CBDim;          // <= BDim - 1: how many of the items in Base are NOT
                       // subgradients (i.e. are constraints)

 HpNum MaxTime;        // maximum time in each call to SolveQP()
 bool AddedOne;        // true if some new item has been added to the bundle
                       // since the last call
 bool LocalOptima;     // true if CalcOptDir() ended in a local optima

 #if( ! LAZY_Q )
  QuRow tempQ;         // Temporary of QuNum's
 #endif

 #if( ! WASTE_MEM )
  HpRow tmpr;          // Extra temporary for MoveSubgradToLastPos()
 #endif

 #if( LOG_MQ )
  std::ostream *MQLog;         // the output stream object

  #if( LOG_MQ > 1 )
   unsigned long int Calls;    // Calls counter
   unsigned long int Step;     // Step counter
   unsigned long int Success;  // Succesfull calls counter
   float SumBDim;              // Sum( all steps i ) BDim_i
   float SumActBDim;           // Sum( all calls i ) ActBDim_i
   float SumGSBase;            // Sum( all steps i ) CBDim_i / BDim_i
   float SumGSbundle;          // Sum( all calls i ) ActCNum_i / ActBDim_i
   unsigned long int Insert;   // Total n. of insertions
   unsigned long int Delete;   // Total n. of deletions
   unsigned long int SGSIncr;  // Total n. of calls to AddSGSpaceDim()
   unsigned long int SGSDecr;  // Total n. of calls to CutSGSpaceDim()
  #endif
 #endif

 OPTtimers *MQt;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class MinQuad )

/*--------------------------------------------------------------------------*/

 }  // end( namespace MinQuad_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MinQuad.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File MinQuad.h -------------------------------*/
/*--------------------------------------------------------------------------*/
