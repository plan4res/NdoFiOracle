/*--------------------------------------------------------------------------*/
/*--------------------------- File BMinQuad.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Definition of the class BMinQuad, implementing the BTT (Boxed Tall and
 * Thin) algorithm for solving the box-constrained Quadratic Problems arising
 * (among others) as tentative descent direction finding subproblem within 
 * box-constrained Bundle algorithms for the minimization of NonDifferentiable
 * convex functions. It derives from the MinQuad class [see MinQuad.h] since
 * it uses the "unconstrained" TT algorithm as a subroutine.
 *
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

#ifndef __BMinQuad
 #define __BMinQuad  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup BMinQuad_MACROS Compile-time switches in BMinQuad.h
 *  These macros control some important details of the class behavior.
 *  Although using macros for activating features of the interface is not
 *  very C++, switching off some unused features may allow some
 *  implementation to be more efficient in running time or memory.
 *  @{ */

/*------------------------------- LOG_BMQ ----------------------------------*/
/** This macro controls how (if any) BMinQuad produces a log of its activities
 * on a ostream object set with the method SetBMQLog() [see below]:
 *
 * 0  =>  no log at all (SetBMQLog() is even removed from the interface);
 *
 * 1  =>  "basic" log: only the ERRORs are reported;
 *
 * 2  =>  as 1, plus the following indices are kept updated and *succinctly*
 *        (in a row, tab-separated) reported when the object is destroyed:
 *	  - the number of variables
 *	  - the number of calls
 *	  - the average number of constrained variables (in each call)
 *	    (+ what is reported by MinQuad, if any);
 *
 * 3  =>  as 2, but performances reports are more verbose and the "Faults" are
 *        also reported: Faults can happen in the normal run of the algorithm,
 *	  but might also depend on erroneous settings of the parameters;
 *
 * 4  =>  a detailed step-by-step log of the algorithm is printed;
 *
 * 5  =>  the log is very verbose, telling the name of all variables that
 *        move in and out of the "base" of the BTT algorithm. */

#define LOG_BMQ 0

/*----------------------------- LAZY_D -------------------------------------*/
/** At any iteration of the algorithm, the value of the *constrained* entries
 * of the solution d[] must be calculated; to do that, it is necessary to
 * know all the entries z[ i ] of the "aggregated subgradient" z for each name
 * `i' that has been declared as a constrained variable. How such entries are
 * calculated depends on LAZY_D: for each value of the switch, a different
 *  form of CalculateZ() [see] is used.
 *
 * LAZY_D == 0  all such entries of z have to be computed in one blow by a
 *              single call to CalculateZ( < all > ): note that some of the
 *		entries may not be actually required in the current iteration;
 *
 * LAZY_D == 1  the entries are divided in two sets, one corresponding to
 *              "active" variables and the other to "inactive" ones, with
 *		calls to CalculateZ( < a set > ); if MaxVar[Add/Rmv] ==
 *		Inf< Index >() [see below] all the variables in the group will
 *		be needed each time the set is required (while those in the
 *		other set may not be needed), but if MaxVar[Add/Rmv] <
 *		Inf< Index >() some of the entries may still be calculated
 *		without a need;
 *
 * LAZY_D == 2  the value of each necessary entry is requested by a call to
 *              CalculateZ( < one > ): in this case, entries corresponding to
 *		*unconstrained* variables are *never* required. */

#define LAZY_D 0

/*------------------------------ TWOSIDED ---------------------------------*/
/** Often, box constraints l[ i ] <= d[ i ] <= u[ i ] are even "too general",
 * since simple lower bounds l[ i ] <= d[ i ] suffice instead; by setting
 * TWOSIDED == 0, we allow the code to deal only with lower bounds, saving
 * space and time. */

#define TWOSIDED 1

/*------------------------------ SIGNAL_XXCHG ------------------------------*/
/** In some pure virtual methods of the class [see GiTG[j](), CalculateZ() and
 * GiTLB()], information is requested from the outside. Two groups of
 * protected fields of the class are used by these methods, respectively
 * (Mult, Base, BDim) by CalculateZ[h]() and (MBase2, MB2Dim, Base2, B2Dim)
 * by GiTG[j]().
 * If SIGNAL_XXCHG is set to one, then a pure virtual method XXHasChgd() is
 * added to the protected interface of the class: this method is called each
 * time the corresponding group of fields changes, giving an "hook" to
 * derived classes to react to the changes if necessary. */

#define SIGNAL_MBCHG 0

#define SIGNAL_B2CHG 0

/*------------------------------ BEXACT ------------------------------------*/
/** Analogous to the EXACT macro of MinQuad.h: if BEXACT == 0, some numbers
 * are kept updated in O( 1 ) per iteration rather than recalculated from
 * scratch in linear time at each iteration. */

#define BEXACT 0

/*----------------------------- CNDVD_TMP ----------------------------------*/
/** In the code, one temporary is (rarely) used to build the primal solution
 * d; it is only used, if ever, within SolveQP() [see below]. If
 * CNDVD_TMP == 0, SDim LMNum's are allocated in the constructor for this
 * temporary, and deallocated in the destructor. If CNDVD_TMP > 0, no memory
 * is allocated by BMinQuad and the caller is required to supply the memory
 * with SettmpD() and reclaim it with GettmpD() [see]. */

#define CNDVD_TMP 1

/**@}  end( group( BMinQuad_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MinQuad.h"

#if( LOG_BMQ )
 #include <iostream>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace MinQuad_di_unipi_it
{
 using namespace OPTtypes_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BMinQuad ------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/
/** BMinQuad solves the modified version of the (QP)
 *
 *   (P)  min   v + 1/(2 ti) || d ||_2^2
 *              s.t.   v * e_i >= G[ i ] * d - Alfa[ i ] , i = 0 .. n - 1
 *                     l[ j ] <= d[ j ] <= u[ j ]     \forall j
 *
 * where some of the lb[ i ] can be - INF and some of the ub[ i ] can be + INF
 * (actually, finite upper bounds are only allowed if TWOSIDED > 0, (see), and
 * the corresponding dual problem. Although these problems are very similar to
 * the standard ones solved by MinQuad, there are some differences, the main
 * one being that it is no longer true that d^* = - ti * z^*, but rather
 *
 *   d[ i ]^* = min( u[ i ] , max( l[ i ] , - ti * z[ i ]^* ) ).
 *
 * The user is assumed to be familiar with the base class MinQuad [see 
 * MinQuad.h], and in particular with the issue of the "virtual items". As a
 * resume, MinQuad never deals directly with items (e.g. vectors of SgNum's),
 * but rather allow the user to add and remove items from the bundle by means
 * of a symbolic "name". The necessary information is then requested by means
 * of the GiTG[j]() method.
 *
 * Analogously, BMinQuad works with "virtual variables: all it knows is that
 * a certain number of variables (<= SDim) have been "declared" by invoking
 * AddVars() or InitialSetUp() [see]. Variables can be either *constrained*
 * or *unconstrained*; in the former case, lower [and upper] bounds have to
 * be given, that can however also be -[+]  infinity. BMinQuad identifies
 * these variables with the "name" `i' (in the range 0 .. SDim - 1) that the
 * user has given them, and that must be unique; then, it may ask for the
 * i-th entry of something (a subgradient, the whole matrix G, the aggregated
 * subgradient z ...), requiring the entry relative to the variable whose
 * "name" is `i'. A constrained variable can be turned into an unconstrained
 * one by invoking MakeVarUnCnstr(), and vice-versa with MakeVarCnstr(); this
 * can be repeated any number of times.
 *
 * BMinQuad "thinks" that the "declared" variables are the *only* ones of the
 * problem; *this is not necessarily true*, however, since other unconstrained
 * variables can be "implicitly" dealt with by "silently" adding their
 * contribution to the results of GiTG[j]() [see] as in the base class
 * MinQuad. In fact, "unnamed" unconstrained variables can be added / removed
 * from the problem by calling the methods [Add/Cut]SGSpaceDim().
 *
 * As its base class MinQuad, BMinQuad is an *abstract class*, since it has
 * *pure virtual methods*; see the "pure virtual methods" part of the
 * protected interface below]. Hence, instances of BMinQuad cannot be built,
 * and a derived class has to be defined where such methods are actually
 * implemented. */

class BMinQuad : public MinQuad {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--									  --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--									  --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor: initialize the object, but, as in the base MinQuad class, it
 * does not allocate memory, as this is done is SetMaxDim() [see below]. */

 BMinQuad( void );

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** Same meaning as in the base MinQuad class; here, however, the third
 * parameter is *not* optional, and 0 is *not* a valid value. */

 void SetMaxDim( Index m , Index n , Index SDim );

/*--------------------------------------------------------------------------*/
/** In the code, a number eD is used as a measure of the relative precision
 * required for the box constraints: l[ i ] <= d[ i ] is satisfied if
 *
 *     d[ i ] >= l[ i ] - t * eD * size( l[ i ] )
 *
 * and analogously for d[ i ] <= u[ i ].
 *
 * eD is automatically increased by SolveQP() to try to face badly conditioned
 * problems: this method allow to set its current value. Passing 0 (clearly,
 * an impossible value) to SetEpsilonD() resets eD to some "default" value. */

 void SetEpsilonD( HpNum NeweD = 0 );

/*--------------------------------------------------------------------------*/
/** At each iteration, variables are added to the "active set", i.e. the set
 * of variables that are fixed to their lower[/upper] bound. This method sets
 * the maximum number of variables that can be added to the active set in one
 * iteration: usually a "large" value helps in avoiding many unnecessary
 * iterations, but this is in principle instance-dependent.
 *
 * If the method is not called, "Inf< Index >()" is assumed, i.e., "no limit".
 */

 void SetMaxVarAdd( cIndex MVA = 1 ) { MaxVarAdd = MVA; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** At each iteration, variables removed from to the "active set", i.e. the
 * set of variables that are fixed to their lower[/upper] bound. This method
 * sets the maximum number of variables that can be removed from the active
 * set in one iteration: usually a "large" value helps in avoiding many
 * unnecessary iterations, but this is in principle instance-dependent.
 *
 * If the method is not called, "Inf< Index >()" is assumed, i.e., "no limit".
 */

 void SetMaxVarRmv( cIndex MVR = 1 ) { MaxVarRmv = MVR; }

/*--------------------------------------------------------------------------*/

#if( LOG_BMQ )
/** The output of the code is directed onto the ostream object pointed by log:
 * if this method is never called, 'clog' (the standard error bufferized) is
 * used as default. Under Unix-like environments, it can be redirected to a
 * file by using ">&" from the command shell. */

 void SetBMQLog( ostream *log = 0 )
 {
  BMQLog = log;

  #if( LOG_BMQ > 2 ) 
   *BMQLog << "BMinQuad: " << SpaceDim << " variables." << endl << endl;
  #endif
  }

#endif

/*--------------------------------------------------------------------------*/
/** SetBMQTime() decides whether or not an OPTtimers object is used to compute
 * the solution time of the BTT algorithm.
 *
 * Note that time accumulates over the calls; calling SetBMQTime(), however,
 * resets the counters, allowing to time specific groups of calls. */

 virtual void SetBMQTime( const bool TimeIt = true )
 {
  if( TimeIt )
   if( BMQt )
    BMQt->ReSet();
   else
    BMQt = new OPTtimers();
  else {
   delete BMQt;
   BMQt = 0;
   }
  }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*--------------------------------------------------------------------------*/

 void AddSubGrad( cIndex n , cHpNum alfan );

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void AddConstr( cIndex n , cHpNum alfan );

/*--------------------------------------------------------------------------*/

 void ChangeAlfa( cIndex i , cHpNum DeltaAlfai )
 {
  if( DeltaAlfai ) {
   RealAlfa[ i ] += DeltaAlfai;
   MinQuad::ChangeAlfa( i , DeltaAlfai );
   Bf = Inf< HpNum >();  // the objective function is allowed to increase
   }
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void ChangeAlfa( cHpNum DeltaAlfa );

/*--------------------------------------------------------------------------*/
/** Has the same effect than that of the base class, but also the bounds on
 * the variables are updated accordingly. Note that steps Tau > ti lead to
 * an *unfeasible* solution if ActiveVars()[] is nonempty (and see the
 * comments to MinQuad::MoveAlongD() for linear constraints).
 *
 * Important note: to update the bounds on the constrained variable named `i',
 * the i-th component of the latest direction d is needed; this is clearly
 * available if i was already constrained at the time when SolveQP() was last
 * invoked, but it is *not* if i has been created in the meantime with
 * AddVars() or MakeVarCnstr( i ).
 *
 * Therefore, in order to make MoveAlongD() work properly, it needs
 *
 * - either that no new constrained variables be added between the end of
 *   SolveQP() and the call to MoveAlongD();
 *
 * - or that SetD( i , ... ) is invoked after the "creation" of the i-th
 *   variable to provide the required information. */

 void MoveAlongD( cHpNum Tau , cHpNum DeltaFi );

/*--------------------------------------------------------------------------*/

 void ReadAlfa( HpRow NewAlfa );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 cHpRow ReadAlfa( void ) { return( RealAlfa ); }

/*--------------------------------------------------------------------------*/

 void SetGTz( cIndex i , cHpNum GTzi );

/*--------------------------------------------------------------------------*/
/** Used to provide the i-th component of the latest direction d (without the
 * (-t) factor) relative to the newely "created" constrained variable i,
 * prior to a call to MoveAlongD(). */

 void SetD( Index i , LMNum Di ) { di[ i ] = Di; }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Like SetD( i , Di ), but the second returns a pointer to a vector such
 * that SetD( i , Di ) is equivalent to SetD()[ i ] = Di. */

 LMRow SetD( void ) { return( di ); }

/*--------------------------------------------------------------------------*/
/** Access the lower bounds on the variables. Return a pointer to a vector
 * containing the current values of the bounds, that changes at each call to
 * MoveAlongD(): LowerBounds()[ i ] holds the lower bound of the variable with
 * name `i'. */

 LMRow LowerBounds( void )
 {
  #if( TWOSIDED )
   return( lb );
  #else
   return( bounds );
  #endif
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 #if( TWOSIDED )

/** If TWOSIDED != 0, access the upper bounds on the variables. Return a
 * pointer to a vector containing the current values of the bounds, that
 * changes at each call to MoveAlongD(): UpperBounds()[ i ] holds the upper
 * bound of the variable with name `i'. */

 LMRow UpperBounds( void )  { return( ub ); }

#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Change the lower [and upper, if any] bounds on the variables. The new
 * values must first be written in the appropriate entries of
 * Lower[Upper]Bounds() (that does *not* return read-only pointers just
 * because of that), and then ChangeBounds() must be called. The same
 * technique is used for passing the bounds to InitialSetUp() and AddVars().
 *
 * Since it is assuned that d = 0 should always be a feasible solution, it is
 * required that l[ i ] <= 0 <= u[ i ] for each constrained variable `i'.
 *
 * If TWOSIDED > 0, it makes sense to allow "mixtures" of variables with only
 * one (upper or lower) bound and variables with two bounds; hence, l[ i ]
 * can be set to - LMINF and u[ i ] can be set to + LMINF, that will be
 * recognized and appropriately handled. However, it is required that at
 * least one of the bounds be finite; this is not restrictive, since a
 * variable i with both infinite bounds is simply unconstrained, and it
 * should be explicitly declared such [see InitialSetUp(), AddVars() and
 * MakeVarUnCnstr()]. Consequently, if TWOSIDED == 0 it is *not* allowed to
 * pass - LMINF as a lower bound for a constrained variable. */

 void ChangeBounds( void );

/*--------------------------------------------------------------------------*/
/** InitialSetUp allows some initializations to be performed more efficiently,
 * and some user-defined knowledge to be passed to the solver. The method
 * reads the first SDim entries of the vector of chars returned by GetVars()
 * and sets the variables according to the values found there.
 *
 * The type of a variable, that is coded in the bits of the char: if
 * vi = GetVars()[ i ], one has
 *
 * - vi & NNVar()  is nonzero <=> `i' is constrained to be nonnegative (this
 *                 is guaranteed to be coded in the first bit);
 *
 * - vi & IsVar()  is nonzero <=> `i' is defined in (P);
 *
 * - vi & AcVar()  is nonzero <=> `i' is a nonnegative variable that is set
 *                 to one of its bounds in the optimal solution of (P);
 *		   if TWOSIDED == 0, the bound is clearly the lower one;
 *
 * - vi & UBVar()  is nonzero <=> `i' is a nonnegative variable that is set
 *                 to its upper bound in the optimal solution of (P)
 *		   (this should be used only if TWOSIDED > 0).
 *
 * For constrained variables, the values written there prior to a call to
 * InitialSetUp() are meant to provide a guess of the variables that will be
 * fixed at their lower[upper] bounds in the optimal solution of the problem:
 * a good guess can considerably decrease the time required to solve the
 * problem. Hence set GetVars()[ i ] to
 *
 * NNVar() | IsVar() | AcVar() | UBVar() if the entry d[ i ] of the optimal
 * solution (direction) corresponding to variable `i' is probably == u[ i ];
 *
 * NNVar() | IsVar() | AcVar() if d[ i ] == l[ i ] (probably);
 *
 * NNVar() | IsVar() if l[ i ] < d[ i ] < u[ i ] (probably);
 *
 * IsVar() if `i' is declared in (P) as unconstrained;
 *
 * NNVar() if `i' is *not* declared in (P), but it will be constrained;
 *
 * 0 if `i' is *not* declared in (P), but it will be unconstrained.
 *
 * The bit-wise coding allow to set each field separately; for instance,
 * vi |= IsVar() sets the variable as declared whatever its "sign" be, while
 * vi &= ~IsVar() sets the variable as *not* declared.
 *
 * The default (i.e. if nothing is written in vi) is that the variable is
 * *not* declared but it is constrained. Note that the "basic type"
 * (constrained or unconstrained) of a variable is maintained even when the
 * variable is removed [see RemoveVars()]; hence, if `i' has been
 * unconstrained previously and it has been removed, GetVar( i ) & NNVar()
 * will be zero instead.
 *
 * If the method is called prior to the first iteration, there are not
 * general linear constraints in the bundle but exactly one subgradient is
 * available, the optimal solution will just be - ti * < the subgradient >
 * "projected" over the box constraints, hence the above sets can be exactly
 * guessed (e.g. - ti * subg[ i ] <= l[ i ] ==> AcVar()).
 *
 * Note that a call to InitialSetUp() will overwrite any previous setting,
 * making *all and only* these variables declared. The values of the bounds
 * for the constrained variables must be provided as in ChangeBounds().
 *
 * If the bundle is empty, the method is very cheap; if there are items, the
 * whole matrix Q of the problem must be recomputed. This is done by invoking
 * ChangeQ() [see MinQuad.h], which in turn *may* call GiTG[j]() [see] for
 * calculating the entries. Hence, it is better that this method be called
 * *before* inserting items in the bundle, if possible. */

 void InitialSetUp( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// mask for variables with lower bound, used in InitialSetUp()

 char NNVar( void ) { return( 1 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// mask for existing variables, used in InitialSetUp()

 char IsVar( void ) { return( 2 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// mask for active variables, used in InitialSetUp()

 char AcVar( void ) { return( 4 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( TWOSIDED )
 /// mask for variables with uppper bound, used in InitialSetUp()

 char UBVar( void ) { return( 12 ); }

#endif

/*--------------------------------------------------------------------------*/
/** Adds `hwmny' new variables to the set of declared variables, those whose
 * names are found in the set whch[], that has to be ordered in increasing
 * sense, without replicated elements and Inf< Index >()-terminated (that is,
 * whch[ hwmny ] == Inf< Index >()). None of the variable can be declared just
 * prior to the call of the metho.
 *
 * The type of the variable must be written in GetVars()[ i ] prior to the
 * call: see InitialSetUp() above for how to choose the value.
 *
 * Note that the "basic type" (constrained or unconstrained) of a variable is
 * maintained even when the variable is removed [see RemoveVars()]; hence, if
 * `i' has been constrained previously, it has been removed and it is created
 * again, it will be taken as constrained unless GetVars()[ i ] is explicitly
 * changed (and the same holds for unconstrained variables).
 *
 * If `i' is a constrained variable, the corresponding lower[upper] bound
 * must be written in Lower[Upper]Bounds()[ i ] prior to the calls as for
 * ChangeBounds(). */

 void AddVars( cIndex_Set whch , cIndex hwmny );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Like AddVars( whch , hwmny ), but adds all those variable with names
 * from `strt' to `strt + hwmny - 1', extremes included, and it can only be
 * called if strt is *larger* than the name of any declared variable. Thus,
 * none of the variable can be declared just prior to the call of the method.
 */

 void AddVars( cIndex strt , cIndex hwmny );

/*--------------------------------------------------------------------------*/
/** Deletes the variable with name `i' from the set of declared variables,
 * regardless to its type (constrained or unconstrained). */

 void RemoveVars( cIndex_Set whch , Index hwmny );


/*--------------------------------------------------------------------------*/
/** Renames the variable with name `i', giving it name `j', but maintaining the
 * same type (constrained or unconstrained), bounds and so on. There must
 * *not* be already a variable with name `j' among the declared ones.
 *
 * In the typical use, `i' is the last variable (the one with largest name);
 * if this is the case, setting iIsLst == true (the default) allows some
 * operations to be performed more efficiently. */

 void MoveVar( cIndex i , cIndex j , const bool iIsLst = true );

/*--------------------------------------------------------------------------*/
/** Renames all the variables to make the set of variable names fit with the
 * elimination of the variables whose names are contained in the set `whch'
 * (ordered in increasing sense, without duplications and
 * Inf< Index >()-terminated).
 *
 * The removal of a variable with RemoveVar() from the set of declared
 * variables can be thought to be a "temporary" operation; in fact, it is
 * arranged in such a way that re-declaring the variable with AddVar() is
 * very simple (for instance, the type of the variable is kept).
 *
 * A more "permanent" removal of a variable can be required, where all the
 * information about the variable is lost. This is what happens to the
 * single variable `j' in MoveVar() [see above]. This method provides means
 * for eliminating a (large) set of variables all at once.
 *
 * The set `whch' must contain names of *undeclared* variables (that is,
 * variables that have never been declared with InitalSetUp() or AddVar(),
 * or finally undeclared with either RemoveVar() or MoveVar()). All the
 * variables, both declared and undeclared ones, are renamed as to fit with
 * the elimination of those variables. For instance, if SDim == 10 [see
 * SetMaxDim()] and whch = { 2, 3, 7 }, after a call to RenameVars() one has
 *
 *    current name  0 1 2 3 4 5 6 7 8 9
 *    previous name 0 1 4 5 6 8 9 n n n
 *
 * where "n" means "this is a new variable, completely unrelated to the
 * previous ones, and currently undeclared". */

 void RenameVars( cIndex_Set whch );

/*--------------------------------------------------------------------------*/
/** Converts the unconstrained variable `i' into a constrained variable,
 * reading its "status" and its lower [and upper] bound like AddVar() does.
 * It is uncorrect to call this method with an i not being the name of a
 * declared variable, and also if `i' is already a constrained variable. */

 void MakeVarCnstr( cIndex i );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Converts the constrained variable `i' into an unconstrained variable. It
 * is uncorrect to call this method with an i not being the name of a
 * declared variable. */

 void MakeVarUnCnstr( cIndex i );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/
/** Solves the problem with the current data, i.e. the current bundle, Alfa,
 * ti and lower[/upper] bounds. The second level of dynamic tolerances
 * adjustment is implemented here inside: hence, eD [see SetEpsilonD() above]
 * can be changed to face increase "requests".
 *
 * Exit codes:
 *
 * - kOK , kQPPrimUnbndd = as in the base class;
 *
 * - kFatal = either MinQuad::SolveQP() returned a kFatal, or even adjusting
 *            eD was not sufficient to handle the problems: that eD has been
 *            increased out of the bound. */

 MQError SolveQP( HpNum ti );

/*--------------------------------------------------------------------------*/

#if( CNDVD_TMP )
/** If CNDVD_TMP > 0, a temporary vector of SDim LMNum's must be provided by
 * the calling code *prior* to a call to SolveQP() by calling SettmpD(). The
 * temporary may become "property" of BMinQuad, but if it happens then the
 * same amount of memory becomes available: hence, if the calling code wants
 * to reclaim that memory, it must require its address by calling GettmpD()
 * *after* SolveQP(). */

 void SettmpD( LMRow td ) { tmpdi = td; }

 LMRow GettmpD( void ) { return( tmpdi ); }

#endif

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR READING RESULTS ---------------------*/
/*--------------------------------------------------------------------------*/
/** Same meaning as in the base class, i.e.,
 *
 *      f() = (1/2) * ti * ReadzNorm() + ReadSigma(),
 *
 * but it is important to remark that now
 *
 *      ReadzNorm() = || d^* / ti ||_2^2 != || z^* ||_2^2.
 *
 * That is, ReadzNorm() reports the norm of the *scaled* direction d^*, that
 * is no longer equal to the "aggregate subgradient" z^*.
 *
 * Important Note: if BEXACT == 0, they may give wrong results if called after
 * MoveAlongD(), ChangeBounds(), RemoveVar() or MakeVarUnCnstr(). */

 HpNum ReadzNorm( void ) { return( Quad + bNorm / ( PrvsTi * PrvsTi ) ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Same meaning as in the base class, i.e.,
 *
 *      f() = (1/2) * ti * ReadzNorm() + ReadSigma(),
 *
 * but it is important to remark that now ReadSigma() contains an extra term
 * concerning the bounds of the "active" variables; since the "active
 * variables" corresponds to constraints in Base[], their contribution can be
 * eliminated (together with that of ordinary constraints) by setting
 * IncldCnst == false. 
 *
 * Important Note: if BEXACT == 0, it may give wrong results if called after
 *  MoveAlongD(), ChangeBounds(), RemoveVar() or MakeVarUnCnstr(). */

 HpNum ReadSigma( const bool IncldCnst = true )
 {
  if( IncldCnst || ( ! ReadCBDim() ) )
   return( Lin - bNorm / PrvsTi );
  else {
   HpNum tS = 0;
   cHpRow tM = Mult;
   cIndex_Set tB = Base;
   for( Index h ; ( h = *(tB++) ) < Inf< Index >() ; tM++ )
    if( ! IsAConst( h ) )
     tS += (*tM) * RealAlfa[ h ];

   return( tS );
   }
  }

/*--------------------------------------------------------------------------*/
/** All the variables "live" in the space of components [ 0 .. SDim - 1 ], so
 * that any object in the primal space (the direction d or a subgradient) is
 * naturally viewed as a vector (of LMnum's) with SDim components. This method
 * return a read-only pointer the "aggregated subgradient" z^*, that is for
 * each declared variable `i' one has
 *
 *     z[ i ] = Sum{ h in Base } G[ Base[ h ] ][ i ] * Mult[ h ].
 *
 * Note that z^* may not be, strictly speaking, an "aggregate subgradient", as
 * constraints might be in Base: however, the constraints are subgradients of
 * the characteristic function of the feasible set.
 *
 * Note that the content of the entries of z[] *not* corresponding on a
 * *declared variable* the depends on LAZY_D and the (user-provided)
 * implementation of CalculateZ(). If LAZY_D > 0 then again only the entries
 * corresponding to constrained variables are significant, but if
 * LAZY_D == 0 then the returned pointer points to the same vector that has
 * been passed to CalculateZ(), and that has not been modified ever since.
 * Hence, if some sort of "correct values" have been written there, they will
 * still be there. */

 cLMRow ReadZ( void ) { return( di ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Like ReadZ( void ), but writes z^* into z. z must be at least as long as
 * the maximum number i that has been used as a name for a variable. Note
 * that z[ i ] is written *if and only if* `i' is the name of a *declared
 * variable*, all the other entries being left untouched. */

 void ReadZ( LMRow z );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Like ReadZ( LMRow ), but gives z in a "sparse" format, i.e. it only uses
 * the first NV positions of z, where NV is the total number of declared
 * variables. The variables are ordered in ascending sense by their name, and
 * the h-th variable (h = 1 .. NV - 1) is written in z[ h ]. The names and
 * the number of variables should be known by the calling program, but they
 * can be queried by a call to ReadVNames(). */

 void ReadZSprs( LMRow z );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Writes the (ordered) names of declared variables in the first NV
 * components of VNames, and returns NV. */

 Index ReadVNames( Index_Set VNames );

/*--------------------------------------------------------------------------*/
/** Analogous ReadZ( LMRow ) for the optimal primal solution d; that is,
 * d[ i ] = min( u[ i ] , max( l[ i ] , - ti * z[ i ] ) ) where z[] is the
 * vector reported by ReadZ[Sprs]().
 *
 * The parameter CpyFrst that all the first CpyFrst entries of the vector d
 * have to be written, even those corresponding to non-declared variables. If
 * CpyFrst > 0, it is assumed that "correct" values for z have been provided
 * during the calculation [see ReadZ()] and they are available in the data
 * structure: - ti * z[ i ] is then written in d[ i ] for all non-declared
 * variables i < CpyFrst. Otherwise, only the entries of d[] corresponding to
 * declared variables are written. */

 void ReadD( LMRow d , cIndex CpyFrst = 0 );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Analogous ReadZSprs( LMRow ) for the optimal primal solution d; that is,
 * d[ i ] = min( u[ i ] , max( l[ i ] , - ti * z[ i ] ) ) where z[] is the
 * vector reported by ReadZSprs(). */

 void ReadDSprs( LMRow d );

/*--------------------------------------------------------------------------*/
/** Returns the set of "active" variables, i.e., where d[ i ] is set to the
 * bound that is violated by - ti * z[ i ]. Such a set is described by the
 * protected field Base2[]; [see the methods in the pure virtual protected
 * interface]. However, to let this information available to non-derived
 * classes, it is also returned by ActiveVars(). */

 cIndex_Set ActiveVars( void ) { return( Base2 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Returns the size of the set of "active" variables, see ActiveVars().

 Index AVDim( void ) { return( B2Dim ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the set of "inactive" variables, i.e., where
 * d[ i ] = - ti * z[ i ]. Such a set is described by the protected field
 * MBase2[]; [see the methods in the pure virtual protected interface].
 * However, to let this information available to non-derived classes, it is
 * also returned by InActiveVars(). */

 cIndex_Set InActiveVars( void ) { return( MBase2 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Returns the size of the set of "inactive" variables, see InActiveVars().

 Index IAVDim( void ) { return( MB2Dim ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Describes the status of each variable;  see InitialSetUp() above for how
 * to decode the returned value.
 *
 * Note that the "basic type" (constrained or unconstrained) of a variable is
 * maintained even when the variable is removed [see RemoveVar() below];
 * hence, GetVar( i ) & NNVar() is nonzero for a variable that is not
 * currently defined if and only if it has been a constrained variable the
 * last time that it has been destroyed. By default, the variables are all
 * constrained initially (even those that are not defined). */

 char GetVar( cIndex i ) { return( GS[ i ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** GetVars()[ i ] == GetVar( i ); the returned pointer is *not* a read-only
 * one, since writing in that memory area is allowed prior to a call to
 * InitialSetUp(), AddVar() or MakeVarCnstr() for passing a "guess" on the
 * initial statuses of the variables. */

 char * GetVars( void ) { return( GS ); }

/*--------------------------------------------------------------------------*/
/** Returns
 *
 *   - ( 1 / ti ) * ( d[] * g[] )
 *
 * i.e. the (scaled) scalar product between the direction d[] and the vector
 * g[]. g[] is assumed to be in the "naive" format, i.e. a vector of SgNum
 * such that g[ i ] is the entry of the variable whose "name" is i.
 * The "- (1 / ti)" scaling factor is used for this result to be immediately
 * passed to the SetGTz() method of the base class [see MinQuad.h and
 *ReadGTz()]. */

 HpNum DPerG( cSgRow g );

/*--------------------------------------------------------------------------*/
/** Set L1 = L2 - ( Tau / ti ) * d^*, a typical "move" in the primal space.
 *
 * The intended "format" of L2 and L1 is the same of Read[Z/D](); that is, L2
 * and L1 are considered at least as long as the maximum current name of a
 * variable, and the result relative to variable `i' is taken from L2[ i ]
 * and written to L1[ i ] (the components not corresponding to declared
 * variables are left untouched).
 *
 * The method recognises the case Tau = t and deals with it explicitly. */

 void AddD( LMRow L1 , cLMRow L2 , cHpNum Tau );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Set L1 = L2 - ( Tau / ti ) * d^*, a typical "move" in the primal space.
 *
 * The intended "format" of L2 and L1 is the same of Read[Z/D]Sprs(); that
 * only the first NV components of L2 and L1 are used.
 *
 * The method recognises the case Tau = t and deals with it explicitly. */

 void AddDSprs( LMRow L1 , cLMRow L2 , cHpNum Tau );

/*--------------------------------------------------------------------------*/
/** Same meaning as in the base class, except that here ReadGTz( i ) returns
 *
 *     - ( d^* / ti ) * G[ i ]
 *
 * i.e., DPerG( G[ i ] ). This is the same as z^* * G[ i ] if there are no
 * "active" variables. */

 HpNum ReadGTz( cIndex i );

/*--------------------------------------------------------------------------*/

 void SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 );

/*--------------------------------------------------------------------------*/
/// Returns the user and sistem time in seconds spent so far by SolveQP()

 void BMQTime( double &t_us , double &t_ss )
 {
  t_us = t_ss = 0;
  if( BMQt )
   BMQt->Read( t_us , t_ss ); 
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Returns the total time in seconds spent so far by SolveQP()

 double BMQTime( void ) { return( BMQt ? BMQt->Read() : 0 ); }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
//// Returns the current value of eD, see SetEpsilonD() above.

 HpNum EpsilonD( void ) { return( eD ); }


/*----------------------------------------------------------------------------

   IMPORTANT NOTE: the methods of the base class MinQuad

   MinQuad::LowQ( i )   and   MinQuad::LowQ( i , j )

   may no more be reliable now to obtain the "full" scalar product
   G[ i ] * G[ j ], since the matrix is modified during the course of the
   algorithm; in fact, LowQ( i , j ) == GiTGj( i , j ) as defined above,
   i.e., the scalar product calculated on a subset only (given by MBase2) of
   the variables.

----------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Memory deallocation. Statistics (if any) are printed.

 virtual ~BMinQuad();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PURE VIRTUAL METHODS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--									  --*/
/*--  The following methods are a part of the protected interface of the  --*/
/*--     class, but they *must* be implemented in the derived classes.    --*/
/*--									  --*/
/*--------------------------------------------------------------------------*/
/** These are the methods used by the base class MinQuad to access the data
 * of the problem. However, within BMinQuad there is an important difference:
 * the scalar products must *not* be performed on *all* the entries of the
 * items, but only on the subset specified by the protected fields MBase2
 * (Index_Set) and MB2Dim (Index), i.e.
 *
 * GiTGj( i , j ) =
 *  Sum{ p = 0 .. MB2Dim - 1 } G[ i ][ MBase2[ p ] ] * G[ j ][ MBase2[ p ] ].
 *
 * In the simplest case (there exist a "naive" form of the items as vectors
 * of SgNum's with SDim components, the variable with "name" `i' corresponds
 * to the i-th entry of those vectors and the item whose "name" is 'h' is the
 * h-th row of the matrix G) the above scalar product can be calculated with
 *
 *     ScalarProduct( MBase2 , G[ i ] , G[ j ] )
 *
 * Note that MBase2 will contain (the names of) a possibly proper subset of
 * the constrained variables, but *all the unconstrained variables*.
 *
 * Hence, assume that an item is a k-array of SgNum's, with k possibly > SDim.
 * Let C be the set of "declared" variables (a subset of { 0 .. SDim - 1 },
 * containing names of both constrained and unconstrained variables) and U be
 * the set of "undeclared" variables (a subset of { SDim .. k - 1 },
 * containing names of only unconstrained variables). Then, a correct GiTGj()
 * would return
 *
 *    ScalarProduct( MBase2 , G[ i ] , G[ j ] ) +
 *         ScalarProduct( U , G[ i ] , G[ j ] )
 *
 * Hence, the objects of class BMinQuad (as their ancestors MinQuad) need not
 * know anything about the "real" implementation of the subgradients: the
 * implementation of GiTGj( i , j ) is *not* provided, and it's due to the
 * derived class, that can access to MBase2 and MB2Dim. The rationale is (as
 * usual) that the items may have some structure (e.g. sparsity) that can be
 * exploited to fasten the calculation of the scalar product.
 *
 * GiTG() is exactly as in the base class, i.e. the pseudo-code
 *
 *   for( j = 0 ; j < iMax ; j++ )
 *    if( IsThere( j ) )
 *     Qi[ j ] = GiTGj( i , j );
 *
 * is still a valid representation of what GiTG( i , Qi ) is required to do
 * provided that "GiTGj()" accomplish the right task.
 *
 * The implementing code can rely on the following invariants, guaranteed by
 * the BMinQuad and MinQuad classes:
 *
 * - MBase2 is ordered, i.e. i < j => MBase2[ i ] < MBase2[ j ];
 *
 * - each element is unique;
 *
 * - MBase2 is "infinity terminated", i.e. MBase2[ MB2Dim ] == Inf< Index >();
 *
 * - GiTGj( i , j ) is always called with i <= j. */

#if( LAZY_Q )

 virtual QuNum GiTGj( cIndex i , cIndex j ) = 0;

#else

 virtual void GiTG( cIndex i , QuRow Qi , cIndex iMax ) = 0;

#endif


/*--------------------------------------------------------------------------*/
/** In this case, the scalar products does not give "enough information" to
 * solve the problem: it is also necessary to access to the matrix G of the
 * subgradients, in a row-wise fashion. This method must return a read-only
 * pointer RG to the i-th row of G, i.e.
 *
 *             / G[ h ][ i ] if h is the "name" of a subgradient in the bundle
 *   RG[ h ] = |
 *             \ anything    otherwise
 *
 * In the simplest case (there exist a "naive" form of items as vectors of
 * SgNum's with SDim components, the variable with "name" `i' corresponds to
 * the i-th entry of those vectors and the item whose "name" is 'h' is the
 * h-th row of the matrix G) the correct RG[] is
 *
 *   forall( < h = index of an item currently in the bundle > )
 *    RG[ h ] = G[ h ][ i ];  */

 virtual cSgRow GiTilde( cIndex i ) = 0;

/*--------------------------------------------------------------------------*/
/** Dealing with constraints on the direction d requires at least checking if
 * they are violated; since d[ i ] = - ti * z[ i ] wherever this does not
 * violate a constraint, the entries of the "aggregate item" z must be
 * computed.
 *
 * In the usual (simple) case where there exist a "naive" form of the items
 * as vectors of SgNum's with SDim components, all of which corresponding to
 * constrained variables, and the subgradient whose "name" is `i' is just the
 * i-th row of the matrix G, CalculateZ( < all > ) (the third form) must do
 *
 *    for( i = 0 ; i < SDim ; i++ )
 *     for( z[ i ] = 0 , h = 0 ; h < BDim ; h++ )
 *      z[ i ] += G[ Base[ h ] ][ i ] * Mult[ i ];
 *
 * where Mult (HpRow), Base (Index_Set) and BDim (Index) are protected fields
 * of class BMinQuad (actually, of base class MinQuad). Note that *the memory
 * of z is provided by the BMinQuad object*. To ease the calculation, Base[]
 * is guaranteed to be Inf< Index >()-terminated, i.e.,
 * Base[ BDim ] == Inf< Index >().
 *
 * If the components of z are requested one by one (the first form), the
 * intended semantic of CalculateZ( < one > ) is that of the following code
 *
 *    cLMRow temp = GiTilde( h );
 *    LMNum zh = 0;
 *
 *    for( Index i = 0 ; i < BDim ; i++ )
 *     zh += temp[ Base[ i ] ] * Mult[ i ];
 *
 *    return( zh );
 *
 * but smarter implementations may be possible in some circumstances.
 *
 * Finally, if the components of z are requested "in slots" (the second form),
 * intended semantic of CalculateZ( < a set > ) is that of
 *
 *   for( Index h = 0 ; Wh[ h ] < Inf< Index >() ; h++ )
 *    z[ Wh[ h ] ] = CalculateZ( Wh[ h ] );
 *
 * Which of these three implementations is better depends on the details of
 * how the subgradients are kept in memory (by rows or by columns, in dense
 * or sparse format) and on the instances to be solved (average number of
 * items vs number of variables, presence of unconstrained variables), so
 * that the final decision should be due to the final user. */

#if( LAZY_D == 2 )

 virtual LMNum CalculateZ( cIndex h ) = 0;

#elif( LAZY_D == 1 )

 virtual void CalculateZ( cIndex_Set Wh , LMRow z ) = 0;

#else

 virtual void CalculateZ( LMRow z ) = 0;

#endif

/*--------------------------------------------------------------------------*/
/** This method must compute the scalar product between the item whose "name"
 * is `i' / all the items and the vector l, limited to the entries whose
 * indices are in the lBd-vector lB (ordered in increasing sense and
 * Inf< Index >()-terminated).
 *
 * As usual, this is better visualized by a fragment of code; for the first
 * form this is
 *
 *    LMNum res = 0;
 *    for( j = 0 ; j < lBd ; j++ )
 *     res += G[ i ][ lB[ j ] ] * l[ lB[ j ] ];
 *
 *    return( res );
 *
 * lB[] is always a set of names of declared *constrained* variables, while
 * l is the vector of the "active" bounds.
 *
 * GiTLB() is called inside Add[SubGrad/Constr]() and ChangeBounds(), so be
 * sure that all the data structures are updated before invoking those
 * methods. */


 virtual LMNum GiTLB( cIndex i , cLMRow l , cIndex_Set lB , cIndex lBd ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Like GiTLB( i , ... ), but the code is
 *
 *    for( j = 0 ; j < iMax ; j++ )
 *     if( < an item of name "j" is in the bundle > )
 *      if( add )
 *       gtlb[ j ] += GiTLB( j , l , lB , lBd );
 *      else
 *       gtlb[ j ] -= GiTLB( j , l , lB , lBd );
 *
 * (but, of course, depending on the implementation different and more
 * efficient ways for modifying gtlb[] may exist). */

 virtual void GiTLB( HpRow gtlb , cLMRow l , cIndex_Set lB , cIndex lBd ,
		     const bool add ) = 0;

/*--------------------------------------------------------------------------*/
/*---------------------- HOOKS FOR DERIVED CLASSES -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods are called by the solver to give the user a   --*/
/*--  better control over the optimization process.                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#if( SIGNAL_MBCHG )
 /// This method is called each time (Mult, Base, BDim) change

 virtual void MBHasChgd( void ) = 0;

#endif

#if( SIGNAL_B2CHG )
 /// This method is called each time (MBase2, MB2Dim, Base2, B2Dim) change

 virtual void B2HasChgd( void ) = 0;

#endif

/*--------------------------------------------------------------------------*/
/*------------------ "REAL" PROTECTED PART OF THE CLASS --------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** In case some data of the problem changes, RealAlfa has to be recomputed
 * from scratch (as opposed to being smartly updated); this method does just
 * that. */

 void RecomputeRealAlfa( void );

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

 Index_Set Base2;       /**< Set of names of defined constrained variables at
			 * their LB or UB in the optimal solution */
 Index B2Dim;           ///< Dimension of Base2

 Index_Set MBase2;      /**< Set of names of defined constrained variables
			 * that are *not* in Base2 or defined unconstrained
			 * variables */
 Index MB2Dim;          ///< Dimension of MBase2 ( <= SpaceDim - B2Dim )

 Index_Set MvdVars;     /**< Set of names of defined variables that are being
                         * moved from Base2 to MBase2 or vice-versa */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--									  --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--									  --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void CalcOptDir( HpNum ti );

/*--------------------------------------------------------------------------*/

  void AddToB2( cIndex_Set MVD , cIndex MVDd );

  void RmvFrmB2( cIndex_Set MVD , cIndex MVDd );

  void AddToMB2( cIndex_Set MVD , cIndex MVDd );

  void RmvFrmMB2( cIndex_Set MVD , cIndex MVDd );

/*--------------------------------------------------------------------------*/

  void CutOffConstrs( cIndex_Set MVD , Index MVDd );

/*--------------------------------------------------------------------------*/

  void PutInConstrs( cIndex_Set MVD , Index MVDd );

/*--------------------------------------------------------------------------*/

  void ClearTabooList( void );

/*--------------------------------------------------------------------------*/

  void MemDealloc( void );

/*--------------------------------------------------------------------------*/

  void CheckB2( void );

  void CheckMB2( void );

  void CheckGS( void );

  void CheckRA( void );

  void CheckDS( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  HpNum eD;           // Relative error for constraints violation tests
  HpNum Bf;           // Objective function value
  HpNum bNorm;        // Norm( bounds , Base2 , B2Dim );

  Index MaxVarAdd;    // How many variables at most can be added ..
  Index MaxVarRmv;    // .. and removed from Base2[] at each iteration

  Index TLDim;        // total number of "taboo" items
  Index NNStop;       // 1 + name of the last NN variable in the bundle

  HpRow RealAlfa;     // The "real" value of Alfa, while Alfa[ i ] contains
                      // RealAlfa[ i ] - G[ i ]{Xsi} * l{Xsi} (or u{Xsi})
		      // where {Xsi} is the set of the active box constr.

  Index SpaceDim;     // Max. number of variables
  char *GS;           // Status of each variable

  LMRow bounds;       // The [active] bounds
  #if( TWOSIDED )
   LMRow lb;          // Lower Bounds
   LMRow ub;          // Upper Bounds
  #endif

  LMRow di;           // The current primal solution
  LMRow tmpdi;        // Temporary for primal solution (with "exact" pricing)

  #if( LOG_BMQ )
   std::ostream *BMQLog;       // the output stream object

   #if( LOG_BMQ > 1 )
    unsigned long int BCalls;  // Calls counter
    unsigned long int BSccss;  // Succesfull calls counter
    float SumAverages;         // Sum{ all the steps i } B2Dim{i}
   #endif
  #endif

  OPTtimers *BMQt;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class BMinQuad )

/*--------------------------------------------------------------------------*/

 }  // end( namespace MinQuad_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BMinQuad.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BMinQuad.h -------------------------------*/
/*--------------------------------------------------------------------------*/
