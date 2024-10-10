/*--------------------------------------------------------------------------*/
/*---------------------------- File Bundle.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the Bundle class, which implements the NDOSolver interface
 * for NonDifferentiable Optimization Solvers, as described NDOSlver.h,
 * using a "Generalized Bundle" algorithm.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni "Generalized Bundle Methods"
 *  SIAM Journal on Optimization 13(1), p. 117 - 156, 2002
 *
 * available at
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#SIOPT02
 * \endlink
 *
 * The class requires that the function to be minimized be available under
 * the form of a FiOracle object, as described in FiOracle.h.
 *
 * The class is parametric over the type of Master Problem used: it just
 * relies over an object of class MPSolver to solve it, see MPSolver.h.
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

#ifndef _Bundle
 #define _Bundle  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Bundle_MACROS Compile-time switches in Bundle.h
    @{ */

/*------------------------------- LOG_BND ----------------------------------*/

#define LOG_BND 1

/* If LOG_BND > 0, the Bundle class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetBLog() [see below]. */

/** @} end( group( Bundle_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NDOSlver.h"

#include "MPSolver.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Bundle_CLASSES Classes in Bundle.h
    @{ */

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The Bundle class implements the NDOSolver interface for NonDifferentiable
    Optimization Solvers [see NDOSlver.h], using a "Generalized Bundle"
    algorithm as described in

     A. Frangioni "Generalized Bundle Methods"
     SIAM Journal on Optimization 13(1), p. 117 - 156, 2002

    This is in fact a "meta" algorithm, in that it is parametric over the
    type of "model" and of "stabilizing term" used, and therefore over the
    exact Master Problem that has to be solved at each iteration: it just
    relies over an object of class MPSolver [see MPSolver.h] to solve it.

    As any other NDOSolver, the Bundle class requires that the function to
    be minimized be available under the form of a FiOracle object [see
    FiOracle.h].

    The interface between the Bundle object and the applications that will
    use it is mostly derived from the interface of the base class NDOSolver
    [see NDOSlvr.h], thus greatly simplifying the task of using different
    NDOSolvers for solving the same NDO problem. */

class Bundle : public NDOSolver
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

/** Public enum which "extends" the enum NDOSolver::NDOParam for handling
    the Bundle-specific algorithmic parameters in (the two overloaded
    versions of) Bundle::SetPar() [see below]. */

   enum BndParam { kBPar1 = kLastNDOParam ,
		   kBPar2 , kBPar3 , kBPar4 , kBPar5 , kBPar6 ,
		   km1 , km3 ,
		   kmxIncr , kmnIncr , kMnSSC , kmxDecr , kmnDecr , kMnNSC , 
		   ktMaior , ktMinor , ktInit , ktCurr , ktSPar1 , ktSPar2 ,
		   kPPar1 , kPPar2 , kPPar3 ,
		   kMPEFsb , kMPEOpt
                   };

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

   Bundle( istream *iStrm = 0 );

/**< Constructor of the class. The parameter `iStrm', if provided, is taken as
   a pointer to a istream from which the algorithmic parameters for the Bundle
   algorithm are sequentially read in the following order. Each parameter
   must be placed at the beginning of a separate line, max 255 characters
   long, with all the rest of the line up to the first newline character '\n'
   (apart from a separating whitespace) being available for comments. Any
   line whose first character is '#' and any blank line is ignored. If 0 is
   passed, the file ends before reaching a given parameter, or some parameter
   is in the wrong format, each non-specified parameter is given a default
   value, shown in [] below.

   Note that `iStrm' is also passed to the constructor of the base class [see
   NDOSolver.h], which reads its algorithmic parameters out of it. It has to
   be remarked that the Bundle class adds its specific "twists" to some of the
   algorithmic parameters of the base class, in particular the ones regarding
   inexact computation of the function:

    HpNum EInit    [1e-2]
    HpNum EFnal    [1e-6]
    HpNum EDcrs    [.95]
    Index EStps    [0]

   For this it has to be remarked that the Bundle has an "emergency mechanism"
   such that if it detects that it is "stalling" because of insufficient
   accuracy in the oracle computation, it will automatically decrease the
   relative precision required to the FiOracle [see SetPrecision() in
   FiOracle.h]. Then, the above parameters now have the following meaning:
 
    HpNum EInit    ABS( EInit ) is the initial, and *maximum*, precision
                   required to the oracle, but the sign tells how the
		   "emergency mechanism" alluded to above interacts with
   the "regular mechanism" controlled by these parameters. In particular,
   if EInit > 0 then the accuracy is monotonically non-increasing: if the
   "emergency mechanism" reduces it, then it will remain "at least as small"
   in all the following iterations, even if the value computed by the
   "regular mechanism" would be larger. If EInit < 0 instead, then each time
   a "regular step" is computed the precision is reset to that dictated by
   the "regular mechanism", even if it is larger than the current one (for
   instance, if EStps == 0, see below, then the precision is set to a "fixed"
   value).

    HpNum EFnal    The other three parameters define a general formula that
                   sets the "regular mechanism" for changing the precision
		   along iterations. Note that EFnal has a completely
   different meaning as the one postulated by the base class (smallest
   precision) because that makes no sense: the "final" precision clearly
   has to be EpsLin. In fact, a value larger than EpsLin would make it
   impossible (in theory) to reach EpsLin-accuracy for the overall
   optimization, and a value smaller than EpsLin is wasteful as a higher
   precision than EpsLin is not required. The idea is that the precision
   should improve along the iterations, and the "speed" at which this
   happens is dictated by EStps and EFnal; however, one can also keep the
   precision "fixed" by setting EStps == 0. In this case, having defined

     EpsU = ( ReadDStart( tStar ) + ReadSigma() ) / ReadFiVal()

   the current estimate of the optimality measure, the formula is

     precision = / ABS( EInit )          if EDcrs >= 0
                 \ ABS( EDcrs ) * EpsU   if EDcrs < 0

   and this is kept fixed along all the iterations (except, EpsU is not
   really fixed); only the "emergency mechanism" will increase it if this
   is absolutely needed. If instead EStps != 0, having defined

     opt = / ABS( EInit )   if EDcrs >= 0
           \ EpsU           if EDcrs < 0

     h = / NrIter()   if EStps > 0    ,     k = ceil( h / ABS( EStps ) )
         \ NrSSs()    if EStps < 0

   the formula is:

     precision = opt * / ABS( EDcrs )^{ EFnal * k }   if EFnal >= 0
                       \ ABS( EDcrs ) * k^{ EFnal }   if EFnal < 0

   (while ensuring precision <= ABS( EInit )).

   Since the Bundle constructor is executed after the one of NDOSolver, the
   following Bundle-specific parameters have to be found in the stream
   *after* those of the base class:

    Index BPar1    [10] If an item has had a zero multiplier [see ReadMult()
                   in NDOSolver.h] for the last BPar1 steps, it is eliminated;
		   if BPar1 is "too small" precious information may be lost,
		   but keeping the "bundle" small obviously makes the Master
		   Problem cheaper.

    Index BPar2    [100] Maximum dimension of the bundle: has more or less the
                   Same "problems" as BPar1, but if the latter is well chosen
		   then BPar2 can be kept big while the "B-strategy" keeps the
    actual number of items low. A small BPar2 can affect the convergence of
    the algorithm, in theory as well as in practice, if aggregation is not
    allowed. However, an unnecessarily large BPar2 may force the MP Solver to
    allocate a large amount of memory without a real need.

    HpNum BPar3    [-1] Maximum ...
    HpNum BPar4    [-1] ... and minimum number of new subgradients/
                   constraints (items) to be fetched from the FiOracle for
		   each function evaluation. Two different ways are given for
    specifying these numbers: positive values are (rounded up and) taken as
    absolute values, while negative numbers are first multiplied by
    FiOracle::GetNrFi()---the number of components of Fi()---(and then
    rounded up); thus, the default vale "-1" stands for "one for each of the
    components of Fi()". Clearly, the "finalized" value of BPar3 has to be
    <= BPar2, and the "finalized" value of BPar4 has to be <= than that.

    HpNum BPar5    [30] These parameters control how the actual number of
    Index BPar6    [0] subgradients/constraints (items) that are requested
                   to the FiOracle varies, between BPar4 and BPar3, as the
		   algorithm proceeds; note that what varies in practice is
    the maximum number, as it is always legal for the FiOracle to refuse
    giving other items, although the Bundle code will complain and stop if
    less than BPar4 are given. In the Bundle code, the number

       EpsU = Sigma + D_{tStar}*( z* ) / max( | ReadFiVal() | , 1 ) ,

    where Sigma = Sum_i Fi[ i ]_{B,Lambda}*( z[ i ]* ) + \sigma_L( w ) and
    z* = - Sum_i z[ i ]* is the optimal solution of the stabilized Dual
    Master Problem [see MPSolver.h], is used as an estimate of the relative
    gap between the current and the optimal solution; that is, IsOptimal()
    returns true if EpsU <= EpsLin. Thus, the number EpsLin / EpsU is always
    smaller than one, and typically increases as the algorithm proceeds.
    Depending on the value of BPar6, the following formulae for the actual
    value of BPar3, aBP3, are used:
      0: aBP3 is set to (the "finalized" value of) BPar3 and never changed;
      1: if BPar5 > 0 then aBP3 is initialized to (...) BPar4 and increased
         every BPar5 iterations, while if BPar5 <= 0 then aBP3 is
	 initialized to BPar3 and decreased every - BPar5 iterations;
      2: aBP3 is set to
          ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 * ( EpsLin / EpsU )
      3: aBP3 is set to
          ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / sqrt( EpsU / EpsLin )
      4: aBP3 is set to
          ( BPar5 > 0 ? BPar4 : BPar3 ) + BPar5 / log10( EpsU / EpsLin )

    HpNum m1       [.1] SS condition: if DeltaFi >= | m1 | * Deltav, then a
                   SS is done. What is taken as Deltav depends on the sign of
		   m1: if m1 > 0 then Deltav = - v* (the decrease predicted
    by the model), while if m1 < 0 then Deltav = - ( v* + D_t( d* ) ), i.e.
    the optimal objective function value the dual Master Problem. Since
    - v* >= - ( v* + D_t( d* ) ), the second condition is weaker and may
    lead to a larger number of SS (and therefore possibly a fater convergence)
    while still ensuring global convergence. The value m1 = 0, i.e., perform a
    SS for whatever small improvement in the objective function, can only be
    used, at least in theory, for some classes of functions (the polyhedral
    ones) and some under assumptions on the MP.

    HpNum m3       [3]  A nevly obtained subgradient is deemed "useless" if
                   Alfa >= m3 * Sigma; in this case, if a NS has to be done,
		   t is decreased. This parameter is mostly critical: if no
    "long-term" t-strategy [see tSPar1 below] is used, values < 2/3 usually
    make t to decrease rather fast to tMinor [see below], possibly making the
    algorithm to perform very short steps and therefore to converge very
    slowly. When a "long-term" t-strategy is used, .9 may be a good value.

    HpNum mxIncr   [10]  each time t grows, the new value of t is chosen in
    HpNum mnIncr   [1.5] the interval [t * mnIncr, t * mxIncr] (t is the
                   previous value)
    Index MnSSC    [0]   minimum number of consecutive SS with the same t that
                         have to be performed before t is allowed to grow

    HpNum mxDecr   [.1]  each time t diminishes, the new value of t is chosen
    HpNum mmDecr   [.66] in the interval [t * mxDecr, t * mnDecr] (t is the
                   previous value)
    Index MnNSC    [0]   minimum number of consecutive NS with the same t that
                         have to be performed before t is allowed to diminish

    HpNum tMaior   [1e+6] Maximum, ...
    HpNum tMinor   [1e-6] ... minimum, and ...
    HpNum tInit    [1] ... initial value of t. These parameters may be
                   critical, but they are not very difficult to set. Usually,
		   there is a "right" order of magnitude for t, that is the
    one that is guessed by the t-heuristics during most of the run, even
    though the starting value is very different. Hence, a good setting for
    tInit is in that order of magnitude, while tMinor should be set small
    enough to never enter into play. Note that t is always kept <= tStar
    [see NDOSolver.h], and that a "good" value for tStar (i.e., one that
    actually ensures that the stopping point is EpsLin-optimal) is usually
    one or two orders of magnitude larger than a "good" tInit.

    Index tSPar1   [0] Select the t-strategy used. This field is coded
                       bit-wise in the following way.
		   The first two bits control which heuristics are used to
		   compute a new value of t when increasing/decreasing it.
    There are two heuristics available, H1 and H2, both based on a quadratic
    interpolation of the function but differing in which derivative is used:
    H1 uses the derivative in the new tentative point, and it is guaranteed
    to produce a value greater than the current one if and only if the scalar
    product between the direction and the newly obtained subgradient is < 0
    (indicating that a longer step along the same direction could have been
    advantageous), while H2 uses the derivative in the current point and it
    does not possess this property. The value of the first two bits of
    tSPar1 has the following meaning:
      bit 0:  which heuristic is used to increase t: 0 = H1, 1 = H2
      bit 1:  which heuristic is used to decrease t: 0 = H2, 1 = H1
    The following bits of tSPar1 tell which long-term t-strategy is used,
    with the following values:
      0 (+ 0):  none, only the heuristics are used
      1 (+ 4):  the "soft" long-term t-strategy is used: an optimality
                estimate EpsU is maintained which estimates how far from the
		optimal value one currently is, and decreases of t are
		inhibited whenever v < tSPar2 * EpsU * | Fi |
      2 (+ 8):  the "hard" long-term t-strategy is used: an optimality
                estimate EpsU is maintained as above and t is increased
		whenever v < tSPar2 * EpsU * | Fi |
      3 (+12):  the "balancing" long-term t-strategy is used, where the two
                terms D*_t( -z* ) and Sigma* are kept of "roughly the same
		size": if D*_1( -z* ) <= tSPar2 * Sigma* then t increases
		are inhibited (increasing t causes a decrease of D*_1( -z* )
		that is already small), if tSPar2 * D*_1( -z* ) >= Sigma*
		then t decreases are inhibited (decreasing t causes an
		increase of D*_1( -z* )	that is already big).
      Still later bits of tSPar1 activate "special cases" t-strategies:
      4 (+16):  the "endgame" t-strategy is used, where if D*_1( -z* ) is
                "small" (~ 1/10 of the current absolute epsilon) t is
                decreased no matter what the other strategies dictated.
                The rationale is that we are "towards the end" of the
                optimization and here t needs decrease. However, note that
                having D*_1( -z* ) "small" is no guarantee that we actually
                are at the end, especially if the FiOracle dynamically
                generates its variables, so use with caution.

    HpNum tSPar2   [.1] Numerical parameter for the long-term t-strategies,
                   see tSPar1 above.

    Index PPar1    [30] Parameters controlling the variables generator: "price
                   in" (discover if new variables have to be added) is done
		   all the first PPar1 iterations ...
    Index PPar2    [10] ... and then every PPar2 iterations; note that the
                   price in is done anyway each time convergence is detected.
		   If PPar2 == 0, all the variables are present from the
		   beginning (PPar1 is ignored if PPar2 == 0).
    Index PPar3    [5] A variable that has been inactive for the last PPar3
                   pricings (this one included) is eliminated: note that the
		   "price out" operation is done every PPar2 iterations, so
    that a variable that is eliminated is likely to have been inactive for
    (about) PPar2 * PPar3 iterations. For PPar3 == 1, a variable is eliminated
    in the very pricing in which it is discovered to be zero (and the
    direction saying that it would stay zero). If PPar3 == 0, variables are
    *never* removed. PPar3 is ignored if PPar2 == 0.

    HpNum MPEFsb   [1e-6] (relative) precision required to the MP Solver as
                   far as constraints satisfaction is concerned.
    HpNum MPEOpt   [1e-6] (relative) precision required to the MP Solver as
                   far as optimality of the solution is concerned.

   The methods SetPar() of the base class are extended to allow the setting
   of the Bundle-specific parameters [see below]. */

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   void SetMPSolver( MPSolver *MPS = 0 );

/**< Gives to the Bundle object a pointer to an object of class MPSolver that
   will be used as Master Problem Solver during the Bundle algorithm.

   The MPSolver can be changed during the life of a Bundle object, but this
   change clearly forces the reset of all the information about the function
   accumulated so far (the bundle). Passing a 0 pointer does exactly this
   job.

   In contrast with the standard rule, SetMPSolver() can be called even if no
   FiOracle is currently set to the NDO solver, that is, SetFiOracle() has not
   already been called or it has been called the last time with 0. */

/*--------------------------------------------------------------------------*/

   void SetFiOracle( FiOracle *Fi = 0 );

/*--------------------------------------------------------------------------*/

   void SetLambda( cLMRow tLambda = 0 );

/*--------------------------------------------------------------------------*/

   void KeepBestLambda( const bool KBL = true );

/*--------------------------------------------------------------------------*/

   inline void SetPar( const int wp , const int value );

/**< Extends NDOSolver::SetPar( , cIndex ) for handling the Bundle-specific
   parameters; the enum BndParam is used (in the obvious way) for selecting
   the parameter to be set.

   Some remarks are needed about setting some of the parameters:

   - the extra parameter Bundle::ktCurr can be used to force "from the
     outside" a change of the current value of t. Note that the change
     only takes effect the next time the Master Problem is solved, and
     could be overridden if the normal t-strategy is executed in the
     meantime. If this is done in GetFiStatus(), a possible trick is to
     have that call to return kFiChgd, so that the next iteration is
     immediately started and the change of t takes effect.

   - changing BPar2 may force the Master Problem solver to restructure its
     internal data structures [see MPSolver::SetDim()], and therefore may be
     expensive; also, diminishing the bundle size forces some subgradients to
     be discarded. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   inline void SetPar( const int wp , cHpNum value );

/**< Extends NDOSolver::SetPar( , cHpNum ) for handling the Bundle-specific
   parameters; the enum BndParam is used (in the obvious way) for selecting
   the parameter to be set.

   Some remarks are needed about setting some of the parameters:

   - if a "variable Fi() precision" scheme is used, the current precision is
     *not* re-initialized to its initial value if Solve() is called more than
     once (it is, however, if the FiOracle changes or ReSetAlg() is called);
     analogously, the counter of Fi() evaluations left before decreasing the
     precision is *not* reset. SetPar( NDOSolver::kEInit , 0 ) forces a
     reset of the scheme. */

/*--------------------------------------------------------------------------*/

   void SetNDOLog( std::ostream *outs = 0 , const char lvl = 0 );

/**< lvl controls the "level of verbosity" of the code. The first four bits
   of lvl have the following meaning:

    0  =>  no log at all (also assumed if log = 0);

    1  =>  "basic" log: only the errors are reported;

    2  =>  a detailed step-by-step log of the algorithm is displayed;

    4 .. 15 unused, available to derived classes;

   The remaining of lvl is coded bit-wise, so that logging of specific
   features can be activated independently:

    bit 4 == 1 (+ 16)  =>  every operation on the set of items (adding,
                           removing, aggregating) is logged;

    bit 5 == 1 (+ 32)  =>  the activity of the variables generator [see PPar*
                           in the comments of Bundle() above] is logged. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   NDOStatus Solve( void );

/**< Tries to minimize the function. If there are constraints the
   minimization is divided into two phases: in Phase 0 a feasible point is
   sought for, and in Phase 1 the function is minimized moving on feasible
   points only.

   Returns   if

   kOK       optimization has been successful: a solution that is "optimal"
             (w.r.t. the current parameters settings) has been found;

   kUnbndd   there has been an error in the FiOracle, i.e. Fi() has returned
             - HpINF, or the function is unbounded below: the latter case can
             be detected only if a lower bound on the min. value of Fi() is
             available [see FiOracle::GetMinusInfinity()];

   kUnfsbl   the polyhedral set defined by the constraints is empty: in this
             case, the primal optimal solution is an unbounded *extreme ray*
	     for the dual problem;

   kStopped  Solve() has been stopped, either by FiOracle::GetFiStatus() or
             by EveryIteration() [see below];

   kStpIter  the max. number of iterations has been exhausted;

   kLwPrcsn  the algorithm cannot proceed because the function cannot be
             computed with enough precision [see FiOracle::SetPrecision()]:
             this means that the function has been minimized up to the
	     maximum extent that is possible due to the limited precision that
             the FiOracle can provide;

   kError    There was an error in the Master Problem solver, and this
             condition has not been corrected by the elimination of items.

   Note that, whatever the exit condition be, the current point is always
   available by calling ReadSol(), and its Fi() value by calling ReadFiVal()
   [see below]. There are constraints and kStpIter has been returned, Phase
   0 may *not* have finished yet: hence, the current point may *not be
   feasible*, so that ReadFiVal() may return + HpINF.

   Solve() is "virtual" in order to allow derived classes to implement
   different "main" strategies: this can be easily done by exploiting the
   methods FormD(), FormLambda1(), FiAndGi(), GotoLambda1(), UpdtCntrs() and
   EveryIteration() in the protected interface of the class [see below]. */

/*--------------------------------------------------------------------------*/

   void ReSetAlg( unsigned char RstLvl = 0 );

/**< Resets the internal state of the Bundle algorithm. Since several
   different things can be reset independently, RstLvl is coded bit-wise:

   - bit 0: if 0, all the algorithmic parameters are reset to the default
     values read by the stream/set by SetPar(), while if 1 they are left
     untouched;

   - bit 1: if 0 the current point is reset to the all-0 vector, while if
     1 it is left untouched;

   - bit 2: if 0, all the subgradients are removed from the bundle, except
     the constant (sub)gradient of the linear 0-th component, while if 1
     the subgradients are left there;

   - bit 3: if 0, all the constraints are removed from the bundle, while
     if 1 the constraints are left there.

   - bit 4: if 0 the value of Fi() in the current point is reset to HpINF
     (i.e., unknown), while if 1 it is left untouched; note that resetting
     the current point [see bit 1] has this as a side-effect, regardless to
     the value of bit 4. */

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

   cLMRow ReadSol( cIndex_Set &I , Index &D );

/**< If Solve() has returned a kOK and the tStar has been properly set, the
   point returned by ReadSol() - and, a fortiori, the one returned by
   ReadBestSol() - is EpsLin-optimal. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   cLMRow ReadBestSol( cIndex_Set &I , Index &D );

/**< If Solve() has returned a kOK and the tStar has been properly set, the
   point returned by ReadSol() - and, a fortiori, the one returned by
   ReadBestSol() - is EpsLin-optimal. */

/*--------------------------------------------------------------------------*/

   inline bool IsOptimal( HpNum eps = 0 ) const;

/*--------------------------------------------------------------------------*/

   inline HpNum ReadSigma( void ) const;

   inline HpNum ReadDStart( cHpNum tt = 1 ) const;

/**< Return the current value of the "two pieces" of the stopping criterion.
   Indeed, IsOptimal() returns true if

     ReadDStart( tStar ) + ReadSigma() <= EpsLin * ReadFiVal()

   where tStar and EpsLin are the general stopping parameters defined in
   NDOSolver [see NDOSolver.h]. For more details on what exactly ReadDStart()
   and ReadSigma() return see the same-named methods in MPSolver.h. */

/*--------------------------------------------------------------------------*/

   inline bool CurrentIsLast( void );

/**< The point returned by ReadSol() [see above], called the `current point',
   may not be the point corresponding to the last call to FiOracle::Fi() [see
   FiOracle.h], because a number of `Null Steps' may have performed done after
   the last `Serious Step'.

   This method returns true if the current point is actually the last point
   for which Fi() was evaluated, and false otherwise. This may be useful in
   some cases. */

/*--------------------------------------------------------------------------*/

   HpNum ReadFiVal( cIndex wFi = Inf< Index >() );

   HpNum ReadBestFiVal( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   HpNum ReadLBMult( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   HpNum Readt( void )

/**< Return the current value of the (tremendous) t parameter. */
   {
    return( t );
    }

/*--------------------------------------------------------------------------*/
/** Returns the current number of Serious Steps; note that analogously to
    the number of iterations [see NrIter()] this is *not* necessarily reset
    each time Solve() is called. */

   inline Index NrSSs( void ) const
   {
    return( ParSS );
    }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   inline void GetPar( const int wp , int &value );

   inline void GetPar( const int wp , HpNum &value );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   void AddVariables( Index NNwVrs , cLMRow IVs = 0 );

   void RemoveVariables( cIndex_Set whch = 0 , Index hwmny = 0 );

   void ChgFiV( cIndex wFi = Inf< Index >() );

   void ChgSbG( cIndex strt = 0 , Index stp = Inf< Index >() ,
		cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   void RemoveItem( cIndex Name );

/**< Remove the item `Name' from the bundle.

   Note that this method is not a part of the NDOSolver interface since it
   is somewhat subsumed by ChgFiV() and ReSetAlg() [see above]; however, it
   may still be of some use, so it is kept. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void RemoveItems( void );

/**< Remove all the items from the bundle, except the (sub)gradient of the
   linear 0-th component of Fi).

   Note that this method is not a part of the NDOSolver interface since it
   is somewhat subsumed by ChgFiV() and ReSetAlg() [see above]; however, it
   may still be of some use, so it is kept. */

/** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

  virtual ~Bundle();

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/** Protected enum describing the possible return values of EveryIteration()
    [see below]. */

   enum EIStatus { kEINorm = 0,
		   kEIAbort ,
		   kEILoopNow ,
		   kEIContAnyway
                   };

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- HOOK FOR DERIVED CLASSES -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Hook for derived classes
    The following method is called by the solver to give derived classes,
    which may implement variants of the Generalized Bundle approach, a
    better control over the optimization process.
    @{ */

   virtual EIStatus EveryIteration( void );

/**< This method is an "hook" for derived classes: it is called at Every
   Iteration, between the computation of the tentative direction and the
   computation of Fi(). It can serve to various purposes, primarily checking
   extra stopping conditions or interfering with the usual stopping conditions
   of the Bundle code: however, any kind of operation can be performed here
   inside, e.g. adding or removing variables from the problem [see
   [Add/Remove]Variable() above].
   More in general, this method can be used to merge the main cycle of the
   Bundle method within any other however complex code: the Bundle gives out
   the control at this time, and resumes its operations when EveryIteration()
   returns. The returned value influences the behaviour of the Bundle for the
   current iteration:

   kEINorm        the current iteration is continued normally;

   kEIAbort       the whole algorithm is aborted, and Solve() is immediately
                  terminated returning kAbort: this is useful for instance
                  to enforce new termination criteria;

   kEILoopNow     the current iteration is aborted, i.e. the stopping
                  condition is *not* checked, and Fi() is *not* called: the
		  next iteration is immediately started, but the iterations
		  count is *not* increased. This is useful e.g. if something
		  has been changed in the data of the problem that suggests
		  to try a new direction, like a new "active" or variable
		  [see AddVariable() above] to be inserted;

   kEIContAnyway  the current iteration is continued normally but for the
                  fact that the stopping condition is *not* checked. */

/** @} ---------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to really extend the functionality of the Bundle class     --*/
/*--  (e.g. by implementing a different subproblem solver) may use these  --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   void FormD( void );

/* When no variables generation is done (PPar2 == 0), FormD() just calls
   SolveMP() once and calculates the direction d: however, it also implements
   some strategies to survive to "fatal" failures in the subproblem solver,
   typically eliminating some of the items in the bundle.

   Set the protected field Result to kOK if (eventually after some "fatal"
   failure) a tentative descent direction could be found, to kUnfsbl if the
   MP is dual unfeasible and to kError if this was returned by SolveMP(): in
   the latter cases, the whole algorithm must abort.

   If variables generation is done (PPar2 > 0), this is where the
   corresponding strategies are implemented: in this case, SolveMP() can be
   called more than once within the same call to FormD(), since the resulting
   direction has to be optimal w.r.t. all the current "active set" of
   variables. */

/*--------------------------------------------------------------------------*/

   void FormLambda1( HpNum Tau );

/* After a (successful) call to FormD(), sets the new tentative point Lambda1
   (a protected field of type LMRow) as Lambda1 = Lambda + ( Tau / t ) * d. */

/*--------------------------------------------------------------------------*/

   bool FiAndGi( void );

/* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
   constraints) in the bundle. Returns true <=> the newly obtained information
   changes the solution of the MP. */

/*--------------------------------------------------------------------------*/

   void GotoLambda1( void );

/* Move the current point to Lambda1. */

/*--------------------------------------------------------------------------*/

   void UpdtCntrs( void );

/* Updates the out-of-base counters for all items in the Bundle. */

/*--------------------------------------------------------------------------*/

   void SimpleBStrat( void );

/* Eliminate outdated items, i.e., these with "large" out-of-base counter. */

/*--------------------------------------------------------------------------*/

   void Log1( void );

   void Log2( void );

/*--------------------------------------------------------------------------*/

   bool CheckAlfa( const bool All = false );

   void StrongCheckAlfa( void );

/*--------------------------------------------------------------------------*/

   void UpdtLowerBound( void );

/*--------------------------------------------------------------------------*/

   void UpdtaBP3( void );

   void CmptaBPX( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  MPSolver *Master;    // (pointer to) the Master Problem Solver

  int BPar1;           // parameter for removal of items (B-strategy)
  int BPar2;           // max Bundle size
  HpNum BPar3;         // max number of items fetched from Fi() at each call
  HpNum BPar4;         // min number of items fetched from Fi() at each call
  HpNum BPar5;         // control how the actual BPar3 changes over time
  int BPar6;           // same as above
 
  HpNum mxIncr;        // max increase/decrease t parameters:
  HpNum mnIncr;        // see the description in the constructor
  int MnSSC;
  HpNum mxDecr;
  HpNum mnDecr;
  int MnNSC;

  HpNum m1;            // parameters for deciding if a SS/NS is done:
  HpNum m3;            // see the description in the constructor

  HpNum tMaior;        // max value for t
  HpNum tMinor;        // min value for t
  HpNum tInit;         // initial value for t

  int tSPar1;          // long-term t-strategy parameters
  HpNum tSPar2;        // see the description in the constructor

  int PPar1;           // pricing related parameters
  int PPar2;           // see the description in the constructor
  int PPar3;

  HpNum MPEFsb;        // precision required to the MP Solver (feasibility)
  HpNum MPEOpt;        // precision required to the MP Solver (optimality)

  Index MaxNumVar;     // maximum number of variables
  Bool_Vec IsEasy;     // tells whether any component of Fi is "easy"

  LMRow Lambda;        // the current point
  LMRow Lambda1;       // the tentative point
  LMRow LmbdBst;       // the best point found so far

  Index_Set LamBase;   // the set of indices of Lambda
  Index_Set Lam1Bse;   // the set of indices of Lambda1
  Index LamDim;        // dimension of LamBase

  bool KpBstL;         // if LmbdBst has to be kept
  bool BHasChgd;       // true if LamBase has changed during the latest
                       // pricing (never set to true if PPar2 == 0, unless
                       // at the very first call to the oracle)
  bool LHasChgd;       // true if Lambda has changed since the latest call
                       // to FiAndGi(): allows repeated calls in the same
                       // Lambda, e.g. with increasing precision
  bool tHasChgd;       // true if t has changed since the last MP

  HpRow FiLambda;      // Fi[ k ]( Lambda )
  HpRow FiBest;        // best value(s) of Fi found so far
  HpRow FiLambda1;     // Fi[ k ]( Lambda1 )
  HpRow RfrncFi;       // the value of Fi[ k ]() where the zero of the Cutting
                       // Plane models are fixed: it is == FiLambda[ k ]() but
                       // when FiLambda[ k ]() == HpINF
  // HpNum b0;         // the constant in the affine 0-th component of Fi
  // we could get to know it, if it was useful (which it is not)

  Index_Set whisZ;     // the position in the bundle where the "aggregate
                       // subgradient" Z[ k ] of "component" k is kept in
                       // whisZ[ k ]; Inf< Index >() means it is not in the
                       // bundle
  Index_Set whisG1;    // name of the "representative subgradient" for each
                       // component of Fi()
  LMRow ScPr1;         // ScalarProduct( dir , G[ WhIsG1[ k ] ] )
  HpRow Alfa1;         // linearization error of G[ WhIsG1[ k ] ] w.r.t. the
                       // current point Lambda
  HpRow DeltaAlfa;     // correction of Fi-values due to inexactness

  HpRow LowerBound;    // Lower Bound over (the various components of) Fi

  HpNum t;             // the (tremendous) t parameter
  HpNum Prevt;         // what t were before being changed for funny reasons

  HpNum Sigma;         // Sigma*: convex combination of the Alfa's
  HpNum DSTS;          // D*_{t*}( -z* ), the other part of the dual objective
  HpNum vStar;         // v*, the predicted improvement
  HpNum Deltav;        // the "desired improvement" in the Fi-value

  HpNum DeltaFi;       // FiLambda - FiLambda1
  HpNum EpsU;          // precision required by the long-term t-strategy

  HpNum EpsCurr;       // the precision currently required to the FiOracle
  HpNum EpsFi;         // the last precision passed to the FiOracle (can be
                       // different from EpsCurr)

  int ParSS;           // number of SS within the present call to Solve()

  int CSSCntr;         // counter of consecutive SS
  int CNSCntr;         // counter of consecutive NS

  Index FreDim;        // number of free positions in the Bundle
  Index_Set FreList;   // list (heap) of free positions in the Bundle

  SIndex_Set OOBase;   // Out-Of-Base counters:
                       // = Inf< SIndex >() means no item is there
                       // = k > 0 means out of base since k iterations
                       // = 0 means in the current base but potentially
                       //   removable
                       // = a *finite* negative value - k means not
                       //   removable for the next k iterations: note that
                       //   some items in base may be such
                       // = - Inf< SIndex >() means unremovable
  Index_Set InctvCtr;  // "out of base" counter for variables
  Index_Set nBase;     // temporary

  bool TrueLB;         // true if LowerBound is a "true" lower bound rather
                       // than just the "minus infinity"
  bool LBHasChgd;      // true some LowerBound has changed
  bool SSDone;         // true if the laste step was a SS

  bool ItemsChgd;      // true if no items have been added to MP

  FiOracle::FiStatus *FiStatus;

  #if( NONMONOTONE )
   HpRow FiVals;       // Fi-values for the last NONMONOTONE SSs
  #endif

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  inline bool DoSS( void );

/*--------------------------------------------------------------------------*/

  inline void Delete( cIndex i );

/*--------------------------------------------------------------------------*/

  inline bool FindNextSG( Index &wFi );

/*--------------------------------------------------------------------------*/

  inline Index BStrategy( cIndex wFi );

/*--------------------------------------------------------------------------*/

  inline Index FindAPlace( cIndex wFi );

/*--------------------------------------------------------------------------*/

  inline HpNum Heuristic1( void );

  inline HpNum Heuristic2( void );

/*--------------------------------------------------------------------------*/

  inline void InitMP( void );

/*--------------------------------------------------------------------------*/

  inline void AggregateZ( cHpRow Mlt , cIndex_Set MBse , Index MBDm ,
			  cIndex wFi , cIndex whr );

/*--------------------------------------------------------------------------*/

  inline HpNum WhichFi( cHpRow FiVect , cIndex wFi );

/*--------------------------------------------------------------------------*/

  inline void MemDealloc( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  Index MBDim;      // number of items in the optimal multiplier base
  Index aBP3;       // current max number of items to be fetched
  Index aBP4;       // min number of items to be fetched

/*--------------------------------------------------------------------------*/

 };  // end( class Bundle )

/** @} end( group( Bundle_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Bundle::SetPar( const int wp , const int value )
{
 switch( wp ) {
  case( kBPar1 ):  BPar1 = value; break;
  case( kBPar2 ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( value != BPar2 ) {
    if( Oracle ) {            // communicate the changes to the Oracle
     Oracle->SetMaxName( value );

     // resize data structures depending on BPar2
     SIndex_Set oldOOBase = OOBase;
     Index_Set oldFreList = FreList;

     OOBase = new SIndex[ value ];
     FreList = new Index[ value ];

     const int copydim = std::min( value , int( BPar2 ) );

     if( oldOOBase ) {
      for( int i = 0 ; i < copydim ; i ++ )
       OOBase[ i ] = oldOOBase[ i ];
      delete[] oldOOBase;
      }
     else
      for( int i = 0 ; i < value ; i ++ )
       OOBase[ i ] = Inf< SIndex >();

     if( oldFreList ) {
      for( int i = 0 ; i < copydim ; i ++ )
       FreList[ i ] = oldFreList[ i ];
      delete[] oldFreList;
      }
     }

    BPar2 = value;          // set the parameter
    if( Master && Oracle )  // communicate the changes to the Master
     Master->SetDim( BPar2 , 0 , PPar2 ? true : false );
    }
   break;
  case( kBPar6 ):  BPar6 = value; CmptaBPX(); break;
  case( kMnSSC ):  MnSSC = value; break;
  case( kMnNSC ):  MnNSC = value; break;
  case( ktSPar1 ): tSPar1 = value; break;
  case( kPPar1 ):  PPar1 = value; break;
  case( kPPar2 ):  PPar2 = value; break;
  case( kPPar3 ):  PPar3 = value; break;
  default:         NDOSolver::SetPar( wp , value );
  }
 }  // end( Bundle::SetPar( Index ) )

/*--------------------------------------------------------------------------*/

inline void Bundle::SetPar( const int wp , cHpNum value )
{
 switch( wp ) {
  case( kBPar3 ):  BPar3 = value; CmptaBPX(); break;
  case( kBPar4 ):  BPar4 = value; CmptaBPX(); break;
  case( kBPar5 ):  BPar5 = value; CmptaBPX(); break;
  case( km1 ):     m1 = value; break;
  case( km3 ):     m3 = value; break;
  case( kmxIncr ): mxIncr = value; if( mxIncr < 1 ) mxIncr = 1; break;
  case( kmnIncr ): mnIncr = value; if( mnIncr < 1 ) mnIncr = 1; break;
  case( kmxDecr ): mxDecr = value; if( mxDecr > 1 ) mxDecr = 1; break;
  case( kmnDecr ): mnDecr = value; if( mnDecr > 1 ) mnDecr = 1; break;
  case( ktMaior ): tMaior = value;
                   if( t > tMaior ) { t = tMaior; tHasChgd = true; }
                   if( tInit > tMaior ) tInit = tMaior;
		   break;
  case( ktMinor ): tMinor = value; break;
                   if( t < tMinor ) { t = tMinor; tHasChgd = true; }
                   if( tInit < tMaior ) tInit = tMinor;
		   break;
  case( ktInit ):  tInit = value;
                   if( t != tInit ) { t = tInit; tHasChgd = true; }
                   break;
  case( ktCurr ):  t = value; tHasChgd = true; break;
  case( ktSPar2 ): tSPar2 = value; break;
  case( kEInit ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( EInit != value ) {
    EInit = value;
    EpsCurr = EInit >= 0 ? EInit : - EInit;
    }
   break;
  case( kMPEFsb ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( Master )
    Master->SetPar( MPSolver::kFsbEps , MPEFsb = value );
   break;
  case( kMPEOpt ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( Master )
    Master->SetPar( MPSolver::kOptEps , MPEOpt = value );
   break;
  default:
   NDOSolver::SetPar( wp , value );
  }
 }  // end( Bundle::SetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

inline bool Bundle::IsOptimal( HpNum eps ) const
{
 HpNum FiL = *FiLambda;

 if( ( FiL == Inf< HpNum >() ) || ( vStar == -Inf< HpNum >() ) )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;

  if( eps <= 0 )
   eps = EpsLin;

  return( DSTS + Sigma <= eps * FiL );
  }
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReadSigma( void ) const
{
 return( Sigma );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReadDStart( cHpNum tt ) const
{
 return( Master->ReadDStart( tt ) );
 }

/*--------------------------------------------------------------------------*/

inline bool Bundle::CurrentIsLast( void )
{
 return( SSDone );
 }

/*--------------------------------------------------------------------------*/

inline void Bundle::GetPar( const int wp , int &value )
{
 switch( wp ) {
  case( kBPar1 ):  value = BPar1; break;
  case( kBPar2 ):  value = BPar2; break;
  case( kBPar6 ):  value = BPar6; break;
  case( kMnSSC ):  value = MnSSC; break;
  case( kMnNSC ):  value = MnNSC; break;
  case( ktSPar1 ): value = tSPar1; break;
  case( kPPar1 ):  value = PPar1; break;
  case( kPPar2 ):  value = PPar2; break;
  case( kPPar3 ):  value = PPar3; break;
  default:         NDOSolver::GetPar( wp , value );
  }
 }  // end( Bundle::GetPar( Index ) )

/*--------------------------------------------------------------------------*/

inline void Bundle::GetPar( const int wp , HpNum &value )
{
 switch( wp ) {
  case( kBPar3 ):  value = BPar3; break;
  case( kBPar4 ):  value = BPar4; break;
  case( kBPar5 ):  value = BPar5; break;
  case( km1 ):     value = m1; break;
  case( km3 ):     value = m3; break;
  case( kmxIncr ): value = mxIncr; break;
  case( kmnIncr ): value = mnIncr; break;
  case( kmxDecr ): value = mxDecr; break;
  case( kmnDecr ): value = mnDecr; break;
  case( ktMaior ): value = tMaior; break;
  case( ktMinor ): value = tMinor; break;
  case( ktInit ):  value = tInit; break;
  case( ktSPar2 ): value = tSPar2; break;
  case( kMPEFsb ): value = MPEFsb; break;
  case( kMPEOpt ): value = MPEOpt; break;
  default:         NDOSolver::GetPar( wp , value );
  }
 }  // end( Bundle::GetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

 }  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Bundle.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Bundle.h ------------------------------*/
/*--------------------------------------------------------------------------*/
