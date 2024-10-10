/*--------------------------------------------------------------------------*/
/*------------------------- File OSIMPSolver.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the OSIMPSolver class, which solves Master Problems for
 * Bundle algorithms using a generic OSISolverInterface object. This class
 * conforms to the interface defined by the class MPSolver [see MPSolver.h].
 *
 * \author Antonio Frangioni (original idea & implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone (implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' della Calabria \n
 *
 * \author Andrea Nerli (implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OSIMPSolver
 #define __OSIMPSolver /* self-identification: #endif at the end of the file*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* As the name implies, OSIMPSolver uses an object derived from
 * OsiSolverInterface to solve the MP. Certain usage patterns cause the MP to
 * be changed, and *after that* the previous optimal primal and dual solution
 * to be accessed. Certain OsiSolvers graciously allow that, by keeping the
 * previous solution even in face of changes until the MP is explicitly
 * re-solved, while others are more strict and delete any previous solution
 * as soon as anything changes in the MP.
 *
 * This macro, if set to 1, dictates that the OSIMPSolver object stores the
 * solution information of the OsiSolver in its own data structures, so that
 * it remains available even after changes in the MP. Doing so is safe, but
 * it may be useless (and therefore better avoided for higher efficiency in
 * both space and time) with certain OsiSolvers. */

#define PRESERVE_OSI_SOLS 1

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MPSolver.h"

#include "OsiSolverInterface.hpp"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The OSIMPSolver class implements a Master Problem Solver for Bundle
 * algorithms, using a generic OSISolverInterface object. This class
 * derives from MPSolver, and therefore conforms to the generic interface
 * defined there for Master Problem Solvers [see MPSolver.h].  */

class OSIMPSolver : public MPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/// Public enum for handling the Bundle stabilization term [see below].

 enum StabFun { unset = 0 , none , boxstep , quadratic };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Public enum which is used to define the algorithm used in Solve() [see
 * below] by the OSISolverInterface object for solving the Master problem,
 * which may be either LP or QP problem:
 *
 *  - kAuto:  automatic
 *  - kPrim:  Primal simplex
 *  - kDual:  Dual simplex
 *  - kNet:   Network simplex
 *  - kBar:   Barrier
 *  - kSif:   Sifting
 *  - kCon:   Concurrent (Dual, Barrier, and Primal) */

 enum OsiAlg { kAuto = 0 , kPrim , kDual , kNet , kBar , kSif , kCon };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Public enum which is used to tell wich reductions (if any) are used by
 * the OSISolverInterface object while solving the Master problem:
 *
 *  - rNo:    no reductions
 *  - rPrim:  primal reductions only
 *  - rDual:  dual reductions only
 *  - rBoth:  both primal and dual reduction */

 enum OsiRed { rNo = 0 , rPrim , rDual , rBoth };

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
 *  @{ */

 OSIMPSolver( std::istream *iStrm = nullptr );

/** @ ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 void SetDim( cIndex MxBSz = 0 , FiOracle * Oracle = nullptr ,
	      const bool UsAvSt = false ) override;

/*--------------------------------------------------------------------------*/

 void Sett( cHpNum tt = 1 ) override;

/*--------------------------------------------------------------------------*/

 void SetPar( const int wp , cHpNum value ) override;

/*--------------------------------------------------------------------------*/

 void SetThreads( int nthreads ) override;

/*--------------------------------------------------------------------------*/

 void SetLowerBound( cHpNum LwBnd = - Inf< HpNum >() ,
		     cIndex wFi = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/
/** lvl controls the "level of verbosity" of the code. The first four bits
 * of lvl have the following meaning:
 *
 *  0  =>  no log at all (also assumed if log = 0);
 *
 *  1  =>  some messages are printed;
 *
 *  2  =>  as, 1 plus the current Master Problem is also saved as a MPS file
 *         just prior to being solved;
 *
 *  3  ==>  as 2, plus very verbose information is also printed
 *
 *  4 .. 15 unused, available to derived classes. */

 void SetMPLog( ostream *outs = 0 , const char lvl = 0 ) override;

/*--------------------------------------------------------------------------*/
/** Change "int" algorithmic parameters of the NDO solver. This method is
 * used to set the method for solving either LP or QP. */

 void SetAlgo( const OsiAlg algo = kAuto , const OsiRed reduc = rNo );

/** @} ---------------------------------------------------------------------*/
/*---------------------- OSIMPSolver-SPECIFIC METHODS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name OSIMPSolver-specific methods
 *  @{ */

/// Provides OSIMPSolver with the actual OsiSolverInterface
/** Provides OSIMPSolver with an object of any class derived from
 * OsiSolverInterface, which is used to actually solve the Master Problem.
 * Note that the object does *not* become "property" of the OSIMPSolver,
 * which will *not* delete it when it is deleted. However, calling SetOsi()
 * does delete the existing solver to replace it with the new one (if any),
 * hence to delete the current OsiSolverInterface one just has to call it
 * with default == nullptr argument. This is a bit awkward, but it'll do
 * for now. */

 void SetOsi( OsiSolverInterface * osi = nullptr );

/*--------------------------------------------------------------------------*/
 /// gives back the current OsiSolverInterface

 OsiSolverInterface * GetOsi( void ) const { return( osiSlvr ); }

/*--------------------------------------------------------------------------*/
/// Sets the stabilizing term in the Master Problem
/** Sets the stabilizing term in the Master Problem. Possible values are:
 *
 * - none: no stabilization;
 *
 * - boxstep: a box in the infinity-norm in the primal, a 1-norm penalty
 *            in the dual;
 *
 * - quadratic: a 2-norm penalty term both in the primal and in the dual. */

 void SetStabType( const StabFun sf = none );

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
 *  @{ */

 MPStatus SolveMP( void ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading results
 *  @{ */

 HpNum ReadFiBLambda( cIndex wFi = Inf< Index >() ) override;

 HpNum ReadDt( cHpNum tt = 1 ) override;

 HpNum ReadSigma( cIndex wFi = Inf< Index >() ) override;

 HpNum ReadDStart( cHpNum tt = 1 ) override;

 cLMRow Readd( bool Fulld = false ) override;

 void ReadZ( LMRow tz , cIndex_Set & I , Index & D ,
	     cIndex wFi = Inf< Index >() ) override;

 cHpRow ReadMult( cIndex_Set & I , Index & D , cIndex wFi = Inf< Index >() ,
		  const bool IncldCnst = true ) override;

 HpNum ReadLBMult( cIndex wFi = Inf< Index >() ) override;

 HpNum ReadGid( cIndex Nm = Inf< Index >() ) override;

 void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau ) override;

 void SensitAnals( HpNum &lp , HpNum &cp ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
 *  @{ */

 Index BSize( cIndex wFi = Inf< Index >() ) override;

 Index BCSize( cIndex wFi = Inf< Index >() ) override;

 Index MaxName( cIndex wFi = Inf< Index >() ) override;

 Index WComponent( cIndex i ) override;

 bool IsSubG( cIndex i ) override;

 Index NumNNVars( void ) override;

 Index NumBxdVars( void ) override;

 bool IsNN( cIndex i ) override;

 void CheckIdentical( const bool Chk = true ) override;

 cHpRow ReadLinErr( void ) override;
  
 HpNum ReadLowerBound( cIndex wFi = Inf< Index >() ) override;

 HpNum EpsilonD( void ) override;

 bool FiBLambdaIsExact( cIndex wFi );

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
 *  @{ */

 SgRow GetItem( cIndex wFi = Inf< Index >() ) override;

 void SetItemBse( cIndex_Set SGBse = nullptr , cIndex SGBDm = 0 ) override;

 Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum & Ai , HpNum & ScPri )
  override;

 Index CheckCnst( HpNum & Ai , HpNum & ScPri , cHpRow CrrPnt ) override;

 bool ChangesMPSol( void ) override;

 void SetItem( cIndex Nm = Inf< Index >() ) override;

 void SubstItem( cIndex Nm ) override;

/*--------------------------------------------------------------------------*/

 void RmvItem( cIndex i ) override;

 void RmvItems( void ) override;

/*--------------------------------------------------------------------------*/

 void SetActvSt( cIndex_Set AVrs = nullptr , cIndex AVDm = 0 ) override;

 void AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs ) override;

 void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs ) override;

/*--------------------------------------------------------------------------*/

 void AddVars( cIndex NNwVrs ) override;

 void RmvVars( cIndex_Set whch = nullptr , Index hwmny = 0 ) override;

/*--------------------------------------------------------------------------*/

 void ChgAlfa( cHpRow DeltaAlfa ) override;

 void ChgAlfa( cHpRow NewAlfa , cIndex wFi ) override;

 void ChgAlfa( cIndex i , cHpNum Ai ) override;

/*--------------------------------------------------------------------------*/

 void ChangeCurrPoint( cLMRow DLambda , cHpRow DFi ) override;

 void ChangeCurrPoint( cHpNum Tau , cHpRow DFi ) override;

/*--------------------------------------------------------------------------*/

 void ChgSubG( cIndex strt , Index stp , cIndex wFi ) override;

/*--------------------------------------------------------------------------*/

 void ChgCosts( Index wFi , cLMRow Lambda ) override;

/*--------------------------------------------------------------------------*/
/** This method allows to change the LHS/RHS of the constraints of the "easy"
 * component wFi. It is illegal if it is called with wFi not an "easy"
 * component. The method calls GetBDesc() (with all parameters 0 save these
 * of the LHS/RHS) to retrieve the new LHS/RHS and change the master
 * problem accordingly.
 *
 * IMPORTANT NOTE: the number of rows in the "easy" component must *not*
 *                 have changed.
 *
 * IMPORTANT NOTE: to make the implementation simple, the method only works
 *                 if the "two sided status" of each constraint is preserved
 * by any change. This (in its simplest terms) means that for each constraint,
 * if a LHS/RHS was 0, finite or -INF/+INF before the change it must have
 * remained 0, finite or -INF/+INF; in other words, only finite nonzero
 * LHS/RHS can change. Also, if the LHS and RHS of a row were equal before the
 * change (an equality constraint) they must be equal (to each other, but
 * possibly different from before) after the change.
 *
 * This is rather restrictive but allowing more flexibility would require a
 * substantial rework of the method that is not possible right now.
 *
 * Failure to comply with these rules may result in an exception being thrown,
 * or, even worse, wrong results being reported with no warning about it. */

 void ChgRLHS( Index wFi ) override;

/*--------------------------------------------------------------------------*/
/** This method allows to change the lower and upper bounds of the variables
 * of the "easy" component wFi. It is illegal if it is called with wFi not an
 * "easy" component. The method calls GetBDesc() (with all parameters 0 save
 * these of lower and upper bounds) to retrieve the new lower and upper
 * bounds and change the master problem accordingly.
 *
 * IMPORTANT NOTE: the number of columns in the "easy" component must *not*
 *                 have changed.
 *
 * IMPORTANT NOTE: to make the implementation simple, the method only works
 *                 if the "finite nonzero" status of the bounds is preserved
 * by any change. This (in its simplest terms) means that if the bound was
 * finite and nonzero before the change it must have remained so after the
 * change. That is, bounds that are either 0 or -INF/+INF must remain in one
 * of the two statuses (but they can swap, i.e., 0 can become -INF/+INF and
 * vice/versa), and bounds that are finite and nonzero must remain finite and
 * nonzero.
 *
 * This is rather restrictive but allowing more flexibility would require a
 * substantial rework of the method that is not possible right now.
 *
 * Failure to comply with these rules may result in an exception being thrown,
 * or, even worse, wrong results being reported with no warning about it. */

 void ChgLUBD( Index wFi ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
 *  @{ */

/** Note: the destructor does *not* delete the OSISolverInterface object,
 * this has to be done by whomever created it in the first place. */

 virtual ~OSIMPSolver();

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 void cleanup( void );

/*--------------------------------------------------------------------------*/
/* Update the prices after either t has changed or if the current point
   has been changed. The last form changes the prices in tmpHP[], assumed to
   hold the whole vector of cost coefficients of the master problem. */

 void tUpdatePrices( cIndex strt = 0 , Index stp = Inf< Index >() );

 void ptUpdatePrices( cIndex strt = 0 , Index stp = Inf< Index >() );

 void ptUpdatePricesInPlace( void );

/*--------------------------------------------------------------------------*/

 void UpdateRhoCol( void );

/*--------------------------------------------------------------------------*/
/* Private method for handling the quadratic stabilization. */

 void switchToQP( void );

/*--------------------------------------------------------------------------*/
//  inline bool isactive( Index i );

 void activate( Index i );

 void deactivate( Index i );

 void resizeHP( Index i );

 void resizeI( Index i );

/*--------------------------------------------------------------------------*/
/* Placeholder method: it should do nothing, but depending on the value of
   the macro CHECK_DS [see OSIMPSolver.C] some extra checking is done in
   there for debugging purposes. */

 void CheckDS( void );

/*--------------------------------------------------------------------------*/

 Index CheckBCopy( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 OsiSolverInterface *osiSlvr;  ///< pointer to the OsiSolver
 FiOracle *FIO;                ///< pointer to the FiOracle

 OsiAlg algorithm;          ///< algorithm used in the solver

 bool useactiveset;	    ///< true if the active set are taken in account
 bool first;		    /**< true if it's the first master problem to
			     * solve */
 bool checkID;
 Bool_Vec weasy;            ///< which Fi-component is easy
 Index NrEasy;              ///< how many easy Fi-components are there

 StabFun stab;	            ///< tell the stabilization term adopted
 HpNum t;		    ///< proximal parameter term

 double MaxTime;            ///< max running time for each call to SolveMP
 double OptEps;             ///< optimality tolerance
 double FsbEps;             ///< feasibility tolerance

 double * RhoCol;           ///< entries of the column of rho
 int  * RhoColBse;          ///< ... in sparse format: indices ...
 Index	RhoColBDm;          ///< ... and size
 Index RhoColSgPos;         ///< position in RhoColBse of the "subgradient"
                            /**< position in RhoColBse of the first entry
			     * corresponding to the "subgradient part", i.e.,
  * >= comp_row[ NrFi ]; if there are no easy component then the first part
  * of RhoCol is dense and therefore RhoColSgPos == comp_row[ NrFi ], but if
  * there are easy components some entries of RhoCol before comp_row[ NrFi ]
  * may be 0, and RhoCol[] is stored in sparse format. */

 SgRow NewItem;	            ///< the new item
 Index NewItemFi;	    ///< the function component relative to new item
 cIndex_Set NewItemBse;     ///< index vector of the new item
 Index NewItemBDm;	    ///< dimension of new item

 Index MaxNZ;		    ///< max dimension of the item

 bool NewItemisSG;	    ///< subgradient vs constraint
 double NewItemprice;	    ///< price of the item
 double NewItemScPri;       ///< Gid product

 Index_Set NSubG;           ///< number of subgradients (per component)
 Index_Set NConst;          ///< number of constraints (per component)

 Index_Set comp_row;	    ///< row vocabulary
 Index_Set dict_item;	    ///< item vocabulary
 Index_Set dict_slack;	    ///< slack vocabulary
 Index_Set dict_stab;       ///< stabilization variables vocabulary

 Index_Set wcomp;           /**< Has several different uses:
			     * - if InINF, it means the item is not there
  * - the value *masked by the two most significant bits* is the component
  *   of the function that the item belongs to
  * - the value of the most significant bit is 1 if the item is a
  *   subgradient, 0 if it is a constraint
  * - the value of the second-most significant bit is 1 if the item was
  *   *not* present at the time that the last Master Problem has been
  *   solved, 0 otherwise */

 Index_Set comp_col;	    ///< column vocabulary

 Index item_maxname;        ///< maximum name among current items + 1

 HpRow tempHP;		    ///< temporary vector of HpNum
 int*  tempI;               ///< temporary vector of ints
 Index tempHP_size;         ///< size of tempHP
 Index tempI_size;          ///< size of tempI

 LMRow Upper;		    ///< upper bounds on the primal variable
 LMRow Lower;		    ///< lower bounds on the primal variable

 Index NNVars;              ///< variables with >= 0 constraints
 Index BxdVars;             ///< variables with >= 0 or <= UB constraints

 cIndex_Set Aset;           ///< active set structure
 Index Asetdim;	            ///< size of active set

 HpRow LwrBnds;             ///< values of lower bounds for each component
  
 HpRow GiPerd;              /**< keep the scalar products between each
			     * item and the direction: while it is not used
 * at all iterations, its computation uses information that may change (e.g.
 * if Alphas change) so it's better to save it */

 #if( PRESERVE_OSI_SOLS )
  double * csol;            ///< column (primal) solution
  int csols;                ///< length of csol
  double * rsol;            ///< row (dual) solution
  int rsols;                ///< length of rsol
  double * rcst;            ///< row reduced cost
 #endif

 CoinMessageHandler *derhand;

 Index NRCall;              ///<  number of calls to SolveMP procedure

/*--------------------------------------------------------------------------*/

 };  // end( class OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* OSIMPSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OSIMPSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/
