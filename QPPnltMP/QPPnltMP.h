/*--------------------------------------------------------------------------*/
/*---------------------------- File QPPnltMP.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class QPPenaltyMP implements the Quadratic-Penalty Master Problem for
 * Bundle algorithms, using the specialized (QP) solvers implemented by the
 * [B]MinQuad class [see [B]MinQuad.h and references therein]. The class
 * conforms to the interface defined by the class MPSolver [see MPSolver.h].
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

#ifndef __QPPnltMP
 #define __QPPnltMP  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup QPPnltMP_MACROS Compile-time switches in QPPnltMP.h
    Some relevant switches are defined in MinQuad.h and BMinQuad.h.
    *mandatory switches* (the code won't compile/run otherwise) are

     LAZY_Q        0   if G_IMPLM > 1               (MinQuad.h)
     SIGNAL_MBCHG  0                                (BMinQuad.h)
     SIGNAL_B2CHG  0                                (BMinQuad.h)

    Setting

     TWOSIDED      1                                (BMinQuad.h)

    is required if upper bounds on the variables are present [see GetUB() in
    FiOracle.h]; if they are not, setting TWOSIDED == 0 saves some memory.

    Obviously, switches in BMinQuad.h have no influence if HV_NNVAR == 0.

    Although the compile-time switches are generally thought to be orthogonal
    to each other, certain combinations make no sense and should be avoided.
    Some *strongly advised switches* are

     CONSTR        0   if ADD_CNST == 0             (MinQuad.h)
    @{ */

/*------------------------------ HV_NNVAR ----------------------------------*/
/** If no Multipliers are constrained to be NonNegative (NN), the QPPenaltyMP
   class can inherit directly from the base MinQuad class rather than from
   the derived class BMinQuad, saving time and space; in fact, if ! HV_NNVAR
   the BMinQuad.* files need not to be compiled and linked with the code.
   Obviously, tring to solve a problem with NN Multipliers and HV_NNVAR == 0
   will result in an exception being thrown. */

#define HV_NNVAR 1

/*------------------------------ G_IMPLM -----------------------------------*/
/** The two computationally expensive tasks of the (QP) algorithm are the
   calculation of the entries of the Hessian matrix Q[][] (each time a new
   item is inserted) and the calculation of the entries of the tentative
   direction d[] (at each iteration, and possibly more than once). How these
   two tasks are accomplished is in large part controlled by the switches
   LAZY_Q in MinQuad.h and LAZY_D in BMinQuad.h: however, another important
   issue is how the set of items (Bundle) is implemented.

   There are essentially two ways: either by columns or by rows. A third
   possibility is clearly that of having both representations. Furthermore,
   in the latter case another degree of freedom is available, that is which
   of the two forms have to be used for each task.

   Each choice has its pros and cons, depending on the instance (max dim. of
   the Bundle vs number of Multipliers), on other algorithmic options (wheter
   or not a LVG is used) and on the target architecture (caching etc.). 
   All these issues are controlled by this switch, that has the following
   four possible values:

   0  =>  the Bundle is implemented by columns only, i.e. as a set of
          NumVar-vectors each one being an item;

   1  =>  the Bundle is implemented by rows *and* columns, that is both forms
          are kept: however, the *column-wise* implementation is used for
          calculating the entries of Q[][] while the *row-wise* implementation
          is used for calculating the entries of d[];

   2  =>  as for 1, the Bundle is implemented by rows *and* columns; however,
          in this case the *row-wise* implementation is used for calculating
          the entries of Q[][] while the *column-wise* implementation is used
          for calculating the entries of d[];

   3  =>  the Bundle is implemented by rows only, i.e. as a set of
          (max Bundle dim.)-vectors each one containing one entry of all the
          items: in this case, LAZY_Q *must* be 0 [see MinQuad.h].

   4  =>  as 3, plus memory is allocated and deallocated on the fly to keep
          the number of allocated rows to the bare minimum.

   Clearly, if only the column-wise implementation is available (G_IMPLM == 0)
   or only the row-wise implementation is available (G_IMPLM >= 3), they are
   used for both the tasks. */

#define G_IMPLM 0

/*------------------------------ ADD_CNST ----------------------------------*/
/** If ADD_CNST == 1, general linear constraints on the Lambda space are
   allowed: automatic generation of constraints is done for cases where the
   feasible set is not know a priori [see Fi() and SetGi() / GetGi() below].

   From a "primal" viewpoint, constraints are usually needed when the set X is
   *not* compact, and therefore the Lagrangean subproblem can be *unbounded*:
   in this case, constraints correspond to *extreme rays* of X. */

#define ADD_CNST 0

/*---------------------------- RECOMPUTE_GTZ -------------------------------*/
/** The scalar products between the subgradients G_i and z* are crucial
 * because they are used to update the linearization errors during a Serious
 * Step. This means that numerical errors accumulate into linearization
 * errors, and there is no way in which they can be checked for their true
 * value and the errors corrected. This may lead to them being very far off
 * their true value, which may impair the correctness of the results of the
 * bundle algorithm.
 *
 * Although the flaw is inherent in the process of continuously updating the
 * linearization errors, it can be minimized by computing the scalar
 * products as accurately as possible. However, in MinQuad they are also
 * computed "indirectly" during the optimality checks, and for efficiency
 * one would typically use these already available values. Since they can
 * be affected to errors, the RECOMPUTE_GTZ switch is available to allow to
 * ignore them and rather use freshly computed scalar products, trading off
 * a higher computational cost for (hopefully) better accuracy.
 *
 * Note that the re-computation is only done when it actually matters, i.e.,
 * right before a Serious Step is performed. */

#define RECOMPUTE_GTZ 0

/** @} end( group( QPPnltMP_MACROS ) ) */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MPSolver.h"

#if( HV_NNVAR )
 #include "BMinQuad.h"
#else
 #include "MinQuad.h"

 #define TWOSIDED 0
#endif

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{
 using namespace MinQuad_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/

class QPPenaltyMP : public MPSolver ,
                    #if( HV_NNVAR )
                     private BMinQuad
                    #else
                     private MinQuad
                    #endif
{

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
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

  QPPenaltyMP( std::istream *iStrm = 0 );

/** Constructor of the class. The parameter `iStrm', if provided, is taken as
   a pointer to a istream from which the algorithmic parameters for the QP
   solver are sequentially read in the following order.  Each parameter
   must be placed at the beginning of a separate line, max 128 characters
   long, with all the rest of the line up to the first newline character
   '\n' (apart from a separating whitespace) being available for comments. 
   If NULL is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

    HpNum CtOff [.01] The "break" value for the pricing in the base class
                MinQuad [see SetPricing() in MinQuad.h]: positive values are
		passed untouched, while any negative value is turned into
		Inf< HpNum >().

    Index MxAdd [0] How many variables can be "moved" at each iteration of the
    Index MxRmv [0] "upper-level method" for solving bound-constrained MPs
                implemented into the base class BMinQuad (if any, see
		HV_NNVAR above, otherwise the parameters are just ignored);
		see SetMaxVar[Add/Rmv]() in BMinQuad.h. If 0 [the default]
		is used, no bound is given. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void SetDim( cIndex MxBSz = 0 , FiOracle *Oracle = 0 , 
		const bool UsAvSt = false ) override;
    
/*--------------------------------------------------------------------------*/
    // "publicize" the MinQuad::SetPricing method
    
    using MinQuad::SetPricing;
    
/*--------------------------------------------------------------------------*/
    
#if( HV_NNVAR )
    // "publicize" the BMinQuad::SetMaxVarAdd method
    using BMinQuad::SetMaxVarAdd;
    
    // "publicize" the BMinQuad::SetMaxVarRmv method
    using BMinQuad::SetMaxVarRmv;
#else
    // provide a dummy implementation of BMinQuad::SetMaxVarAdd
    inline void SetMaxVarAdd( cIndex MVA = 1 ) { }
    
    // provide a dummy implementation of BMinQuad::SetMaxVarRmv
    inline void SetMaxVarRmv( cIndex MVR = 1 ) { }
#endif

/*--------------------------------------------------------------------------*/

   void Sett( cHpNum tt = 1 ) override { tCurr = tt; }

/*--------------------------------------------------------------------------*/

   void SetPar( const int wp , cHpNum value ) override
   {
    switch( wp ) {
     case( kMaxTme ): SetMaxTime( value ); break;
     case( kOptEps ): SetEpsilonR( value ); break;
     case( kFsbEps ):
      #if( HV_NNVAR )
       SetEpsilonD( value );
      #else
       EpsD = value;
      #endif
      break;
      // case( kZero ):  SetMinFVal( value ); break;
     default: throw( NDOException( "SetPar( HpNum ): unknown parameter" ) );
     }
    }  // end( QPPenaltyMP::SetPar( HpNum ) )

/*--------------------------------------------------------------------------*/
/** Setting an individual lower bound on the component 1 is not supported
 * unless RHSp == nullptr, i.e., there is no component 0; in this case
 * clearly the lower bound on component 1 is a global lower bound for all
 * the model, which can be set. */

   void SetLowerBound( cHpNum LwBnd = - Inf< HpNum >() ,
		       cIndex wFi = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

   void CheckIdentical( const bool Chk = true ) override { ChkIde = Chk; }

/*--------------------------------------------------------------------------*/
/** The meaning of lvl is:

    0  =>  no log at all (also assumed if log = 0), except if errors
           are found in the Check****() routines;

    1  =>  "basic" log: only the errors are reported;

    2  =>  a more detailed log where also non-necessarily erroneous but
           "strange" situations are reported;

   Note that most of the log is demanded to the ancestor classes [B]MinQuad,
   that actually solve the problem. They will use outs as the log stream,
   but their "verbosity" is rather set by the compile-time switches LOG_MQ
   (in MinQuad.h) and LOG_BMQ (in BMinQuad.h). */

   void SetMPLog( std::ostream *outs = 0 , const char lvl = 0 ) override
   {
    if( ( MPLog = outs ) )
     MPLLvl = lvl;
    else
     MPLLvl = 0;

    #if( LOG_MQ )
     if( outs )
      SetMQLog( outs );
    #endif
    #if( HV_NNVAR )
     #if( LOG_BMQ )
      if( outs )
       SetBMQLog( outs );
     #endif
    #endif
    }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   MPStatus SolveMP( void ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   HpNum ReadFiBLambda( cIndex wFi = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

   HpNum ReadDt( cHpNum tt = 1 ) override
   {
    return( ( ( tt == tCurr ? tt : ( tt * tt ) / tCurr ) / 2 ) * ReadzNorm() );
    }

/*--------------------------------------------------------------------------*/

   HpNum ReadSigma( cIndex wFi = Inf< Index >() ) override
   {
    #if( HV_NNVAR )
     return( BMinQuad::ReadSigma() );
    #else
     return( MinQuad::ReadSigma() );
    #endif
    }

/*--------------------------------------------------------------------------*/

   HpNum ReadDStart( cHpNum tt = 1 ) override
   {
    return( ( tt / 2 ) * ReadzNorm() );
    }

/*--------------------------------------------------------------------------*/

   cLMRow Readd( bool Fulld = false ) override;

/*--------------------------------------------------------------------------*/

   void ReadZ( LMRow tz , cIndex_Set &I , Index &D ,
	       cIndex wFi = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

   cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf< Index >() ,
		    const bool IncldCnst = false ) override;

/*--------------------------------------------------------------------------*/

   HpNum ReadLBMult( cIndex wFi = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

   HpNum ReadGid( cIndex Nm = Inf< Index >() ) override;

/*--------------------------------------------------------------------------*/

   void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau ) override;

/*--------------------------------------------------------------------------*/

   void SensitAnals( HpNum &lp , HpNum &cp ) override;

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   Index BSize( cIndex wFi = Inf< Index >() ) override { return( ActBDim ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index BCSize( cIndex wFi = Inf< Index >() ) override { return( ActCNum ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index MaxName( cIndex wFi = Inf< Index >() ) override
   {
    return( MaxItemN() );
    }

/*--------------------------------------------------------------------------*/

   Index WComponent( cIndex i ) override
   {
    return( ( ( i < MaxItemN() ) && IsThere( i ) ) ? 1 : Inf< Index >() );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   bool IsSubG( cIndex i ) override
   {
    return( ( i < MaxItemN() ? IsASubG( i ) : false ) );
    }

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   Index NumNNVars( void ) override { return( NNVars ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index NumBxdVars( void ) override
   {
    #if( TWOSIDED )
     return( BxdVars );
    #else
     return( NNVars );
    #endif
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   bool IsNN( cIndex i ) override
   {
    return( LowerBounds()[ i ] > - Inf< HpNum >() );
    }

#endif

/*--------------------------------------------------------------------------*/

   cHpRow ReadLinErr( void ) override { return( ReadAlfa() ); }

/*--------------------------------------------------------------------------*/

   HpNum ReadLowerBound( cIndex wFi = Inf< Index >() ) override {
    if( ! wFi )
     throw( NDOException( "QPPenaltyMP::ReadLowerBound: wFi == 0" ) );

    if( ( wFi > 1 ) || ( ! RHSp ) )
     return( ReadLB() );

    return( - Inf< HpNum >() );
    }
    
/*--------------------------------------------------------------------------*/

   Index ReadGi( cIndex Nm , SgRow Gi , cIndex_Set &GB );

/*--------------------------------------------------------------------------*/

   HpNum EpsilonD( void ) override
   {
    #if( HV_NNVAR )
     return( BMinQuad::EpsilonD() );
    #else
     return( EpsD * ( BDim ? BDim : 1 ) );
    #endif
    }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   SgRow GetItem( cIndex wFi = Inf< Index >() ) override;

   void SetItemBse( cIndex_Set SGBse = 0 , cIndex SGBDm = 0 ) override;

   Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai , HpNum &ScPri )
    override;

   Index CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt ) override;

   bool ChangesMPSol( void ) override;

   void SetItem( cIndex Nm = Inf< Index >() ) override;

   void SubstItem( cIndex Nm ) override;

/*--------------------------------------------------------------------------*/

   void RmvItem( cIndex i ) override;

   void RmvItems( void ) override;

/*--------------------------------------------------------------------------*/

   void SetActvSt( cIndex_Set AVrs = 0 , cIndex AVDm = 0 ) override;

/*--------------------------------------------------------------------------*/

   void AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs ) override;

   void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs ) override;

/*--------------------------------------------------------------------------*/

   void AddVars( cIndex NNwVrs ) override;

   void RmvVars( cIndex_Set whch = 0 , Index hwmny = 0 ) override;

/*--------------------------------------------------------------------------*/

   void ChgAlfa( cHpRow DeltaAlfa ) override;

   void ChgAlfa( cHpRow NewAlfa , cIndex wFi ) override;

   void ChgAlfa( cIndex i , cHpNum Ai ) override;

/*--------------------------------------------------------------------------*/

   void ChangeCurrPoint( cLMRow DLambda ,  cHpRow DFi ) override;

   void ChangeCurrPoint( cHpNum Tau ,  cHpRow DFi ) override;

/*--------------------------------------------------------------------------*/

   void ChgSubG( cIndex strt , Index stp , cIndex wFi ) override;

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~QPPenaltyMP();   /// just the destructor

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

/*--------------------------------------------------------------------------*/
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

   void ComputeRHSxd( void );

/*--------------------------------------------------------------------------*/

   void Computedir( bool Fulld );

/*--------------------------------------------------------------------------*/

   Index CheckBCopy( void );

/*--------------------------------------------------------------------------*/

#if( LAZY_Q )

   QuNum GiTGj( cIndex i , cIndex j ) override;

#else

   void GiTG( cIndex i , QuRow Qi , cIndex iMax ) override;

#endif

/* The "standard" interface with the (QP) solver, common to both the MinQuad
   and BMinQuad classes. */

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   cSgRow GiTilde( cIndex i ) override;

/* This method of the BMinQuad-specific part of the interface with the (QP)
   solver returns a pointer to the i-th row of the matrix of subgradients:
   if such a row is not available for free, GiTilde() constructs it in some
   temporary memory (the vector TSGk) and returns a pointer to it. */

#endif

/*--------------------------------------------------------------------------*/

   LMNum CalculateZ( cIndex h );

   void CalculateZ( cIndex_Set Wh , LMRow z );

#if( HV_NNVAR )

   void CalculateZ( LMRow z ) override;

#else

   void CalculateZ( LMRow z );

#endif

/* These methods are used to construct the tentative direction. If HV_NNVAR,
   one of the CalculateZ()'s (depending on LAZY_D) is called by the BMinQuad
   class at each iteration; otherwise they are called by the methods of the
   QPPenaltyMP class when the direction needs to be computed.

   If ActvVrs != NULL, CalculateZ( z ) only calculates the entries of z
   corresponding to *existing* variables, those whose names are in ActvVrs[].
   The other two variants can also be used *only* to calculate entries of z
   corresponding to existing variables. */

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   LMNum GiTLB( cIndex h , cLMRow l , cIndex_Set lB , cIndex lBd ) override;

   void GiTLB( HpRow gtlb , cLMRow l , cIndex_Set lB , cIndex lBd ,
	       const bool add ) override;

/* These methods of the BMinQuad-specific part of the interface with the (QP)
   solver needs to be defined only if HV_NNVAR > 0. */

#endif

/*--------------------------------------------------------------------------*/

#if( ( G_IMPLM > 0 ) && ( G_IMPLM < 3 ) )

   void GetNewGTRow( cIndex i );

   void SetGTRow( cIndex i );

   void PutGTRow( cIndex i );

/* GetNewGTRow() and PutGTRow() can *only* be called if ActvVrs != NULL, i.e.,
   Buffer[] is defined. */

#endif

/*--------------------------------------------------------------------------*/

#if( G_IMPLM < 1 )

   SgRow TSGi( cIndex i );

#endif

/*--------------------------------------------------------------------------*/

   void FullZ( LMRow z , cIndex_Set tB , cHpRow tM );

/* Construct the full Z into z, using tB and tM as the Base and Mult. */

/*--------------------------------------------------------------------------*/

   void CompleteZ( bool Fullz = true );

/* If necessary, complete a partly calculated z*. If  Fullz == true a "full"
   z* is computed, otherwise only the *existing* entries of z* are computed.
   The thing is set in the internal data structures of BMinQuad if HV_NNVAR,
   and in the field dir[] of QPPenaltyMP otherwise. */

/*--------------------------------------------------------------------------*/

#if( ADD_CNST )

   void ComputeZBase( void );

/* Extract from the current base the part relative to constraints and store
   it in ZBase/ZMult. */

#endif

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   void CheckBounds( void );

#endif

/*--------------------------------------------------------------------------*/

   void ChgRHS( cIndex strt , cIndex stp );

/* Changes the current RHS to the new one contained in tmpG1. Only the
   entries with indices between strt (included) and stp (excluded) have
   changed w.r.t. the previous value; it is assumed that strt >= 0 and
   stp <= CrrSGLen.

   tmpG1 is assumed to contain the new values of the RHS in "dense" format,
   i.e., tmpG1[ i ] contains the RHS of variable strt + i for 0 <= i <
   stp - strt.

   Also, G1Perz is assumed to be zero if the new values of the RHS are all
   zero (this is done automatically SetItemBse()). */

/*--------------------------------------------------------------------------*/

   void CheckSG( void );

/*--------------------------------------------------------------------------*/

   void CheckGT( void );

/*--------------------------------------------------------------------------*/

   void MemDealloc( void );

/* Deallocates all the memory of the class. */

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  FiOracle *CrrOrcl;   // the current oracle

  Index DimMinQuad;    // dimension of the Bundle "seen" by [B]MinQuad

  HpNum MinNewAlfa;    // the min. among the Alfa's of "new" items
  HpNum tCurr;         // current value of t
  HpNum Alfa1;         // linearization error of the newly inserted item:
                       // Alfa1 == Inf< HpNum >() signals that the item is
                       // the RHS (0-th component)
  HpNum G1Perz;        // scalar product between the newly inserted item and z
  #if( ADD_CNST )
   bool G1IsSubg;      // true if the newly inserted item is a subgradient
  #endif

  LMRow dir;           // when HV_NNVAR > 0, z* is contained into the internal
                       // data structure of the BMinQuad class, and dir[] is
                       // used as the temporary for SettmpD() and to hold the
                       // actual direction; when HV_NNVAR == 0 instead, dir[]
                       // holds either z* or d*, since only the scaling factor
                       // "- t" distinguishes the two

  char ZCmptd;         // 0 = the (QP) has not been (correctly) solved
                       // 1 = no entries of z* have been calculated
                       // 2 = only the "defined" entries for NN variables
                       //     (only happens with HV_NNVAR && LAZY_D == 2)
                       // 3 = only the "defined" entries for *all* variables
                       // 4 = *all* the entries
                       // note that if HV_NNVAR > 0 then the entries are
                       // written and held into the internal data structure
                       // of the BMinQuad class, otherwise into dir[]

  char dirCmptd;       // tells the status of dir[]:
                       // 0 = no entries of dir[] have been calculated
                       // 1 = only the "defined" entries (in ActvVrs[])
                       // 2 = *all* the entries have been computed

  SgRow RHSp;          // the Right Hand Side of the variables
  HpNum RHSxd;         // contains ( - 1 / t ) * scalar product between the
                       // RHS and d*; if there are no "active bounds" (which
                       // is true in particular if HV_NNVAR == 0) this is
                       // the scalar product between the RHS and z*,
                       // otherwise the two numbers are different

  cIndex_Set ActvVrs;  // the "active set" of variables
  Index ActvVDim;      // size of the active set

  #if( HV_NNVAR )
   Index NNVars;       // number of the variables that are constrained in sign
   #if( TWOSIDED )
    Index BxdVars;     // number of variables that have either an upper or
                       // a lower bound
   #endif
  #else
   HpNum EpsD;         // precision in checking feasibility
  #endif

  #if( G_IMPLM < 3 )
   SgMat SubG;         // the bundle (itself)

   #if( G_IMPLM < 1 )
    SgRow TSGk;        // temporary: a row in the bundle
   #endif
  #else
   Index Insrtd;       // the "name" of G1 (see below)
  #endif

  SgRow tmpG1;         // temporary: the current subgradient
  
  #if( G_IMPLM > 0 )
   SgMat TSubG;        // Transpose of G (GT)
   Index NrAllRws;     // how many rows of GT have been allocated so far

   #if( ( G_IMPLM > 0 ) && ( G_IMPLM < 3 ) )
    SgMat Buffer;      // Temporary buffer of rows of G
    Index BFCntr;      // current size of Buffer
   #endif
  #endif

  #if( ADD_CNST )
   Index_Set ZBase;    // the "base" for doing aggregation
   HpRow ZMult;        // its multipliers
   Index ZBDim;        // its size
   Index CBZDim;       // number of constraints in the last optimal Base
  #endif

  bool ChkIde;         // true if checks for discovering and discarding
                       // identical items are performed

/*--------------------------------------------------------------------------*/

 };  // end( class QPPenaltyMP )

/*--------------------------------------------------------------------------*/

 }  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* QPPnltMP.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File QPPnltMP.h ----------------------------*/
/*--------------------------------------------------------------------------*/
